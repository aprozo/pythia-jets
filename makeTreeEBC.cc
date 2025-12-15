// main.cc - MC Dijet Embedding in Pythia Background
#include "Pythia8/Pythia.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

using namespace Pythia8;

// ============================================================================
// Configuration Parameters (all in one place)
// ============================================================================
struct Config {
  // Pythia background settings
  double pythia_ptHatMin = 3;
  double pythia_ptHatMax = -1; // low-pT soft events

  // Particle-level cuts
  double part_ptMin = 0.15;
  double part_etaMax = 1.0;

  // Jet reconstruction
  double jet_R = 0.4;
  double jet_ptMin = 1.0;
  double jet_etaMax = part_etaMax - jet_R; // for reconstructed jets

  // MC dijet generation
  double mc_jet_ptMin = 1.0;
  double mc_jet_ptMax = 50.0;
  double mc_jet_ptPower = -5.0; // power law: dN/dpT ~ pT^power
  double mc_jet_etaMax =
      part_etaMax - jet_R; // truth jet eta range (leave room for particles)
  double mc_jet_deltaPhiSmear = 0.1; // back-to-back smearing (radians)

  int mc_nPart_min = 2;
  int mc_nPart_max = 20;
  double mc_particle_deltaR = 0.15; // spread of particles around jet axis

  // Background transverse regions
  // double bkg_deltaPhi_min = M_PI / 3.0;       // π/3
  // double bkg_deltaPhi_max = 2.0 * M_PI / 3.0; // 2π/3

  // Event generation
  int nEvents = 10000;

  std::string outputPrefix = "embedded_dijets";
};

// ============================================================================
// Utility Functions
// ============================================================================

inline double wrapPhi(double phi) {
  while (phi <= -M_PI)
    phi += 2.0 * M_PI;
  while (phi > M_PI)
    phi -= 2.0 * M_PI;
  return phi;
}
double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = wrapPhi(phi1 - phi2);
  return std::sqrt(deta * deta + dphi * dphi);
}

// ============================================================================
// MC Particle Structure
// ============================================================================
struct MCParticle {
  double px, py, pz, E;
  double pt, eta, phi;
  int jetID; // 1 or 2 for MC jets, 0 for Pythia background

  MCParticle() : px(0), py(0), pz(0), E(0), pt(0), eta(0), phi(0), jetID(0) {}

  MCParticle(double pt_, double eta_, double phi_)
      : pt(pt_), eta(eta_), phi(phi_), jetID(0) {
    px = pt * std::cos(phi);
    py = pt * std::sin(phi);
    pz = pt * std::sinh(eta);
    // use pion mass
    E = std::sqrt(px * px + py * py + pz * pz);
  }

  fastjet::PseudoJet toPseudoJet() const {
    return fastjet::PseudoJet(px, py, pz, E);
  }
};

// ============================================================================
// MC Dijet Structure
// ============================================================================
struct MCDijet {
  // Truth jet kinematics
  double jet1_pt, jet1_eta, jet1_phi;
  int jet1_nPart;
  double jet2_pt, jet2_eta, jet2_phi;
  int jet2_nPart;

  // Truth particles
  std::vector<MCParticle> jet1_particles;
  std::vector<MCParticle> jet2_particles;

  MCDijet()
      : jet1_pt(0), jet1_eta(0), jet1_phi(0), jet1_nPart(0), jet2_pt(0),
        jet2_eta(0), jet2_phi(0), jet2_nPart(0) {}
};

// ============================================================================
// MC Dijet Generator
// ============================================================================
class MCDijetGenerator {
private:
  const Config &cfg;
  TF1 ptSpectrum;
  TH1D *qaParticlePt;
  TH1D *qaParticlePhi;
  TH1D *qaParticleEta;
  TH2D *qaJetPt1VsJetPt2;
  TH2D *qaJetPhi1VsJetPhi2;
  TH2D *qaJetEta1VsJetEta2;
  TH2D *qaJetN1VsJetN2;

public:
  MCDijetGenerator(const Config &config)
      : cfg(config), ptSpectrum("ptSpec", "pow(x, [0])", config.mc_jet_ptMin,
                                config.mc_jet_ptMax) {
    ptSpectrum.SetParameter(0, config.mc_jet_ptPower);
    qaParticlePt =
        new TH1D("qaParticlePt", "Particle pT; p_{T} (GeV/c); Counts", 100, 0,
                 config.mc_jet_ptMax);
    qaParticlePhi = new TH1D("qaParticlePhi", "Particle Phi; #phi;Counts", 100,
                             -M_PI, M_PI);
    qaParticleEta = new TH1D("qaParticleEta", "Particle Eta; #eta;Counts", 100,
                             -config.part_etaMax, config.part_etaMax);
    qaJetPt1VsJetPt2 =
        new TH2D("qaJetPt1VsJetPt2",
                 "Jet1 pT vs Jet2 pT; p_{T,1} (GeV/c); p_{T,2} (GeV/c)", 100, 0,
                 config.mc_jet_ptMax, 100, 0, config.mc_jet_ptMax);
    qaJetPhi1VsJetPhi2 = new TH2D("qaJetPhi1VsJetPhi2",
                                  "Jet1 Phi vs Jet2 Phi; #phi_{1}; #phi_{2}",
                                  100, -M_PI, M_PI, 100, -M_PI, M_PI);
    qaJetEta1VsJetEta2 = new TH2D("qaJetEta1VsJetEta2",
                                  "Jet1 Eta vs Jet2 Eta; #eta_{1}; #eta_{2}",
                                  100, -config.part_etaMax, config.part_etaMax,
                                  100, -config.part_etaMax, config.part_etaMax);
    qaJetN1VsJetN2 = new TH2D(
        "qaJetN1VsJetN2", "Jet1 N vs Jet2 N; N_{1}; N_{2}", config.mc_nPart_max,
        0, config.mc_nPart_max, config.mc_nPart_max, 0, config.mc_nPart_max);
  }

  void writeQaHists(TFile *outFile) {
    outFile->cd();
    qaParticlePt->Write();
    qaParticlePhi->Write();
    qaParticleEta->Write();
    qaJetPt1VsJetPt2->Write();
    qaJetPhi1VsJetPhi2->Write();
    qaJetEta1VsJetEta2->Write();
    qaJetN1VsJetN2->Write();
  }

  MCDijet generate() {
    MCDijet dijet;

    // Sample jet1 kinematics
    dijet.jet1_pt = ptSpectrum.GetRandom();
    dijet.jet1_eta = gRandom->Uniform(-cfg.mc_jet_etaMax, cfg.mc_jet_etaMax);
    dijet.jet1_phi = gRandom->Uniform(-M_PI, M_PI);

    // Sample jet2 kinematics (back-to-back with smearing)

    const double relSmear = abs(gRandom->Gaus(0.0, 0.1));
    dijet.jet2_pt = dijet.jet1_pt * std::max(0.1, 1.0 - relSmear);

    dijet.jet2_eta = gRandom->Uniform(-cfg.mc_jet_etaMax, cfg.mc_jet_etaMax);
    dijet.jet2_phi =
        dijet.jet1_phi + M_PI + gRandom->Gaus(0, cfg.mc_jet_deltaPhiSmear);
    dijet.jet2_phi = wrapPhi(dijet.jet2_phi);
    // KEY CORRELATION: Same number of particles in both jets
    int nPart = gRandom->Integer(cfg.mc_nPart_max - cfg.mc_nPart_min + 1) +
                cfg.mc_nPart_min;
    dijet.jet1_nPart = nPart;
    dijet.jet2_nPart = nPart;

    qaJetPt1VsJetPt2->Fill(dijet.jet1_pt, dijet.jet2_pt);
    qaJetPhi1VsJetPhi2->Fill(dijet.jet1_phi, dijet.jet2_phi);
    qaJetEta1VsJetEta2->Fill(dijet.jet1_eta, dijet.jet2_eta);
    qaJetN1VsJetN2->Fill(dijet.jet1_nPart, dijet.jet2_nPart);

    // Generate particles for jet 1
    dijet.jet1_particles = generateJetParticles(dijet.jet1_pt, dijet.jet1_eta,
                                                dijet.jet1_phi, nPart, 1);

    // Generate particles for jet 2
    dijet.jet2_particles = generateJetParticles(dijet.jet2_pt, dijet.jet2_eta,
                                                dijet.jet2_phi, nPart, 2);

    return dijet;
  }

private:
  std::vector<MCParticle> generateJetParticles(double jetPt, double jetEta,
                                               double jetPhi, int nPart,
                                               int jetID) {

    std::vector<MCParticle> particles;
    std::vector<double> particle_pts;
    double sumW = 0;
    // Sample particle pTs (will rescale to match jet pT)
    for (int i = 0; i < nPart; ++i) {
      double w = gRandom->Exp(1.0); // exponential distribution
      particle_pts.push_back(w);
      sumW += w;
    }
    // Rescale so sum matches jet pT (approximately)
    for (int i = 0; i < nPart; ++i) {
      particle_pts[i] *= (jetPt / sumW);
    }

    // Generate particle positions around jet axis
    for (int i = 0; i < nPart; ++i) {
      double dr = 0.0;
      do {
        dr = std::fabs(gRandom->Gaus(0.0, cfg.mc_particle_deltaR));
      } while (dr > cfg.jet_R);

      const double alpha = gRandom->Uniform(0.0, 2.0 * M_PI);

      const double dEta = dr * std::cos(alpha);
      const double dPhi = dr * std::sin(alpha);

      double part_eta = jetEta + dEta;
      double part_phi = wrapPhi(jetPhi + dPhi);

      MCParticle p(particle_pts[i], part_eta, part_phi);
      p.jetID = jetID;
      particles.push_back(p);

      qaParticlePt->Fill(p.pt);
      qaParticleEta->Fill(p.eta);
      qaParticlePhi->Fill(p.phi);
    }

    return particles;
  }
};

// ============================================================================
// Tracking efficiency
// ============================================================================
bool isAcceptedTrack(double pt, TF1 &eff) {
  // if (pt > 30.0)
  //   return false;
  return gRandom->Rndm() < eff.Eval(pt);
}

// ============================================================================
// Background Analysis
// ============================================================================
struct BackgroundResult {
  int multA; // multiplicity in transverse region A
  int multB; // multiplicity in transverse region B

  BackgroundResult() : multA(0), multB(0) {}
};

BackgroundResult
analyzeBackground(const std::vector<fastjet::PseudoJet> &jets,
                  const std::vector<fastjet::PseudoJet> &particles,
                  const Config &cfg) {

  BackgroundResult result;

  if (jets.size() < 2)
    return result;

  // Define dijet axis as average of two leading jets

  fastjet::PseudoJet dijet = jets[0] + jets[1];
  double phi_dijet = dijet.phi(); // in [-pi, pi)
  double eta_dijet = dijet.eta();

  double phiA = wrapPhi(phi_dijet + M_PI / 2);
  double phiB = wrapPhi(phi_dijet - M_PI / 2);

  for (const auto &p : particles) {
    const double pt = p.pt();
    if (pt < cfg.part_ptMin)
      continue;

    const double eta = p.eta();
    if (std::abs(eta) > cfg.part_etaMax)
      continue;

    const double phi = p.phi();
    if (deltaR(eta, phi, eta_dijet, phiA) < cfg.jet_R)
      result.multA++;
    else if (deltaR(eta, phi, eta_dijet, phiB) < cfg.jet_R)
      result.multB++;
  }
  return result;
}

// ============================================================================
// Main Program
// ============================================================================
int main(int argc, char *argv[]) {

  Config cfg;

  // Parse command line arguments
  if (argc > 1)
    cfg.nEvents = std::atoi(argv[1]);
  if (argc > 2)
    cfg.outputPrefix = argv[2];

  std::cout << "\n=== MC Dijet Embedding Configuration ===\n";
  std::cout << "Events: " << cfg.nEvents << "\n";
  std::cout << "Pythia pThat: [" << cfg.pythia_ptHatMin << ", "
            << cfg.pythia_ptHatMax << "] GeV\n";
  std::cout << "MC jet pT: [" << cfg.mc_jet_ptMin << ", " << cfg.mc_jet_ptMax
            << "] GeV\n";
  std::cout << "MC jet eta: ±" << cfg.mc_jet_etaMax << "\n";
  std::cout << "MC particles per jet: [" << cfg.mc_nPart_min << ", "
            << cfg.mc_nPart_max << "]\n";
  std::cout << "Jet radius R: " << cfg.jet_R << "\n";
  std::cout << "========================================\n\n";

  // Tracking efficiency function
  TF1 eff("eff", "[0]*(1-exp(-pow(x/[1],[2])))", 0, 30);
  eff.SetParameters(0.88, 0.25, 1.2);

  // =========================================================================
  // Setup Pythia for soft background
  // =========================================================================
  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200.");
  // pythia.readString("HardQCD:all = on");
  pythia.readString("SoftQCD:nonDiffractive = on");

  pythia.readString("PhaseSpace:pTHatMin = " +
                    std::to_string(cfg.pythia_ptHatMin));
  pythia.readString("PhaseSpace:pTHatMax = " +
                    std::to_string(cfg.pythia_ptHatMax));

  if (!pythia.init()) {
    std::cerr << "[ERROR] Pythia initialization failed.\n";
    return 1;
  }

  // =========================================================================
  // Setup MC Dijet Generator
  // =========================================================================
  MCDijetGenerator mcGen(cfg);

  // =========================================================================
  // Setup FastJet
  // =========================================================================
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfg.jet_R);

  // =========================================================================
  // Setup ROOT output
  // =========================================================================
  std::string outFile = cfg.outputPrefix + ".root";
  TFile *fout = new TFile(outFile.c_str(), "RECREATE");

  // Tree 1: Truth MC dijet
  TTree *tMC = new TTree("mcTruth", "MC Truth Dijets");
  int eventID;
  double mc_jet1_pt, mc_jet1_eta, mc_jet1_phi;
  int mc_jet1_nPart;
  double mc_jet2_pt, mc_jet2_eta, mc_jet2_phi;
  int mc_jet2_nPart;

  tMC->Branch("eventID", &eventID);
  tMC->Branch("mc_jet1_pt", &mc_jet1_pt);
  tMC->Branch("mc_jet1_eta", &mc_jet1_eta);
  tMC->Branch("mc_jet1_phi", &mc_jet1_phi);
  tMC->Branch("mc_jet1_nPart", &mc_jet1_nPart);
  tMC->Branch("mc_jet2_pt", &mc_jet2_pt);
  tMC->Branch("mc_jet2_eta", &mc_jet2_eta);
  tMC->Branch("mc_jet2_phi", &mc_jet2_phi);
  tMC->Branch("mc_jet2_nPart", &mc_jet2_nPart);

  // Tree 2: Reconstructed dijets
  TTree *tReco = new TTree("recoJets", "Reconstructed Jets with Matching");
  int reco_eventID;
  bool matched1, matched2;
  double reco_jet1_pt, reco_jet1_eta, reco_jet1_phi;
  int reco_jet1_nPart;
  double reco_jet2_pt, reco_jet2_eta, reco_jet2_phi;
  int reco_jet2_nPart;
  double dR_mc1_reco1, dR_mc2_reco2;
  int bkg_multA, bkg_multB;

  tReco->Branch("eventID", &reco_eventID);
  tReco->Branch("matched1", &matched1);
  tReco->Branch("matched2", &matched2);
  tReco->Branch("mc_jet1_pt", &mc_jet1_pt);
  tReco->Branch("mc_jet1_eta", &mc_jet1_eta);
  tReco->Branch("mc_jet1_phi", &mc_jet1_phi);
  tReco->Branch("mc_jet1_nPart", &mc_jet1_nPart);
  tReco->Branch("reco_jet1_pt", &reco_jet1_pt);
  tReco->Branch("reco_jet1_eta", &reco_jet1_eta);
  tReco->Branch("reco_jet1_phi", &reco_jet1_phi);
  tReco->Branch("reco_jet1_nPart", &reco_jet1_nPart);
  tReco->Branch("dR_mc1_reco1", &dR_mc1_reco1);
  tReco->Branch("mc_jet2_pt", &mc_jet2_pt);
  tReco->Branch("mc_jet2_eta", &mc_jet2_eta);
  tReco->Branch("mc_jet2_phi", &mc_jet2_phi);
  tReco->Branch("mc_jet2_nPart", &mc_jet2_nPart);
  tReco->Branch("reco_jet2_pt", &reco_jet2_pt);
  tReco->Branch("reco_jet2_eta", &reco_jet2_eta);
  tReco->Branch("reco_jet2_phi", &reco_jet2_phi);
  tReco->Branch("reco_jet2_nPart", &reco_jet2_nPart);
  tReco->Branch("dR_mc2_reco2", &dR_mc2_reco2);

  tReco->Branch("multA", &bkg_multA);
  tReco->Branch("multB", &bkg_multB);

  // =========================================================================
  // Event Loop
  // =========================================================================
  int nProcessed = 0;
  int nMatched = 0;

  std::cout << "Generating " << cfg.nEvents << " events...\n";

  for (int iEvent = 0; iEvent < cfg.nEvents; ++iEvent) {

    if ((iEvent + 1) % 1000 == 0) {
      std::cout << "  Processed " << (iEvent + 1) << " events...\r"
                << std::flush;
    }

    eventID = iEvent;

    // -----------------------------------------------------------------------
    // Step 1: Generate MC dijet (independent of Pythia)
    // -----------------------------------------------------------------------
    MCDijet mcDijet = mcGen.generate();

    // Store truth info
    mc_jet1_pt = mcDijet.jet1_pt;
    mc_jet1_eta = mcDijet.jet1_eta;
    mc_jet1_phi = mcDijet.jet1_phi;
    mc_jet1_nPart = mcDijet.jet1_nPart;

    mc_jet2_pt = mcDijet.jet2_pt;
    mc_jet2_eta = mcDijet.jet2_eta;
    mc_jet2_phi = mcDijet.jet2_phi;
    mc_jet2_nPart = mcDijet.jet2_nPart;

    // -----------------------------------------------------------------------
    // Step 2: Generate Pythia soft background
    // -----------------------------------------------------------------------
    if (!pythia.next())
      continue;

    tMC->Fill(); // fill after pythia generated event

    std::vector<fastjet::PseudoJet> allParticles;

    // Extract Pythia particles
    for (int i = 0; i < pythia.event.size(); ++i) {
      const auto &p = pythia.event[i];
      if (!p.isFinal() || !p.isVisible() || !p.isCharged())
        continue;
      if (std::abs(p.eta()) > cfg.part_etaMax)
        continue;
      if (p.pT() < cfg.part_ptMin)
        continue;
      // if (!isAcceptedTrack(p.pT(), eff))
      //   continue;

      fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
      pj.set_user_index(0); // 0 = background
      allParticles.push_back(pj);
    }

    // -----------------------------------------------------------------------
    // Step 3: Embed MC particles into the event
    // -----------------------------------------------------------------------
    for (const auto &p : mcDijet.jet1_particles) {
      // if (!isAcceptedTrack(p.pt, eff))
      // continue;
      fastjet::PseudoJet pj = p.toPseudoJet();
      pj.set_user_index(1); // 1 = MC jet 1
      allParticles.push_back(pj);
    }

    for (const auto &p : mcDijet.jet2_particles) {
      // if (!isAcceptedTrack(p.pt, eff))
      //   continue;
      fastjet::PseudoJet pj = p.toPseudoJet();
      pj.set_user_index(2); // 2 = MC jet 2
      allParticles.push_back(pj);
    }

    // -----------------------------------------------------------------------
    // Step 4: Reconstruct jets from combined system
    // -----------------------------------------------------------------------
    fastjet::ClusterSequence cs(allParticles, jetDef);
    auto all_jets = fastjet::sorted_by_pt(cs.inclusive_jets(cfg.jet_ptMin));

    std::vector<fastjet::PseudoJet> selectedJets;
    for (const auto &jet : all_jets) {
      if (std::abs(jet.eta()) < cfg.jet_etaMax) {
        selectedJets.push_back(jet);
      }
    }

    // -----------------------------------------------------------------------
    // Step 5: Match reconstructed jets to MC truth
    // -----------------------------------------------------------------------
    reco_eventID = iEvent;
    matched1 = false;
    matched2 = false;
    dR_mc1_reco1 = 999.0;
    dR_mc2_reco2 = 999.0;

    reco_jet1_pt = -999;
    reco_jet1_eta = -999;
    reco_jet1_phi = -999;
    reco_jet1_nPart = -999;

    reco_jet2_pt = -999;
    reco_jet2_eta = -999;
    reco_jet2_phi = -999;
    reco_jet2_nPart = -999;

    // Find best match for jet 1
    for (const auto &jet : selectedJets) {
      double dr = deltaR(mc_jet1_eta, mc_jet1_phi, jet.eta(), jet.phi());
      if (dr < dR_mc1_reco1 && dr < cfg.jet_R) {
        dR_mc1_reco1 = dr;
        matched1 = true;
        reco_jet1_pt = jet.pt();
        reco_jet1_eta = jet.eta();
        reco_jet1_phi = jet.phi();
        reco_jet1_nPart = jet.constituents().size();
      }
    }

    // Find best match for jet 2
    for (const auto &jet : selectedJets) {
      double dr = deltaR(mc_jet2_eta, mc_jet2_phi, jet.eta(), jet.phi());
      if (dr < dR_mc2_reco2 && dr < cfg.jet_R) {
        dR_mc2_reco2 = dr;
        matched2 = true;
        reco_jet2_pt = jet.pt();
        reco_jet2_eta = jet.eta();
        reco_jet2_phi = jet.phi();
        reco_jet2_nPart = jet.constituents().size();
      }
    }

    if (matched1 && matched2)
      nMatched++;

    // -----------------------------------------------------------------------
    // Step 6: Analyze background transverse system
    // -----------------------------------------------------------------------

    BackgroundResult bkg = analyzeBackground(selectedJets, allParticles, cfg);
    bkg_multA = bkg.multA;
    bkg_multB = bkg.multB;

    tReco->Fill();

    nProcessed++;
  }

  std::cout << "\n\n=== Summary ===\n";
  std::cout << "Total events processed: " << nProcessed << "\n";
  std::cout << "Both jets matched: " << nMatched << " ("
            << (100.0 * nMatched / nProcessed) << "%)\n";

  // =========================================================================
  // Save and cleanup
  // =========================================================================
  fout->Write();

  mcGen.writeQaHists(fout);

  fout->Close();
  delete fout;

  pythia.stat();

  std::cout << "\nOutput written to: " << outFile << "\n";
  return 0;
}