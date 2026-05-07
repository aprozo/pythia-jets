#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace Pythia8;

static double wrapPhi(double phi) {
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

// ---------------------------- config ----------------------------
struct Config {

  // acceptance (charged final only)
  double part_ptMin;
  double part_etaMax;

  // jet finder
  double jet_R;
  double jet_ptMin;
  double jet_etaMax;

  // generation safety cone (keep constituents well inside R to avoid splitting)
  double gen_Rmax; // <= ~R/2 recommended

  int nEvents;

  Config()
      : part_ptMin(0.15), part_etaMax(1.0), jet_R(0.4), jet_ptMin(8.0),
        jet_etaMax(1.0 - 0.4), gen_Rmax(0.20), nEvents(10000) {}
};

struct JetConstituents {
  float z[100];
  float jt[100];
  float px[100];
  float py[100];
  float pz[100];
  float e[100];
};

// ============================================================================
// MC Dijet Structure (store fastjet objects directly; no MCParticle needed)
// ============================================================================
struct McDijet {
  double jet1_pt, jet1_eta, jet1_phi;
  int jet1_nPart;

  double jet2_pt, jet2_eta, jet2_phi;
  int jet2_nPart;

  JetConstituents jet1_constituents;
  JetConstituents jet2_constituents;

  McDijet()
      : jet1_pt(0), jet1_eta(0), jet1_phi(0), jet1_nPart(0), jet2_pt(0),
        jet2_eta(0), jet2_phi(0), jet2_nPart(0) {}
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

// ---------------------------- main embedding loop
// ----------------------------
int main(int argc, char **argv) {
  Config cfg;

  std::string outFile = "embeddedSpecialDijets.root";

  if (argc > 1)
    cfg.nEvents = std::atoi(argv[1]);
  if (argc > 2)
    outFile = argv[2];

  cfg.jet_etaMax = cfg.part_etaMax - cfg.jet_R;

  // Tracking efficiency function
  TF1 eff("eff", "[0]*(1-exp(-pow(x/[1],[2])))", 0, 30);
  eff.SetParameters(0.88, 0.25, 1.2);

  // background Pythia
  Pythia pythiaBackground;
  pythiaBackground.readString("Beams:idA = 2212");
  pythiaBackground.readString("Beams:idB = 2212");
  pythiaBackground.readString("Beams:eCM = 200.");
  pythiaBackground.readString("SoftQCD:nonDiffractive = on");

  if (!pythiaBackground.init()) {
    std::cerr << "[ERROR] Pythia background init failed\n";
    return 1;
  }

  fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm, cfg.jet_R);

  // Input
  TFile *inputFile = TFile::Open("specialDijets.root", "READ");
  TTree *inputTree = (TTree *)inputFile->Get("tree");
  McDijet mcDijet;

  inputTree->SetBranchAddress("jet1_pt", &mcDijet.jet1_pt);
  inputTree->SetBranchAddress("jet1_eta", &mcDijet.jet1_eta);
  inputTree->SetBranchAddress("jet1_phi", &mcDijet.jet1_phi);

  inputTree->SetBranchAddress("jet2_pt", &mcDijet.jet2_pt);
  inputTree->SetBranchAddress("jet2_eta", &mcDijet.jet2_eta);
  inputTree->SetBranchAddress("jet2_phi", &mcDijet.jet2_phi);

  inputTree->SetBranchAddress("jet1_nNch", &mcDijet.jet1_nPart);
  inputTree->SetBranchAddress("jet2_nNch", &mcDijet.jet2_nPart);

  inputTree->SetBranchAddress("px1_constituents",
                              &mcDijet.jet1_constituents.px);
  inputTree->SetBranchAddress("py1_constituents",
                              &mcDijet.jet1_constituents.py);
  inputTree->SetBranchAddress("pz1_constituents",
                              &mcDijet.jet1_constituents.pz);
  inputTree->SetBranchAddress("e1_constituents", &mcDijet.jet1_constituents.e);
  inputTree->SetBranchAddress("jt1_constituents",
                              &mcDijet.jet1_constituents.jt);
  inputTree->SetBranchAddress("z1_constituents", &mcDijet.jet1_constituents.z);
  inputTree->SetBranchAddress("px2_constituents",
                              &mcDijet.jet2_constituents.px);
  inputTree->SetBranchAddress("py2_constituents",
                              &mcDijet.jet2_constituents.py);
  inputTree->SetBranchAddress("pz2_constituents",
                              &mcDijet.jet2_constituents.pz);
  inputTree->SetBranchAddress("e2_constituents", &mcDijet.jet2_constituents.e);
  inputTree->SetBranchAddress("jt2_constituents",
                              &mcDijet.jet2_constituents.jt);
  inputTree->SetBranchAddress("z2_constituents", &mcDijet.jet2_constituents.z);

  // Output
  TFile *fout = new TFile(outFile.c_str(), "RECREATE");

  TTree *outTree =
      new TTree("events", "Reco jets in embedded SoftQCD background");

  // Variables requested
  int eventID = 0;

  double mc_jet1_pt = 0, mc_jet1_eta = -99, mc_jet1_phi = -99;
  int mc_jet1_nPart = 0;

  double mc_jet2_pt = 0, mc_jet2_eta = -999, mc_jet2_phi = -99;
  int mc_jet2_nPart = 0;

  int matched1, matched2;
  double reco_jet1_pt, reco_jet1_eta, reco_jet1_phi;
  int reco_jet1_nPart;
  double reco_jet2_pt, reco_jet2_eta, reco_jet2_phi;
  int reco_jet2_nPart;
  double dR_mc1_reco1, dR_mc2_reco2;
  int bkg_multA, bkg_multB;

  outTree->Branch("eventID", &eventID);
  outTree->Branch("matched1", &matched1);
  outTree->Branch("matched2", &matched2);
  outTree->Branch("mc_jet1_pt", &mc_jet1_pt);
  outTree->Branch("mc_jet1_eta", &mc_jet1_eta);
  outTree->Branch("mc_jet1_phi", &mc_jet1_phi);
  outTree->Branch("mc_jet1_nPart", &mc_jet1_nPart);
  outTree->Branch("reco_jet1_pt", &reco_jet1_pt);
  outTree->Branch("reco_jet1_eta", &reco_jet1_eta);
  outTree->Branch("reco_jet1_phi", &reco_jet1_phi);
  outTree->Branch("reco_jet1_nPart", &reco_jet1_nPart);
  outTree->Branch("dR_mc1_reco1", &dR_mc1_reco1);

  outTree->Branch("mc_jet2_pt", &mc_jet2_pt);
  outTree->Branch("mc_jet2_eta", &mc_jet2_eta);
  outTree->Branch("mc_jet2_phi", &mc_jet2_phi);
  outTree->Branch("mc_jet2_nPart", &mc_jet2_nPart);
  outTree->Branch("reco_jet2_pt", &reco_jet2_pt);
  outTree->Branch("reco_jet2_eta", &reco_jet2_eta);
  outTree->Branch("reco_jet2_phi", &reco_jet2_phi);
  outTree->Branch("reco_jet2_nPart", &reco_jet2_nPart);
  outTree->Branch("dR_mc2_reco2", &dR_mc2_reco2);

  outTree->Branch("multA", &bkg_multA);
  outTree->Branch("multB", &bkg_multB);

  int maxEvents = std::min(cfg.nEvents, (int)inputTree->GetEntries());
  // Event loop
  for (int iev = 0; iev < maxEvents; ++iev) {
    eventID = iev;
    if (!pythiaBackground.next())
      continue;

    inputTree->GetEntry(iev);

    // Store truth info
    mc_jet1_pt = mcDijet.jet1_pt;
    mc_jet1_eta = mcDijet.jet1_eta;
    mc_jet1_phi = mcDijet.jet1_phi;
    mc_jet1_nPart = mcDijet.jet1_nPart;

    mc_jet2_pt = mcDijet.jet2_pt;
    mc_jet2_eta = mcDijet.jet2_eta;
    mc_jet2_phi = mcDijet.jet2_phi;
    mc_jet2_nPart = mcDijet.jet2_nPart;

    // ------------------ build reco particles by smearing truthParticles
    // ------------------

    // ------------------ background truth particles ------------------
    std::vector<fastjet::PseudoJet> truthParticles;
    truthParticles.reserve(4000);

    for (int i = 0; i < pythiaBackground.event.size(); ++i) {
      const Particle &p = pythiaBackground.event[i];
      if (!p.isFinal() || !p.isVisible() || !p.isCharged())
        continue;
      if (std::fabs(p.eta()) > cfg.part_etaMax)
        continue;
      if (p.pT() < cfg.part_ptMin)
        continue;

      fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.e());
      pj.set_user_index(0); // background
      truthParticles.push_back(pj);
    }

    // -----------------------------------------------------------------------
    // Step 3: Embed MC dijet particles into the event and smear whole event
    // -----------------------------------------------------------------------
    // add jet 1 constituents
    for (int i = 0; i < mcDijet.jet1_nPart; ++i) {
      fastjet::PseudoJet pj(
          mcDijet.jet1_constituents.px[i], mcDijet.jet1_constituents.py[i],
          mcDijet.jet1_constituents.pz[i], mcDijet.jet1_constituents.e[i]);
      pj.set_user_index(1); // signal
      truthParticles.push_back(pj);
    }
    // add jet 2 constituents
    for (int i = 0; i < mcDijet.jet2_nPart; ++i) {
      fastjet::PseudoJet pj(
          mcDijet.jet2_constituents.px[i], mcDijet.jet2_constituents.py[i],
          mcDijet.jet2_constituents.pz[i], mcDijet.jet2_constituents.e[i]);
      pj.set_user_index(1); // signal
      truthParticles.push_back(pj);
    }

    // ------------------ recoParticles (smearing hook) ------------------
    std::vector<fastjet::PseudoJet> recoParticles;
    recoParticles.reserve(truthParticles.size());
    for (size_t i = 0; i < truthParticles.size(); ++i) {
      fastjet::PseudoJet smeared = truthParticles[i];
      recoParticles.push_back(smeared);
    }
    // -----------------------------------------------------------------------
    // Step 4: Reconstruct jets from combined system
    // -----------------------------------------------------------------------
    fastjet::ClusterSequence cs(recoParticles, jetDefinition);
    std::vector<fastjet::PseudoJet> allRecoJets =
        fastjet::sorted_by_pt(cs.inclusive_jets(cfg.jet_ptMin));

    std::vector<fastjet::PseudoJet> selectedJets;
    for (const auto &jet : allRecoJets) {
      if (std::abs(jet.eta()) < cfg.jet_etaMax) {
        selectedJets.push_back(jet);
      }
    }

    // Reset reco outputs
    matched1 = matched2 = 0;
    reco_jet1_pt = reco_jet1_eta = reco_jet1_phi = -999;
    reco_jet2_pt = reco_jet2_eta = reco_jet2_phi = -999;
    reco_jet1_nPart = reco_jet2_nPart = -1;
    dR_mc1_reco1 = dR_mc2_reco2 = 1e3;
    bkg_multA = bkg_multB = 0;

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

    // ------------------ background multiplicities (transverse regions)
    // ------------------
    BackgroundResult bkg = analyzeBackground(selectedJets, recoParticles, cfg);
    bkg_multA = bkg.multA;
    bkg_multB = bkg.multB;

    outTree->Fill();
  }

  outTree->Write();
  fout->Close();
  delete fout;

  pythiaBackground.stat();
  return 0;
}
