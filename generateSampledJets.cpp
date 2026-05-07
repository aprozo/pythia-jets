// embed_controlled_dijets.cc
// Embed controlled back-to-back dijets (sampled from your templates.root)
// into a soft Pythia background event, while enforcing:
//   sum(p_i_vec) == target jet momentum vector  (=> jet pT/eta/phi reproduced
//   if clustered)
// using constrained sampling of z and jT.

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

// ---------------------------- templates ----------------------------
struct JetDistributionSamples {
  TH1D *hLeadPt;
  TH1D *hZ;
  TH1D *hJT;
  TH1D *hXJ;
  TH1D *hDphi;
  TH1D *hDeta;
  TH1D *hNch;

  JetDistributionSamples()
      : hLeadPt(0), hZ(0), hJT(0), hXJ(0), hDphi(0), hDeta(0), hNch(0) {}
};

// ---------------------------- templates (FIXED: clone histograms)
// ----------------------------
static TH1D *getH1Clone(TFile *file, const char *name) {
  TH1D *h = dynamic_cast<TH1D *>(file->Get(name));
  if (!h)
    return 0;
  TH1D *c = dynamic_cast<TH1D *>(h->Clone());
  if (!c)
    return 0;
  c->SetDirectory(0);
  return c;
}

static bool loadTemplates(const std::string &templateFile,
                          JetDistributionSamples &t) {
  TFile *f = TFile::Open(templateFile.c_str(), "READ");
  if (!f || f->IsZombie())
    return false;

  t.hLeadPt = getH1Clone(f, "hLeadPt");
  t.hZ = getH1Clone(f, "hZ");
  t.hJT = getH1Clone(f, "hJT");
  t.hXJ = getH1Clone(f, "hXJ");
  t.hDphi = getH1Clone(f, "hDphi");
  t.hDeta = getH1Clone(f, "hDeta");
  t.hNch = getH1Clone(f, "hNch");

  f->Close();
  delete f;

  return (t.hLeadPt && t.hZ && t.hJT && t.hXJ && t.hDphi && t.hDeta && t.hNch);
}

// ============================================================================
// MC Dijet Structure (store fastjet objects directly; no MCParticle needed)
// ============================================================================
struct McDijet {
  double jet1_pt, jet1_eta, jet1_phi;
  int jet1_nPart;

  double jet2_pt, jet2_eta, jet2_phi;
  int jet2_nPart;

  std::vector<fastjet::PseudoJet> jet1_particles; // user_index = 1
  std::vector<fastjet::PseudoJet> jet2_particles; // user_index = 2

  McDijet()
      : jet1_pt(0), jet1_eta(0), jet1_phi(0), jet1_nPart(0), jet2_pt(0),
        jet2_eta(0), jet2_phi(0), jet2_nPart(0) {}
};

// Helper: sample histogram value with hard bounds (uses global gRandom set in
// main)
static double sampleHistogramWithBounds(TH1D *histogram, double minValue,
                                        double maxValue, int maxTries = 5000) {
  if (!histogram)
    return 0.0;

  double value = 0.0;
  for (int it = 0; it < maxTries; ++it) {
    value = histogram->GetRandom();
    if (value >= minValue && value <= maxValue)
      return value;
  }
  // last resort: clamp
  if (value < minValue)
    value = minValue;
  if (value > maxValue)
    value = maxValue;
  return value;
}

// Helper: build two unit vectors transverse to the jet axis (right-handed)
static void buildTransverseBasis(const TVector3 &jetAxisUnitVector,
                                 TVector3 &transverseUnitVector1,
                                 TVector3 &transverseUnitVector2) {
  TVector3 referenceVector(0.0, 0.0, 1.0);
  if (std::fabs(jetAxisUnitVector.Dot(referenceVector)) > 0.9)
    referenceVector.SetXYZ(0.0, 1.0, 0.0);

  transverseUnitVector1 = jetAxisUnitVector.Cross(referenceVector);
  if (transverseUnitVector1.Mag() <= 0.0) {
    referenceVector.SetXYZ(1.0, 0.0, 0.0);
    transverseUnitVector1 = jetAxisUnitVector.Cross(referenceVector);
  }
  transverseUnitVector1 = transverseUnitVector1.Unit();
  transverseUnitVector2 = jetAxisUnitVector.Cross(transverseUnitVector1).Unit();
}

// ============================================================================
// Tracking efficiency
// ============================================================================
bool isAcceptedTrack(double pt, TF1 &eff) {
  // if (pt > 30.0)
  //   return false;
  return gRandom->Rndm() < eff.Eval(pt);
}

// ============================================================================
// McDijetGenerator
// ============================================================================
struct McDijetGenerator {
  const Config &cfg;
  const JetDistributionSamples &templ;
  TRandom3 &rng; // single RNG shared with gRandom (TH1::GetRandom uses gRandom)

  // optional QA histograms (write them from main if you want)
  TH1D *hQaJetPt;
  TH1D *hQaNch;
  TH1D *hQaXJ;
  TH1D *hQaDphi;
  TH1D *hQaDeta;
  TH1D *hQaJt;
  TH1D *hQaZ;

  TH1D *hQaParticlePt;
  TH1D *hQaParticlePhi;
  TH1D *hQaParticleEta;
  TH2D *hQaJetPt1VsJetPt2;
  TH2D *hQaJetPhi1VsJetPhi2;
  TH2D *hQaJetEta1VsJetEta2;
  TH2D *hQaJetN1VsJetN2;

  McDijetGenerator(const Config &config,
                   const JetDistributionSamples &templates, TRandom3 &random)
      : cfg(config), templ(templates), rng(random) {
    hQaJetPt =
        new TH1D("hQaJetPt", "jet pT; p_{T,1} [GeV]; counts", 240, 0, 120);
    hQaNch = new TH1D("hQaNch", "N_{ch}; N_{ch}; counts", 120, 0, 120);
    hQaXJ = new TH1D("hQaXJ", "x_{J}; x_{J}; counts", 240, 0, 1.2);
    hQaDphi = new TH1D("hQaDphi", "|#Delta#phi|; |#Delta#phi|; counts", 200,
                       2.0, 3.25);
    hQaDeta =
        new TH1D("hQaDeta", "#Delta#eta; #Delta#eta; counts", 200, -2.0, 2.0);
    hQaJt = new TH1D("hQaJt", "j_{T}; j_{T} [GeV]; counts", 1600, 0, 8);
    hQaZ = new TH1D("hQaZ", "z; z; counts", 1200, 0, 1.2);

    hQaParticlePt = new TH1D("hQaParticlePt",
                             "Particle pT; p_{T} (GeV/c); Counts", 500, 0, 50);
    hQaParticlePhi = new TH1D("hQaParticlePhi", "Particle Phi; #phi;Counts",
                              100, -M_PI, M_PI);
    hQaParticleEta = new TH1D("hQaParticleEta", "Particle Eta; #eta;Counts",
                              100, -config.part_etaMax, config.part_etaMax);
    hQaJetPt1VsJetPt2 =
        new TH2D("hQaJetPt1VsJetPt2",
                 "Jet1 pT vs Jet2 pT; p_{T,1} (GeV/c); p_{T,2} (GeV/c)", 100, 0,
                 100, 100, 0, 100);
    hQaJetPhi1VsJetPhi2 = new TH2D("hQaJetPhi1VsJetPhi2",
                                   "Jet1 Phi vs Jet2 Phi; #phi_{1}; #phi_{2}",
                                   100, -M_PI, M_PI, 100, -M_PI, M_PI);
    hQaJetEta1VsJetEta2 = new TH2D(
        "hQaJetEta1VsJetEta2", "Jet1 Eta vs Jet2 Eta; #eta_{1}; #eta_{2}", 100,
        -config.part_etaMax, config.part_etaMax, 100, -config.part_etaMax,
        config.part_etaMax);
    hQaJetN1VsJetN2 =
        new TH2D("hQaJetN1VsJetN2", "Jet1 N vs Jet2 N; N_{1}; N_{2}", 30, 0, 30,
                 30, 0, 30);

    hQaJetPt->SetDirectory(0);
    hQaNch->SetDirectory(0);
    hQaXJ->SetDirectory(0);
    hQaDphi->SetDirectory(0);
    hQaDeta->SetDirectory(0);
    hQaJt->SetDirectory(0);
    hQaZ->SetDirectory(0);
    hQaParticlePt->SetDirectory(0);
    hQaParticlePhi->SetDirectory(0);
    hQaParticleEta->SetDirectory(0);
    hQaJetPt1VsJetPt2->SetDirectory(0);
    hQaJetPhi1VsJetPhi2->SetDirectory(0);
    hQaJetEta1VsJetEta2->SetDirectory(0);
    hQaJetN1VsJetN2->SetDirectory(0);
  }

  // Build massless constituents so that sum(p_i_vec) == jetMomentumVector
  std::vector<fastjet::PseudoJet>
  generateJetConstituents(TRandom3 &rng, const JetDistributionSamples &templ,
                          const Config &cfg, double targetJetPt,
                          double targetJetEta, double targetJetPhi,
                          int nParticles, int jetUserIndex) {
    std::vector<fastjet::PseudoJet> outParticles;

    if (nParticles <= 0) {
      cout << "Warning: nParticles <= 0 in generateJetConstituents" << endl;
      return outParticles;
    }
    // ------------------ target jet 3-momentum ------------------
    const double jetPx = targetJetPt * std::cos(targetJetPhi);
    const double jetPy = targetJetPt * std::sin(targetJetPhi);
    const double jetPz = targetJetPt * std::sinh(targetJetEta);

    const TVector3 jetMomentumVector(jetPx, jetPy, jetPz);
    const double jetMomentumMagnitude = jetMomentumVector.Mag();
    if (jetMomentumMagnitude <= 0.0)
      return outParticles;

    const TVector3 jetAxisUnitVector = jetMomentumVector.Unit();

    // Special case: N = 1 (single particle = the jet)
    if (nParticles == 1) {
      hQaZ->Fill(1.0);
      const TVector3 particleMomentumVector = jetMomentumVector;
      const double particleEnergy = particleMomentumVector.Mag();

      fastjet::PseudoJet particlePseudoJet(
          particleMomentumVector.X(), particleMomentumVector.Y(),
          particleMomentumVector.Z(), particleEnergy);
      particlePseudoJet.set_user_index(jetUserIndex);
      outParticles.push_back(particlePseudoJet);
      return outParticles;
    }

    // ------------------ transverse basis around jet axis ------------------
    TVector3 transverseUnitVector1, transverseUnitVector2;
    buildTransverseBasis(jetAxisUnitVector, transverseUnitVector1,
                         transverseUnitVector2);

    // Keep constituents inside a narrow cone so they cluster into one jet
    const double maxTransverseOverLongitudinal =
        std::tan(cfg.gen_Rmax); // jT / p_parallel < tan(Rmax)

    // z lower bound: keep >0, and try to avoid ultra-soft fragments
    const double zMin = 1e-3;

    // ------------------ constrained sampling (retry loop) ------------------
    const int maxAttempts = 300;
    for (int attempt = 0; attempt < maxAttempts; ++attempt) {

      // ---------- 1) sample z_i with exact sum(z_i) = 1 ----------
      std::vector<double> zFraction(nParticles, 0.0);

      double remainingZ = 1.0;
      bool zSamplingOk = true;

      for (int i = 0; i < nParticles - 1; ++i) {
        const int remainingParticlesAfterThis = (nParticles - 1 - i);
        const double minReservedForRest = zMin * remainingParticlesAfterThis;
        const double maxThisZ = remainingZ - minReservedForRest;

        if (maxThisZ < zMin) {
          zSamplingOk = false;
          break;
        }

        const double thisZ =
            sampleHistogramWithBounds(templ.hZ, zMin, maxThisZ, 10000);
        zFraction[i] = thisZ;
        hQaZ->Fill(thisZ);
        remainingZ -= thisZ;
      }

      if (!zSamplingOk) {
        cout << "Warning: z sampling failed in generateJetConstituents" << endl;
        continue;
      }

      zFraction[nParticles - 1] = remainingZ;
      hQaZ->Fill(remainingZ);

      // Shuffle so the "last particle" is not special
      for (int i = nParticles - 1; i > 0; --i) {
        const int j = rng.Integer(i + 1);
        std::swap(zFraction[i], zFraction[j]);
      }

      // ---------- 2) sample transverse kicks with exact sum(kT_vec) = 0
      // ----------
      std::vector<TVector3> transverseKickVector(nParticles, TVector3(0, 0, 0));

      TVector3 sumTransverseKick(0.0, 0.0, 0.0);
      bool kSamplingOk = true;

      for (int i = 0; i < nParticles - 1; ++i) {
        const double particleLongitudinalMomentum =
            zFraction[i] * jetMomentumMagnitude;
        if (particleLongitudinalMomentum <= 0.0) {
          kSamplingOk = false;
          break;
        }

        const double maxAllowedTransverseKick =
            particleLongitudinalMomentum * maxTransverseOverLongitudinal;

        const double sampledTransverseKickMagnitude =
            sampleHistogramWithBounds(templ.hJT, 0.0, maxAllowedTransverseKick);

        const double sampledAngle = rng.Uniform(0.0, 2.0 * M_PI);

        const double kickComponentAlongE1 =
            sampledTransverseKickMagnitude * std::cos(sampledAngle);
        const double kickComponentAlongE2 =
            sampledTransverseKickMagnitude * std::sin(sampledAngle);

        const TVector3 kickVector =
            transverseUnitVector1 * kickComponentAlongE1 +
            transverseUnitVector2 * kickComponentAlongE2;

        transverseKickVector[i] = kickVector;
        sumTransverseKick += kickVector;
        hQaJt->Fill(sampledTransverseKickMagnitude);
      }

      if (!kSamplingOk) {
        cout << "Warning: kT sampling failed in generateJetConstituents"
             << endl;
        continue;
      }

      // Last transverse kick closes the vector sum to zero
      hQaJt->Fill(sumTransverseKick.Mag());
      transverseKickVector[nParticles - 1] = -sumTransverseKick;

      // Check last kick also respects the cone constraint
      {
        const double particleLongitudinalMomentum =
            zFraction[nParticles - 1] * jetMomentumMagnitude;
        if (particleLongitudinalMomentum <= 0.0)
          continue;

        const double maxAllowedTransverseKick =
            particleLongitudinalMomentum * maxTransverseOverLongitudinal;

        const double lastKickMagnitude =
            transverseKickVector[nParticles - 1].Mag();
        if (lastKickMagnitude > maxAllowedTransverseKick)
          continue;
      }

      // ---------- 3) build particle momenta and validate acceptance ----------
      outParticles.clear();
      outParticles.reserve(nParticles);

      TVector3 sumParticleMomentum(0.0, 0.0, 0.0);
      bool particleAcceptanceOk = true;

      for (int i = 0; i < nParticles; ++i) {
        const double particleLongitudinalMomentum =
            zFraction[i] * jetMomentumMagnitude;

        const TVector3 particleMomentumVector =
            jetAxisUnitVector * particleLongitudinalMomentum +
            transverseKickVector[i];

        const double particleEnergy = particleMomentumVector.Mag(); // massless

        const double particlePt = particleMomentumVector.Pt();
        if (particlePt < cfg.part_ptMin) {
          particleAcceptanceOk = false;
          break;
        }
        hQaParticlePt->Fill(particlePt);

        const double particleEta = particleMomentumVector.Eta();
        if (std::fabs(particleEta) > cfg.part_etaMax) {
          particleAcceptanceOk = false;
          break;
        }
        hQaParticleEta->Fill(particleEta);

        const double particlePhi = particleMomentumVector.Phi();
        hQaParticlePhi->Fill(particlePhi);

        fastjet::PseudoJet particlePseudoJet(
            particleMomentumVector.X(), particleMomentumVector.Y(),
            particleMomentumVector.Z(), particleEnergy);
        particlePseudoJet.set_user_index(jetUserIndex);
        outParticles.push_back(particlePseudoJet);

        sumParticleMomentum += particleMomentumVector;
      }

      if (!particleAcceptanceOk)
        continue;

      // ---------- 4) final closure check: sum(p_i_vec) == jetMomentumVector
      // ----------
      const TVector3 momentumDifference =
          sumParticleMomentum - jetMomentumVector;
      const double relativeMismatch =
          momentumDifference.Mag() / std::max(1e-12, jetMomentumMagnitude);

      if (relativeMismatch > 1e-4)
        continue;
    }

    return outParticles;
  }

  bool generate(McDijet &dijet) {
    dijet = McDijet();

    // --- sample leading jet kinematics ---
    dijet.jet1_pt = templ.hLeadPt->GetRandom();
    dijet.jet1_eta = rng.Uniform(-cfg.jet_etaMax, cfg.jet_etaMax);
    dijet.jet1_phi = rng.Uniform(-M_PI, M_PI);

    // --- sample dijet correlations ---
    const double xJ = templ.hXJ->GetRandom();
    const double dphi = templ.hDphi->GetRandom(); // |Δphi|
    const double deta = templ.hDeta->GetRandom();

    dijet.jet2_pt = xJ * dijet.jet1_pt;
    dijet.jet2_eta = rng.Uniform(-cfg.jet_etaMax, cfg.jet_etaMax);

    // int nEtaMaxTries = 1000;
    // for (int i = 0; i < nEtaMaxTries; ++i) {
    //   dijet.jet2_eta = dijet.jet1_eta + deta;
    //   if (dijet.jet2_eta > cfg.jet_etaMax || dijet.jet2_eta <
    //   -cfg.jet_etaMax)
    //     continue;
    // }

    const double sign = (rng.Uniform() < 0.5) ? -1.0 : +1.0;
    dijet.jet2_phi = wrapPhi(dijet.jet1_phi + sign * dphi);

    // --- multiplicities (from template; replace here if you want NBD or forced
    // N) ---
    dijet.jet1_nPart = std::max(1, (int)templ.hNch->GetRandom());
    dijet.jet2_nPart = dijet.jet1_nPart;

    hQaJetPt1VsJetPt2->Fill(dijet.jet1_pt, dijet.jet2_pt);
    hQaJetPhi1VsJetPhi2->Fill(dijet.jet1_phi, dijet.jet2_phi);
    hQaJetEta1VsJetEta2->Fill(dijet.jet1_eta, dijet.jet2_eta);
    hQaJetN1VsJetN2->Fill(dijet.jet1_nPart, dijet.jet2_nPart);

    // --- generate constituents (enforces sum(p_i_vec) == target jet momentum)
    // ---
    dijet.jet1_particles =
        generateJetConstituents(rng, templ, cfg, dijet.jet1_pt, dijet.jet1_eta,
                                dijet.jet1_phi, dijet.jet1_nPart, 1);
    if (dijet.jet1_particles.empty())
      return false;

    dijet.jet2_particles =
        generateJetConstituents(rng, templ, cfg, dijet.jet2_pt, dijet.jet2_eta,
                                dijet.jet2_phi, dijet.jet2_nPart, 2);
    if (dijet.jet2_particles.empty())
      return false;

    // --- QA ---
    hQaJetPt->Fill(dijet.jet1_pt);
    hQaNch->Fill(dijet.jet1_nPart);

    double qaDphi = std::fabs(wrapPhi(dijet.jet1_phi - dijet.jet2_phi));
    double qaXJ = dijet.jet2_pt / dijet.jet1_pt;
    double qaDeta = dijet.jet1_eta - dijet.jet2_eta;
    hQaXJ->Fill(qaXJ);
    hQaDphi->Fill(qaDphi);
    hQaDeta->Fill(qaDeta);

    return true;
  }

  void writeQA(TFile *outFile) {
    outFile->cd();
    hQaJetPt->Write();
    hQaNch->Write();
    hQaXJ->Write();
    hQaDphi->Write();
    hQaDeta->Write();
    hQaJt->Write();
    hQaZ->Write();
    hQaParticlePt->Write();
    hQaParticlePhi->Write();
    hQaParticleEta->Write();
    hQaJetPt1VsJetPt2->Write();
    hQaJetPhi1VsJetPhi2->Write();
    hQaJetEta1VsJetEta2->Write();
    hQaJetN1VsJetN2->Write();
  }
};

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
  std::string templateFile = "templates.root";
  std::string outFile = "out.root";

  if (argc > 1)
    cfg.nEvents = std::atoi(argv[1]);
  if (argc > 2)
    templateFile = argv[2];
  if (argc > 3)
    outFile = argv[3];

  cfg.jet_etaMax = cfg.part_etaMax - cfg.jet_R;

  JetDistributionSamples templ;
  if (!loadTemplates(templateFile, templ)) {
    std::cerr << "[ERROR] failed to load templates from: " << templateFile
              << "\n";
    return 1;
  }

  // ROOT histogram sampling uses global gRandom
  TRandom3 rng(0);
  gRandom = &rng;

  McDijetGenerator mcGen(cfg, templ, rng);

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

  // Output
  TFile *fout = new TFile(outFile.c_str(), "RECREATE");

  TTree *tree = new TTree("tree", "Reco jets in embedded SoftQCD background");

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

  tree->Branch("eventID", &eventID);
  tree->Branch("matched1", &matched1);
  tree->Branch("matched2", &matched2);
  tree->Branch("mc_jet1_pt", &mc_jet1_pt);
  tree->Branch("mc_jet1_eta", &mc_jet1_eta);
  tree->Branch("mc_jet1_phi", &mc_jet1_phi);
  tree->Branch("mc_jet1_nPart", &mc_jet1_nPart);
  tree->Branch("reco_jet1_pt", &reco_jet1_pt);
  tree->Branch("reco_jet1_eta", &reco_jet1_eta);
  tree->Branch("reco_jet1_phi", &reco_jet1_phi);
  tree->Branch("reco_jet1_nPart", &reco_jet1_nPart);
  tree->Branch("dR_mc1_reco1", &dR_mc1_reco1);

  tree->Branch("mc_jet2_pt", &mc_jet2_pt);
  tree->Branch("mc_jet2_eta", &mc_jet2_eta);
  tree->Branch("mc_jet2_phi", &mc_jet2_phi);
  tree->Branch("mc_jet2_nPart", &mc_jet2_nPart);
  tree->Branch("reco_jet2_pt", &reco_jet2_pt);
  tree->Branch("reco_jet2_eta", &reco_jet2_eta);
  tree->Branch("reco_jet2_phi", &reco_jet2_phi);
  tree->Branch("reco_jet2_nPart", &reco_jet2_nPart);
  tree->Branch("dR_mc2_reco2", &dR_mc2_reco2);

  tree->Branch("multA", &bkg_multA);
  tree->Branch("multB", &bkg_multB);

  // Event loop
  for (int iev = 0; iev < cfg.nEvents; ++iev) {
    eventID = iev;
    if (!pythiaBackground.next())
      continue;

    McDijet mcDijet;
    if (!mcGen.generate(mcDijet))
      continue;

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
    // Step 3: Embed MC particles into the event and smear whole event
    // -----------------------------------------------------------------------
    truthParticles.insert(truthParticles.end(), mcDijet.jet1_particles.begin(),
                          mcDijet.jet1_particles.end());
    truthParticles.insert(truthParticles.end(), mcDijet.jet2_particles.begin(),
                          mcDijet.jet2_particles.end());

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

    tree->Fill();
  }

  tree->Write();
  mcGen.writeQA(fout);
  fout->Close();
  delete fout;

  pythiaBackground.stat();
  return 0;
}
