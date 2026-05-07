// make_pythia_templates_1bin.cc
// One-bin (global) sampling templates from Pythia dijets:
//  - hLeadPt : leading jet pT
//  - hZ      : z = (p_i · n̂) / |p_jet|
//  - hJT     : jT = |p_i - (p_i·n̂)n̂|
//  - hNch    : charged constituent multiplicity per jet
//  - hXJ     : xJ = pT2/pT1
//  - hDphi   : |Δphi(j2,j1)| (away-side selected)
//  - hDeta   : Δη(j2,j1)
//  - nbd_mu, nbd_k saved as TParameter<double> estimated from hNch
//  (method-of-moments)
//
// Compile:
//   g++ -O2 -std=c++98 make_pythia_templates_1bin.cc -o
//   make_pythia_templates_1bin \
//     $(pythia8-config --cxxflags --libs) $(fastjet-config --cxxflags --libs) \
//     $(root-config --cflags --libs)
//
// Run:
//   ./make_pythia_templates_1bin 200000 templates.root 12 0.4
//
// Args: nEvents outFile pTHatMin R

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TParameter.h"
#include "TTree.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace Pythia8;

// ---------------------------- helpers ----------------------------
static double wrapPhi(double phi) {
  while (phi <= -M_PI)
    phi += 2.0 * M_PI;
  while (phi > M_PI)
    phi -= 2.0 * M_PI;
  return phi;
}

struct JetConstituents {
  float z[100];
  float jT[100];
  float px[100];
  float py[100];
  float pz[100];
  float e[100];
};

// z, jT, and Nch for one jet
JetConstituents fillFrag(const fastjet::PseudoJet &jet, TH1D *hZ, TH1D *hJT,
                         TH1D *hNch) {
  if (!hZ || !hJT || !hNch)
    return JetConstituents();
  JetConstituents jetConstituents;

  const double pxJ = jet.px();
  const double pyJ = jet.py();
  const double pzJ = jet.pz();
  const double pJ = std::sqrt(pxJ * pxJ + pyJ * pyJ + pzJ * pzJ);
  if (pJ <= 0.0)
    return JetConstituents();

  const double nx = pxJ / pJ;
  const double ny = pyJ / pJ;
  const double nz = pzJ / pJ;

  const std::vector<fastjet::PseudoJet> c = jet.constituents();
  hNch->Fill((double)c.size());

  for (size_t i = 0; i < c.size(); ++i) {
    const fastjet::PseudoJet &p = c[i];
    const double px = p.px();
    const double py = p.py();
    const double pz = p.pz();

    const double ppar = px * nx + py * ny + pz * nz; // longitudinal component
    if (ppar <= 0.0)
      continue;

    const double z = ppar / pJ;
    hZ->Fill(z);

    const double ptx = px - ppar * nx;
    const double pty = py - ppar * ny;
    const double ptz = pz - ppar * nz;
    const double jt = std::sqrt(ptx * ptx + pty * pty + ptz * ptz);
    hJT->Fill(jt);
    jetConstituents.z[i] = (float)z;
    jetConstituents.jT[i] = (float)jt;
    jetConstituents.px[i] = (float)px;
    jetConstituents.py[i] = (float)py;
    jetConstituents.pz[i] = (float)pz;
    jetConstituents.e[i] = (float)p.e();
  }
  return jetConstituents;
}

// NBD method-of-moments estimate from multiplicity histogram:
// Var = mu + mu^2/k  => k = mu^2 / (Var - mu)
static void estimateNBD(const TH1D *hNch, double &mu, double &k) {
  mu = 0.0;
  k = 1e9; // ~Poisson
  if (!hNch)
    return;

  const double m = hNch->GetMean();
  const double v = hNch->GetRMS();
  const double var = v * v;

  mu = m;
  if (var > mu + 1e-9)
    k = mu * mu / (var - mu);
}

// ---------------------------- config ----------------------------
struct Config {
  int nEvents;
  std::string outFile;

  double eCM;
  double pTHatMin;
  double pTHatMax;

  double part_ptMin;
  double part_etaMax;

  double R;
  double jet_ptMin;
  double jet_etaMax;

  double awayDphiMin;

  Config()
      : nEvents(2000000), outFile("specialDijets.root"), eCM(200.0),
        pTHatMin(30.0), pTHatMax(-1.0), part_ptMin(0.15), part_etaMax(1.0),
        R(0.4), jet_ptMin(3.0), jet_etaMax(1.0 - 0.4),
        awayDphiMin(3 * TMath::Pi() / 4) {}
};

static std::string toStr(double x) {
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

// ---------------------------- main ----------------------------
int main(int argc, char **argv) {
  Config cfg;

  if (argc > 1)
    cfg.nEvents = std::atoi(argv[1]);
  if (argc > 2)
    cfg.outFile = argv[2];
  if (argc > 3)
    cfg.pTHatMin = std::atof(argv[3]);
  if (argc > 4)
    cfg.R = std::atof(argv[4]);

  cfg.jet_etaMax = cfg.part_etaMax - cfg.R;

  // --- histograms (one bin set) ---
  TH1D *hLeadPt = new TH1D("hLeadPt", "Leading jet pT; p_{T,1} [GeV]; counts",
                           1200, 0, 120);
  TH1D *hZ = new TH1D("hZ", "z; z; counts", 1200, 0, 1.2);
  TH1D *hJT = new TH1D("hJT", "j_{T} [GeV]; j_{T}; counts", 1600, 0, 8);
  TH1D *hNch = new TH1D("hNch", "N_{ch} in jet; N_{ch}; counts", 120, 0, 120);
  TH1D *hXJ = new TH1D("hXJ", "x_{J}=pT2/pT1; x_{J}; counts", 1200, 0, 1.2);
  TH1D *hDphi =
      new TH1D("hDphi", "|#Delta#phi|; |#Delta#phi|; counts", 1200, 2.0, 3.25);
  TH1D *hDeta =
      new TH1D("hDeta", "#Delta#eta; #Delta#eta; counts", 1200, -2.0, 2.0);

  TTree *tree = new TTree("tree", "Reco jets in embedded SoftQCD background");
  // Variables requested
  int eventID = 0;

  double jet1_pt = 0, jet1_eta = -99, jet1_phi = -99;
  double jet2_pt = 0, jet2_eta = -999, jet2_phi = -99;

  int jet1_nNch = 0;
  JetConstituents jet1Constituents;

  int jet2_nNch = 0;
  JetConstituents jet2Constituents;

  tree->Branch("eventID", &eventID);

  tree->Branch("jet1_pt", &jet1_pt);
  tree->Branch("jet1_eta", &jet1_eta);
  tree->Branch("jet1_phi", &jet1_phi);

  tree->Branch("jet2_pt", &jet2_pt);
  tree->Branch("jet2_eta", &jet2_eta);
  tree->Branch("jet2_phi", &jet2_phi);

  tree->Branch("jet1_nNch", &jet1_nNch);
  tree->Branch("jet2_nNch", &jet2_nNch);

  tree->Branch("px1_constituents", &jet1Constituents.px,
               "px1_constituents[jet1_nNch]/F");
  tree->Branch("py1_constituents", &jet1Constituents.py,
               "py1_constituents[jet1_nNch]/F");
  tree->Branch("pz1_constituents", &jet1Constituents.pz,
               "pz1_constituents[jet1_nNch]/F");
  tree->Branch("e1_constituents", &jet1Constituents.e,
               "e1_constituents[jet1_nNch]/F");

  tree->Branch("jt1_constituents", &jet1Constituents.jT,
               "jt1_constituents[jet1_nNch]/F");
  tree->Branch("z1_constituents", &jet1Constituents.z,
               "z1_constituents[jet1_nNch]/F");

  tree->Branch("px2_constituents", &jet2Constituents.px,
               "px2_constituents[jet2_nNch]/F");
  tree->Branch("py2_constituents", &jet2Constituents.py,
               "py2_constituents[jet2_nNch]/F");
  tree->Branch("pz2_constituents", &jet2Constituents.pz,
               "pz2_constituents[jet2_nNch]/F");
  tree->Branch("e2_constituents", &jet2Constituents.e,
               "e2_constituents[jet2_nNch]/F");
  tree->Branch("jt2_constituents", &jet2Constituents.jT,
               "jt2_constituents[jet2_nNch]/F");
  tree->Branch("z2_constituents", &jet2Constituents.z,
               "z2_constituents[jet2_nNch]/F");
  // keep in memory until write
  hLeadPt->SetDirectory(0);
  hZ->SetDirectory(0);
  hJT->SetDirectory(0);
  hNch->SetDirectory(0);
  hXJ->SetDirectory(0);
  hDphi->SetDirectory(0);
  hDeta->SetDirectory(0);

  // --- Pythia hard QCD for dijets ---
  Pythia p;
  p.readString("Beams:idA = 2212");
  p.readString("Beams:idB = 2212");
  p.readString("Beams:eCM = " + toStr(cfg.eCM));
  p.readString("HardQCD:all = on");
  p.readString("PhaseSpace:pTHatMin = " + toStr(cfg.pTHatMin));
  if (cfg.pTHatMax > 0.0)
    p.readString("PhaseSpace:pTHatMax = " + toStr(cfg.pTHatMax));

  if (!p.init()) {
    std::cerr << "[ERROR] Pythia init failed\n";
    return 1;
  }

  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfg.R);

  int nDijetsKept = 0;

  for (int iev = 0; iev < cfg.nEvents; ++iev) {
    if (!p.next())
      continue;

    // collect charged final visible particles within acceptance
    std::vector<fastjet::PseudoJet> parts;
    parts.reserve(2000);

    for (int i = 0; i < p.event.size(); ++i) {
      const Particle &pp = p.event[i];
      if (!pp.isFinal() || !pp.isVisible() || !pp.isCharged())
        continue;
      if (std::fabs(pp.eta()) > cfg.part_etaMax)
        continue;
      if (pp.pT() < cfg.part_ptMin)
        continue;
      parts.push_back(fastjet::PseudoJet(pp.px(), pp.py(), pp.pz(), pp.e()));
    }
    if (parts.size() < 2)
      continue;

    // jets
    fastjet::ClusterSequence cs(parts, jetDef);
    std::vector<fastjet::PseudoJet> jets =
        fastjet::sorted_by_pt(cs.inclusive_jets(cfg.jet_ptMin));
    if (jets.size() < 2)
      continue;

    // leading jet in acceptance
    const fastjet::PseudoJet &j1 = jets[0];
    if (std::fabs(j1.eta()) > cfg.jet_etaMax)
      continue;

    // away-side partner: highest pT satisfying Δphi + acceptance
    int idx2 = -1;
    double bestPt = -1.0;

    for (size_t j = 1; j < jets.size(); ++j) {
      if (std::fabs(jets[j].eta()) > cfg.jet_etaMax)
        continue;
      const double dphi = std::fabs(wrapPhi(jets[j].phi() - j1.phi()));
      if (dphi < cfg.awayDphiMin)
        continue;
      if (jets[j].pt() > bestPt) {
        bestPt = jets[j].pt();
        idx2 = (int)j;
      }
    }
    if (idx2 < 0)
      continue;

    const fastjet::PseudoJet &j2 = jets[idx2];

    // compare the number of constituents
    if (j1.constituents().size() != j2.constituents().size())
      continue;

    // store tree variables
    eventID = iev;
    jet1_pt = j1.pt();
    jet1_eta = j1.eta();
    jet1_phi = j1.phi();
    jet1_nNch = j1.constituents().size();

    jet2_pt = j2.pt();
    jet2_eta = j2.eta();
    jet2_phi = j2.phi();
    jet2_nNch = j2.constituents().size();

    // fill dijet-level templates
    hLeadPt->Fill(j1.pt());
    hXJ->Fill(j2.pt() / j1.pt());
    hDphi->Fill(std::fabs(wrapPhi(j2.phi() - j1.phi())));
    hDeta->Fill(j2.eta() - j1.eta());

    // fill fragmentation from both jets (one-bin, global)
    JetConstituents jet1Constituents = fillFrag(j1, hZ, hJT, hNch);
    JetConstituents jet2Constituents = fillFrag(j2, hZ, hJT, hNch);

    ++nDijetsKept;
    tree->Fill();
  }

  // NBD parameters from hNch
  double nbd_mu = 0.0, nbd_k = 1e9;
  estimateNBD(hNch, nbd_mu, nbd_k);

  // write output
  TFile *fout = new TFile(cfg.outFile.c_str(), "RECREATE");

  hLeadPt->Write();
  hZ->Write();
  hJT->Write();
  hNch->Write();
  hXJ->Write();
  hDphi->Write();
  hDeta->Write();

  tree->Write();
  fout->Close();
  delete fout;

  std::cerr << "[templates] dijets kept: " << nDijetsKept << "\n";
  p.stat();
  return 0;
}
