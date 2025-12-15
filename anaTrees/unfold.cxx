// analyzeEmbedding_unfoldN_vsPt.cxx
// Embedded dijet analysis + joint (N1,N2) unfolding vs pT
// Requires RooUnfold (RooUnfoldBayes, RooUnfoldResponse)

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include <cmath>
#include <iostream>
#include <vector>

void drawLabel(TPad *pad, float x = 0.57, TString extra = "") {
  pad->cd();
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->DrawLatex(x, 0.85, "PYTHIA 8.315");
  latex->DrawLatex(x, 0.80, "#it{pp} at #sqrt{#it{s}} = 200 GeV");
  latex->DrawLatex(x, 0.75, "|#eta|<0.6, anti-k_{T}, #it{R}=0.4");
  latex->DrawLatex(x, 0.70, "jet p_{t} > 3 GeV/c");
  latex->DrawLatex(x, 0.65, "p_{t,track}>0.15 GeV/c, |#eta_{track}|<1");
  latex->DrawLatex(x, 0.60, "|#phi_{1} - #phi_{2}| > 3#pi/4");
  latex->DrawLatex(x, 0.55, extra);
}

// ============================================================================
// Entropy and "covariance" (mutual information) helpers
// ============================================================================

double getEntropy(TH1D *h) {
  double entropy = 0.0;
  int nBins = h->GetNbinsX();
  double total = h->Integral();
  if (total <= 0.0)
    return 0.0;

  for (int i = 1; i <= nBins; ++i) {
    double p = h->GetBinContent(i) / total;
    if (p > 0.0)
      entropy -= p * std::log(p);
  }
  return entropy;
}

double getEntropy(TH2D *h) {
  double entropy = 0.0;
  int nBinsX = h->GetNbinsX();
  int nBinsY = h->GetNbinsY();
  double total = h->Integral();
  if (total <= 0.0)
    return 0.0;

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      double p = h->GetBinContent(i, j) / total;
      if (p > 0.0)
        entropy -= p * std::log(p);
    }
  }
  return entropy;
}

// "Covariance" = S1 + S2 - S12 (actually mutual information)
double getCovariance(TH2D *h) {
  TH1D *px = (TH1D *)h->ProjectionX("_px", 1, -1, "e");
  px->SetDirectory(0);
  double S1 = getEntropy(px);
  TH1D *py = (TH1D *)h->ProjectionY("_py", 1, -1, "e");
  py->SetDirectory(0);
  double S2 = getEntropy(py);
  double S12 = getEntropy(h);
  delete px;
  delete py;
  return S1 + S2 - S12;
}

// For TH3D: loop over pT bins and compute covariance slice-by-slice
TH1D *getCovariance(TH3D *h, TString title = "", TString name = "") {
  if (name == "")
    name = TString(h->GetName()) + "_cov";

  TString z_title = h->GetZaxis()->GetTitle();
  if (z_title.Index(";") >= 0)
    z_title = z_title(0, z_title.Index(";"));

  TH1D *cov =
      new TH1D(name, title + ";" + z_title + ";Covariance", h->GetNbinsZ(),
               h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax());

  for (int i = 1; i <= h->GetNbinsZ(); ++i) {
    h->GetZaxis()->SetRange(i, i);
    TH2D *h2 = (TH2D *)h->Project3D("xy");
    h2->SetDirectory(0);
    h2->SetName(Form("%s_slice%d", h->GetName(), i));

    double c = getCovariance(h2);
    cov->SetBinContent(i, c);

    if (h2->GetEntries() > 0) {
      double err = 1.0 / std::sqrt(h2->GetEntries());
      cov->SetBinError(i, err);
    }

    delete h2;
  }

  h->GetZaxis()->SetRange(0, 0); // reset
  return cov;
}

// Label helper (if you want to draw later)
void drawLabel(TPad *pad, float x = 0.57, TString extra = "") {
  pad->cd();
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->DrawLatex(x, 0.85, "PYTHIA 8.315");
  latex->DrawLatex(x, 0.80, "#it{pp} at #sqrt{#it{s}} = 200 GeV");
  latex->DrawLatex(x, 0.75, "|#eta|<0.6, anti-k_{T}, #it{R}=0.4");
  latex->DrawLatex(x, 0.70, "jet p_{T} > 3 GeV/c");
  latex->DrawLatex(x, 0.65, "p_{T,track}>0.15 GeV/c, |#eta_{track}|<1");
  latex->DrawLatex(x, 0.60, "|#phi_{1} - #phi_{2}| > 3#pi/4");
  latex->DrawLatex(x, 0.55, extra);
}

// ============================================================================
// Main analysis + joint (N1,N2) unfolding vs pT
// ============================================================================

void unfold(const char *inputFile = "embedded_dijets.root",
            const char *outputFile = "analysis_unfoldN_vsPt.root") {

  std::cout
      << "\n=== Embedded Dijet Analysis with joint N1,N2 unfolding vs pT ===\n";
  std::cout << "Input:  " << inputFile << "\n";
  std::cout << "Output: " << outputFile << "\n\n";

  // Open input
  TFile *fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: Cannot open input file!\n";
    return;
  }

  TTree *tMC = (TTree *)fin->Get("mcTruth");
  TTree *tReco = (TTree *)fin->Get("recoJets");
  if (!tMC || !tReco) {
    std::cerr << "ERROR: Cannot find required trees!\n";
    fin->Close();
    return;
  }

  Long64_t nEvents = tMC->GetEntries();
  std::cout << "Total MC events: " << nEvents << "\n\n";

  // Output file
  TFile *fout = TFile::Open(outputFile, "RECREATE");

  // Binning
  const int nMultBins = 25;
  const double multMin = 0.0;
  const double multMax = 25.0;

  const int nPtBins = 200;
  const double ptMin = 0.0;
  const double ptMax = 60.0;

  const int nEtaBins = 20;
  const double etaMin = -1.0;
  const double etaMax = 1.0;

  const int nPhiBins = 32;
  const double phiMin = -M_PI;
  const double phiMax = M_PI;

  const int nDeltaBins = 40;
  const double deltaMin = -20.0;
  const double deltaMax = 20.0;

  const int nRatioBins = 50;
  const double ratioMin = 0.0;
  const double ratioMax = 2.0;

  const int nDRBins = 50;
  const double dRMin = 0.0;
  const double dRMax = 0.5;

  // =======================================================================
  // MC truth histograms
  // =======================================================================
  fout->mkdir("mcTruth");
  fout->cd("mcTruth");

  TH1D *hMC_jet1_pt =
      new TH1D("hMC_jet1_pt", "MC Truth Jet 1 p_{T};p_{T} [GeV/c];Events",
               nPtBins, ptMin, ptMax);
  TH1D *hMC_jet2_pt =
      new TH1D("hMC_jet2_pt", "MC Truth Jet 2 p_{T};p_{T} [GeV/c];Events",
               nPtBins, ptMin, ptMax);

  TH1D *hMC_jet1_eta =
      new TH1D("hMC_jet1_eta", "MC Truth Jet 1 #eta;#eta;Events", nEtaBins,
               etaMin, etaMax);
  TH1D *hMC_jet2_eta =
      new TH1D("hMC_jet2_eta", "MC Truth Jet 2 #eta;#eta;Events", nEtaBins,
               etaMin, etaMax);

  TH1D *hMC_jet1_phi =
      new TH1D("hMC_jet1_phi", "MC Truth Jet 1 #phi;#phi;Events", nPhiBins,
               phiMin, phiMax);
  TH1D *hMC_jet2_phi =
      new TH1D("hMC_jet2_phi", "MC Truth Jet 2 #phi;#phi;Events", nPhiBins,
               phiMin, phiMax);

  TH1D *hMC_jet1_nPart =
      new TH1D("hMC_jet1_nPart", "MC Truth Jet 1 N_{part};N_{part};Events",
               nMultBins, multMin, multMax);
  TH1D *hMC_jet2_nPart =
      new TH1D("hMC_jet2_nPart", "MC Truth Jet 2 N_{part};N_{part};Events",
               nMultBins, multMin, multMax);

  TH2D *hMC_deltaPhi = new TH2D(
      "hMC_deltaPhi", "MC Truth #Delta#phi Distribution;#phi_{1};#phi_{2}",
      nPhiBins, phiMin, phiMax, nPhiBins, phiMin, phiMax);

  TH1D *hMC_ptBalance = new TH1D(
      "hMC_ptBalance", "MC Truth p_{T} Balance;p_{T,2} / p_{T,1};Events",
      nRatioBins, ratioMin, ratioMax);

  TH2D *hMC_N1vsN2 =
      new TH2D("hMC_N1vsN2",
               "MC Truth: N_{part}^{jet1} vs N_{part}^{jet2};"
               "N_{part}^{jet1};N_{part}^{jet2}",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  // TH3D: N1 vs N2 vs pT (pT = truth jet1 pT)
  TH3D *hMC_N1vsN2vsPt =
      new TH3D("hMC_N1vsN2vsPt",
               "MC Truth: N1 vs N2 vs p_{T};"
               "N_{part}^{jet1};N_{part}^{jet2};p_{T}^{jet1} [GeV/c]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax,
               nPtBins, ptMin, ptMax);

  // -----------------------------------------------------------------------
  // Fill MC truth histograms
  // -----------------------------------------------------------------------
  std::cout << "Filling MC truth histograms...\n";

  int eventID;
  double mc_jet1_pt, mc_jet1_eta, mc_jet1_phi;
  int mc_jet1_nPart;
  double mc_jet2_pt, mc_jet2_eta, mc_jet2_phi;
  int mc_jet2_nPart;

  tMC->SetBranchAddress("eventID", &eventID);
  tMC->SetBranchAddress("mc_jet1_pt", &mc_jet1_pt);
  tMC->SetBranchAddress("mc_jet1_eta", &mc_jet1_eta);
  tMC->SetBranchAddress("mc_jet1_phi", &mc_jet1_phi);
  tMC->SetBranchAddress("mc_jet1_nPart", &mc_jet1_nPart);
  tMC->SetBranchAddress("mc_jet2_pt", &mc_jet2_pt);
  tMC->SetBranchAddress("mc_jet2_eta", &mc_jet2_eta);
  tMC->SetBranchAddress("mc_jet2_phi", &mc_jet2_phi);
  tMC->SetBranchAddress("mc_jet2_nPart", &mc_jet2_nPart);

  for (Long64_t i = 0; i < tMC->GetEntries(); ++i) {
    tMC->GetEntry(i);

    hMC_jet1_pt->Fill(mc_jet1_pt);
    hMC_jet2_pt->Fill(mc_jet2_pt);
    hMC_jet1_eta->Fill(mc_jet1_eta);
    hMC_jet2_eta->Fill(mc_jet2_eta);
    hMC_jet1_phi->Fill(mc_jet1_phi);
    hMC_jet2_phi->Fill(mc_jet2_phi);
    hMC_jet1_nPart->Fill(mc_jet1_nPart);
    hMC_jet2_nPart->Fill(mc_jet2_nPart);

    hMC_deltaPhi->Fill(mc_jet1_phi, mc_jet2_phi);
    if (mc_jet1_pt > 0.0)
      hMC_ptBalance->Fill(mc_jet2_pt / mc_jet1_pt);

    hMC_N1vsN2->Fill(mc_jet1_nPart, mc_jet2_nPart);
    hMC_N1vsN2vsPt->Fill(mc_jet1_nPart, mc_jet2_nPart, mc_jet1_pt);
  }

  // =======================================================================
  // Reconstructed-level histograms
  // =======================================================================
  fout->mkdir("recoJets");
  fout->cd("recoJets");

  TH1D *hReco_matchedFraction =
      new TH1D("hReco_matchedFraction", "Jet Matching;Status;Events", 4, 0, 4);
  hReco_matchedFraction->GetXaxis()->SetBinLabel(1, "None");
  hReco_matchedFraction->GetXaxis()->SetBinLabel(2, "Jet1 only");
  hReco_matchedFraction->GetXaxis()->SetBinLabel(3, "Jet2 only");
  hReco_matchedFraction->GetXaxis()->SetBinLabel(4, "Both");

  TH1D *hReco_dR1 = new TH1D(
      "hReco_dR1", "Matching #DeltaR for Jet 1;#DeltaR(MC, Reco);Events",
      nDRBins, dRMin, dRMax);
  TH1D *hReco_dR2 = new TH1D(
      "hReco_dR2", "Matching #DeltaR for Jet 2;#DeltaR(MC, Reco);Events",
      nDRBins, dRMin, dRMax);

  TH1D *hReco_jet1_pt =
      new TH1D("hReco_jet1_pt", "Reco Jet 1 p_{T};p_{T} [GeV/c];Events",
               nPtBins, ptMin, ptMax);
  TH1D *hReco_jet2_pt =
      new TH1D("hReco_jet2_pt", "Reco Jet 2 p_{T};p_{T} [GeV/c];Events",
               nPtBins, ptMin, ptMax);

  TH1D *hReco_jet1_nPart =
      new TH1D("hReco_jet1_nPart", "Reco Jet 1 N_{part};N_{part};Events",
               nMultBins, multMin, multMax);
  TH1D *hReco_jet2_nPart =
      new TH1D("hReco_jet2_nPart", "Reco Jet 2 N_{part};N_{part};Events",
               nMultBins, multMin, multMax);

  TH2D *hReco_N1vsN2 =
      new TH2D("hReco_N1vsN2",
               "Reco: N_{part}^{jet1} vs N_{part}^{jet2};"
               "N_{part}^{jet1};N_{part}^{jet2}",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  // TH3D reco: N1 vs N2 vs pT (pT = truth jet1 pT for consistency)
  TH3D *hReco_N1vsN2vsPt =
      new TH3D("hReco_N1vsN2vsPt",
               "Reco: N1 vs N2 vs p_{T};"
               "N_{part}^{jet1};N_{part}^{jet2};p_{T}^{jet1} [GeV/c]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax,
               nPtBins, ptMin, ptMax);

  // =======================================================================
  // Response and closure histograms (1D pT and 1D N)
  // =======================================================================
  fout->mkdir("response");
  fout->cd("response");

  TH2D *hResp_pt1 =
      new TH2D("hResp_pt1",
               "Jet 1 p_{T} Response;p_{T}^{MC} [GeV/c];p_{T}^{Reco} [GeV/c]",
               nPtBins, ptMin, ptMax, nPtBins, ptMin, ptMax);
  TH2D *hResp_pt2 =
      new TH2D("hResp_pt2",
               "Jet 2 p_{T} Response;p_{T}^{MC} [GeV/c];p_{T}^{Reco} [GeV/c]",
               nPtBins, ptMin, ptMax, nPtBins, ptMin, ptMax);

  TH2D *hResp_N1 = new TH2D(
      "hResp_N1", "Jet 1 N_{part} Response;N_{part}^{MC};N_{part}^{Reco}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);
  TH2D *hResp_N2 = new TH2D(
      "hResp_N2", "Jet 2 N_{part} Response;N_{part}^{MC};N_{part}^{Reco}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  TH1D *hDelta_pt1 =
      new TH1D("hDelta_pt1",
               "Jet 1 #Deltap_{T};p_{T}^{Reco} - p_{T}^{MC} [GeV/c];Events",
               nDeltaBins, deltaMin, deltaMax);
  TH1D *hDelta_pt2 =
      new TH1D("hDelta_pt2",
               "Jet 2 #Deltap_{T};p_{T}^{Reco} - p_{T}^{MC} [GeV/c];Events",
               nDeltaBins, deltaMin, deltaMax);

  TH1D *hDelta_N1 = new TH1D(
      "hDelta_N1", "Jet 1 #DeltaN;N_{part}^{Reco} - N_{part}^{MC};Events",
      nDeltaBins, deltaMin, deltaMax);
  TH1D *hDelta_N2 = new TH1D(
      "hDelta_N2", "Jet 2 #DeltaN;N_{part}^{Reco} - N_{part}^{MC};Events",
      nDeltaBins, deltaMin, deltaMax);

  TH2D *hResp_N1vsN2_MC =
      new TH2D("hResp_N1vsN2_MC", "MC: N_{1} vs N_{2};N_{1}^{MC};N_{2}^{MC}",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax);
  TH2D *hResp_N1vsN2_Reco = new TH2D(
      "hResp_N1vsN2_Reco", "Reco: N_{1} vs N_{2};N_{1}^{Reco};N_{2}^{Reco}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  // =======================================================================
  // Background histograms
  // =======================================================================
  fout->mkdir("background");
  fout->cd("background");

  TH1D *hBkg_multA =
      new TH1D("hBkg_multA", "Background Multiplicity A;N_{ch}^{A};Events",
               nMultBins, multMin, multMax);
  TH1D *hBkg_multB =
      new TH1D("hBkg_multB", "Background Multiplicity B;N_{ch}^{B};Events",
               nMultBins, multMin, multMax);

  TH2D *hBkg_multAvsB = new TH2D(
      "hBkg_multAvsB", "Background: N_{A} vs N_{B};N_{ch}^{A};N_{ch}^{B}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  TH3D *hBkg_N1vsN2vsPt = new TH3D("hBkg_N1vsN2vsPt",
                                   "Background: N_{A} vs N_{B} vs p_{T};"
                                   "N_{ch}^{A};N_{ch}^{B};p_{T}^{jet1} [GeV/c]",
                                   nMultBins, multMin, multMax, nMultBins,
                                   multMin, multMax, nPtBins, ptMin, ptMax);

  // =======================================================================
  // Additional objects for joint (N1,N2) unfolding vs pT
  // =======================================================================
  fout->cd("response");

  const int nStates = nMultBins * nMultBins;

  // Global state response (integrated over pT)
  TH1D *hTruth_state_global = new TH1D(
      "hTruth_state_global", "Truth state index (N1,N2);state_{MC};Events",
      nStates, 0.5, nStates + 0.5);

  TH1D *hReco_state_global = new TH1D(
      "hReco_state_global", "Reco state index (N1,N2);state_{Reco};Events",
      nStates, 0.5, nStates + 0.5);

  TH2D *hResp_state_global = new TH2D(
      "hResp_state_global", "Response in state space;state_{Reco};state_{MC}",
      nStates, 0.5, nStates + 0.5, nStates, 0.5, nStates + 0.5);

  // State distributions per truth pT bin
  std::vector<TH1D *> hTruth_state_pt(nPtBins, (TH1D *)nullptr);
  std::vector<TH1D *> hReco_state_pt(nPtBins, (TH1D *)nullptr);

  for (int ipt = 0; ipt < nPtBins; ++ipt) {
    hTruth_state_pt[ipt] = new TH1D(
        Form("hTruth_state_ptBin%d", ipt + 1),
        Form("Truth state index in pT bin %d;state_{MC};Events", ipt + 1),
        nStates, 0.5, nStates + 0.5);

    hReco_state_pt[ipt] = new TH1D(
        Form("hReco_state_ptBin%d", ipt + 1),
        Form("Reco state index in pT bin %d;state_{Reco};Events", ipt + 1),
        nStates, 0.5, nStates + 0.5);
  }

  // =======================================================================
  // Fill reconstructed + response + state histograms
  // =======================================================================
  std::cout << "Filling reconstruction, response and state histograms...\n";

  int reco_eventID;
  bool matched1, matched2;
  double reco_jet1_pt, reco_jet1_eta, reco_jet1_phi;
  int reco_jet1_nPart;
  double reco_jet2_pt, reco_jet2_eta, reco_jet2_phi;
  int reco_jet2_nPart;
  double dR_mc1_reco1, dR_mc2_reco2;
  int bkg_multA, bkg_multB;

  tReco->SetBranchAddress("eventID", &reco_eventID);
  tReco->SetBranchAddress("matched1", &matched1);
  tReco->SetBranchAddress("matched2", &matched2);
  tReco->SetBranchAddress("mc_jet1_pt", &mc_jet1_pt);
  tReco->SetBranchAddress("mc_jet1_eta", &mc_jet1_eta);
  tReco->SetBranchAddress("mc_jet1_phi", &mc_jet1_phi);
  tReco->SetBranchAddress("mc_jet1_nPart", &mc_jet1_nPart);
  tReco->SetBranchAddress("reco_jet1_pt", &reco_jet1_pt);
  tReco->SetBranchAddress("reco_jet1_eta", &reco_jet1_eta);
  tReco->SetBranchAddress("reco_jet1_phi", &reco_jet1_phi);
  tReco->SetBranchAddress("reco_jet1_nPart", &reco_jet1_nPart);
  tReco->SetBranchAddress("dR_mc1_reco1", &dR_mc1_reco1);

  tReco->SetBranchAddress("mc_jet2_pt", &mc_jet2_pt);
  tReco->SetBranchAddress("mc_jet2_eta", &mc_jet2_eta);
  tReco->SetBranchAddress("mc_jet2_phi", &mc_jet2_phi);
  tReco->SetBranchAddress("mc_jet2_nPart", &mc_jet2_nPart);
  tReco->SetBranchAddress("reco_jet2_pt", &reco_jet2_pt);
  tReco->SetBranchAddress("reco_jet2_eta", &reco_jet2_eta);
  tReco->SetBranchAddress("reco_jet2_phi", &reco_jet2_phi);
  tReco->SetBranchAddress("reco_jet2_nPart", &reco_jet2_nPart);
  tReco->SetBranchAddress("dR_mc2_reco2", &dR_mc2_reco2);

  tReco->SetBranchAddress("multA", &bkg_multA);
  tReco->SetBranchAddress("multB", &bkg_multB);

  int nBothMatched = 0;
  int nOneMatched = 0;
  int nNoneMatched = 0;

  for (Long64_t i = 0; i < tReco->GetEntries(); ++i) {
    tReco->GetEntry(i);

    int matchStatus = 0;
    if (matched1 && matched2) {
      matchStatus = 3;
      ++nBothMatched;
    } else if (matched1 || matched2) {
      matchStatus = matched1 ? 1 : 2;
      ++nOneMatched;
    } else {
      ++nNoneMatched;
    }
    hReco_matchedFraction->Fill(matchStatus);

    if (matched1) {
      hReco_dR1->Fill(dR_mc1_reco1);
      hReco_jet1_pt->Fill(reco_jet1_pt);
      hReco_jet1_nPart->Fill(reco_jet1_nPart);

      hResp_pt1->Fill(mc_jet1_pt, reco_jet1_pt);
      hResp_N1->Fill(mc_jet1_nPart, reco_jet1_nPart);

      hDelta_pt1->Fill(reco_jet1_pt - mc_jet1_pt);
      hDelta_N1->Fill(reco_jet1_nPart - mc_jet1_nPart);
    }

    if (matched2) {
      hReco_dR2->Fill(dR_mc2_reco2);
      hReco_jet2_pt->Fill(reco_jet2_pt);
      hReco_jet2_nPart->Fill(reco_jet2_nPart);

      hResp_pt2->Fill(mc_jet2_pt, reco_jet2_pt);
      hResp_N2->Fill(mc_jet2_nPart, reco_jet2_nPart);

      hDelta_pt2->Fill(reco_jet2_pt - mc_jet2_pt);
      hDelta_N2->Fill(reco_jet2_nPart - mc_jet2_nPart);
    }

    if (matched1 && matched2) {
      hReco_N1vsN2->Fill(reco_jet1_nPart, reco_jet2_nPart);
      // use truth jet1 pT for z-axis
      hReco_N1vsN2vsPt->Fill(reco_jet1_nPart, reco_jet2_nPart, mc_jet1_pt);

      hResp_N1vsN2_MC->Fill(mc_jet1_nPart, mc_jet2_nPart);
      hResp_N1vsN2_Reco->Fill(reco_jet1_nPart, reco_jet2_nPart);

      // ---------- Joint (N1,N2) state mapping for unfolding ----------
      int binN1_MC = hMC_jet1_nPart->GetXaxis()->FindBin(mc_jet1_nPart);
      int binN2_MC = hMC_jet2_nPart->GetXaxis()->FindBin(mc_jet2_nPart);
      int binN1_Reco = hReco_jet1_nPart->GetXaxis()->FindBin(reco_jet1_nPart);
      int binN2_Reco = hReco_jet2_nPart->GetXaxis()->FindBin(reco_jet2_nPart);

      if (binN1_MC < 1 || binN1_MC > nMultBins)
        continue;
      if (binN2_MC < 1 || binN2_MC > nMultBins)
        continue;
      if (binN1_Reco < 1 || binN1_Reco > nMultBins)
        continue;
      if (binN2_Reco < 1 || binN2_Reco > nMultBins)
        continue;

      int state_MC = binN1_MC + (binN2_MC - 1) * nMultBins;
      int state_Reco = binN1_Reco + (binN2_Reco - 1) * nMultBins;

      hTruth_state_global->Fill(state_MC);
      hReco_state_global->Fill(state_Reco);
      hResp_state_global->Fill(state_Reco, state_MC);

      // assign state to truth pT bin
      int iptTruth = hMC_jet1_pt->GetXaxis()->FindBin(mc_jet1_pt);
      if (iptTruth >= 1 && iptTruth <= nPtBins) {
        hTruth_state_pt[iptTruth - 1]->Fill(state_MC);
        hReco_state_pt[iptTruth - 1]->Fill(state_Reco);
      }

      // background TH3D (vs reco jet1 pt)
      hBkg_N1vsN2vsPt->Fill(bkg_multA, bkg_multB, reco_jet1_pt);
    }

    hBkg_multA->Fill(bkg_multA);
    hBkg_multB->Fill(bkg_multB);
    hBkg_multAvsB->Fill(bkg_multA, bkg_multB);
  }

  std::cout << "Matched: both=" << nBothMatched << ", one=" << nOneMatched
            << ", none=" << nNoneMatched << "\n";

  // =======================================================================
  // Joint (N1,N2) unfolding vs truth pT using flattened state index
  // =======================================================================
  fout->cd("response");

  RooUnfoldResponse resp_state(hReco_state_global, hTruth_state_global,
                               hResp_state_global);

  const int nIter_ = 4; // iterations for Bayes

  TH3D *hUnfold_N1vsN2vsPt =
      new TH3D("hUnfold_N1vsN2vsPt",
               "Unfolded N_{1} vs N_{2} vs p_{T};"
               "N_{part}^{jet1};N_{part}^{jet2};p_{T}^{jet1} [GeV/c]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax,
               nPtBins, ptMin, ptMax);

  for (int ipt = 0; ipt < nPtBins; ++ipt) {
    TH1D *hMeas_state = hReco_state_pt[ipt];
    if (!hMeas_state || hMeas_state->GetEntries() <= 0)
      continue;

    RooUnfoldBayes unfold_state_pt(&resp_state, hMeas_state, nIter);
    TH1D *hUnfold_state_pt = (TH1D *)unfold_state_pt.Hunfold();
    hUnfold_state_pt->SetDirectory(0);

    for (int k = 1; k <= nStates; ++k) {
      double content = hUnfold_state_pt->GetBinContent(k);
      if (content <= 0.0)
        continue;

      int idx = k - 1;
      int binN1 = (idx % nMultBins) + 1;
      int binN2 = (idx / nMultBins) + 1;

      double old = hUnfold_N1vsN2vsPt->GetBinContent(binN1, binN2, ipt + 1);
      hUnfold_N1vsN2vsPt->SetBinContent(binN1, binN2, ipt + 1, old + content);
    }

    delete hUnfold_state_pt;
  }

  // =======================================================================
  // Covariance vs pT from TH3D: MC, Reco, Unfolded, Background
  // =======================================================================
  fout->cd();

  TH1D *hCov_MC = getCovariance(
      hMC_N1vsN2vsPt, "MC Truth COV(N_{1}, N_{2}) vs p_{T}", "hCov_MC");
  TH1D *hCov_Reco = getCovariance(
      hReco_N1vsN2vsPt, "Reco COV(N_{1}, N_{2}) vs p_{T}", "hCov_Reco");
  TH1D *hCov_Unfold = getCovariance(
      hUnfold_N1vsN2vsPt, "Unfolded COV(N_{1}, N_{2}) vs p_{T}", "hCov_Unfold");
  TH1D *hCov_Bkg = getCovariance(
      hBkg_N1vsN2vsPt, "Background COV(N_{A}, N_{B}) vs p_{T}", "hCov_Bkg");

  // Canvas 3: Covariance comparison
  TCanvas *can = new TCanvas("can", "Covariance vs pT", 800, 600);

  hCov_MC->GetYaxis()->SetRangeUser(-0.2, hCov_MC->GetMaximum() * 2);
  hCov_MC->SetLineWidth(2);
  hCov_Reco->SetLineWidth(2);
  hCov_Bkg->SetLineWidth(2);

  hCov_MC->Draw("E1");
  hCov_Reco->Draw("E1 same");
  hCov_Bkg->Draw("E1 same");

  {
    TLegend *leg = new TLegend(0.6, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hCov_MC, "MC Truth", "lep");
    leg->AddEntry(hCov_Reco, "Reconstructed", "lep");
    leg->AddEntry(hCov_Unfold, "Unfolded", "lep");
    leg->AddEntry(hCov_Bkg, "Background", "lep");
    leg->Draw();
  }
  drawLabel((TPad *)gPad, 0.15, "Covariance vs p_{T}");

  TCanvas *can2 =
      new TCanvas("covRation", "Covariance True/Unfolded", 800, 600);
  TH1D *hCov_Ratio = (TH1D *)hCov_Unfold->Clone("hCov_Ratio");
  hCov_Ratio->SetTitle(
      "Covariance Ratio: Unfolded / MC Truth;p_{T} [GeV/c];Ratio");
  hCov_Ratio->Divide(hCov_MC);
  hCov_Ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hCov_Ratio->SetLineWidth(2);
  hCov_Ratio->Draw("E1");
  drawLabel((TPad *)gPad, 0.15, "Covariance Ratio vs p_{T}");

  // =======================================================================
  // Write everything
  // =======================================================================

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. Results saved to " << outputFile << "\n";
}
