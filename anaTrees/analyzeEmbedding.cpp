// analyzeEmbedding.cc
// Analysis macro for embedded MC dijet closure tests
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TTree.h"
#include <TColor.h>
#include <iostream>
#include <vector>

// ============================================================================
// Utility Functions
// ============================================================================

double getEntropy(TH1D *h) {
  double entropy = 0;
  int nBins = h->GetNbinsX();
  double total = h->Integral();
  if (total <= 0)
    return 0;

  for (int i = 1; i <= nBins; ++i) {
    double p = h->GetBinContent(i) / total;
    if (p > 0) {
      entropy -= p * std::log(p);
    }
  }
  return entropy;
}

double getEntropy(TH2D *h) {
  double entropy = 0;
  int nBinsX = h->GetNbinsX();
  int nBinsY = h->GetNbinsY();
  double total = h->Integral();
  if (total <= 0)
    return 0;

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      double p = h->GetBinContent(i, j) / total;
      if (p > 0) {
        entropy -= p * std::log(p);
      }
    }
  }
  return entropy;
}

double getCovariance(TH2D *h) {
  double S1 = getEntropy((TH1D *)h->ProjectionX("_px", 1, -1, "e"));
  double S2 = getEntropy((TH1D *)h->ProjectionY("_py", 1, -1, "e"));
  double S12 = getEntropy(h);
  return S1 + S2 - S12;
}

TH1D *getCovariance(TH3D *h, TString title = "", TString name = "") {
  if (name == "")
    name = TString(h->GetName()) + "_cov";

  TString z_title = h->GetZaxis()->GetTitle();
  // Strip off everything after ; in z_title
  if (z_title.Index(";") >= 0)
    z_title = z_title(0, z_title.Index(";"));

  TH1D *cov =
      new TH1D(name, title + ";" + z_title + ";Covariance", h->GetNbinsZ(),
               h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax());

  for (int i = 1; i <= h->GetNbinsZ(); ++i) {
    h->GetZaxis()->SetRange(i, i);
    TH2D *h2 = (TH2D *)h->Project3D("xy");
    h2->SetName(Form("%s_slice%d", h->GetName(), i));

    double c = getCovariance(h2);
    cov->SetBinContent(i, c);

    if (h2->GetEntries() > 0) {
      double err = 1.0 / std::sqrt(h2->GetEntries());
      cov->SetBinError(i, err);
    }

    delete h2;
  }

  h->GetZaxis()->SetRange(0, 0); // Reset range
  return cov;
}

// ============================================================================
// Main Analysis
// ============================================================================

void analyzeEmbedding(const char *inputFile = "embedded_dijets.root",
                      const char *outputFile = "analysis_output.root") {

  std::cout << "\n=== Embedded Dijet Analysis ===\n";
  std::cout << "Input:  " << inputFile << "\n";
  std::cout << "Output: " << outputFile << "\n\n";

  // Open input file
  TFile *fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: Cannot open input file!\n";
    return;
  }

  // Get trees
  TTree *tMC = (TTree *)fin->Get("mcTruth");
  TTree *tReco = (TTree *)fin->Get("recoJets");

  if (!tMC || !tReco) {
    std::cerr << "ERROR: Cannot find required trees!\n";
    fin->Close();
    return;
  }

  Long64_t nEvents = tMC->GetEntries();
  std::cout << "Total events: " << nEvents << "\n\n";

  // Set up styles
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  TH3::SetDefaultSumw2(true);

  // Create output file
  TFile *fout = TFile::Open(outputFile, "RECREATE");

  // =========================================================================
  // Histogram Binning
  // =========================================================================
  const int nMultBins = 25;
  const double multMin = 0;
  const double multMax = 25;

  const int nPtBins = 20;
  const double ptMin = 0;
  const double ptMax = 60;

  const int nEtaBins = 20;
  const double etaMin = -1.0;
  const double etaMax = 1.0;

  const int nPhiBins = 32;
  const double phiMin = 0;
  const double phiMax = 2.0 * M_PI;

  const int nDeltaBins = 40;
  const double deltaMin = -20;
  const double deltaMax = 20;

  const int nRatioBins = 50;
  const double ratioMin = 0;
  const double ratioMax = 2.0;

  const int nDRBins = 50;
  const double dRMin = 0;
  const double dRMax = 0.5;

  // =========================================================================
  // MC Truth Level Histograms
  // =========================================================================
  fout->mkdir("mcTruth");
  fout->cd("mcTruth");

  // Jet kinematics
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

  // Multiplicity
  TH1D *hMC_jet1_nPart =
      new TH1D("hMC_jet1_nPart", "MC Truth Jet 1 N_{part};N_{part};Events",
               nMultBins, multMin, multMax);
  TH1D *hMC_jet2_nPart =
      new TH1D("hMC_jet2_nPart", "MC Truth Jet 2 N_{part};N_{part};Events",
               nMultBins, multMin, multMax);

  // Dijet correlations
  TH2D *hMC_deltaPhi = new TH2D(
      "hMC_deltaPhi", "MC Truth #Delta#phi Distribution;#phi_{1};#phi_{2}",
      nPhiBins, phiMin, phiMax, nPhiBins, phiMin, phiMax);

  TH1D *hMC_ptBalance = new TH1D(
      "hMC_ptBalance", "MC Truth p_{T} Balance;p_{T,2} / p_{T,1};Events",
      nRatioBins, ratioMin, ratioMax);

  // KEY PLOT: N1 vs N2 correlation (should be diagonal)
  TH2D *hMC_N1vsN2 =
      new TH2D("hMC_N1vsN2",
               "MC Truth: N_{part}^{jet1} vs "
               "N_{part}^{jet2};N_{part}^{jet1};N_{part}^{jet2}",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  // N1 vs N2 vs pT (for covariance analysis)
  TH3D *hMC_N1vsN2vsPt =
      new TH3D("hMC_N1vsN2vsPt",
               "MC Truth: N1 vs N2 vs "
               "p_{T};N_{part}^{jet1};N_{part}^{jet2};p_{T}^{jet1} [GeV/c]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax,
               nPtBins, ptMin, ptMax);

  // =========================================================================
  // Reconstructed Level Histograms
  // =========================================================================
  fout->mkdir("recoJets");
  fout->cd("recoJets");

  // Matching quality
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

  // Reconstructed kinematics (matched jets only)
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

  // Reco N1 vs N2 correlation
  TH2D *hReco_N1vsN2 =
      new TH2D("hReco_N1vsN2",
               "Reco: N_{part}^{jet1} vs "
               "N_{part}^{jet2};N_{part}^{jet1};N_{part}^{jet2}",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  TH3D *hReco_N1vsN2vsPt =
      new TH3D("hReco_N1vsN2vsPt",
               "Reco: N1 vs N2 vs "
               "p_{T};N_{part}^{jet1};N_{part}^{jet2};p_{T}^{jet1} [GeV/c]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax,
               nPtBins, ptMin, ptMax);

  // =========================================================================
  // Response and Closure Histograms
  // =========================================================================
  fout->mkdir("response");
  fout->cd("response");

  // pT response
  TH2D *hResp_pt1 =
      new TH2D("hResp_pt1",
               "Jet 1 p_{T} Response;p_{T}^{MC} [GeV/c];p_{T}^{Reco} [GeV/c]",
               nPtBins, ptMin, ptMax, nPtBins, ptMin, ptMax);
  TH2D *hResp_pt2 =
      new TH2D("hResp_pt2",
               "Jet 2 p_{T} Response;p_{T}^{MC} [GeV/c];p_{T}^{Reco} [GeV/c]",
               nPtBins, ptMin, ptMax, nPtBins, ptMin, ptMax);

  TProfile *hProf_pt1 = new TProfile(
      "hProf_pt1",
      "Jet 1 p_{T} Response Profile;p_{T}^{MC} [GeV/c];<p_{T}^{Reco}> [GeV/c]",
      nPtBins, ptMin, ptMax);
  TProfile *hProf_pt2 = new TProfile(
      "hProf_pt2",
      "Jet 2 p_{T} Response Profile;p_{T}^{MC} [GeV/c];<p_{T}^{Reco}> [GeV/c]",
      nPtBins, ptMin, ptMax);

  // Multiplicity response (closure test)
  TH2D *hResp_N1 = new TH2D(
      "hResp_N1", "Jet 1 N_{part} Response;N_{part}^{MC};N_{part}^{Reco}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);
  TH2D *hResp_N2 = new TH2D(
      "hResp_N2", "Jet 2 N_{part} Response;N_{part}^{MC};N_{part}^{Reco}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  TProfile *hProf_N1 = new TProfile(
      "hProf_N1",
      "Jet 1 N_{part} Response Profile;N_{part}^{MC};<N_{part}^{Reco}>",
      nMultBins, multMin, multMax);
  TProfile *hProf_N2 = new TProfile(
      "hProf_N2",
      "Jet 2 N_{part} Response Profile;N_{part}^{MC};<N_{part}^{Reco}>",
      nMultBins, multMin, multMax);

  // Differences (for resolution studies)
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

  // 2D response for correlation studies
  TH2D *hResp_N1vsN2_MC =
      new TH2D("hResp_N1vsN2_MC", "MC: N_{1} vs N_{2};N_{1}^{MC};N_{2}^{MC}",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax);
  TH2D *hResp_N1vsN2_Reco = new TH2D(
      "hResp_N1vsN2_Reco", "Reco: N_{1} vs N_{2};N_{1}^{Reco};N_{2}^{Reco}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  // =========================================================================
  // Background Histograms
  // =========================================================================

  TH1D *hBkg_multA =
      new TH1D("hBkg_multA", "Background Multiplicity A;N_{ch}^{A};Events",
               nMultBins, multMin, multMax);
  TH1D *hBkg_multB =
      new TH1D("hBkg_multB", "Background Multiplicity B;N_{ch}^{B};Events",
               nMultBins, multMin, multMax);

  TH2D *hBkg_multAvsB = new TH2D(
      "hBkg_multAvsB", "Background: N_{A} vs N_{B};N_{ch}^{A};N_{ch}^{B}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  TH3D *hBkg_N1vsN2vsPt =
      new TH3D("hBkg_N1vsN2vsPt",
               "Background: N_{A} vs N_{B} vs "
               "p_{T};N_{ch}^{A};N_{ch}^{B};p_{T}^{jet1} [GeV/c]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax,
               nPtBins, ptMin, ptMax);

  // =========================================================================
  // Fill MC Truth Histograms
  // =========================================================================
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
    hMC_ptBalance->Fill(mc_jet2_pt / mc_jet1_pt);

    // KEY: N1 vs N2 correlation
    hMC_N1vsN2->Fill(mc_jet1_nPart, mc_jet2_nPart);
    hMC_N1vsN2vsPt->Fill(mc_jet1_nPart, mc_jet2_nPart, mc_jet1_pt);
  }

  // =========================================================================
  // Fill Reconstructed and Response Histograms
  // =========================================================================
  std::cout << "Filling reconstruction and response histograms...\n";

  int reco_eventID;
  bool matched1, matched2;
  double reco_jet1_pt, reco_jet1_eta, reco_jet1_phi;
  int reco_jet1_nPart;
  double reco_jet2_pt, reco_jet2_eta, reco_jet2_phi;
  int reco_jet2_nPart;
  double dR_mc1_reco1, dR_mc2_reco2;

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

  int bkg_eventID;
  int bkg_multA, bkg_multB;
  double bkg_ptA, bkg_ptB;

  tReco->SetBranchAddress("eventID", &bkg_eventID);
  tReco->SetBranchAddress("multA", &bkg_multA);
  tReco->SetBranchAddress("multB", &bkg_multB);
  tReco->SetBranchAddress("ptA", &bkg_ptA);
  tReco->SetBranchAddress("ptB", &bkg_ptB);

  int nBothMatched = 0;
  int nOneMatched = 0;
  int nNoneMatched = 0;

  for (Long64_t i = 0; i < tReco->GetEntries(); ++i) {
    tReco->GetEntry(i);

    // Matching statistics
    int matchStatus = 0;
    if (matched1 && matched2) {
      matchStatus = 3;
      nBothMatched++;
    } else if (matched1 || matched2) {
      matchStatus = matched1 ? 1 : 2;
      nOneMatched++;
    } else {
      nNoneMatched++;
    }
    hReco_matchedFraction->Fill(matchStatus);

    if (matched1) {
      hReco_dR1->Fill(dR_mc1_reco1);
      hReco_jet1_pt->Fill(reco_jet1_pt);
      hReco_jet1_nPart->Fill(reco_jet1_nPart);

      // Response
      hResp_pt1->Fill(mc_jet1_pt, reco_jet1_pt);
      hProf_pt1->Fill(mc_jet1_pt, reco_jet1_pt);
      hResp_N1->Fill(mc_jet1_nPart, reco_jet1_nPart);
      hProf_N1->Fill(mc_jet1_nPart, reco_jet1_nPart);

      // Differences
      hDelta_pt1->Fill(reco_jet1_pt - mc_jet1_pt);
      hDelta_N1->Fill(reco_jet1_nPart - mc_jet1_nPart);
    }

    if (matched2) {
      hReco_dR2->Fill(dR_mc2_reco2);
      hReco_jet2_pt->Fill(reco_jet2_pt);
      hReco_jet2_nPart->Fill(reco_jet2_nPart);

      // Response
      hResp_pt2->Fill(mc_jet2_pt, reco_jet2_pt);
      hProf_pt2->Fill(mc_jet2_pt, reco_jet2_pt);
      hResp_N2->Fill(mc_jet2_nPart, reco_jet2_nPart);
      hProf_N2->Fill(mc_jet2_nPart, reco_jet2_nPart);

      // Differences
      hDelta_pt2->Fill(reco_jet2_pt - mc_jet2_pt);
      hDelta_N2->Fill(reco_jet2_nPart - mc_jet2_nPart);
    }

    // Correlation (both matched)
    if (matched1 && matched2) {
      hReco_N1vsN2->Fill(reco_jet1_nPart, reco_jet2_nPart);
      hReco_N1vsN2vsPt->Fill(reco_jet1_nPart, reco_jet2_nPart, reco_jet1_pt);

      hResp_N1vsN2_MC->Fill(mc_jet1_nPart, mc_jet2_nPart);
      hResp_N1vsN2_Reco->Fill(reco_jet1_nPart, reco_jet2_nPart);
    }

    // =========================================================================
    // Fill Background Histograms
    // =========================================================================

    hBkg_multA->Fill(bkg_multA);
    hBkg_multB->Fill(bkg_multB);
    hBkg_multAvsB->Fill(bkg_multA, bkg_multB);

    if (matched1 && matched2) {
      hBkg_N1vsN2vsPt->Fill(bkg_multA, bkg_multB, reco_jet1_pt);
    }
  }

  // =========================================================================
  // Compute Covariances
  // =========================================================================
  std::cout << "Computing covariances...\n";

  fout->cd();

  TH1D *hCov_MC =
      getCovariance(hMC_N1vsN2vsPt, "MC Truth COV(N_{1}, N_{2})", "hCov_MC");
  hCov_MC->SetLineColor(kBlue);
  hCov_MC->SetMarkerColor(kBlue);
  hCov_MC->SetMarkerStyle(20);

  TH1D *hCov_Reco =
      getCovariance(hReco_N1vsN2vsPt, "Reco COV(N_{1}, N_{2})", "hCov_Reco");
  hCov_Reco->SetLineColor(kRed);
  hCov_Reco->SetMarkerColor(kRed);
  hCov_Reco->SetMarkerStyle(21);

  TH1D *hCov_Bkg = getCovariance(hBkg_N1vsN2vsPt,
                                 "Background COV(N_{A}, N_{B})", "hCov_Bkg");
  hCov_Bkg->SetLineColor(kGreen + 2);
  hCov_Bkg->SetMarkerColor(kGreen + 2);
  hCov_Bkg->SetMarkerStyle(22);

  // =========================================================================
  // Print Statistics
  // =========================================================================
  std::cout << "\n=== Analysis Summary ===\n";
  std::cout << "Total events: " << nEvents << "\n";
  std::cout << "Both jets matched: " << nBothMatched << " ("
            << (100.0 * nBothMatched / nEvents) << "%)\n";
  std::cout << "One jet matched: " << nOneMatched << " ("
            << (100.0 * nOneMatched / nEvents) << "%)\n";
  std::cout << "No jets matched: " << nNoneMatched << " ("
            << (100.0 * nNoneMatched / nEvents) << "%)\n\n";

  std::cout << "MC Truth covariance: " << getCovariance(hMC_N1vsN2) << "\n";
  std::cout << "Reco covariance: " << getCovariance(hReco_N1vsN2) << "\n";
  std::cout << "Background covariance: " << getCovariance(hBkg_multAvsB)
            << "\n\n";

  // =========================================================================
  // Create Summary Plots
  // =========================================================================
  std::cout << "Creating summary plots...\n";

  fout->cd();

  // Canvas 1: MC Truth vs Reco multiplicity correlation
  TCanvas *c1 = new TCanvas("c1", "Multiplicity Correlations", 1600, 800);
  c1->Divide(2, 1);

  c1->cd(1);
  gPad->SetRightMargin(0.15);
  hMC_N1vsN2->Draw("COLZ");
  TLine *line1 = new TLine(multMin, multMin, multMax, multMax);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(2);
  line1->SetLineStyle(2);
  line1->Draw("same");

  c1->cd(2);
  gPad->SetRightMargin(0.15);
  hReco_N1vsN2->Draw("COLZ");
  TLine *line2 = new TLine(multMin, multMin, multMax, multMax);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(2);
  line2->SetLineStyle(2);
  line2->Draw("same");

  c1->Write();

  // Canvas 2: Response matrices
  TCanvas *c2 = new TCanvas("c2", "Multiplicity Response", 1600, 800);
  c2->Divide(2, 1);

  c2->cd(1);
  gPad->SetRightMargin(0.15);
  hResp_N1->Draw("COLZ");
  TLine *line3 = new TLine(multMin, multMin, multMax, multMax);
  line3->SetLineColor(kRed);
  line3->SetLineWidth(2);
  line3->Draw("same");

  c2->cd(2);
  gPad->SetRightMargin(0.15);
  hResp_N2->Draw("COLZ");
  TLine *line4 = new TLine(multMin, multMin, multMax, multMax);
  line4->SetLineColor(kRed);
  line4->SetLineWidth(2);
  line4->Draw("same");

  c2->Write();

  // Canvas 3: Covariance comparison
  TCanvas *c3 = new TCanvas("c3", "Covariance vs pT", 800, 600);

  hCov_MC->GetYaxis()->SetRangeUser(0, hCov_MC->GetMaximum() * 1.3);
  hCov_MC->Draw("E1");
  hCov_Reco->Draw("E1 same");
  hCov_Bkg->Draw("E1 same");

  TLegend *leg = new TLegend(0.6, 0.65, 0.88, 0.88);
  leg->AddEntry(hCov_MC, "MC Truth", "lep");
  leg->AddEntry(hCov_Reco, "Reconstructed", "lep");
  leg->AddEntry(hCov_Bkg, "Background", "lep");
  leg->Draw();

  c3->Write();

  // Canvas 4: Response profiles
  TCanvas *c4 = new TCanvas("c4", "Response Profiles", 1600, 800);
  c4->Divide(2, 1);

  c4->cd(1);
  hProf_N1->Draw("E1");
  hProf_N1->SetLineColor(kBlue);
  hProf_N1->SetMarkerColor(kBlue);
  hProf_N1->SetMarkerStyle(20);
  TLine *line5 = new TLine(multMin, multMin, multMax, multMax);
  line5->SetLineColor(kRed);
  line5->SetLineStyle(2);
  line5->Draw("same");

  c4->cd(2);
  hProf_N2->Draw("E1");
  hProf_N2->SetLineColor(kBlue);
  hProf_N2->SetMarkerColor(kBlue);
  hProf_N2->SetMarkerStyle(20);
  TLine *line6 = new TLine(multMin, multMin, multMax, multMax);
  line6->SetLineColor(kRed);
  line6->SetLineStyle(2);
  line6->Draw("same");

  c4->Write();

  // =========================================================================
  // Save and Close
  // =========================================================================
  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "\nAnalysis complete!\n";
  std::cout << "Output saved to: " << outputFile << "\n";
  std::cout << "=========================\n\n";
}