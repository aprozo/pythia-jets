// unfold_Nmultiplicity.C
// Bayesian unfolding of jet multiplicity N_reco vs N_mc (jet1 and jet2)
// using RooUnfoldBayes. Designed to work with the same trees/branches
// as in analyzeEmbedding() above.

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

#include <iostream>
#include <vector>


// Main macro
void unfold(const char *inputFile = "embedded_dijets.root",
            const char *outputFile = "unfold_N.root") {
  std::cout
      << "\n=== Bayesian unfolding of jet multiplicity N_reco vs N_mc ===\n";
  std::cout << "Input:  " << inputFile << "\n";
  std::cout << "Output: " << outputFile << "\n\n";

  // ------------------------------------------------------------
  // Open input and get reco tree
  // ------------------------------------------------------------
  TFile *fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open input file " << inputFile << "\n";
    return;
  }

  TTree *tReco = (TTree *)fin->Get("recoJets");
  if (!tReco) {
    std::cerr << "ERROR: TTree 'recoJets' not found in file\n";
    fin->Close();
    return;
  }

  Long64_t nEntries = tReco->GetEntries();
  std::cout << "Entries in recoJets: " << nEntries << "\n";

  // ------------------------------------------------------------
  // Branches (same as in analyzeEmbedding)
  // ------------------------------------------------------------
  int reco_eventID;
  bool matched1, matched2;
  double mc_jet1_pt, mc_jet1_eta, mc_jet1_phi;
  int mc_jet1_nPart;
  double reco_jet1_pt, reco_jet1_eta, reco_jet1_phi;
  int reco_jet1_nPart;
  double dR_mc1_reco1;

  double mc_jet2_pt, mc_jet2_eta, mc_jet2_phi;
  int mc_jet2_nPart;
  double reco_jet2_pt, reco_jet2_eta, reco_jet2_phi;
  int reco_jet2_nPart;
  double dR_mc2_reco2;

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

  // ------------------------------------------------------------
  // Binning for N_part (same as in analyzeEmbedding)
  // ------------------------------------------------------------
  const int nMultBins = 25;
  const double multMin = 0.0;
  const double multMax = 25.0;

  // ------------------------------------------------------------
  // Training / test histograms for Jet 1 and Jet 2
  // ------------------------------------------------------------
  // Jet 1
  TH1D *hTruthN1_train =
      new TH1D("hTruthN1_train", "Jet1 N_{part}^{MC} (train);N_{part};Events",
               nMultBins, multMin, multMax);
  TH1D *hRecoN1_train =
      new TH1D("hRecoN1_train", "Jet1 N_{part}^{Reco} (train);N_{part};Events",
               nMultBins, multMin, multMax);

  TH1D *hTruthN1_test = (TH1D *)hTruthN1_train->Clone("hTruthN1_test");
  hTruthN1_test->SetTitle("Jet1 N_{part}^{MC} (test)");
  TH1D *hRecoN1_test = (TH1D *)hRecoN1_train->Clone("hRecoN1_test");
  hRecoN1_test->SetTitle("Jet1 N_{part}^{Reco} (test)");

  // Jet 2
  TH1D *hTruthN2_train =
      new TH1D("hTruthN2_train", "Jet2 N_{part}^{MC} (train);N_{part};Events",
               nMultBins, multMin, multMax);
  TH1D *hRecoN2_train =
      new TH1D("hRecoN2_train", "Jet2 N_{part}^{Reco} (train);N_{part};Events",
               nMultBins, multMin, multMax);

  TH1D *hTruthN2_test = (TH1D *)hTruthN2_train->Clone("hTruthN2_test");
  hTruthN2_test->SetTitle("Jet2 N_{part}^{MC} (test)");
  TH1D *hRecoN2_test = (TH1D *)hRecoN2_train->Clone("hRecoN2_test");
  hRecoN2_test->SetTitle("Jet2 N_{part}^{Reco} (test)");

  // Response matrices: X = measured (Reco), Y = truth (MC)
  TH2D *hRespN1 = new TH2D(
      "hRespN1", "Jet1 N_{part} response;N_{part}^{Reco};N_{part}^{MC}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  TH2D *hRespN2 = new TH2D(
      "hRespN2", "Jet2 N_{part} response;N_{part}^{Reco};N_{part}^{MC}",
      nMultBins, multMin, multMax, nMultBins, multMin, multMax);

  // ------------------------------------------------------------
  // Fill: train/test split
  // ------------------------------------------------------------
  TRandom3 splitRng(0);
  const double testFraction = 0.20;

  for (Long64_t i = 0; i < nEntries; ++i) {
    if (i % 100000 == 0)
      std::cout << "Entry " << i << "/" << nEntries << "\r" << std::flush;

    tReco->GetEntry(i);

    const bool useForTraining = (splitRng.Uniform() > testFraction);

    // Jet 1
    if (matched1) {
      if (useForTraining) {
        hRespN1->Fill(reco_jet1_nPart, mc_jet1_nPart);
        hRecoN1_train->Fill(reco_jet1_nPart);
        hTruthN1_train->Fill(mc_jet1_nPart);
      } else {
        hRecoN1_test->Fill(reco_jet1_nPart);
        hTruthN1_test->Fill(mc_jet1_nPart);
      }
    }

    // Jet 2
    if (matched2) {
      if (useForTraining) {
        hRespN2->Fill(reco_jet2_nPart, mc_jet2_nPart);
        hRecoN2_train->Fill(reco_jet2_nPart);
        hTruthN2_train->Fill(mc_jet2_nPart);
      } else {
        hRecoN2_test->Fill(reco_jet2_nPart);
        hTruthN2_test->Fill(mc_jet2_nPart);
      }
    }
  }
  std::cout << "\nFinished filling response and train/test histograms.\n";

  // ------------------------------------------------------------
  // Build RooUnfoldResponse objects
  // ------------------------------------------------------------
  RooUnfoldResponse respN1(hRecoN1_train, hTruthN1_train, hRespN1);
  respN1.SetName("respN1");
  respN1.SetTitle("Jet1 N_{part} response (Reco vs MC)");

  RooUnfoldResponse respN2(hRecoN2_train, hTruthN2_train, hRespN2);
  respN2.SetName("respN2");
  respN2.SetTitle("Jet2 N_{part} response (Reco vs MC)");

  // ------------------------------------------------------------
  // Bayes unfolding for several iterations
  // ------------------------------------------------------------
  const int bayesIterations[] = {1, 2, 3, 4};
  const int nBayesIterations =
      sizeof(bayesIterations) / sizeof(bayesIterations[0]);

  std::vector<TH1D *> unfoldedN1(nBayesIterations, nullptr);
  std::vector<TH1D *> unfoldedN2(nBayesIterations, nullptr);

  for (int i = 0; i < nBayesIterations; ++i) {
    int nIter = bayesIterations[i];

    // Jet 1
    {
      RooUnfoldBayes unfoldN1(&respN1, hRecoN1_test, nIter);
      TH1D *hUnfN1 = (TH1D *)unfoldN1.Hunfold();
      hUnfN1->SetDirectory(0);
      hUnfN1->SetName(Form("hUnfold_N1_iter%d", nIter));
      hUnfN1->SetTitle(Form("Jet1 unfolded N_{part} (Bayes, %d it.)", nIter));
      unfoldedN1[i] = hUnfN1;
    }

    // Jet 2
    {
      RooUnfoldBayes unfoldN2(&respN2, hRecoN2_test, nIter);
      TH1D *hUnfN2 = (TH1D *)unfoldN2.Hunfold();
      hUnfN2->SetDirectory(0);
      hUnfN2->SetName(Form("hUnfold_N2_iter%d", nIter));
      hUnfN2->SetTitle(Form("Jet2 unfolded N_{part} (Bayes, %d it.)", nIter));
      unfoldedN2[i] = hUnfN2;
    }
  }

  // ------------------------------------------------------------
  // Simple closure canvas for Jet 1 (optional visual check)
  // ------------------------------------------------------------
  TCanvas *cN1 = new TCanvas("cN1", "Jet1 N_{part} unfolding", 800, 800);
  cN1->Divide(1, 2);

  cN1->cd(1);
  gPad->SetLogy();

  // normalize by bin width (discrete bins, but keep same logic as pT code)
  TH1D *hTruthN1_testW = (TH1D *)hTruthN1_test->Clone("hTruthN1_testW");
  TH1D *hRecoN1_testW = (TH1D *)hRecoN1_test->Clone("hRecoN1_testW");
  hTruthN1_testW->Scale(1.0, "width");
  hRecoN1_testW->Scale(1.0, "width");

  hTruthN1_testW->SetMarkerStyle(20);
  hTruthN1_testW->SetLineColor(kBlack);
  hRecoN1_testW->SetMarkerStyle(24);
  hRecoN1_testW->SetLineColor(kBlue + 1);

  hTruthN1_testW->Draw("E1");
  hRecoN1_testW->Draw("E1 SAME");

  TLegend *leg1 = new TLegend(0.55, 0.60, 0.88, 0.88);
  leg1->AddEntry(hTruthN1_testW, "Truth N_{1} (test)", "lp");
  leg1->AddEntry(hRecoN1_testW, "Measured N_{1} (test)", "lp");

  for (int i = 0; i < nBayesIterations; ++i) {
    TH1D *hUnfW = (TH1D *)unfoldedN1[i]->Clone(
        Form("%s_width", unfoldedN1[i]->GetName()));
    hUnfW->Scale(1.0, "width");
    hUnfW->SetMarkerStyle(20);
    hUnfW->SetMarkerColor(kMagenta + i);
    hUnfW->SetLineColor(kMagenta + i);
    hUnfW->Draw("E1 SAME");
    leg1->AddEntry(hUnfW, Form("Bayes %d it.", bayesIterations[i]), "lp");
  }
  leg1->Draw();

  // Ratio unfolded / truth
  cN1->cd(2);
  TH1D *frame = (TH1D *)hTruthN1_testW->Clone("ratioFrameN1");
  frame->Reset();
  frame->GetYaxis()->SetTitle("Unfolded / Truth");
  frame->GetYaxis()->SetRangeUser(0.5, 1.5);
  frame->Draw("");

  for (int i = 0; i < nBayesIterations; ++i) {
    TH1D *hUnfW = (TH1D *)unfoldedN1[i]->Clone(
        Form("%s_ratio", unfoldedN1[i]->GetName()));
    hUnfW->Scale(1.0, "width");
    hUnfW->Divide(hTruthN1_testW);
    hUnfW->SetMarkerStyle(20);
    hUnfW->SetMarkerColor(kMagenta + i);
    hUnfW->SetLineColor(kMagenta + i);
    hUnfW->Draw(i == 0 ? "E1" : "E1 SAME");
  }

  TLine *lUp = new TLine(frame->GetXaxis()->GetXmin(), 1.1,
                         frame->GetXaxis()->GetXmax(), 1.1);
  lUp->SetLineStyle(2);
  lUp->SetLineColor(kGray + 1);
  lUp->Draw();

  TLine *lDown = new TLine(frame->GetXaxis()->GetXmin(), 0.9,
                           frame->GetXaxis()->GetXmax(), 0.9);
  lDown->SetLineStyle(2);
  lDown->SetLineColor(kGray + 1);
  lDown->Draw();

  // ------------------------------------------------------------
  // Write everything out
  // ------------------------------------------------------------
  TFile *fout = TFile::Open(outputFile, "RECREATE");

  // Response and train/test histograms
  fout->mkdir("responseN");
  fout->cd("responseN");
  hRespN1->Write();
  hRespN2->Write();
  hTruthN1_train->Write();
  hRecoN1_train->Write();
  hTruthN1_test->Write();
  hRecoN1_test->Write();
  hTruthN2_train->Write();
  hRecoN2_train->Write();
  hTruthN2_test->Write();
  hRecoN2_test->Write();

  // Unfolded spectra and covariance matrices
  fout->mkdir("unfoldedN");
  fout->cd("unfoldedN");
  for (int i = 0; i < nBayesIterations; ++i) {
    if (unfoldedN1[i])
      unfoldedN1[i]->Write();
    if (unfoldedN2[i])
      unfoldedN2[i]->Write();
  }

  fout->mkdir("plots");
  fout->cd("plots");
  cN1->Write("cN1_unfold_multiplicity");

  fout->Close();
  fin->Close();

  std::cout << "Unfolding done. Results written to: " << outputFile << "\n";
}
