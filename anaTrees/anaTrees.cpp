#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>

#include <vector>

double getEntropy(TH1D *h)
{
   double entropy = 0;
   int nBins = h->GetNbinsX();
   double total = h->Integral();
   for (int i = 1; i <= nBins; ++i) {
      double p = h->GetBinContent(i) / total;
      if (p > 0) {
         entropy -= p * std::log(p);
      }
   }
   return entropy;
}

double getEntropy(TH2D *h)
{
   double entropy = 0;
   int nBinsX = h->GetNbinsX();
   int nBinsY = h->GetNbinsY();
   double total = h->Integral();
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

void drawLabel(TPad *pad)
{
   pad->cd();
   TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.04);
   latex->SetTextFont(42);
   float x = 0.57;
   // PYTHIA 8.315
   // pp at \sqrt{s} = 200 GeV, |eta|<0.6, anti-k_{T}, R=0.4
   // min jet pt =3 GeV/c, pT_track>0.15 GeV/c, |eta_track|<1
   latex->DrawLatex(x, 0.85, "PYTHIA 8.315");
   latex->DrawLatex(x, 0.80, "#it{pp} at #sqrt{#it{s}} = 200 GeV");
   latex->DrawLatex(x, 0.75, "|#eta|<0.6, anti-k_{T}, #it{R}=0.4");
   latex->DrawLatex(x, 0.70, "jet p_{T} > 3 GeV/c");
   latex->DrawLatex(x, 0.65, "p_{T,track}>0.15 GeV/c, |#eta_{track}|<1");
   latex->DrawLatex(x, 0.60, ">1 Dijet per Event");
}

// pp200_pThat_11_15.root  pp200_pThat_2_3.root    pp200_pThat_35_45.root  pp200_pThat_55_inf.root pp200_pThat_9_11.root
// pp200_pThat_15_20.root  pp200_pThat_25_35.root  pp200_pThat_45_55.root  pp200_pThat_5_7.root
// pp200_pThat_20_25.root  pp200_pThat_3_4.root    pp200_pThat_4_5.root    pp200_pThat_7_9.root

void anaTrees()
{
   TH1::SetDefaultSumw2(true); // proper errors when scaling
   // only show n entries in stats
   gStyle->SetOptStat(11);
   TString prefix = "output/sum_pp200_ptHat_";
   vector<TString> ptHatBins = {"2_3",   "3_4",   "4_5",   "5_7",   "7_9",   "9_11", "11_15",
                                "15_20", "20_25", "25_35", "35_45", "45_55", "55_-1"};

   int lead_n_charged, sub_n_charged;
   float lead_pt, sub_pt, lead_eta, sub_eta, lead_phi, sub_phi, closeness;

   TFile *outFile = TFile::Open("anaTrees.root", "RECREATE");
   if (!outFile || outFile->IsZombie()) {
      std::cerr << "Error: could not create output file anaTrees.root" << std::endl;
      return;
   }

   TH1D *hMultLead = new TH1D("hMultLead", ";N_{ch}^{lead}", 30, 0, 30);
   TH1D *hMultSublead = new TH1D("hMultSublead", ";N_{ch}^{sublead}", 30, 0, 30);
   TH2D *hMultLeadVsSub = new TH2D(
      "hMultLeadVsSub", "Dijet multiplicity; N_{ch}^{lead};N_{ch}^{sublead}; d#sigma/dN [mb]", 30, 0, 30, 30, 0, 30);

   TH1D *hPtAll = new TH1D("hPtAll", "All Jet p_{T}; p_{T} (GeV/c); d^{2}#sigma/(d#eta dp_{T}) [mb]", 100, 0, 100);
   TH1D *hPtLead =
      new TH1D("hPtLead", "Leading Jet p_{T}; p_{T} (GeV/c); d^{2}#sigma/(d#eta dp_{T}) [mb]", 100, 0, 100);
   TH1D *hPtSub =
      new TH1D("hPtSub", "Subleading Jet p_{T}; p_{T} (GeV/c); d^{2}#sigma/(d#eta dp_{T}) [mb]", 100, 0, 100);

   TH1D *hPtDiff = new TH1D("hPtDiff",
                            "p_{T}^{lead} - p_{T}^{sublead}; p_{T}^{lead} - p_{T}^{sublead} (GeV/c); "
                            "d#sigma/d(p_{T}^{lead} - p_{T}^{sublead}) [mb]",
                            100, 0, 100);

   TH1D *stats;
   vector<TString> statNames = {"nEvents", "nAccepted", "ptHatMin", "ptHatMax", "sigmaGen_mb", "sigmaErr_mb"};

   for (const auto &bin : ptHatBins) {
      TString fileName = prefix + bin + ".root";

      TFile *f = TFile::Open(fileName);
      if (!f || f->IsZombie()) {
         cerr << "Error: could not open file " << fileName << endl;
         continue;
      }

      stats = (TH1D *)f->Get("stats");
      if (!stats) {
         cerr << "Error: stats histogram not found in file " << fileName << endl;
         f->Close();
         continue;
      }
      int nEvents = stats->GetBinContent(1);
      double xsec = stats->GetBinContent(5); // mb

      // get it from name in bin  ptmin_ptmax using TString operations
      double ptHatMin = -1;
      double ptHatMax = -1;

      ptHatMin = TString(bin(0, bin.Index("_"))).Atof();
      ptHatMax = TString(bin(bin.Index("_") + 1, bin.Length())).Atof();

      double weight = xsec / nEvents;
      cout << "nEvents = " << nEvents / 1e6 << "M, xsec = " << xsec << " mb, ptHat " << ptHatMin << "-" << ptHatMax
           << endl;
      TTree *t = (TTree *)f->Get("events");
      t->SetBranchAddress("lead_pt", &lead_pt);
      t->SetBranchAddress("sub_pt", &sub_pt);
      t->SetBranchAddress("lead_eta", &lead_eta);
      t->SetBranchAddress("sub_eta", &sub_eta);
      t->SetBranchAddress("lead_phi", &lead_phi);
      t->SetBranchAddress("sub_phi", &sub_phi);
      t->SetBranchAddress("lead_n_charged", &lead_n_charged);
      t->SetBranchAddress("sub_n_charged", &sub_n_charged);
      t->SetBranchAddress("closeness", &closeness);

      Long64_t nEntries = t->GetEntries();
      for (Long64_t i = 0; i < nEntries; ++i) {
         t->GetEntry(i);
         hMultLead->Fill(lead_n_charged, weight);
         hMultSublead->Fill(sub_n_charged, weight);
         hMultLeadVsSub->Fill(lead_n_charged, sub_n_charged, weight);
         hPtAll->Fill(lead_pt, weight);
         hPtAll->Fill(sub_pt, weight);
         hPtLead->Fill(lead_pt, weight);
         hPtSub->Fill(sub_pt, weight);
         hPtDiff->Fill(lead_pt - sub_pt, weight);
      }
   }

   outFile->cd();
   TCanvas *can = new TCanvas("can", "can", 800, 600);
   can->SetLogy();
   hPtAll->Draw();
   drawLabel(can);
   can->Write();

   can->Clear();
   hPtDiff->Draw();
   drawLabel(can);
   can->Write();

   can->Clear();
   can->SetLogy(0);
   can->SetLogz();
   hMultLeadVsSub->Draw("COLZ");
   drawLabel(can);
   can->Write();

   // add text with jet parameters as latex to Histograms

   double S1 = getEntropy(hMultLead);
   double S2 = getEntropy(hMultSublead);
   double S12 = getEntropy(hMultLeadVsSub);

   cout << "S(lead) = " << S1 << ", S(sublead) = " << S2 << ", S(lead, sublead) = " << S12 << endl;

   double covariance = S1 + S2 - S12;
   cout << "Covariance: " << covariance << endl;

   double correlation = covariance / sqrt(S1 * S2);
   cout << "Correlation: " << correlation << endl;

   outFile->Write();
   outFile->Close();
}