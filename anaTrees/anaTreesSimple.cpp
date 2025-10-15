#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <TColor.h>
#include <iostream>

#include <vector>

static const double balanceCut = 0.2;

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

double getCovariance(TH2D *h)
{
   double S1 = getEntropy((TH1D *)h->ProjectionX());
   double S2 = getEntropy((TH1D *)h->ProjectionY());
   double S12 = getEntropy(h);
   return S1 + S2 - S12;
}

TH1D *getCovariance(TH3D *h, TString title = "")
{
   TString name = TString(h->GetName()) + "_cov";
   TString z_title = h->GetZaxis()->GetTitle();
   // strip off everything after ; in z_title
   if (z_title.Index(";") >= 0)
      z_title = z_title(0, z_title.Index(";"));

   TH1D *cov = new TH1D(name, title + ";" + z_title + ";" + title, h->GetNbinsZ(), h->GetZaxis()->GetXmin(),
                        h->GetZaxis()->GetXmax());
   for (int i = 1; i <= h->GetNbinsZ(); ++i) {
      h->GetZaxis()->SetRange(i, i);
      TH2D *h2 = (TH2D *)h->Project3D("xy");
      double c = getCovariance(h2);
      cov->SetBinContent(i, c);
      if (h2->GetEntries() != 0) {
         double err = 1 / sqrt(h2->GetEntries());
         cov->SetBinError(i, err);
      }
   }
   return cov;
}

void anaTreesSimple()
{
   TH3::SetDefaultSumw2(true);
   // only show n entries in stats
   gStyle->SetOptStat(0);
   TString prefix = "output/sum_pp200_ptHat_";
   vector<TString> ptHatBins = {"2_3",   "3_4",   "4_5",   "5_7",   "7_9",   "9_11", "11_15",
                                "15_20", "20_25", "25_35", "35_45", "45_55", "55_-1"};

   int lead_n_charged, sub_n_charged;
   double lead_pt, sub_pt, lead_eta, sub_eta, lead_phi, sub_phi, closeness, background_mult_A, background_mult_B;

   TFile *outFile = TFile::Open("anaTrees.root", "RECREATE");
   // add text with jet parameters as latex to Histograms

   const int nPtBins = 20;
   const double ptMin = 0;
   const double ptMax = 100;

   const int nMultBins = 30;
   const double multMin = 0;
   const double multMax = 30;

   TH3D *hMult3D = new TH3D("hMult3D",
                            "Dijet multiplicity; N_{ch}^{lead};N_{ch}^{sublead};p_{t}^{lead} (GeV/c); d#sigma/dN "
                            "[mb]",
                            nMultBins, multMin, multMax, nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);

   TH3D *hBackgroundMultAVsMultBVsPt =
      new TH3D("hBackgroundMultAVsMultBVsPt",
               "Background multiplicity; N_{ch}^{A};N_{ch}^{B};p_{t}^{lead} (GeV/c); d#sigma/dN [mb]", nMultBins,
               multMin, multMax, nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);

   TH1D *stats;
   vector<TString> statNames = {"nEvents", "nAccepted", "ptHatMin", "ptHatMax", "sigmaGen_mb", "sigmaErr_mb"};

   for (const auto &bin : ptHatBins) {
      TString fileName = prefix + bin + ".root";

      TFile *f = TFile::Open(fileName);

      stats = (TH1D *)f->Get("stats");
      int nEvents = stats->GetBinContent(1);
      double xsec = stats->GetBinContent(5); // mb

      // get it from name in bin  ptmin_ptmax using TString operations
      double ptHatMin = -1;
      double ptHatMax = -1;

      ptHatMin = TString(bin(0, bin.Index("_"))).Atof();
      ptHatMax = TString(bin(bin.Index("_") + 1, bin.Length())).Atof();

      double weight = xsec / nEvents;

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
      t->SetBranchAddress("background_mult_A", &background_mult_A);
      t->SetBranchAddress("background_mult_B", &background_mult_B);

      Long64_t nEntries = t->GetEntries();

      for (Long64_t i = 0; i < nEntries; ++i) {
         t->GetEntry(i);
         double balance = sub_pt / lead_pt;

         if (balance < balanceCut)
            continue; // remove unbalanced dijets

         hMult3D->Fill(lead_n_charged, sub_n_charged, lead_pt, weight);
         hBackgroundMultAVsMultBVsPt->Fill(background_mult_A, background_mult_B, lead_pt, weight);
      }
   }
   TH1D *covVsPt = getCovariance(hMult3D, "COV(N_{ch}^{lead},N_{ch}^{sublead})");
   TH1D *backgroundCovVsPt = getCovariance(hBackgroundMultAVsMultBVsPt, "COV(UE_{A},UE_{B})");

   outFile->Write();
   outFile->Close();
}