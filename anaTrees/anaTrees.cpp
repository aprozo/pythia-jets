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

void drawLabel(TPad *pad, float x = 0.57, TString extra = "")
{
   pad->cd();
   TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.04);
   latex->SetTextFont(42);
   // PYTHIA 8.315
   // pp at \sqrt{s} = 200 GeV, |eta|<0.6, anti-k_{T}, R=0.4
   // min jet pt =3 GeV/c, pT_track>0.15 GeV/c, |eta_track|<1
   latex->DrawLatex(x, 0.85, "PYTHIA 8.315");
   latex->DrawLatex(x, 0.80, "#it{pp} at #sqrt{#it{s}} = 200 GeV");
   latex->DrawLatex(x, 0.75, "|#eta|<0.6, anti-k_{T}, #it{R}=0.4");
   latex->DrawLatex(x, 0.70, "jet p_{t} > 3 GeV/c");
   latex->DrawLatex(x, 0.65, "p_{t,track}>0.15 GeV/c, |#eta_{track}|<1");
   latex->DrawLatex(x, 0.60, "1 dijet per event");
   latex->DrawLatex(x, 0.55, Form("p_{t}^{sublead}/p_{t}^{lead} > %.1f", balanceCut));
   latex->DrawLatex(x, 0.50, "|#phi_{lead} - #phi_{sublead}| > 3#pi/4");
   latex->DrawLatex(x, 0.45, extra);
}

void drawProjectionsPt(TCanvas *can, TH2D *h)
{
   vector<int> colors = {2000, 2001, 2002, 2004, 2005, 2006, 2007};
   TLegend *leg2 = new TLegend(0.75, 0.3, 0.95, 0.7);
   leg2->SetBorderSize(0);
   can->Clear();
   can->SetName(TString("draw_projections_") + h->GetName());
   can->SetLogy(1);
   int counter = 0;
   int step = h->GetNbinsX() / 5;
   for (int i = 1; i <= h->GetNbinsX() - step; i += step) {

      int binMax = i + step - 1;
      if (binMax > h->GetNbinsX())
         binMax = h->GetNbinsX();
      TH1D *proj = h->ProjectionY(Form("proj_%d_%d", i, binMax), i, binMax);
      proj->GetYaxis()->SetTitle("Normalized entries");
      proj->SetTitle("");
      proj->SetLineColor(colors[counter]);
      proj->SetMarkerColor(colors[counter]);
      proj->SetMarkerStyle(20);
      proj->SetMarkerSize(0.6);
      proj->GetYaxis()->SetRangeUser(1e-4, 10);
      proj->Rebin(4);
      proj->DrawNormalized(i == 1 ? "E1" : "same E1");

      leg2->AddEntry(proj,
                     TString(h->GetXaxis()->GetTitle()) +
                        Form("%2.f - %2.f GeV/c", h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(binMax)),
                     "lp");
      counter++;
   }
   leg2->Draw();
}

// pp200_pThat_11_15.root  pp200_pThat_2_3.root    pp200_pThat_35_45.root  pp200_pThat_55_inf.root pp200_pThat_9_11.root
// pp200_pThat_15_20.root  pp200_pThat_25_35.root  pp200_pThat_45_55.root  pp200_pThat_5_7.root
// pp200_pThat_20_25.root  pp200_pThat_3_4.root    pp200_pThat_4_5.root    pp200_pThat_7_9.root

void anaTrees()
{
   TH1::SetDefaultSumw2(true); // proper errors when scaling
   TH3::SetDefaultSumw2(true);
   TH2::SetDefaultSumw2(true);
   // only show n entries in stats
   gStyle->SetOptStat(0);
   TString prefix = "output/sum_pp200_ptHat_";
   vector<TString> ptHatBins = {"2_3",   "3_4",   "4_5",   "5_7",   "7_9",   "9_11", "11_15",
                                "15_20", "20_25", "25_35", "35_45", "45_55", "55_-1"};
   // vector<TString> ptHatBins = {"2_3", "55_-1"};

   int lead_n_charged, sub_n_charged;
   double lead_pt, sub_pt, lead_eta, sub_eta, lead_phi, sub_phi, closeness, background_mult_A, background_mult_B;

   TFile *outFile = TFile::Open("anaTrees.root", "RECREATE");
   if (!outFile || outFile->IsZombie()) {
      std::cerr << "Error: could not create output file anaTrees.root" << std::endl;
      return;
   }
   // add text with jet parameters as latex to Histograms

   TH1D *hStatistics = new TH1D("hStatistics", "nEvents; ptHat range;nEvents", ptHatBins.size(), 0, ptHatBins.size());
   for (size_t i = 0; i < ptHatBins.size(); ++i) {
      hStatistics->GetXaxis()->SetBinLabel(i + 1, ptHatBins[i]);
   }

   TH1D *hAcceptedEvents = (TH1D *)hStatistics->Clone("hAcceptedEvents");
   hAcceptedEvents->SetTitle("nAcceptedEvents; ptHat range;nAcceptedEvents");

   const int nPtBins = 20;
   const double ptMin = 0;
   const double ptMax = 100;

   const int nMultBins = 30;
   const double multMin = 0;
   const double multMax = 30;

   TH1D *hMultLead = new TH1D("hMultLead", ";N_{ch}^{lead}", nMultBins, multMin, multMax);
   TH1D *hMultSublead = new TH1D("hMultSublead", ";N_{ch}^{sublead}", nMultBins, multMin, multMax);
   TH2D *hMultLeadVsSub =
      new TH2D("hMultLeadVsSub", "Dijet multiplicity; N_{ch}^{lead};N_{ch}^{sublead}; d#sigma/dN [mb]", nMultBins,
               multMin, multMax, nMultBins, multMin, multMax);

   TH3D *hMult3D = new TH3D("hMult3D",
                            "Dijet multiplicity; N_{ch}^{lead};N_{ch}^{sublead};p_{t}^{lead} (GeV/c); d#sigma/dN "
                            "[mb]",
                            nMultBins, multMin, multMax, nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);
   TH3D *hMult3DLeSub =
      new TH3D("hMult3DLeSub",
               "Dijet multiplicity; N_{ch}^{lead};N_{ch}^{sublead};p_{t}^{lead} - p_{t}^{sublead} (GeV/c);  d#sigma/dN "
               "[mb]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);

   TH2D *hBackgroundAverageMult =
      new TH2D("hBackgroundAverageMult",
               "Background avg multiplicity; (N_{ch}^{A}+N_{ch}^{B})/2;p_{t}^{lead} (GeV/c); d#sigma/dN [mb]",
               nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);

   TH3D *hBackgroundMultAVsMultBVsPt =
      new TH3D("hBackgroundMultAVsMultBVsPt",
               "Background multiplicity; N_{ch}^{A};N_{ch}^{B};p_{t}^{lead} (GeV/c); d#sigma/dN [mb]", nMultBins,
               multMin, multMax, nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);

   TH3D *hBackgroundMultAVsMultBVsLeSub =
      new TH3D("hBackgroundMultAVsMultBVsLeSub",
               "Background multiplicity; N_{ch}^{A};N_{ch}^{B};p_{t}^{lead} - p_{t}^{sublead} (GeV/c); d#sigma/dN [mb]",
               nMultBins, multMin, multMax, nMultBins, multMin, multMax, nPtBins, ptMin, ptMax);

   // QA histograms
   TH1D *hPtLeSub = new TH1D("hPtLeSub",
                             "p_{t}^{lead} - p_{t}^{sublead}; p_{t}^{lead} - p_{t}^{sublead} (GeV/c); "
                             "d#sigma/dp_{t} [mb]",
                             nPtBins, ptMin, ptMax);
   TH1D *hPtAll =
      new TH1D("hPtAll", "All Jet p_{t}; p_{t} (GeV/c); d^{2}#sigma/(d#eta dp_{t}) [mb]", nPtBins, ptMin, ptMax);
   TH1D *hPtLead =
      new TH1D("hPtLead", "Leading Jet p_{t}; p_{t} (GeV/c); d^{2}#sigma/(d#eta dp_{t}) [mb]", nPtBins, ptMin, ptMax);
   TH1D *hPtSub =
      new TH1D("hPtSub", "Subleading Jet p_{t}; p_{t} (GeV/c); d^{2}#sigma/(d#eta dp_{t}) [mb]", nPtBins, ptMin, ptMax);
   TH1D *hBalance = new TH1D("hBalance", "p_{t}^{lead}/p_{t}^{sublead}; p_{t}^{lead}/p_{t}^{sublead}", 50, 0, 1);
   TH2D *hBalanceVsPt =
      new TH2D("hBalanceVsPt",
               "p_{t}^{sublead}/p_{t}^{lead} vs p_{t}^{lead}; p_{t}^{lead} (GeV/c); p_{t}^{sublead}/p_{t}^{lead}", 100,
               ptMin, ptMax, 50, 0, 1);
   TH2D *hBalanceVsLeSub = new TH2D(
      "hBalanceVsLeSub",
      "p_{t}^{sublead}/p_{t}^{lead} vs p_{t}^{lead} - p_{t}^{sublead}; p_{t}^{lead} - p_{t}^{sublead} (GeV/c); "
      "p_{t}^{sublead}/p_{t}^{lead}",
      100, ptMin, ptMax, 50, 0, 1);

   TH1D *hCloseness = new TH1D("hCloseness", "Closeness; |#phi_{lead} - #phi_{sublead} - #pi/2|", 50, 0, 1);
   TH2D *hClosenessVsPt = new TH2D(
      "hClosenessVsPt", "Closeness vs p_{t}^{lead}; p_{t}^{lead} (GeV/c); |#phi_{lead} - #phi_{sublead}- #pi/2|", 100,
      ptMin, ptMax, 100, 0, 1);
   TH2D *hClosenessVsLeSub = new TH2D("hClosenessVsLeSub",
                                      "Closeness vs p_{t}^{lead} - p_{t}^{sublead}; p_{t}^{lead} - p_{t}^{sublead} "
                                      "(GeV/c); |#phi_{lead} - #phi_{sublead}- #pi/2|",
                                      100, ptMin, ptMax, 100, 0, 1);

   TH1D *stats;
   vector<TString> statNames = {"nEvents", "nAccepted", "ptHatMin", "ptHatMax", "sigmaGen_mb", "sigmaErr_mb"};

   double totalXsec = 0;
   double totalNevents = 0;
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
      totalNevents += stats->GetBinContent(1);
      totalXsec += stats->GetBinContent(5); // mb
      hStatistics->SetBinContent(hStatistics->GetXaxis()->FindBin(bin), stats->GetBinContent(1));
      hAcceptedEvents->SetBinContent(hAcceptedEvents->GetXaxis()->FindBin(bin), stats->GetBinContent(2));
      f->Close();
   }
   cout << "Total cross section from all ptHat bins: " << totalXsec << " mb" << endl;

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
      cout << "nEvents = " << nEvents / 1e6 << "M, accepted = " << stats->GetBinContent(2) / 1e6 << "M, xsec = " << xsec
           << " mb, ptHat " << ptHatMin << "-" << ptHatMax << endl;
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
         hBalance->Fill(balance, weight);
         hBalanceVsPt->Fill(lead_pt, balance, weight);
         hBalanceVsLeSub->Fill(lead_pt - sub_pt, balance, weight);

         hCloseness->Fill(closeness, weight);
         hClosenessVsPt->Fill(lead_pt, closeness, weight);
         hClosenessVsLeSub->Fill(lead_pt - sub_pt, closeness, weight);

         if (balance < balanceCut)
            continue; // remove unbalanced dijets

         hMult3D->Fill(lead_n_charged, sub_n_charged, lead_pt, weight);
         hMult3DLeSub->Fill(lead_n_charged, sub_n_charged, lead_pt - sub_pt, weight);
         if (lead_pt > 70) {
            hMultLead->Fill(lead_n_charged, weight);
            hMultSublead->Fill(sub_n_charged, weight);
            hMultLeadVsSub->Fill(lead_n_charged, sub_n_charged, weight);
         }
         hPtAll->Fill(lead_pt, weight);
         hPtAll->Fill(sub_pt, weight);
         hPtLead->Fill(lead_pt, weight);
         hPtSub->Fill(sub_pt, weight);
         hPtLeSub->Fill(lead_pt - sub_pt, weight);
         hBackgroundMultAVsMultBVsPt->Fill(background_mult_A, background_mult_B, lead_pt, weight);
         hBackgroundMultAVsMultBVsLeSub->Fill(background_mult_A, background_mult_B, lead_pt - sub_pt, weight);

         double avgBackgroundMult = (background_mult_A + background_mult_B) / 2.0;

         hBackgroundAverageMult->Fill(avgBackgroundMult, lead_pt, weight);
      }
   }

   TH1D *covVsPt = getCovariance(hMult3D, "COV(N_{ch}^{lead},N_{ch}^{sublead})");
   TH1D *covVsLeSub = getCovariance(hMult3DLeSub, "COV(N_{ch}^{lead},N_{ch}^{sublead})");
   TH1D *backgroundCovVsPt = getCovariance(hBackgroundMultAVsMultBVsPt, "COV(UE_{A},UE_{B})");
   TH1D *backgroundCovVsLeSub = getCovariance(hBackgroundMultAVsMultBVsLeSub, "COV(UE_{A},UE_{B})");

   TCanvas *can = new TCanvas("can", "can", 800, 600);
   TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.04);
   latex->SetTextFont(42);

   new TColor(2000, (255. / 255.), (89. / 255.), (74. / 255.));   // red-ish
   new TColor(2001, (25. / 255.), (170. / 255.), (25. / 255.));   // green-ish
   new TColor(2002, (66. / 255.), (98. / 255.), (255. / 255.));   // blue-ish
   new TColor(2001, (153. / 255.), (0. / 255.), (153. / 255.));   // magenta-ish
   new TColor(2004, (255. / 255.), (166. / 255.), (33. / 255.));  // yellow-ish
   new TColor(2005, (0. / 255.), (170. / 255.), (255. / 255.));   // azur-ish
   new TColor(2006, (204. / 255.), (153. / 255.), (255. / 255.)); // violet-ish
   new TColor(2007, (107. / 255.), (142. / 255.), (35. / 255.));  // olive
   new TColor(2008, (100. / 255.), (149. / 255.),
              (237. / 255.));                                   // corn flower blue
   new TColor(2009, (255. / 255.), (69. / 255.), (0. / 255.));  // orange red
   new TColor(2010, (0. / 255.), (128. / 255.), (128. / 255.)); // teal
   new TColor(2011, (176. / 255.), (196. / 255.),
              (222. / 255.));                                   // light steel blue
   new TColor(2012, (255. / 255.), (215. / 255.), (0. / 255.)); // gold-ish

   new TColor(3000, (251. / 255.), (228. / 255.), (216. / 255.));
   new TColor(3001, (223. / 255.), (182. / 255.), (178. / 255.));
   new TColor(3002, (133. / 255.), (79. / 255.), (108. / 255.));
   new TColor(3003, (82. / 255.), (43. / 255.), (91. / 255.));
   new TColor(3004, (43. / 255.), (18. / 255.), (76. / 255.));
   new TColor(3005, (25. / 255.), (0. / 255.), (25. / 255.));
   new TColor(3006, (49. / 255.), (61. / 255.), (90. / 255.));

   can->SaveAs("anaTrees.pdf[");

   outFile->cd();
   can->SetLogy();
   can->SetName("draw_ptAll");
   hPtAll->SetLineColor(2002);
   hPtAll->SetMarkerColor(2002);
   hPtAll->Draw("E1");
   drawLabel(can);
   can->Write();
   can->SaveAs("anaTrees.pdf");

   can->SetLogy();
   can->SetName("draw_ptLeSub");
   hPtLeSub->SetLineColor(2002);
   hPtLeSub->SetMarkerColor(2002);
   hPtLeSub->Draw("E1");
   drawLabel(can);
   can->Write();
   can->SaveAs("anaTrees.pdf");

   can->SetLogy(0);

   can->SetLogz(1);

   can->SetName("draw_MultLeadVsSub");
   hMultLeadVsSub->DrawNormalized("colz");
   drawLabel(can, 0.53, "p_{t}^{lead} > 70 GeV/c");
   can->Write();
   can->SaveAs("anaTrees.pdf");

   can->SetName("draw_MultLead");
   hMultLead->SetLineColor(2002);
   hMultLead->SetMarkerColor(2002);
   hMultLead->DrawNormalized("E1");
   drawLabel(can, 0.53, "p_{t}^{lead} > 70 GeV/c");
   can->Write();
   can->SaveAs("anaTrees.pdf");

   can->SetName("draw_MultSublead");
   hMultSublead->SetLineColor(2002);
   hMultSublead->SetMarkerColor(2002);
   hMultSublead->DrawNormalized("E1");
   drawLabel(can, 0.53, "p_{t}^{lead} > 70 GeV/c");
   can->Write();
   can->SaveAs("anaTrees.pdf");

   /// COVARIANCE
   covVsPt->SetLineWidth(2);
   covVsPt->SetLineColor(2002);
   covVsPt->SetMarkerColor(2002);
   covVsPt->GetYaxis()->SetRangeUser(-0.1, 1.6);
   covVsPt->Draw("E1");
   backgroundCovVsPt->SetLineWidth(2);
   backgroundCovVsPt->SetLineColor(2000);
   backgroundCovVsPt->SetMarkerColor(2000);

   backgroundCovVsPt->Draw("same E1");

   TH1D *subtractedCovVsPt = (TH1D *)covVsPt->Clone("subtractedCovVsPt");
   subtractedCovVsPt->Add(backgroundCovVsPt, -1);
   subtractedCovVsPt->GetZaxis()->SetTitle("COV(N_{ch}^{lead},N_{ch}^{sublead}) - COV(UE_{A},UE_{B})");
   subtractedCovVsPt->SetLineWidth(2);
   subtractedCovVsPt->SetLineColor(kBlack);
   subtractedCovVsPt->SetMarkerColor(kBlack);
   subtractedCovVsPt->Draw("same E1");
   drawLabel(can, 0.15, "|#phi_{lead} - #phi_{A/B}| = #pi/2");

   TLegend *leg = new TLegend(0.15, 0.2, 0.5, 0.4);
   leg->SetBorderSize(0);
   leg->AddEntry(covVsPt, "COV(N_{ch}^{lead},N_{ch}^{sublead})", "lp");
   leg->AddEntry(backgroundCovVsPt, " COV(UE_{A},UE_{B})", "lp");
   leg->AddEntry(subtractedCovVsPt, "COV(N_{ch}^{lead},N_{ch}^{sublead}) -  COV(UE_{A},UE_{B})", "lp");
   leg->Draw();

   can->SetName("draw_subtractedCovVsPt");
   can->Write();
   can->SaveAs("anaTrees.pdf");

   /// leSub
   can->Clear();
   covVsLeSub->SetLineWidth(2);
   covVsLeSub->SetLineColor(2002);
   covVsLeSub->SetMarkerColor(2002);
   covVsLeSub->GetYaxis()->SetRangeUser(-0.1, 1.6);
   covVsLeSub->Draw("E1");

   backgroundCovVsLeSub->SetLineColor(2000);
   backgroundCovVsLeSub->SetMarkerColor(2000);
   backgroundCovVsLeSub->Draw("same E1");

   TH1D *subtractedCovVsLeSub = (TH1D *)covVsLeSub->Clone("subtractedCovVsLeSub");
   subtractedCovVsLeSub->Add(backgroundCovVsLeSub, -1);
   subtractedCovVsLeSub->SetTitle("COV(N_{ch}^{lead},N_{ch}^{sublead}) - COV(UE_{A},UE_{B})");

   subtractedCovVsLeSub->SetLineColor(kBlack);
   subtractedCovVsLeSub->SetMarkerColor(kBlack);
   subtractedCovVsLeSub->Draw("same E1");
   drawLabel(can, 0.53, "|#phi_{lead} - #phi_{A/B}| = #pi/2");
   // set position of legend
   leg->SetX1NDC(0.53);
   leg->SetX2NDC(0.9);

   leg->Draw();

   can->SetName("draw_subtractedCovVsLeSub");
   can->Write();
   can->SaveAs("anaTrees.pdf");

   // make projections
   drawProjectionsPt(can, hBalanceVsPt);
   can->Write();
   can->SaveAs("anaTrees.pdf");

   drawProjectionsPt(can, hBalanceVsLeSub);
   can->Write();
   can->SaveAs("anaTrees.pdf");

   drawProjectionsPt(can, hClosenessVsPt);
   can->Write();
   can->SaveAs("anaTrees.pdf");

   drawProjectionsPt(can, hClosenessVsLeSub);
   can->Write();
   can->SaveAs("anaTrees.pdf");

   can->SaveAs("anaTrees.pdf]");

   outFile->Write();
   outFile->Close();
}