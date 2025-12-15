// Define custom color palette once
void initCustomColors() {
  static bool initialized = false;
  if (initialized)
    return;
  initialized = true;

  new TColor(2000, (255. / 255.), (89. / 255.), (74. / 255.));   // red-ish
  new TColor(2001, (25. / 255.), (170. / 255.), (25. / 255.));   // green-ish
  new TColor(2002, (66. / 255.), (98. / 255.), (255. / 255.));   // blue-ish
  new TColor(2003, (153. / 255.), (0. / 255.), (153. / 255.));   // magenta-ish
  new TColor(2004, (255. / 255.), (166. / 255.), (33. / 255.));  // yellow-ish
  new TColor(2005, (0. / 255.), (170. / 255.), (255. / 255.));   // azur-ish
  new TColor(2006, (204. / 255.), (153. / 255.), (255. / 255.)); // violet-ish
  new TColor(2007, (107. / 255.), (142. / 255.), (35. / 255.));  // olive
  new TColor(2008, (100. / 255.), (149. / 255.),
             (237. / 255.)); // corn flower blue
  new TColor(2009, (255. / 255.), (69. / 255.), (0. / 255.));  // orange red
  new TColor(2010, (0. / 255.), (128. / 255.), (128. / 255.)); // teal
  new TColor(2011, (176. / 255.), (196. / 255.),
             (222. / 255.)); // light steel blue
  new TColor(2012, (255. / 255.), (215. / 255.), (0. / 255.)); // gold-ish

  new TColor(3000, (251. / 255.), (228. / 255.), (216. / 255.));
  new TColor(3001, (223. / 255.), (182. / 255.), (178. / 255.));
  new TColor(3002, (133. / 255.), (79. / 255.), (108. / 255.));
  new TColor(3003, (82. / 255.), (43. / 255.), (91. / 255.));
  new TColor(3004, (43. / 255.), (18. / 255.), (76. / 255.));
  new TColor(3005, (25. / 255.), (0. / 255.), (25. / 255.));
  new TColor(3006, (49. / 255.), (61. / 255.), (90. / 255.));
}

// Label block (kept from old code; adjust text if you change selections)
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

void plot() {
  TFile *f = TFile::Open("analysis_unfoldN_vsPt.root");
  TH1D *hCov_MC = (TH1D *)f->Get("hCov_MC");
  hCov_MC->Rebin(4);
  TH1D *hCov_Reco = (TH1D *)f->Get("hCov_Reco");
  hCov_Reco->Rebin(4);
  TH1D *hCov_Unfold = (TH1D *)f->Get("hCov_Unfold");
  hCov_Unfold->Rebin(4);

  hCov_MC->SetLineColor(2002);
  hCov_MC->SetMarkerColor(2002);

  hCov_Reco->SetLineColor(2000);
  hCov_Reco->SetMarkerColor(2000);

  hCov_Unfold->SetLineColor(2004);
  hCov_Unfold->SetMarkerColor(2004);

  TH1D *hRatio = (TH1D *)hCov_Unfold->Clone("hRatio");
  hRatio->SetTitle("Covariance Ratio: Unfolded / MC Truth;p_{T} "
                   "[GeV/c];Unfolded / MC Truth");
  hRatio->Divide(hCov_MC);
  hRatio->GetYaxis()->SetRangeUser(0.0, 2.0);

  // Canvas 3: Covariance comparison
  TCanvas *can = new TCanvas("c3", "Covariance vs pT", 800, 600);
  can->cd();
  hCov_MC->GetYaxis()->SetRangeUser(-0.2, hCov_MC->GetMaximum() * 2);
  hCov_MC->SetLineWidth(2);
  hCov_Reco->SetLineWidth(2);
  hCov_Unfold->SetLineWidth(2);

  hCov_MC->Draw("E1");
  hCov_Reco->Draw("E1 same");
  hCov_Unfold->Draw("E1 same");

  {
    TLegend *leg = new TLegend(0.6, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hCov_MC, "MC Truth", "lep");
    leg->AddEntry(hCov_Reco, "Reconstructed", "lep");
    leg->AddEntry(hCov_Unfold, "Unfolded", "lep");

    leg->Draw();
  }
  drawLabel((TPad *)gPad, 0.15);
  can->SaveAs("unfolded.pdf");

  can->Clear();
  hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hRatio->GetYaxis()->SetTitleOffset(0.9);
  hRatio->SetLineWidth(2);
  hRatio->Draw("E1");
  TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                          hRatio->GetXaxis()->GetXmax(), 1.0);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("same");

  drawLabel((TPad *)gPad, 0.15);
  can->SaveAs("covariance_ratio.pdf");
}