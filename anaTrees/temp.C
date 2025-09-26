void temp()
{
   // ε(pT​)=eff_max[1−exp(−(pT​/p0​)^n)],
   double eff_max = 0.88;
   double p0 = 0.25;
   double n = 1.2;
   TF1 *eff = new TF1("eff", "[0]*(1-exp(-pow(x/[1],[2])))", 0, 30);
   eff->SetParameters(eff_max, p0, n);

   TCanvas *can = new TCanvas("can", "can", 800, 600);
   can->cd();
   eff->SetTitle("Track reconstruction efficiency; p_{t} (GeV/c); efficiency");
   eff->GetXaxis()->SetRangeUser(0, 3);
   eff->Draw();

   TLatex *latex = new TLatex();
   latex->SetNDC();
   latex->SetTextSize(0.08);
   latex->SetTextFont(42);
   latex->DrawLatex(0.25, 0.25, "eff = 0.88 *(1-e^{-(#frac{p_{t}}{0.25 GeV/c})^{1.2}})");

   can->SaveAs("track_efficiency.pdf");
}