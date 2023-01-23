//____________________________________________________________________________________________________
void combine_histos()
{

  // This function plots 2D Pmiss vs. theta_rq histos for the 3 different recoil angle settings for the 22 GeV Jlab workshop plots
  // as well as the projections of each of these ar a given central missing momentum setting 
  // And finally, the projections are combined to get an overall count rate
  
  int font_type = 132;

  gStyle->SetOptStat("ei");
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetTitleFontSize(.08);

  // Create TCanvas
  TCanvas *c1 = new TCanvas("c1", "", 1500, 500);
  c1->Divide(3,2);
  
  // Create TCanvas
  TCanvas *c2 = new TCanvas("c2", "", 500, 500);

  // Create TCanvas
  // TCanvas *c3 = new TCanvas("c3", "", 500, 500);

  //Open  ROOT files;
  TFile *file1 = NULL;
  TFile *file2 = NULL;
  TFile *file3 = NULL;

  TString file1_path = "Q2_4p5/withSpecBoundary/d2_Eb11_Pr1_thrq27_norad_output_5degbins.root";
  TString file2_path = "Q2_4p5/withSpecBoundary/d2_Eb11_Pr1_thrq49_norad_output_5degbins.root";
  TString file3_path = "Q2_4p5/withSpecBoundary/d2_Eb11_Pr1_thrq65_norad_output_5degbins.root";
  
  file1 = new TFile(file1_path.Data());
  file2 = new TFile(file2_path.Data());
  file3 = new TFile(file3_path.Data());

  // declare 1D histos
  TH2F *H_hist1 = 0;  
  TH2F *H_hist2 = 0;
  TH2F *H_hist3 = 0;



  double axis_label_size=0.05;
  // get histogram objects
  file1->cd();
  H_hist1 = (TH2F*)file1->Get("H_Pm_vs_thrq");
  H_hist1->GetXaxis()->SetLabelSize(axis_label_size);
  H_hist1->GetYaxis()->SetLabelSize(axis_label_size);

  file2->cd();
  H_hist2 = (TH2F*)file2->Get("H_Pm_vs_thrq");
  H_hist2->GetXaxis()->SetLabelSize(axis_label_size);
  H_hist2->GetYaxis()->SetLabelSize(axis_label_size);

  file3->cd();
  H_hist3 = (TH2F*)file3->Get("H_Pm_vs_thrq");
  H_hist3->GetXaxis()->SetLabelSize(axis_label_size);
  H_hist3->GetYaxis()->SetLabelSize(axis_label_size);

  // Draw 2d histos
  c1->cd(1);
  H_hist1->Draw("colz");
  c1->cd(2);
  H_hist2->Draw("colz");
  c1->cd(3);
  H_hist3->Draw("colz");

  // Get projections for a certaint missing momentum bin range
  Double_t Pr_c = 1;
  
  TH1D * projh2X_h1 = H_hist1->ProjectionX("H1",  H_hist1->GetYaxis()->FindBin(Pr_c), H_hist1->GetYaxis()->FindBin(Pr_c));
  projh2X_h1->SetFillColorAlpha(kGreen, 0.25);
  projh2X_h1->GetXaxis()->SetLabelSize(axis_label_size);
  projh2X_h1->GetYaxis()->SetLabelSize(axis_label_size);
  c1->cd(4);
  projh2X_h1->Draw("histE0");
  c2->cd();
  projh2X_h1->Draw("histE0");
  
  TH1D * projh2X_h2 = H_hist2->ProjectionX("H2", H_hist2->GetYaxis()->FindBin(Pr_c), H_hist2->GetYaxis()->FindBin(Pr_c));
  projh2X_h2->SetFillColorAlpha(kRed+1, 0.25);
  projh2X_h2->GetXaxis()->SetLabelSize(axis_label_size);
  projh2X_h2->GetYaxis()->SetLabelSize(axis_label_size);
  c1->cd(5);
  projh2X_h2->Draw("histE0");
  c2->cd();
  projh2X_h2->Draw("histE0sames");
   
  TH1D * projh2X_h3 = H_hist3->ProjectionX("", H_hist3->GetYaxis()->FindBin(Pr_c), H_hist3->GetYaxis()->FindBin(Pr_c));
  projh2X_h3->SetFillColorAlpha(kBlue+2, 0.25);
  projh2X_h3->GetXaxis()->SetLabelSize(axis_label_size);
  projh2X_h3->GetYaxis()->SetLabelSize(axis_label_size);
  c1->cd(6);
  projh2X_h3->Draw("histE0");
  c2->cd();
  projh2X_h3->Draw("histE0sames");


  /*
  TH1D *h1 = (TH1D*) projh2X_h1->Clone();
  TH1D *h2 = (TH1D*) projh2X_h2->Clone();
  TH1D *h3 = (TH1D*) projh2X_h3->Clone();

  TH1D *h1_s = (TH1D*) projh2X_h1->Clone();
  TH1D *h2_s = (TH1D*) projh2X_h2->Clone();
  TH1D *h3_s = (TH1D*) projh2X_h3->Clone();

  
  h1->SetFillColor(14);
  h2->SetFillColor(14);
  h2->SetFillColor(14);
  h1->SetFillStyle(3025);
  h2->SetFillStyle(3025);
  h3->SetFillStyle(3025);
  
  // scale by 160/80 uA (double current)
  h1_s->Scale(2.);
  h2_s->Scale(2.);
  h3_s->Scale(2.);

  h1_s->SetFillColorAlpha(kGreen, 0.25);
  h2_s->SetFillColorAlpha(kRed+1, 0.25);
  h3_s->SetFillColorAlpha(kBlue+2, 0.25);
  
  c3->cd();
  
  h1->Draw("histE0same");
  h2->Draw("histE0same");
  h3->Draw("histE0same");
  
  h1_s->Draw("histE0sames");
  h2_s->Draw("histE0sames");
  h3_s->Draw("histE0sames");
  */


  
  /*
  projh2X_h1->Draw("histE0");
  projh2X_h2->Draw("histE0same");
  projh2X_h3->Draw("histE0same");
  */
  

  /*
  h1->Sumw2();
  h1->SetFillColorAlpha(40, 0.9);

  c3->cd();
  h1->Add(h2);
  h1->Add(h3);

  h1->Draw();
  */
  
  /*
  // Add all projections
  TH1D * projh2X_comb =  H_hist1->ProjectionX("H1",  H_hist1->GetYaxis()->FindBin(Pr_c), H_hist1->GetYaxis()->FindBin(Pr_c));
  projh2X_comb->Add(projh2X_h2);
  projh2X_comb->Add(projh2X_h3);

  // draw combined projections to canvas
  c3->cd();
  projh2X_comb->SetFillColorAlpha(kMagenta+2, 0.25);
  projh2X_comb->Draw("");
  
  */

}
