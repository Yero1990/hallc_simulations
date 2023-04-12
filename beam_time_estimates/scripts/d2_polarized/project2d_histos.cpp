//#ifndef PROJECT_2D_H
//#define PROJECT_2D_H

/*
  Author: C. Yero
  Date: Feb 15, 2023 
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;



// helper function : check if a string only contain numebers or not
int is_digits(string& str)
{
    for (char ch : str) {
        int v = ch; // ASCII Val converted
        if (!(ch >= 48 && ch <= 57)) {   // The ASCII value of '0' is 48 , '9' is 57.
            
	  return 0;  // characters in string have letters
        }
	
	
    }
 
    return 1; // no charactes have letter (therefore true number)
}

void project2d_deut( TH2F *hist2d=0, TString setting="", Bool_t display_plots=0 ){

  //avoid display
  gROOT->SetBatch(kTRUE);
  
  /*
    Brief: Projects 2d histograms onto 1d slices in X or Y
    and plots it on a canvas subplot, and also plots relative
    errors on a separate canvas subplot

    For now is specific for deuteron, but can easily be modified for other use
   */
  
  TString basename=Form("deut_stats_monitoring_setting_%s_", setting.Data());
  
  // Set basefilename to save projections to root file
  TString ofile="DEUT_OUTPUT/ROOT/" + basename + "output.root";
  
  TFile *fout = new TFile(ofile.Data(), "RECREATE");

  //set output names of projections to be saved to divided canvas
  TString fout_2dHist      = "./" + basename + "2Dhisto.pdf";
  TString fout_projHist    = "./" + basename + "projY.pdf";
  TString fout_projHistErr = "./" + basename + "projY_relError.pdf";
  
  // set global title/label sizes
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetLabelSize(.1, "XY");
  gStyle->SetTitleY(1.01); // offset title vertically

  //TH2F *hist2d = new TH2F("hist2d", "", 30,-15,15, 30, -15,15);


  // get total number of x,y bins of hist2d
  int nxbins = hist2d->GetXaxis()->GetNbins();
  int nybins = hist2d->GetYaxis()->GetNbins();

  hist2d->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");

  // plot 2d histogra
  TCanvas *c0 = new TCanvas("c0", "", 1200,800);
  c0->cd();
  hist2d->Draw("contz");
  hist2d->Write();
  
  // --- define variables for calculative/plotting of relative stats. error on Pmiss (need to reset vector per projection bin) ---
  vector<double> y_val;       // this will be set to 0 (as reference)
  vector<double> y_err;       // relative error on pmiss bin counts
  vector<double> x_val;       // this is the central value of x (pmiss)
  vector<double> x_err;       // this is actually width on x-axis 
  
  double pm_counts;
  double relative_err;
  double pm_center;
  double thrq_center;
  double thrq_width;
  
  // ----------------------------

  
  // --- define canvas subplots (x,y) based on number of bins in hist2d ---
  int yc, xc;
  // if projection is along y, should make a canvas square out of nxbins
  float remainder = sqrt(nxbins) - int(sqrt(nxbins));

  if(remainder>0 && remainder<0.5){
    yc = round(sqrt(nxbins));
    xc = round(sqrt(nxbins))+1;
  }
  else{
    yc = round(sqrt(nxbins));
    xc = round(sqrt(nxbins));
  }
  
  //cout << Form("rounded area = (%d, %d) ", yc, xc) << endl;
  TCanvas *c1 = new TCanvas("c1", "", 1400,1100);
  TCanvas *c2 = new TCanvas("c2", "", 1400,1100);
  
  c1->Divide(yc, xc);
  c2->Divide(yc, xc);
  //----------------------------

  TH1D *h1 = 0;
  //loop over xbins of hist2d 
  for(int i=1; i<=nxbins; i++){

    // get xbin center value and width
    float thrq_center  = hist2d->GetXaxis()->GetBinCenter(i);
    float thrq_width   = hist2d->GetXaxis()->GetBinWidth(i); 

    //cout << "bin: " << i << ", x-val: " << thrq_center << endl;

    // project hist2d along y-axis (different bins in x)
    h1 = hist2d->ProjectionY(Form("proj_Pm_thrq%.1f", thrq_center), i, i);

    // define integrated counts on projected bin
    float counts = h1->Integral();

    // reduce divisions of histos (for de-cluttering)
    h1->GetYaxis()->SetNdivisions(5);
    h1->GetXaxis()->SetNdivisions(10);
    h1->GetXaxis()->SetLabelSize(0.1);
    //cout << "thrq_center, counts (v1) = " << thrq_center << ", " << counts << endl;
    h1->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    h1->SetTitleSize(10);

   
    //cout << Form("bin#: %d,  x-center: %.1f, counts: %.3f", i, thrq_center, counts) << endl;

    
    //-----------------------------------
    // Compute Relative Error on Counts
    //-----------------------------------
    
    // For each h1 hsitogram, get bin content to calculate its error and plot relative error
    y_val.clear();
    y_err.clear();
    x_val.clear();
    x_err.clear();
    
    pm_counts = 0;
    relative_err = 0;
    
    // set statistical lower (inner) limit for guidance during online data-taking
    int inner_stats = 15;  // +/- 15 %
  
    // loop over all bins of Pmiss (h1 histogram) 
    for(int i =1; i<=h1->GetNbinsX(); i++){

      // get bin content
      pm_counts = h1->GetBinContent(i);

      // get bin center
      pm_center = h1->GetBinCenter(i);
      
      // calculate relative error of bin
      relative_err = (sqrt(pm_counts) / pm_counts ) * 100.; // in %

      // push values to vector
      y_val.push_back(0);
      y_err.push_back( relative_err );
      x_val.push_back(pm_center);
      x_err.push_back(0); // no need to set this width
      
      
    }

    // at the end, should have vector of length N for plotting
    int n=h1->GetNbinsX();
    
    TGraphErrors *gr = new TGraphErrors(n, &x_val[0], &y_val[0], &x_err[0], &y_err[0]);
    gr->SetTitle(Form("proj_Pm_thrq%.1f_relErr", thrq_center));
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(0.);
    gr->SetMarkerStyle(21);

    

    TLine *lo_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), -inner_stats, *(max_element(x_val.begin(), x_val.end())) , -inner_stats);
    TLine *up_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), inner_stats, *(max_element(x_val.begin(), x_val.end())) , inner_stats);

    lo_limit->SetLineColor(kRed);
    lo_limit->SetLineStyle(1);
    lo_limit->SetLineWidth(1);

    up_limit->SetLineColor(kRed);
    up_limit->SetLineStyle(1);
    up_limit->SetLineWidth(1);
    
  
    
    //cout << "thrq_center, counts (v2) = " << thrq_center << ", " << counts << endl;
    gr->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    
    // reduce divisions of histos (for de-cluttering)
    gr->GetYaxis()->SetNdivisions(5);
    gr->GetXaxis()->SetNdivisions(10);


    //---------------------------------------------------
    
    // change to canvas subplot of ith bin projection
    c1->cd(i);
    gPad->Modified(); gPad->Update();
    h1->DrawClone("histE0");
    h1->Write();
      
    c2->cd(i);                                                                                                                                           
    gPad->Modified(); gPad->Update();
    // draw to graph
    gr->Draw("AP");
    lo_limit->Draw();
    up_limit->Draw();
     
    // add legend
    if(i==1){
      auto legend = new TLegend(0.5,0.6,0.9,0.8);
      legend->AddEntry("gr",Form("#pm %d %%",inner_stats),"%d");
      legend->SetBorderSize(0);
      legend->SetTextSize(0.15);
      legend->SetTextColor(kRed);
      legend->Draw();
      /*
      auto txt = new TPaveLabel(0.3, 0.3, 0.95, 0.95, "#pm 15 %", "br");
      txt -> SetFillColor(0);
      txt -> SetTextColor(kRed);
      txt -> Draw();
      */
    }
    

    gr->Write();

  } // end loop over 2D xbins [th_rq]

  // save canvas
  gStyle->SetOptStat(0);
  c0->SaveAs( fout_2dHist.Data()      );
  c1->SaveAs( fout_projHist.Data()    );
  c2->SaveAs( fout_projHistErr.Data() );

  if(display_plots){
    
    //open plots with evince or any other viewer
    //gSystem->Exec(Form("evince %s", fout_projHistErr.Data() )); 
    gSystem->Exec(Form("open %s", fout_2dHist.Data() ));
    gSystem->Exec(Form("open %s", fout_projHist.Data() ));
    gSystem->Exec(Form("open %s", fout_projHistErr.Data() )); 

  }
  
}


void project2d_histos() {


  TString file;

  file = "d2_pm300_Q2_3p5_rad_output.root";
  
 
  
  
  // for each root file, get the desired histogram to be combined
  TFile *data_file =  new TFile(file.Data(), "READ");

  TH2F *myhist2d = 0;
  data_file->GetObject("kin_plots/H_Pm_vs_thrq", myhist2d);
  

  // func2: projectes the 2d histo (slices of x-bins along y-axis) onto 1d bins, the pm_setting is just for histogram naming purposes
  // and should be consistent with the histogram range chosen
  project2d_deut( myhist2d, "300_Q2_3.5", true ); // the bool flag is to display the plots (otherwise, they will be saved)

}

//#endif
