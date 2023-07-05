#ifndef HIST_UTILS_H
#define HIST_UTILS_H

using namespace std;

void extract_2d_hist(TH2F *h2, TString xlabel, TString ylabel, TString out_fname)
{

  
  //Brief: This function extracts 2D Histogram Bin Information
  


  //Output File Stream to write datafile
  ofstream ofile;

  //----------------------------------
  //----Histogram Bin Information-----
  //----------------------------------
  
  //Get the bin content (and error) corresponding to the 2d bin (y,x)
  Double_t z_cont, z_cont_err;
  
  //Get Number of Bins in Y and X-axis
  Int_t y_Nbins = h2->GetNbinsY();
  Int_t x_Nbins = h2->GetNbinsX();

  //Get Bin Number for x, y and a global bin numver for each 2d bin
  Int_t xb, yb, ib;
  
  //Get Bin Width (Assuming same width for all bins)
  Double_t ybin_width = h2->GetYaxis()->GetBinWidth(1);
  Double_t xbin_width = h2->GetXaxis()->GetBinWidth(1);
  
  //Get Central Value / Low Edge / UpEdge for each bin
  Double_t y0, ylow, yup;
  Double_t x0, xlow, xup;


  ofile.open(out_fname.Data());

  TString header = Form("# 2D Histogram Bin Extraction \n"
			"#                   \n"
			"# histogram parameters:  \n"
			"# ybins      = %d   | ylabel: %s \n"
			"# xbins      = %d   | xlabel: %s \n"
			"# ybin width = %.3f \n"
			"# xbin width = %.3f \n"
			"#                   \n"
			"# header definitions:\n"
			"# ib:         (x,y)  bin number \n"
			"# xb:         x-axis bin number \n"
			"# yb:         y-axis bin number \n"
			"# x0:         x-axis central bin value \n"
			"# xlow:       x-axis low-edge bin value \n"
			"# xup:        x-axis up-edge bin value \n"
			"# y0:         y-axis central bin value \n"
			"# ylow:       y-axis low-edge bin value \n"
			"# yup:        y-axis up-edge bin value \n"
			"# zcont:      bin content (z-axis) \n"
			"# zcont_err:  bin content error (z-axis) \n"
			"#                                        \n"
			"ib,xb,yb,x0,xlow,xup,y0,ylow,yup,zcont,zcont_err"
			,y_Nbins, ylabel.Data(), x_Nbins, xlabel.Data(), ybin_width, xbin_width  );

  //Write header to data file
  ofile << header.Data() << endl;
  
  //loop over y-bins
  for(int yb=1; yb<=y_Nbins; yb++){

    //loop over x-bins
    for(int xb=1; xb<=x_Nbins; xb++){

      //Get 2d bin number
      ib = h2->GetBin(xb, yb);

      //Get bin content and error for each (y, x) bin
      z_cont = h2->GetBinContent(xb, yb);
      z_cont_err = h2->GetBinError(xb, yb);
      
      //Get value at bin center / low / up edges of the bin
      ylow = h2->GetYaxis()->GetBinLowEdge(yb);
      yup  = h2->GetYaxis()->GetBinUpEdge(yb);
      y0   = (ylow + yup)/2.; //h2->GetYaxis()->GetBinCenter(yb);
      
      xlow = h2->GetXaxis()->GetBinLowEdge(xb);
      xup  = h2->GetXaxis()->GetBinUpEdge(xb);
      x0   = (xlow + xup)/2.; //h2->GetXaxis()->GetBinCenter(xb);

      
      ofile  << std::setw(7) << ib << std::setw(10) << xb << std::setw(10) << yb << std::setw(14) << x0 << std::setw(10) << xlow << std::setw(10) << xup << std::setw(12) << y0 << std::setw(10) << ylow << std::setw(12) << yup << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;

    }// end x-bins loop
    
  } // end y-bins loop

  ofile.close();

    
}

void extract_1d_hist(TH1F *h1, TString xlabel, TString ylabel, TString out_fname)
{

  
  //Brief: This function extracts 2D Histogram Bin Information
  


  //Output File Stream to write datafile
  ofstream ofile;

  //----------------------------------
  //----Histogram Bin Information-----
  //----------------------------------
  
  //Get the bin content (and error) corresponding to the 1d bin (y,y_err)
  Double_t y_cont, y_cont_err;
  
  //Get Number of Bins in  X-axis
  Int_t x_Nbins = h1->GetNbinsX();

  //Get Bin Number for x and a global bin numver for each 2d bin
  Int_t xb, ib;
  
  //Get Bin Width (Assuming same width for all bins)
  Double_t xbin_width = h1->GetXaxis()->GetBinWidth(1);
  
  //Get Central Value / Low Edge / UpEdge for each bin
  Double_t x0, xlow, xup;


  ofile.open(out_fname.Data());

  TString header = Form("# 1D Histogram Bin Extraction \n"
			"#                   \n"
			"# histogram parameters:  \n"
			"# xbins      = %d   | xlabel: %s \n"
			"# xbin width = %.3f \n"
			"#                   \n"
			"# header definitions:\n"
			"# xb:         x-axis bin number \n"
			"# x0:         x-axis central bin value \n"
			"# xlow:       x-axis low-edge bin value \n"
			"# xup:        x-axis up-edge bin value \n"
			"# ycont:      bin content (y-axis) \n"
			"# ycont_err:  bin content error (y-axis) \n"
			"#                                        \n"
			"xb,x0,xlow,xup,ycont,ycont_err"
			, x_Nbins, xlabel.Data(),xbin_width  );

  //Write header to data file
  ofile << header.Data() << endl;
  
  //loop over x-bins
  for(int xb=1; xb<=x_Nbins; xb++){
    
    //Get bin content and error for each (y, x) bin
    y_cont = h1->GetBinContent(xb);
    y_cont_err = h1->GetBinError(xb);
      
    xlow = h1->GetXaxis()->GetBinLowEdge(xb);
    xup  = h1->GetXaxis()->GetBinUpEdge(xb);
    x0   = (xlow + xup)/2.; //h2->GetXaxis()->GetBinCenter(xb);

    
    ofile  << Form( "%.1d,%.3f,%.3f,%.3f,%.3f,%.3f" , xb, x0, xlow, xup, y_cont, y_cont_err) << endl;
    
  }// end x-bins loop
  
  ofile.close();

    
}



#endif
