#ifndef HIST_UTILS_H
#define HIST_UTILS_H

using namespace std;

//Brief: This function extracts 1D Histogram Bin Information
// Useful for extracting bin contents, and errors for plotting in python
void extract_1d_hist(TH1F *h1, TString xlabel, TString out_fname)
{

  /*
    User Arguments:
    h1: 1D Histogram Object
    xlabel, ylabel : axes label
    out_fname: output file name to write 1D Histo information

  */
 
  //Output File Stream to write datafile
  ofstream ofile;

  //----------------------------------
  //----Histogram Bin Information-----
  //----------------------------------
  
  //Get the bin content (and error) corresponding to the 2d bin (y,x)
  Double_t z_cont, z_cont_err;
  
  //Get Number of Bins in X-axis
  Int_t x_Nbins = h1->GetNbinsX();

  //Get Bin Number for x and a global bin number for each 1d bin
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
			"# ib:          x-bin number \n"
			"# xb:          x-axis bin number \n"
			"# x0:          x-axis central bin value \n"
			"# xlow:        x-axis low-edge bin value \n"
			"# xup:         x-axis up-edge bin value \n"
			"# bin_cnt:     bin content (y-axis) \n"
			"# bin_cnt_err: bin content error (y-axis) \n"
			"#                                        \n"
			"#! ib[i,0]/  xb[i,1]/  x0[f,2]/  xlow[f,3]/  xup[f,4]/  bin_cnt[f,5]/  bin_cnt_err[f,6]/  "
			,x_Nbins, xlabel.Data(), xbin_width  );

  //Write header to data file
  ofile << header.Data() << endl;
  
  //loop over x-bins
  for(int xb=1; xb<=x_Nbins; xb++){
    
    //Get 1d bin number
    ib = h1->GetBin(xb);

    //Get bin content and error for each (y, x) bin
    z_cont = h1->GetBinContent(xb);
    z_cont_err = h1->GetBinError(xb);
    
    //Get value at bin center / low / up edges of the bin
    xlow = h1->GetXaxis()->GetBinLowEdge(xb);
    xup  = h1->GetXaxis()->GetBinUpEdge(xb);
    x0   = (xlow + xup)/2.; //h2->GetXaxis()->GetBinCenter(xb);
    
    
    ofile  << std::setw(7) << ib << std::setw(10) << xb << std::setw(14) << x0 << std::setw(10) << xlow << std::setw(10) << xup << std::setw(12)  << z_cont << std::setw(14) << z_cont_err << endl;
    
  }// end x-bins loop
  

  ofile.close();

    
}

//Brief: This function extracts 2D Histogram Bin Information
// Useful for extracting yields or cross sections binned in 2D kinematics
void extract_2d_hist(TH2F *h2, TString xlabel, TString ylabel, TString out_fname, Bool_t csv_flg=false)
{

  /*
    User Arguments:
    h2: 2D Histogram Object
    xlabel, ylabel : axes label
    out_fname: output file name to write 2D Histo information

  */
 
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
  TString header;
  if(csv_flg){
    header = Form("# 2D Histogram Bin Extraction \n"
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
		  "ib, xb, yb,  x0,  xlow,  xup,  y0,  ylow,  yup,  zcont,  zcont_err,  "
		  ,y_Nbins, ylabel.Data(), x_Nbins, xlabel.Data(), ybin_width, xbin_width  );
  }
  else{
    header = Form("# 2D Histogram Bin Extraction \n"
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
		  "#! ib[i,0]/  xb[i,1]/  yb[i,2]/  x0[f,3]/  xlow[f,4]/  xup[f,5]/  y0[f,6]/  ylow[f,7]/  yup[f,8]/  zcont[f,9]/  zcont_err[f,10]/  "
		  ,y_Nbins, ylabel.Data(), x_Nbins, xlabel.Data(), ybin_width, xbin_width  );
  }
  
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

      if(csv_flg){
	ofile  << Form("%i, %i, %i, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.3E, %.3E \n", ib, xb, yb, x0, xlow, xup, y0, ylow, yup, z_cont, z_cont_err);
      }
      else{
	ofile  << std::setw(7) << ib << std::setw(10) << xb << std::setw(10) << yb << std::setw(14) << x0 << std::setw(10) << xlow << std::setw(10) << xup << std::setw(12) << y0 << std::setw(10) << ylow << std::setw(12) << yup << std::setw(12) << z_cont << std::setw(14) << z_cont_err << endl;
      }
    }// end x-bins loop
    
  } // end y-bins loop

  ofile.close();

    
}



#endif
