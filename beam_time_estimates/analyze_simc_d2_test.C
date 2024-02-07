#include "utils/parse_utils.h"
#include "utils/hist_utils.h"


TRotation       fToLabRot;              //Rotation matrix from TRANSPORT to lab
Double_t        fThetaGeo;              //In-plane geographic central angle (rad)
Double_t        fPhiGeo;                //Out-of-plane geographic central angle (rad)
Double_t        fThetaSph, fPhiSph;     //Central angles in spherical coords. (rad)
Double_t        fSinThGeo, fCosThGeo;   //Sine and cosine of central angles
Double_t        fSinPhGeo, fCosPhGeo;   // in geographical coordinates
Double_t        fSinThSph, fCosThSph;   //Sine and cosine of central angles in 
Double_t        fSinPhSph, fCosPhSph;   // spherical coordinates


//_____________________________________________________
void GeoToSph( Double_t  th_geo, Double_t  ph_geo, Double_t& th_sph, Double_t& ph_sph)
{
  
  // Convert geographical to spherical angles. Units are rad.
  
  static const Double_t twopi = 2.0*TMath::Pi();
  Double_t ct = cos(th_geo), cp = cos(ph_geo);
  Double_t tmp = ct*cp;
  th_sph = acos( tmp );
  tmp = sqrt(1.0 - tmp*tmp);
  ph_sph = (fabs(tmp) < 1e-6 ) ? 0.0 : acos( sqrt(1.0-ct*ct)*cp/tmp );
  if( th_geo/twopi-floor(th_geo/twopi) > 0.5 ) ph_sph = TMath::Pi() - ph_sph;
  if( ph_geo/twopi-floor(ph_geo/twopi) > 0.5 ) ph_sph = -ph_sph;
}

//_______________________________________________________________
void SetCentralAngles(Double_t th_cent=0, Double_t ph_cent=0)
{
  
    
  fThetaGeo = TMath::DegToRad()*th_cent; fPhiGeo = TMath::DegToRad()*ph_cent;
  GeoToSph( fThetaGeo, fPhiGeo, fThetaSph, fPhiSph );
  fSinThGeo = TMath::Sin( fThetaGeo ); fCosThGeo = TMath::Cos( fThetaGeo );
  fSinPhGeo = TMath::Sin( fPhiGeo );   fCosPhGeo = TMath::Cos( fPhiGeo );
  Double_t st, ct, sp, cp;
  st = fSinThSph = TMath::Sin( fThetaSph ); ct = fCosThSph = TMath::Cos( fThetaSph );
  sp = fSinPhSph = TMath::Sin( fPhiSph );   cp = fCosPhSph = TMath::Cos( fPhiSph );
  
  Double_t norm = TMath::Sqrt(ct*ct + st*st*cp*cp);
  TVector3 nx( st*st*sp*cp/norm, -norm, st*ct*sp/norm );
  TVector3 ny( ct/norm,          0.0,   -st*cp/norm   );
  TVector3 nz( st*cp,            st*sp, ct            );
  
  fToLabRot.SetToIdentity().RotateAxes( nx, ny, nz );
  }

//____________________________________________________________________________________
void TransportToLab( Double_t p, Double_t xptar, Double_t yptar, TVector3& pvect) 
{
  
  TVector3 v( xptar, yptar, 1.0 );
  v *= p/TMath::Sqrt( 1.0+xptar*xptar+yptar*yptar );
  pvect = fToLabRot * v;
}



void analyze_simc_d2_test(TString basename="", Bool_t heep_check=false){
  
  
  /* 
     User Input:
     basename: generic file name used in input and simulated root file
     e.g. basename: d2_pm400_Q2_4.0_fsi_rad.root (polarized deuteron proposal generic filename)
     e.g. basename: c12_pm300_Q2_3p5_fsi_rad (polarized deuteron proposal bkg estimates generic filename)
     e.g. basename: d2_pm800_thrq55_fsi_rad.root (FSI deuteron proposal generic filename
  
     Instructions:
     1) need to select which analysis to do via 'analysis_flag' option: polarized deuteron ("d2pol") or deuteron fsi studies ("d2fsi")
     2) check the input and output generic filenames to be read have the proper path and structure
     search for lines containing 'input_',  'simc_' and 'output_', which will have the path 
  */
  
  using namespace std;
  
  TString analysis_flag="d2pol";  // "d2pol" or "d2fsi" 
  
  Bool_t debug = true;
  
  int pm_set=0;
  int thrq_set=0;
  float Q2_set=0.;

  TString tgt_name;
  TString thrq_str;
  TString pm_str;
  TString Q2_str;
  
  TString model="";


  


  
  // extract the kinematics setting from the filename
  if( analysis_flag == "d2fsi") {
    
    // split base into kin setting values
    thrq_str = split(split(split(basename.Data(), '_')[0], '_')[0], '_')[1];
    pm_str = split(split(split(split(basename.Data(), '_')[0], '_')[0], '_')[0], '_')[1];

    // strip non-numeric characters and conver to int
    thrq_set = stoi(std::regex_replace(thrq_str.Data(), std::regex(R"([^\d])"), ""));
    pm_set = stoi(std::regex_replace(pm_str.Data(), std::regex(R"([^\d])"), ""));
    model = split(split(basename.Data(), '_')[0], '_')[1];
  }
  
  
  if( analysis_flag == "d2pol") {

    // split base into kin setting values
    Q2_str = split(split(split(basename.Data(), '_')[0], '_')[0], '_')[1];
    
    pm_str =split( split(split(split(split(basename.Data(), '_')[0], '_')[0], '_')[0], '_')[0], '_')[1];
    
    // strip non-numeric characters and conver to int
    Q2_str = std::regex_replace(Q2_str.Data(), std::regex(R"([^\d])"), "."); // replace the "p" with decimal point "."
    Q2_set = stof(Q2_str.Data());
    pm_set = stoi(std::regex_replace(pm_str.Data(), std::regex(R"([^\d])"), ""));
    model = split(split(basename.Data(), '_')[0], '_')[1];
    tgt_name = split(split(split(split(split(basename.Data(), '_')[0], '_')[0], '_')[0], '_')[0], '_')[0];
      
    if(debug) {
      cout << "Settings Read: " << endl;
      cout << Form("target: %s", tgt_name.Data()) << endl;
      cout << Form("Q2_set: %.1f", Q2_set) << endl;
      cout << Form("Pm_set: %d ", pm_set) << endl;
      cout << "Model: " << model.Data() << endl;
    }
  }
  
  
  TString h_arm_name = "HMS";
  TString e_arm_name = "SHMS";
  
  Double_t N_av = 6.023e23;  // avogadro's number
  Double_t gram2GeV = 5.6095883571872e32 / 1e9 ;  //(5.6e32 eV * 1 GeV / 1e9 eV) = 1 gram
  Double_t pi = 3.141592654;
  Double_t dtr = pi/180.;
  Double_t MP = 0.938272; //GeV
  Double_t MD = 1.87561; //GeV
  Double_t MN = 0.939566; //GeV
  Double_t me = 0.000510998; //GeV


  
  
  //-------------------
  // READ FILENAMES
  //-------------------
  
  //Input parameter controls filenames
  TString input_CutFileName;
  TString input_HBinFileName;

  //Declare TFile Pointers (reading/writing ROOTfiles)
  TFile *inROOT;
  TFile *outROOT;

  //SIMC Kin. Input Filename
  TString simc_infile;

  //Input/Outpu ROOTfile Name (to be read)
  TString simc_InputFileName;
  TString simc_OutputFileName;
  
  TString temp; //temporary string placeholder

  //Output .txt files
  TString output_file;
  TString output_hist_data;

      
  if( analysis_flag == "d2fsi") {
    
    //---Read In File Names with cuts and histogram binning information
    input_CutFileName  = "inp/d2_fsi/set_basic_cuts_d2fsi.inp";
    input_HBinFileName = "inp/d2_fsi/set_basic_histos_d2fsi.inp";
    
    //Define File Name Patterns
    simc_infile         = Form("infiles/deuteron/fsi_deuteron/Q2_4p5/%s.data",    basename.Data());
    simc_InputFileName  = Form("worksim/d2_fsi/raw/%s.root",                      basename.Data());
    simc_OutputFileName = Form("worksim/d2_fsi/analyzed/%s_output.root",          basename.Data());
    
    // define output file to write the rates
    output_file = "yield_estimates/d2_fsi/output_rates_d2fsi.txt";
    
    // define output directory where numerical histogram .txt will be placed
    output_hist_data= Form("yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s", pm_set, thrq_set, model.Data());
  }
  
  if( analysis_flag == "d2pol") {
    
    //---Read In File Names with cuts and histogram binning information
    input_CutFileName  = "inp/d2_pol/set_basic_cuts_d2pol.inp";
    input_HBinFileName = "inp/d2_pol/set_basic_histos_d2pol.inp";
    
    //Define File Name Patterns
    simc_infile         = Form("infiles/deuteron/d2_polarized/smallFSI/%s.data",    basename.Data());
    simc_InputFileName  = Form("worksim/d2_pol/smallFSI/optimized/raw/%s.root",                      basename.Data());
    simc_OutputFileName = Form("worksim/d2_pol/smallFSI/optimized/analyzed/%s_output.root",          basename.Data());
    
    // define output file to write the rates
    output_file = "yield_estimates/d2_pol/smallFSI/optimized/output_rates_d2pol_optim.txt";
    
    // define output directory where numerical histogram .txt will be placed
    output_hist_data= Form("yield_estimates/d2_pol/smallFSI/optimized/histogram_data/pm%d_Q2_%.1f_%s", pm_set, Q2_set, model.Data());
    
    
    if (debug) {
      cout << "---- Set Filenames ----" << endl;
      cout << Form("simc_infile: %s", simc_infile.Data()) << endl;
      cout << Form("raw ROOTfile: %s", simc_InputFileName.Data()) << endl;
	cout << Form("analyzed ROOTfile: %s", simc_OutputFileName.Data()) << endl;
	
	
    }
    
  }
    
  //---------------------------------------------------------------------------------------------------------
  
    //----------------------------
    // READ CENTRAL KIN. SETTINGS
    //----------------------------
  
    // READ TARGET MASS in amu (and convert amu to GeV)  1 amu = 1 gram/mol = 1 gram x 1 mol / N_avogadro x 5.62e23 GeV / 1 gram
    Double_t tgt_mass = ( stod(split(split(FindString("targ%mass_amu", simc_infile.Data())[0], '!')[0], '=')[1]) ) * gram2GeV / N_av;
    //beam energy (GeV)
    Double_t beam_e = (stod(split(split(FindString("Ebeam", simc_infile.Data())[0], '!')[0], '=')[1]))/1000.; 
    //e- arm central momentum setting (GeV/c)
    Double_t e_Pcen = (stod(split(split(FindString("spec%e%P", simc_infile.Data())[0], '!')[0], '=')[1]))/1000.; 
    //proton arm central momentum setting (GeV/c)
    Double_t h_Pcen = (stod(split(split(FindString("spec%p%P", simc_infile.Data())[0], '!')[0], '=')[1]))/1000.; 
    
    //e- arm angle (deg)
    Double_t the_central = (stod(split(split(FindString("spec%e%theta", simc_infile.Data())[0], '!')[0], '=')[1])); 
    //p arm angle (deg)
    Double_t thp_central = (stod(split(split(FindString("spec%p%theta", simc_infile.Data())[0], '!')[0], '=')[1])); 
    
    if(debug) {
  cout << "SIMC KIN" << endl;
  cout << Form("tgt_mass [GeV]: %.3f", tgt_mass) << endl;
  cout << Form("beam_e: %.3f", beam_e) << endl;
  cout << Form("e_Pcen: %.3f", e_Pcen) << endl; 
  cout << Form("h_Pcen: %.3f", h_Pcen) << endl; 
  cout << Form("the_central: %.3f", the_central) << endl; 
  cout << Form("thp_central: %.3f", thp_central) << endl; 
  cout << "-----------------" << endl;
    }
    
    // spectrometers are in-plane
    Double_t phe_central=0;  
    Double_t php_central=0;
    
    
    //--------------------
    // READ ANALYSIS CUTS
    //-------------------

    //-------Collimator Study-------
    Bool_t hmsCollCut_flag;      //flag to enable/disable collimator cut
    Bool_t shmsCollCut_flag;
    
    Bool_t hmsColl_Cut;
    Bool_t shmsColl_Cut;
    
    TCutG *hms_Coll_gCut;   //HMS Collimator Graphical Cut
    TCutG *shms_Coll_gCut;  //SHMS Collimator Graphical Cut
    
    //HMS Octagonal Collimator Size (Each of the octagonal points is a multiple of 1 or 1/2 of these values)
    Double_t hms_hsize = 4.575;  //cm
    Double_t hms_vsize = 11.646;
    
    //SHMS Octagonal Collimator Size (Each of the octagonal points is a multiple of 1 or 1/2 of these values)
    Double_t shms_hsize = 8.5;  //cm
    Double_t shms_vsize = 12.5;
    
    
    //------------------------------
    
    //---------Kinematics Cuts---------- 
    //4-Momentum Transfers
    Bool_t Q2_cut_flag = stoi(split(FindString("Q2_cut_flag", input_CutFileName.Data())[0], '=')[1]);
    Bool_t c_Q2;
    Double_t c_Q2_min = stod(split(FindString("c_Q2_min", input_CutFileName.Data())[0], '=')[1]);;
    Double_t c_Q2_max = stod(split(FindString("c_Q2_max", input_CutFileName.Data())[0], '=')[1]);
    
    //Missing Energy
    Bool_t Em_cut_flag = stoi(split(FindString("Em_cut_flag", input_CutFileName.Data())[0], '=')[1]);
    Bool_t c_Em;
    Double_t c_Em_min = stod(split(FindString("c_Em_min", input_CutFileName.Data())[0], '=')[1]);
    Double_t c_Em_max = stod(split(FindString("c_Em_max", input_CutFileName.Data())[0], '=')[1]);
    
    
    //----------Acceptance Cuts------------  
    // Hadron Arm Momentum Acceptance, Delta [%] |
    Bool_t hdelta_cut_flag = stoi(split(FindString("hdelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
    Bool_t c_hdelta;
    Double_t c_hdelta_min = stod(split(FindString("c_hdelta_min", input_CutFileName.Data())[0], '=')[1]);
    Double_t c_hdelta_max = stod(split(FindString("c_hdelta_max", input_CutFileName.Data())[0], '=')[1]);
    
    // Electron Arm Momentum Acceptance, Delta [%] 
    Bool_t edelta_cut_flag = stoi(split(FindString("edelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
    Bool_t c_edelta;
    Double_t c_edelta_min = stod(split(FindString("c_edelta_min", input_CutFileName.Data())[0], '=')[1]);
    Double_t c_edelta_max = stod(split(FindString("c_edelta_max", input_CutFileName.Data())[0], '=')[1]);
    
    // Angular Acceptance Cut
    
    //Scaling factor to scale collimator cuts from original size cut
    Double_t hms_scale=1.;   //Default
    Double_t shms_scale=1.;
    
    //==Collimator==
    hmsCollCut_flag = stoi(split(FindString("hmsCollCut_flag", input_CutFileName.Data())[0], '=')[1]);
    shmsCollCut_flag = stoi(split(FindString("shmsCollCut_flag", input_CutFileName.Data())[0], '=')[1]);
    
    //==Read Collimator Cut Scale Factor==
    hms_scale  =  stod(split(FindString("hms_scale", input_CutFileName.Data())[0], '=')[1]);
    shms_scale =  stod(split(FindString("shms_scale", input_CutFileName.Data())[0], '=')[1]);
    
    
    
    // Z-Reaction Vertex Difference Cut
    Bool_t ztarDiff_cut_flag = stoi(split(FindString("ztarDiff_cut_flag", input_CutFileName.Data())[0], '=')[1]);
    Bool_t c_ztarDiff;
    Double_t c_ztarDiff_min = stod(split(FindString("c_ztarDiff_min", input_CutFileName.Data())[0], '=')[1]);
    Double_t c_ztarDiff_max = stod(split(FindString("c_ztarDiff_max", input_CutFileName.Data())[0], '=')[1]);
    
    //COMBINE ALL CUTS
    Bool_t c_allCuts;    //combined cuts
    Bool_t c_accpCuts;  //acceptance cuts
    Bool_t c_kinCuts;   //kinematics cuts
    
    //---------------------------------------------------------------------------------------------------------
    
    //--------------------
    // READ HISTO BINS
    //-------------------
    
    //---------------------------------
    // Kinematics Histograms Binning
    //---------------------------------
    
    //Primary Kinematics
    Double_t kf_nbins	= stod(split(FindString("kf_nbins",  	input_HBinFileName.Data())[0], '=')[1]);   //final electron momentum
    Double_t kf_xmin	= stod(split(FindString("kf_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t kf_xmax	= stod(split(FindString("kf_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t the_nbins 	= stod(split(FindString("the_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t the_xmin  	= stod(split(FindString("the_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t the_xmax  	= stod(split(FindString("the_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t Q2_nbins  	= stod(split(FindString("Q2_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t Q2_xmin   	= stod(split(FindString("Q2_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t Q2_xmax   	= stod(split(FindString("Q2_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t X_nbins   	= stod(split(FindString("X_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t X_xmin    	= stod(split(FindString("X_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t X_xmax    	= stod(split(FindString("X_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t nu_nbins  	= stod(split(FindString("nu_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t nu_xmin   	= stod(split(FindString("nu_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t nu_xmax   	= stod(split(FindString("nu_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t q_nbins   	= stod(split(FindString("q_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t q_xmin    	= stod(split(FindString("q_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t q_xmax    	= stod(split(FindString("q_xmax",  input_HBinFileName.Data())[0], '=')[1]);              				          
    
    Double_t thq_nbins 	= stod(split(FindString("thq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t thq_xmin  	= stod(split(FindString("thq_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t thq_xmax  	= stod(split(FindString("thq_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t W_nbins  =  stod(split(FindString("W_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t W_xmin   =  stod(split(FindString("W_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t W_xmax   =  stod(split(FindString("W_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    //Secondary Kinematics
    Double_t Pf_nbins    = stod(split(FindString("Pf_nbins",  	input_HBinFileName.Data())[0], '=')[1]);   //final proton momentum
    Double_t Pf_xmin     = stod(split(FindString("Pf_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pf_xmax     = stod(split(FindString("Pf_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t thp_nbins   = stod(split(FindString("thx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  //proton(hadron) angle
    Double_t thp_xmin    = stod(split(FindString("thx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t thp_xmax    = stod(split(FindString("thx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t Em_nbins   = stod(split(FindString("Em_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Em_xmin    = stod(split(FindString("Em_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Em_xmax    = stod(split(FindString("Em_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t Pm_nbins   = stod(split(FindString("Pm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pm_xmin    = stod(split(FindString("Pm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pm_xmax    = stod(split(FindString("Pm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);               				          	       
    
    Double_t Pmx_nbins   = stod(split(FindString("Pmx_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pmx_xmin    = stod(split(FindString("Pmx_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pmx_xmax    = stod(split(FindString("Pmx_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);               				        	  	       
    
    Double_t Pmy_nbins   = stod(split(FindString("Pmy_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pmy_xmin    = stod(split(FindString("Pmy_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pmy_xmax    = stod(split(FindString("Pmy_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);               				                 	       				  	       
    
    Double_t Pmz_nbins   = stod(split(FindString("Pmz_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pmz_xmin    = stod(split(FindString("Pmz_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t Pmz_xmax    = stod(split(FindString("Pmz_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);               				                 	       				  	       
    
    
    Double_t MM_nbins   = stod(split(FindString("MM_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t MM_xmin    = stod(split(FindString("MM_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t MM_xmax    = stod(split(FindString("MM_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t MM2_nbins  = stod(split(FindString("MM2_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
    Double_t MM2_xmin   = stod(split(FindString("MM2_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t MM2_xmax   = stod(split(FindString("MM2_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t thpq_nbins  = stod(split(FindString("thxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t thpq_xmin   = stod(split(FindString("thxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t thpq_xmax   = stod(split(FindString("thxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t phpq_nbins  = stod(split(FindString("phxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t phpq_xmin   = stod(split(FindString("phxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t phpq_xmax   = stod(split(FindString("phxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t thrq_nbins  = stod(split(FindString("thrq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t thrq_xmin   = stod(split(FindString("thrq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
    Double_t thrq_xmax   = stod(split(FindString("thrq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);               				              
    
    
    //---------------------------------
    // Acceptance Histograms Binning
    //---------------------------------
    
    //----Electron Arm Focal Plane-----
    Double_t exfp_nbins 	= stod(split(FindString("exfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t exfp_xmin  	= stod(split(FindString("exfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t exfp_xmax  	= stod(split(FindString("exfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t expfp_nbins	= stod(split(FindString("expfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t expfp_xmin 	= stod(split(FindString("expfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t expfp_xmax 	= stod(split(FindString("expfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t eyfp_nbins 	= stod(split(FindString("eyfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eyfp_xmin  	= stod(split(FindString("eyfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eyfp_xmax  	= stod(split(FindString("eyfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t eypfp_nbins	= stod(split(FindString("eypfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eypfp_xmin 	= stod(split(FindString("eypfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eypfp_xmax 	= stod(split(FindString("eypfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    //----Electron Arm Reconstructed-----
    Double_t eytar_nbins 	= stod(split(FindString("eytar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eytar_xmin  	= stod(split(FindString("eytar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eytar_xmax  	= stod(split(FindString("eytar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t eyptar_nbins	= stod(split(FindString("eyptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eyptar_xmin 	= stod(split(FindString("eyptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eyptar_xmax 	= stod(split(FindString("eyptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t exptar_nbins	= stod(split(FindString("exptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t exptar_xmin 	= stod(split(FindString("exptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t exptar_xmax 	= stod(split(FindString("exptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t edelta_nbins	= stod(split(FindString("edelta_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t edelta_xmin 	= stod(split(FindString("edelta_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
    Double_t edelta_xmax 	= stod(split(FindString("edelta_xmax",  input_HBinFileName.Data())[0], '=')[1]);   
    
    
    //----Hadron Arm Focal Plane-----
    Double_t hxfp_nbins 	= stod(split(FindString("hxfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hxfp_xmin  	= stod(split(FindString("hxfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hxfp_xmax  	= stod(split(FindString("hxfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t hxpfp_nbins	= stod(split(FindString("hxpfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hxpfp_xmin 	= stod(split(FindString("hxpfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hxpfp_xmax 	= stod(split(FindString("hxpfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t hyfp_nbins 	= stod(split(FindString("hyfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hyfp_xmin  	= stod(split(FindString("hyfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hyfp_xmax  	= stod(split(FindString("hyfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t hypfp_nbins	= stod(split(FindString("hypfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hypfp_xmin 	= stod(split(FindString("hypfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hypfp_xmax 	= stod(split(FindString("hypfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    //----Hadron Arm Reconstructed-----
    Double_t hytar_nbins 	= stod(split(FindString("hytar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hytar_xmin  	= stod(split(FindString("hytar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hytar_xmax  	= stod(split(FindString("hytar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t hyptar_nbins	= stod(split(FindString("hyptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hyptar_xmin 	= stod(split(FindString("hyptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hyptar_xmax 	= stod(split(FindString("hyptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t hxptar_nbins	= stod(split(FindString("hxptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hxptar_xmin 	= stod(split(FindString("hxptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hxptar_xmax 	= stod(split(FindString("hxptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t hdelta_nbins	= stod(split(FindString("hdelta_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hdelta_xmin 	= stod(split(FindString("hdelta_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
    Double_t hdelta_xmax 	= stod(split(FindString("hdelta_xmax",  input_HBinFileName.Data())[0], '=')[1]);   
    
    
    //----Target Quantities----
    //(Use same binning for hadron/electron reconstructed at target)
    Double_t tarx_nbins	= stod(split(FindString("tarx_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t tarx_xmin 	= stod(split(FindString("tarx_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t tarx_xmax 	= stod(split(FindString("tarx_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t tary_nbins	= stod(split(FindString("tary_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t tary_xmin 	= stod(split(FindString("tary_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t tary_xmax 	= stod(split(FindString("tary_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t tarz_nbins	= stod(split(FindString("tarz_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t tarz_xmin 	= stod(split(FindString("tarz_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t tarz_xmax 	= stod(split(FindString("tarz_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t ztar_diff_nbins	= stod(split(FindString("ztar_diff_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t ztar_diff_xmin 	= stod(split(FindString("ztar_diff_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t ztar_diff_xmax 	= stod(split(FindString("ztar_diff_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    
    //----Collimator Quantities----
    Double_t hXColl_nbins	= stod(split(FindString("hXColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t hXColl_xmin 	= stod(split(FindString("hXColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
    Double_t hXColl_xmax 	= stod(split(FindString("hXColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);   
    
    Double_t hYColl_nbins	= stod(split(FindString("hYColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);                                           
    Double_t hYColl_xmin 	= stod(split(FindString("hYColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);                                                     
    Double_t hYColl_xmax 	= stod(split(FindString("hYColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t eXColl_nbins	= stod(split(FindString("eXColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eXColl_xmin 	= stod(split(FindString("eXColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);
    Double_t eXColl_xmax 	= stod(split(FindString("eXColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    Double_t eYColl_nbins	= stod(split(FindString("eYColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);      
    Double_t eYColl_xmin 	= stod(split(FindString("eYColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);                                                                      
    Double_t eYColl_xmax 	= stod(split(FindString("eYColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);
    
    
    //---------------------------------------------------------------------------------------------------------
    
    //--------------------
    // DECLARE HISTOGRAMS
    //-------------------
    
    //Create TLists to store categorical histograms
    TList *kin_HList  = new TList();
    TList *accp_HList = new TList();
    
    
    //--------------------------------------------------------
    //---------HISTOGRAM CATEGORY: Kinematics  (KIN)----------
    //--------------------------------------------------------
    
    //Primary (electron) Kinematics (14 histos)
    TH1F *H_kf      = new TH1F("H_kf", "Final e^{-} Momentum; e^{-} momentum, k_{f} [GeV/c]; Counts", kf_nbins, kf_xmin, kf_xmax);
    H_kf->Sumw2();
    H_kf->SetDefaultSumw2(); 
    TH1F *H_the     = new TH1F("H_the",   "Electron Scattering Angle, #theta_{e}; e^{-} scattering angle, #theta_{e}  [deg]; Counts", the_nbins, the_xmin, the_xmax);
    TH1F *H_Q2      = new TH1F("H_Q2",    "4-Momentum Transfer, Q^{2}; 4-momentum transfer, Q^{2} [GeV^{2}]; Counts", Q2_nbins, Q2_xmin, Q2_xmax);
    TH1F *H_Q2_nsc  = new TH1F("H_Q2_nsc","4-Momentum Transfer, Q^{2}; 4-momentum transfer, Q^{2} [GeV^{2}]; Counts", Q2_nbins, Q2_xmin, Q2_xmax);  //nsc stands for no self-cut (i.e, all cuts except on itself)
    
    TH1F *H_xbj     = new TH1F("H_xbj",     "x-Bjorken; x-bjorken, x_{Bj}; Counts", X_nbins, X_xmin, X_xmax);  
    TH1F *H_nu      = new TH1F("H_nu",      "Energy Transfer, #nu; energy transfer, #nu [GeV]; Counts", nu_nbins, nu_xmin, nu_xmax); 
    TH1F *H_q       = new TH1F("H_q",       "3-Momentum Transfer, |#vec{q}|; 3-momentum transfer, |#vec{q}| [GeV]; Counts", q_nbins, q_xmin, q_xmax);
    TH1F *H_thq     = new TH1F("H_thq",     "In-Plane Angle w.r.t +z(lab), #theta_{q}; in-plane recoil angle, #theta_{q}; Counts", thq_nbins, thq_xmin, thq_xmax); 
    TH1F *H_W       = new TH1F("H_W",       "Invariant Mass, W; invariant mass, W [GeV]; Counts", W_nbins, W_xmin, W_xmax);  
    TH1F *H_W_noCut = new TH1F("H_W_noCut", "Invariant Mass, W (no cuts, DAQ rates); invariant mass, W [Gev]; Counts", W_nbins, W_xmin, W_xmax);  
    
    //Secondary (Hadron) Kinematics (recoil and missing are used interchageably) ()
    TH1F *H_Pf      = new TH1F("H_Pf", "Final Hadron Momentum (detected), p_{f}; final hadron momentum, p_{f} [GeV/c]; Counts", Pf_nbins, Pf_xmin, Pf_xmax);
    TH1F *H_thp     = new TH1F("H_thp", "Hadron Scattering Angle (detected), #theta_{p}; hadron scattering angle, #theta_{p} [deg]; Counts", thp_nbins, thp_xmin, thp_xmax);
    TH1F *H_Em      = new TH1F("H_Em","Missing Energy; missing energy, E_{m} [GeV]; Counts", Em_nbins, Em_xmin, Em_xmax); 
    TH1F *H_Em_nuc      = new TH1F("H_Em_nuc","Nuclear Missing Energy; nuclear missing energy, E_{m, nuc} [GeV]; Counts", Em_nbins, Em_xmin, Em_xmax); 
    TH1F *H_Em_nuc_nsc      = new TH1F("H_Em_nuc_nsc","Nuclear Missing Energy; nuclear missing energy, E_{m} [GeV]; Counts", Em_nbins, Em_xmin, Em_xmax); 
    
    TH1F *H_Pm      = new TH1F("H_Pm","Missing Momentum, P_{miss}; missing momentum, P_{m} [GeV/c]; Counts", Pm_nbins, Pm_xmin, Pm_xmax);
    TH1F *H_Pm_noCut      = new TH1F("H_Pm_noCut","Missing Momentum, P_{miss} (no cuts, DAQ rates); missing momentum, P_{m} [GeV/c]; Counts", Pm_nbins, Pm_xmin, Pm_xmax); 
    
    TH1F *H_MM      = new TH1F("H_MM","Missing Mass, M_{miss}; missing mass, M_{miss} [GeV]", MM_nbins, MM_xmin, MM_xmax);        
    TH1F *H_MM2     = new TH1F("H_MM2","Missing Mass Squared, M^{2}_{miss}; missing mass, MM^{2} [GeV^{2}]; Counts", MM2_nbins, MM2_xmin, MM2_xmax); 
    TH1F *H_thpq    = new TH1F("H_thpq", "In-Plane Angle, #theta_{pq}; in-plane angle, #theta_{pq} [deg]; Counts", thpq_nbins, thpq_xmin, thpq_xmax);
    TH1F *H_thrq    = new TH1F("H_thrq", "In-Plane Angle, #theta_{rq}; in-plane angle, #theta_{rq} [deg]; Counts", thrq_nbins, thrq_xmin, thrq_xmax);
    
    TH1F *H_phi_pq  = new TH1F("H_phi_pq", "Out-of-Plane Angle, #phi_{pq}; out-of-plane angle, #phi_{pq} [deg]; Counts", phpq_nbins, phpq_xmin, phpq_xmax);
    TH1F *H_cphi_pq  = new TH1F("H_cphi_pq", "Out-of-Plane Angle, cos(#phi_{pq}); out-of-plane angle, cos(#phi_{pq}); Counts", 100, -1.2, 1.2);
    
    
    //2D Pm vs. thrq (for cross section calculation)
    TH2F *H_Pm_vs_thrq       = new TH2F("H_Pm_vs_thrq",      "Pm vs. #theta_{rq} (yield); in-plane angle, #theta_{rq} [deg]; missing momentum, P_{m} [GeV/c]", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_Pm_vs_thrq_ps    = new TH2F("H_Pm_vs_thrq_ps",   "Pm vs. #theta_{rq} (phase space); in-plane angle, #theta_{rq} [deg]; missing momentum, P_{m} [GeV/c]", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_Pm_vs_thrq_xsec  = new TH2F("H_Pm_vs_thrq_xsec", "Pm vs. #theta_{rq} (xsec); in-plane angle, #theta_{rq} [deg]; missing momentum, P_{m} [GeV/c]", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    
    //SIMC 2D Average Kinematics Histograms (Pmiss vs. th_rq averaged over different kinematics) 
    TH2F *H_Pm_vs_thrq_v   = new TH2F("H_Pm_vs_thrq_v", "Pm vs. #theta_{rq} (vertex); #theta_{rq} [deg]; P_{m} [GeV/c]", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_Ein_2Davg      = new TH2F("H_Ein_2Davg", "Ein (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_kf_2Davg       = new TH2F("H_kf_2Davg", "Final e^{-} Momentum (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_the_2Davg      = new TH2F("H_the_2Davg", "Electron Scattering Angle (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax); 
    TH2F *H_Pf_2Davg       = new TH2F("H_Pf_2Davg", "Final Proton Momentum (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_thp_2Davg      = new TH2F("H_thp_2Davg", "Proton Scattering Angle (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax); 
    TH2F *H_q_2Davg           = new TH2F("H_q_2Davg", "q-vector, |q| (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_theta_q_2Davg     = new TH2F("H_theta_q_2Davg", "#theta_{q} (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax); 
    TH2F *H_Q2_2Davg          = new TH2F("H_Q2_2Davg","Q2 (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax); 
    TH2F *H_nu_2Davg          = new TH2F("H_nu_2Davg","Energy Transfer, #nu (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax); 
    TH2F *H_xbj_2Davg         = new TH2F("H_xbj_2Davg", "x-Bjorken (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);  
    TH2F *H_theta_pq_2Davg    = new TH2F("H_theta_pq_2Davg", "#theta_{pq} (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_phi_pq_2Davg    = new TH2F("H_phi_pq_2Davg", "#phi_{pq} (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_cphi_pq_2Davg     = new TH2F("H_cphi_pq_2Davg", "cos(#phi_{pq}) (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_sphi_pq_2Davg     = new TH2F("H_sphi_pq_2Davg", "sin(#phi_{pq}) (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    TH2F *H_Pm_2Davg     = new TH2F("H_Pm_2Davg","Missing Momentum (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax); 
    TH2F *H_thrq_2Davg   = new TH2F("H_thrq_2Davg", "#theta_{rq} (2D Average); #theta_{rq} [deg]; P_{m} [GeV/c]",thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    
    // define some vertex histos for comparing with recon. histos
    TH1F *H_Em_nuc_v      = new TH1F("H_Em_nuc_v","Nuclear Missing Energy; nuclear missing energy, E_{m, nuc} [GeV]; Counts", Em_nbins, Em_xmin, Em_xmax); 
    TH1F *H_Pm_v      = new TH1F("H_Pm_v","Missing Momentum, P_{miss} (vertex); missing momentum, P_{m} [GeV/c]; Counts", Pm_nbins, Pm_xmin, Pm_xmax);
    TH1F *H_Q2_v      = new TH1F("H_Q2_v",    "4-Momentum Transfer, Q^{2} (vertex); 4-momentum transfer, Q^{2} [GeV^{2}]; Counts", Q2_nbins, Q2_xmin, Q2_xmax);
    TH1F *H_the_v     = new TH1F("H_the_v",   "Electron Scattering Angle, #theta_{e} (vertex); e^{-} scattering angle, #theta_{e}  [deg]; Counts", the_nbins, the_xmin, the_xmax);
    TH1F *H_thpq_v    = new TH1F("H_thpq_v", "In-Plane Angle, #theta_{pq} (vertex); in-plane angle, #theta_{pq} [deg]; Counts", thpq_nbins, thpq_xmin, thpq_xmax);
    TH1F *H_thrq_v    = new TH1F("H_thrq_v", "In-Plane Angle, #theta_{rq} (vertex); in-plane angle, #theta_{rq} [deg]; Counts", thrq_nbins, thrq_xmin, thrq_xmax);
    TH1F *H_phi_pq_v  = new TH1F("H_phi_pq_v", "Out-of-Plane Angle, #phi_{pq} (vertex); out-of-plane angle, #phi_{pq} [deg]; Counts", phpq_nbins, phpq_xmin, phpq_xmax);
    TH1F *H_cphi_pq_v  = new TH1F("H_cphi_pq_v", "Out-of-Plane Angle, cos(#phi_{pq}); out-of-plane angle, cos(#phi_{pq}); Counts", 100, -1.2, 1.2);

    
    //Add Kin Histos to TList
    
    //Add Primary Kin Histos
    kin_HList->Add( H_kf     );
    kin_HList->Add( H_the    );
    kin_HList->Add( H_Q2     );
    kin_HList->Add( H_Q2_nsc );
    kin_HList->Add( H_xbj    );
    kin_HList->Add( H_nu     );
    kin_HList->Add( H_q      );
    kin_HList->Add( H_thq    );
    kin_HList->Add( H_W      );
    kin_HList->Add( H_W_noCut      );
    
    //Add Secondary Kin Histos
    kin_HList->Add( H_Pf       );
    kin_HList->Add( H_thp      );
    kin_HList->Add( H_Em       );
    kin_HList->Add( H_Em_nuc       );
    kin_HList->Add( H_Em_nuc_nsc   );
    kin_HList->Add( H_Pm       );
    kin_HList->Add( H_Pm_noCut       );
    
    kin_HList->Add( H_MM       );
    kin_HList->Add( H_MM2      );
    kin_HList->Add( H_thpq     );
    kin_HList->Add( H_thrq     );
    kin_HList->Add( H_phi_pq   );
    kin_HList->Add( H_cphi_pq  );
    
    kin_HList->Add( H_Pm_vs_thrq );
    kin_HList->Add( H_Pm_vs_thrq_ps );
    kin_HList->Add( H_Pm_vs_thrq_xsec );
    
    
    // Add averaged kin. histos
    kin_HList->Add( H_Pm_vs_thrq_v );
    kin_HList->Add( H_Ein_2Davg    );
    kin_HList->Add( H_kf_2Davg     );
    kin_HList->Add( H_the_2Davg    );
    kin_HList->Add( H_thp_2Davg    );
    kin_HList->Add( H_Pf_2Davg     );
    kin_HList->Add( H_Pm_2Davg     );
    kin_HList->Add( H_thrq_2Davg   );
    
    kin_HList->Add( H_q_2Davg          );
    kin_HList->Add( H_theta_q_2Davg    );
    kin_HList->Add( H_Q2_2Davg         );
    kin_HList->Add( H_nu_2Davg         );
    kin_HList->Add( H_xbj_2Davg        );
    kin_HList->Add( H_theta_pq_2Davg   );
    kin_HList->Add( H_phi_pq_2Davg     );
    kin_HList->Add( H_cphi_pq_2Davg    );
    kin_HList->Add( H_sphi_pq_2Davg    );
    

    kin_HList->Add( H_Em_nuc_v    );
    kin_HList->Add( H_Pm_v    );
    kin_HList->Add( H_Q2_v    );
    kin_HList->Add( H_the_v   );
    kin_HList->Add( H_thpq_v  );
    kin_HList->Add( H_thrq_v  );
    kin_HList->Add( H_phi_pq_v);
    kin_HList->Add( H_cphi_pq_v);
    
  
  //----------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
  //----------------------------------------------------------------------


  //Electron Arm Focal Plane Quantities
  TH1F *H_exfp = new TH1F("H_exfp", Form("%s X_{fp}; X_{fp} [cm]; Counts ", e_arm_name.Data()), exfp_nbins, exfp_xmin, exfp_xmax);
  TH1F *H_eyfp = new TH1F("H_eyfp", Form("%s Y_{fp}; Y_{fp} [cm]; Counts ", e_arm_name.Data()), eyfp_nbins, eyfp_xmin, eyfp_xmax);
  TH1F *H_expfp = new TH1F("H_expfp", Form("%s X'_{fp}; X'_{fp} [rad]; Counts ", e_arm_name.Data()), expfp_nbins, expfp_xmin, expfp_xmax);
  TH1F *H_eypfp = new TH1F("H_eypfp", Form("%s Y'_{fp}; Y'_{fp} [rad]; Counts ", e_arm_name.Data()), eypfp_nbins, eypfp_xmin, eypfp_xmax);
  
  //Electron Arm Reconstructed Quantities 
  TH1F *H_eytar = new TH1F("H_eytar", Form("%s Y_{tar}; Y_{tar} [cm]; Counts ", e_arm_name.Data()), eytar_nbins, eytar_xmin, eytar_xmax);
  TH1F *H_exptar = new TH1F("H_exptar", Form("%s X'_{tar}; X'_{tar} [rad]; Counts ", e_arm_name.Data()), exptar_nbins, exptar_xmin, exptar_xmax);
  TH1F *H_eyptar = new TH1F("H_eyptar", Form("%s Y'_{tar}; Y'_{tar} [rad]; Counts ", e_arm_name.Data()), eyptar_nbins, eyptar_xmin, eyptar_xmax);
  TH1F *H_edelta = new TH1F("H_edelta", Form("%s Momentum Acceptance, #delta; #delta [%%]; Counts ", e_arm_name.Data()), edelta_nbins, edelta_xmin, edelta_xmax);
  TH1F *H_edelta_nsc = new TH1F("H_edelta_nsc", Form("%s Momentum Acceptance, #delta; #delta [%%]; Counts ", e_arm_name.Data()), edelta_nbins, edelta_xmin, edelta_xmax);
  
  //Hadron arm Focal Plane Quantities
  TH1F *H_hxfp = new TH1F("H_hxfp", Form("%s  X_{fp}; X_{fp} [cm]; Counts ", h_arm_name.Data()), hxfp_nbins, hxfp_xmin, hxfp_xmax);
  TH1F *H_hyfp = new TH1F("H_hyfp", Form("%s  Y_{fp}; Y_{fp} [cm]; Counts ", h_arm_name.Data()), hyfp_nbins, hyfp_xmin, hyfp_xmax);
  TH1F *H_hxpfp = new TH1F("H_hxpfp", Form("%s  X'_{fp}; X'_{fp} [rad]; Counts ", h_arm_name.Data()), hxpfp_nbins, hxpfp_xmin, hxpfp_xmax );
  TH1F *H_hypfp = new TH1F("H_hypfp", Form("%s  Y'_{fp}; Y'_{fp} [rad]; Counts ", h_arm_name.Data()), hypfp_nbins, hypfp_xmin, hypfp_xmax);

  //Hadron arm Reconstructed Quantities 
  TH1F *H_hytar = new TH1F("H_hytar", Form("%s  Y_{tar}; Y_{tar} [cm]; Counts ", h_arm_name.Data()), hytar_nbins, hytar_xmin, hytar_xmax);
  TH1F *H_hxptar = new TH1F("H_hxptar", Form("%s  X'_{tar}; X'_{tar} [rad]; Counts ", h_arm_name.Data()), hxptar_nbins, hxptar_xmin, hxptar_xmax);
  TH1F *H_hyptar = new TH1F("H_hyptar", Form("%s  Y'_{tar}; Y'_{tar} [rad]; Counts ", h_arm_name.Data()), hyptar_nbins, hyptar_xmin, hyptar_xmax );
  TH1F *H_hdelta = new TH1F("H_hdelta", Form("%s  Momentum Acceptance, #delta; #delta [%%]; Counts ", h_arm_name.Data()), hdelta_nbins, hdelta_xmin, hdelta_xmax);
  TH1F *H_hdelta_nsc = new TH1F("H_hdelta_nsc", Form("%s  Momentum Acceptance, #delta; #delta [%%]; Counts ", h_arm_name.Data()), hdelta_nbins, hdelta_xmin, hdelta_xmax);
  

  //Target Reconstruction (Hall Coord. System) 
  TH1F *H_htar_x = new TH1F("H_htar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts ", h_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  TH1F *H_htar_y = new TH1F("H_htar_y", Form("%s y_Target (Lab); y-Target [cm]; Counts ", h_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  TH1F *H_htar_z = new TH1F("H_htar_z", Form("%s z_Target (Lab); z-Target [cm]; Counts ", h_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);
  TH1F *H_etar_x = new TH1F("H_etar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts ", e_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  TH1F *H_etar_y = new TH1F("H_etar_y", Form("%s y-Target (Lab); y-Target [cm]; Counts ", e_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  TH1F *H_etar_z = new TH1F("H_etar_z", Form("%s z-Target (Lab); z-Target [cm]; Counts ", e_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);

  //difference in reaction vertex z (user-defined)
  TH1F *H_ztar_diff = new TH1F("H_ztar_diff", "Ztar Difference; z-Target Difference [cm]; Counts ", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);
  TH1F *H_ztar_diff_nsc = new TH1F("H_ztar_diff_nsc", "Ztar Difference; z-Target Difference [cm]; Counts ", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);


  //2D Collimator Histos
  TH2F *H_hXColl_vs_hYColl = new TH2F("H_hXColl_vs_hYColl", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", h_arm_name.Data(), h_arm_name.Data(), h_arm_name.Data()), hYColl_nbins, hYColl_xmin, hYColl_xmax,  hXColl_nbins, hXColl_xmin, hXColl_xmax);
  TH2F *H_eXColl_vs_eYColl = new TH2F("H_eXColl_vs_eYColl", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", e_arm_name.Data(), e_arm_name.Data(), e_arm_name.Data()), eYColl_nbins, eYColl_xmin, eYColl_xmax, eXColl_nbins, eXColl_xmin, eXColl_xmax); 
  
  TH2F *H_hXColl_vs_hYColl_nsc = new TH2F("H_hXColl_vs_hYColl_nsc", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", h_arm_name.Data(), h_arm_name.Data(), h_arm_name.Data()), hYColl_nbins, hYColl_xmin, hYColl_xmax,  hXColl_nbins, hXColl_xmin, hXColl_xmax);
  TH2F *H_eXColl_vs_eYColl_nsc = new TH2F("H_eXColl_vs_eYColl_nsc", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", e_arm_name.Data(), e_arm_name.Data(), e_arm_name.Data()), eYColl_nbins, eYColl_xmin, eYColl_xmax, eXColl_nbins, eXColl_xmin, eXColl_xmax); 

  
  //2D Hour Glass Histos
  TH2F *H_hxfp_vs_hyfp  = new TH2F("H_hxfp_vs_hyfp", Form("%s  X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", h_arm_name.Data()),  hyfp_nbins, hyfp_xmin, hyfp_xmax, hxfp_nbins, hxfp_xmin, hxfp_xmax);
  TH2F *H_exfp_vs_eyfp  = new TH2F("H_exfp_vs_eyfp", Form("%s  X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", e_arm_name.Data()),  eyfp_nbins, eyfp_xmin, eyfp_xmax, exfp_nbins, exfp_xmin, exfp_xmax);

  //2D HMS v. SHMS Acceptance Correlations
  TH2F *H_hxptar_vs_exptar = new TH2F("H_hxptar_vs_exptar", "HMS vs. SHMS, X'_{tar}", exptar_nbins, exptar_xmin, exptar_xmax, hxptar_nbins, hxptar_xmin, hxptar_xmax);
  TH2F *H_hyptar_vs_eyptar = new TH2F("H_hyptar_vs_eyptar", "HMS vs. SHMS, Y'_{tar}", eyptar_nbins, eyptar_xmin, eyptar_xmax, hyptar_nbins, hyptar_xmin, hyptar_xmax);
  TH2F *H_hdelta_vs_edelta = new TH2F("H_hdelta_vs_edelta", "HMS vs. SHMS, #delta",   edelta_nbins, edelta_xmin, edelta_xmax, hdelta_nbins, hdelta_xmin, hdelta_xmax);
  
  

  
  //Add ACCP Histos to TList
  accp_HList->Add( H_exfp       );
  accp_HList->Add( H_eyfp       );
  accp_HList->Add( H_expfp      );
  accp_HList->Add( H_eypfp      );

  accp_HList->Add( H_eytar       );
  accp_HList->Add( H_exptar      );
  accp_HList->Add( H_eyptar      );
  accp_HList->Add( H_edelta      );
  accp_HList->Add( H_edelta_nsc  );
  
  accp_HList->Add( H_hxfp       );
  accp_HList->Add( H_hyfp       );
  accp_HList->Add( H_hxpfp      );
  accp_HList->Add( H_hypfp      );
  
  accp_HList->Add( H_hytar       );
  accp_HList->Add( H_hxptar      );
  accp_HList->Add( H_hyptar      );
  accp_HList->Add( H_hdelta      );
  accp_HList->Add( H_hdelta_nsc  );

  accp_HList->Add( H_htar_x       );
  accp_HList->Add( H_htar_y       );
  accp_HList->Add( H_htar_z       );
  accp_HList->Add( H_etar_x       );
  accp_HList->Add( H_etar_y       );
  accp_HList->Add( H_etar_z       );
  accp_HList->Add( H_ztar_diff    );
  accp_HList->Add( H_ztar_diff_nsc);

  accp_HList->Add( H_hXColl_vs_hYColl  );
  accp_HList->Add( H_hXColl_vs_hYColl_nsc );
  accp_HList->Add( H_eXColl_vs_eYColl  );
  accp_HList->Add( H_eXColl_vs_eYColl_nsc );
  
  accp_HList->Add( H_hxfp_vs_hyfp  );
  accp_HList->Add( H_exfp_vs_eyfp  );

  accp_HList->Add( H_hxptar_vs_exptar );
  accp_HList->Add( H_hyptar_vs_eyptar );
  accp_HList->Add( H_hdelta_vs_edelta );
  
  //---------------------------------------------------------------------------------------------------------
  
  //--------------------------------
  // READ TREE / SET BRANCH ADDRESS
  //--------------------------------

  TTree *tree;
  Long64_t nentries;
  //Read ROOTfile
  inROOT = new TFile(simc_InputFileName.Data(), "READ");

  //Get the data tree
  tree = (TTree*)inROOT->Get("SNT");
  nentries = tree->GetEntries();

  
  //--- Define Variables for Calculation ----
  
  // 4-vectors 
  TLorentzVector fP0;
  TLorentzVector fP1;
  TLorentzVector fA;
  TLorentzVector fMp; 
  TLorentzVector fA1;
  TLorentzVector fMp1;
  TLorentzVector fQ;
  TLorentzVector fX;
  TLorentzVector fB;

  // 3-vectors
  TVector3 Pf_vec;             
  TVector3 kf_vec;             
  TVector3 bq;   
  TVector3 xq;   
  TVector3 p_miss_q;  
  
  // rotaion from +z to +q
  TRotation rot_to_q;
 

  //Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)
  Double_t Ei;
  Double_t ki;                    //initial electron momentum  
  Double_t kf;                    //final electron momentum
  Double_t the;               //Central electron arm angle relative to +z (hall coord. system)
  Double_t Q2;                   //Four-momentum trasfer
  Double_t nu;                   //Energy Transfer
  Double_t q;                  //magnitude of the 3-vector q 
  Double_t X;                    //B-jorken X  scaling variable
  Double_t th_q;                 //angle between q and +z (hall coord. system)
  Double_t ph_q;                 //out-of-plane angle between q and +z (hall coord. system)
  Double_t W;                    // invariant mass
  Double_t W2;                    // invariant mass squared
  
  //Secondary Kinematics (detected hadron /  recoil system)
  Double_t Ep;                      //final proton energy (needs to be calculated)
  Double_t Pf;                     //final proton momentum
  Double_t fXangle;                    // opening angle between HMS/SHMS
  Double_t thp;                      // proton angle
  Double_t Erecoil;               //Total energy of recoil system (GeV)
  Double_t Em;                     //Standard Missing Energy 
  Double_t Em_nuc;                     //nuclear definition of Missing Energy (nu - det kin E - recoil kin E)
  Double_t Pm;                     //Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
  Double_t Pmx_lab;                //x-comp. Missing Momentum 
  Double_t Pmy_lab;                 //y-comp. Missing Momentum 
  Double_t Pmz_lab;                 //z-comp. Missing Momentum
  Double_t Pmx_q;                //x-comp. Missing Momentum in q-frame  (perpendicular)?
  Double_t Pmy_q;                 //y-comp. Missing Momentum in q-frame (out-of-plane) ?
  Double_t Pmz_q;                 //z-comp. Missing Momentum in q-frame (parallel)
  Double_t Tx;                     // detected particle X kinetic energy
  Double_t Tr;                    // recoil system kinetic energy  
  Double_t MM;                   //Missing Mass (neutron Mass)
  Double_t MM2;                   //Missing Mass Squared
  Double_t th_pq;                  //detected particle in-plane angle w.r.to q-vector
  Double_t th_rq;                  //recoil particle in-plane angle w.r.to q-vector
  Double_t ph_pq;                   // out-of-plane angle between scattering and reaction plane (Oop of proton or detected particle X)
  Double_t ph_rq;                   // out-of-plane angle between scattering and reaction plane (Oop of proton or detected particle X)
  
  //Electron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t e_xfp;
  Double_t e_xpfp;
  Double_t e_yfp;
  Double_t e_ypfp;
  
  Double_t e_ytar;
  Double_t e_yptar;
  Double_t e_xptar;
  Double_t e_delta;

  //Hadron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t h_xfp;
  Double_t h_xpfp;
  Double_t h_yfp;
  Double_t h_ypfp;
  
  Double_t h_ytar;
  Double_t h_yptar;
  Double_t h_xptar;
  Double_t h_delta;
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System (USED BY DATA AND SIMC)
  Double_t tar_x;    //For SIMC ONLY (It is the same for HMS/SHMS)
  Double_t  htar_y;
  Double_t  htar_z;
  Double_t  etar_y;
  Double_t  etar_z;

  Double_t ztar_diff;

  //Collimators
  Double_t hXColl, hYColl, eXColl, eYColl;
 
  //SIMC Specific TTree Variable Names
  Double_t Normfac;               //normalization factor, defined as : normfac = luminosity * accepted / generated, luminosity = EXPER charge / targetfac
  Double_t Weight;               //This Weight has the cross section in it
  Double_t Jacobian_corr;
  Double_t prob_abs;  // Probability of absorption of particle in the HMS Collimator
                      //(Must be multiplies by the weight. If particle interation is
                      //NOT simulated, it is set to 1.)  
  //SIMC Collimator
  Double_t htarx_corr;
  Double_t etarx_corr;

  //--------------------------
  // --- Vertex Quantities ---
  //--------------------------
  
  // 4-vectors (@ vertex)
  TLorentzVector fP0_v;           // Beam 4-momentum
  TLorentzVector fP1_v;           // Scattered electron 4-momentum
  TLorentzVector fA_v;            // Target 4-momentum
  TLorentzVector fMp_v;           // 4-momentum of target (assuming proton)
  TLorentzVector fA1_v;           // Final system 4-momentum
  TLorentzVector fMp1_v;
  TLorentzVector fQ_v;            // Momentum transfer 4-vector
  TLorentzVector fX_v;            // Detected secondary particle 4-momentum (GeV)
  TLorentzVector fB_v;            // Recoil system 4-momentum (GeV)

  // 3-vectors (@ vertex)
  TVector3 Pf_vec_v;             
  TVector3 kf_vec_v;             
  TVector3 bq_v;  
  TVector3 xq_v;   
  TVector3 p_miss_q_v;  
  
  // rotaion from +z to +q
  TRotation rot_to_q_v;

  //Primary (electron) Kinematics (at vertex)
  Double_t Ein_v;               //incident beam energy at vertex (simulates external rad. has rad. tail) ???
  Double_t Ef_v;                //final electron energy at the vertex
  Double_t ki_v;
  Double_t kf_v;
  Double_t the_v;
  Double_t Q2_v;                //Q2 (vertex)  
  Double_t nu_v;                //energy transfer = Ein_v - Ef_v
  Double_t q_v;             //magintude of 3-vector q
  Double_t X_v;
  Double_t thq_v;
  Double_t phq_v;
  Double_t W_v;
  Double_t W2_v;

  
  //Secondary Kinematics (detected hadron /  recoil system)
  Double_t Ep_v;                
  Double_t Pf_v;               
  Double_t fXangle_v;
  Double_t thp_v;
  Double_t Erecoil_v;          
  Double_t Em_v;                
  Double_t Em_nuc_v;            
  Double_t Pm_v;                
  Double_t Pmx_lab_v;
  Double_t Pmy_lab_v;
  Double_t Pmz_lab_v;
  Double_t Pmx_q_v;
  Double_t Pmy_q_v;
  Double_t Pmz_q_v;
  Double_t Tx_v;  
  Double_t Tr_v;   
  Double_t MM_v;  
  Double_t MM2_v;
  Double_t th_pq_v;     
  Double_t th_rq_v;    
  Double_t ph_pq_v;     
  Double_t ph_rq_v;      

 //Electron Arm Focal Plane / Reconstructed Quantities (@ vertex)
  Double_t e_xptar_v;           
  Double_t e_yptar_v;
  Double_t h_xptar_v;
  Double_t h_yptar_v;

  //Light-Cone Momentum Variables (at vertex)
  Double_t PmPar_v;     //parallel component of recoil momentum relative to q-vector
  Double_t PmPerp_v;    //transverse component of recoil momentum relative to q-vector
  Double_t alpha_n_v;   //light-cone momentum fraction of the recoil neutron
  Double_t alpha_v;     //momentum fraction of struck nucleon (normalized such that: alpha + alpha_n = 2)

  
  //----- Set Branch Address ------

  //Primary Kinematics (electron kinematics)
  //ki needs to be calculated (initial e- momentum)
  tree->SetBranchAddress("Ein", &Ei); 
  tree->SetBranchAddress("e_pf", &kf);
  tree->SetBranchAddress("theta_e", &the);
  tree->SetBranchAddress("Q2", &Q2);  
  //Xbj needs to be calculated in the event loop
  tree->SetBranchAddress("nu", &nu);
  tree->SetBranchAddress("q", &q);
  //th_q needs to be calculated in the event loop
  tree->SetBranchAddress("W", &W);

  //Secondary Kinematics (hadron kinematics)
  tree->SetBranchAddress("h_pf",    &Pf);
  tree->SetBranchAddress("theta_p", &thp);
  tree->SetBranchAddress("Em", &Em_nuc);
  tree->SetBranchAddress("Pm", &Pm);
  tree->SetBranchAddress("Pmx", &Pmx_lab);
  tree->SetBranchAddress("Pmy", &Pmy_lab);
  tree->SetBranchAddress("Pmz", &Pmz_lab);

  //Missing Mass (MM) and MM2 will be defined in entry loop
  tree->SetBranchAddress("theta_pq", &th_pq); 
  tree->SetBranchAddress("theta_rq", &th_rq);  
  tree->SetBranchAddress("phi_pq", &ph_pq);  
 
  
  //Electron Arm Focal Plane / Reconstructed Quantities 
  tree->SetBranchAddress("e_xfp",  &e_xfp);
  tree->SetBranchAddress("e_xpfp", &e_xpfp);
  tree->SetBranchAddress("e_yfp",  &e_yfp);
  tree->SetBranchAddress("e_ypfp", &e_ypfp);
  
  tree->SetBranchAddress("e_ytar",  &e_ytar);
  tree->SetBranchAddress("e_yptar", &e_yptar);
  tree->SetBranchAddress("e_xptar", &e_xptar);
  tree->SetBranchAddress("e_delta", &e_delta);
  
  //Hadron Arm Focal Plane / Reconstructed Quantities 
  tree->SetBranchAddress("h_xfp",  &h_xfp);
  tree->SetBranchAddress("h_xpfp", &h_xpfp);
  tree->SetBranchAddress("h_yfp",  &h_yfp);
  tree->SetBranchAddress("h_ypfp", &h_ypfp);
        
  tree->SetBranchAddress("h_ytar",  &h_ytar);
  tree->SetBranchAddress("h_yptar", &h_yptar);
  tree->SetBranchAddress("h_xptar", &h_xptar);
  tree->SetBranchAddress("h_delta", &h_delta);
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System
  tree->SetBranchAddress("tar_x", &tar_x);
  tree->SetBranchAddress("h_yv",  &htar_y);
  tree->SetBranchAddress("h_zv",  &htar_z);
  tree->SetBranchAddress("e_yv",  &etar_y);
  tree->SetBranchAddress("e_zv",  &etar_z);
  
  //SIMC-SPECIFIC LEAF VARIABLES 
  tree->SetBranchAddress("Normfac",  &Normfac);
  tree->SetBranchAddress("Weight",   &Weight);
  tree->SetBranchAddress("Jacobian_corr", &Jacobian_corr);
  tree->SetBranchAddress("probabs", &prob_abs);

  // ---- SIMC Variables at the vertex (used for calcualting averaged kinematics) ----

  tree->SetBranchAddress("Ein_v", &Ein_v);
  tree->SetBranchAddress("Q2_v", &Q2_v);
  tree->SetBranchAddress("nu_v", &nu_v);
  tree->SetBranchAddress("q_lab_v", &q_v);
  tree->SetBranchAddress("pm_v", &Pm_v);
  tree->SetBranchAddress("pf_v", &Pf_v);
  tree->SetBranchAddress("Ep_v", &Ep_v);
  tree->SetBranchAddress("Ef_v", &Ef_v);
  tree->SetBranchAddress("e_xptar_v", &e_xptar_v);
  tree->SetBranchAddress("e_yptar_v", &e_yptar_v);
  tree->SetBranchAddress("h_xptar_v", &h_xptar_v);
  tree->SetBranchAddress("h_yptar_v", &h_yptar_v);
 
  
  //-----------------------------------
  // DEFINE CUT ON COLLIMATOR GEOMETRY
  //-----------------------------------
  
  //Scaling the HMS/SHMS Collimator Cuts
  hms_hsize = hms_scale*hms_hsize;  //The scale factor is read from set_heep_cuts.inp
  hms_vsize = hms_scale*hms_vsize;
  
  shms_hsize = shms_scale*shms_hsize;
  shms_vsize = shms_scale*shms_vsize;  

  //Define HMS Collimator Shape
  hms_Coll_gCut = new TCutG("hmsCollCut", 8 );
  hms_Coll_gCut->SetVarX("X");
  hms_Coll_gCut->SetVarY("Y");
 
  hms_Coll_gCut->SetPoint(0,  hms_hsize,     hms_vsize/2.);
  hms_Coll_gCut->SetPoint(1,  hms_hsize/2.,  hms_vsize   );
  hms_Coll_gCut->SetPoint(2, -hms_hsize/2.,  hms_vsize   );
  hms_Coll_gCut->SetPoint(3, -hms_hsize,     hms_vsize/2.);
  hms_Coll_gCut->SetPoint(4, -hms_hsize,    -hms_vsize/2.);
  hms_Coll_gCut->SetPoint(5, -hms_hsize/2., -hms_vsize   );
  hms_Coll_gCut->SetPoint(6,  hms_hsize/2., -hms_vsize   );
  hms_Coll_gCut->SetPoint(7,  hms_hsize,    -hms_vsize/2.);
  hms_Coll_gCut->SetPoint(8,  hms_hsize,     hms_vsize/2.);

  //Define SHMS Collimator Shape
  shms_Coll_gCut = new TCutG("shmsCollCut", 8 );
  shms_Coll_gCut->SetVarX("X");
  shms_Coll_gCut->SetVarY("Y");
 
  shms_Coll_gCut->SetPoint(0,  shms_hsize,     shms_vsize/2.);
  shms_Coll_gCut->SetPoint(1,  shms_hsize/2.,  shms_vsize   );
  shms_Coll_gCut->SetPoint(2, -shms_hsize/2.,  shms_vsize   );
  shms_Coll_gCut->SetPoint(3, -shms_hsize,     shms_vsize/2.);
  shms_Coll_gCut->SetPoint(4, -shms_hsize,    -shms_vsize/2.);
  shms_Coll_gCut->SetPoint(5, -shms_hsize/2., -shms_vsize   );
  shms_Coll_gCut->SetPoint(6,  shms_hsize/2., -shms_vsize   );
  shms_Coll_gCut->SetPoint(7,  shms_hsize,    -shms_vsize/2.);
  shms_Coll_gCut->SetPoint(8,  shms_hsize,     shms_vsize/2.);

  //-------------------------------
  //  DEFINE FULL WEIGHT VARIABLE
  //-------------------------------
  // STEP1: Determine the charge factor:
  // definition: total charge deposited on target over a time period
  // SIMC input files are set to 'events / 1mC'
   
  // Charge factor is the total integrated charge assuming a beam current and run time
  Double_t Ib;       //beam current 
  Double_t time;    // estimated beam-on-target time at the kinematics setting
  Double_t charge_factor;

  
  Double_t e_trk;
  Double_t h_trk;
  Double_t daq_lt;   
  Double_t tgt_boil;
  Double_t proton_abs;
  Double_t eff_factor;

    
  Double_t FullWeight;
  Double_t FullWeight_forRates; //this is the full weight for true rate estimates (so no inefficiencies accounted for)  
  Double_t PhaseSpace;

  // define additional scaling for carbon-to-nitrogen and carbon-to-he4 (for background estimates)

  // nuclear transparencies
  // Transparency function: T = c * A ** alpha (Q2), where alpha ~ -0.24 for Q2 >= 2 GeV^2, and c=1, A -> mass number
  // reference: https://arxiv.org/abs/1211.2826  "Color Transparency: past, present and future"
  Double_t alpha=-0.24;
  Double_t T_C12 = pow(12, 
  Double_t T_N   = pow(14, alpha);  
  Double_t T_He4 = pow(4, alpha);  

  

  if( analysis_flag == "d2fsi") {
    Ib = 80;       //beam current in (uA) microAmps (micro-Coulombs / sec),   1 mC = 1000 uC
    time = 168.0;     //estimated time (in hours) a run takes (start - end) of run
    charge_factor = Ib * time * 3600. / 1000.; // in mC

    // efficiencies (assume 1 for now)
    e_trk      = 1.;   //0.964;
    h_trk      = 1.;  //0.988;
    daq_lt     = 1.;  //0.98;
    tgt_boil   = 1.;  // boiling factor correction (assume 1 for now)
    proton_abs = 1.;  // 0.95

    eff_factor = 1; // e_trk * h_trk * daq_lt * tgt_boil * proton_abs;

  }

  if( analysis_flag == "d2pol") {
    Ib = 0.1;         // beam current in (uA) microAmps (micro-Coulombs / sec),   1 mC = 1000 uC (0.1 uA -> 100 nA)
    time = 168.0;     //estimated time (in hours) a run takes (start - end) of run
    charge_factor = Ib * time * 3600. / 1000.; // in mC

    // efficiencies (assume 1 for now)
    e_trk      = 1.;   //0.964;
    h_trk      = 1.;  //0.988;
    daq_lt     = 1.;  //0.98;
    tgt_boil   = 1.;  // boiling factor correction (assume 1 for now)
    proton_abs = 1.;  // 0.95

    eff_factor = 1; // e_trk * h_trk * daq_lt * tgt_boil * proton_abs;



    

    
  }



  
  /* ----------------------------------------------------------
     
     NOTE on FullWeight and Phase Space definitions : 
     
     FullWeight = Weight * Normfac / n_acc   (assuming 1 mC of charge and detector efficiencies are 100 % (in fraction is 1.0 ) )
     
     1a) main%sig = main%SF_weight*main%jacobian*main%sigcc
	   main%weight = main%SF_weight*main%jacobian*main%jacobian_corr*main%gen_weight*main%sigcc*main%coul_corr
                       = main%sig * (main%jacobian_corr*main%gen_weight*main%coul_corr)

     1b)  Weight = (sigma_theory*Jacobian) * main%gen_weight * jacobian_corr * coulomb_corr (see event.f) for 'LAGET_DEUT'
          --> The Jacobian is explicitly included in file LagetXsec.f of SIMC
     main%gen_weight = main%gen_weight*(Emax-Emin)/(gen%e%E%max-gen%e%E%min) ----> using the definitions below:  main%gen_weight * (min(gen%sumEgen%max,gen%sumEgen%max) - max(gen%sumEgen%min,gen%sumEgen%min)) / (gen%sumEgen%max - gen%sumEgen%min)

     Emax = =gen%e%E%max = gen%sumEgen%max,  Emin = gen%e%E%min = gen%sumEgen%min  (see init.f and event.f)
     Emin = max(Emin,gen%sumEgen%min)
     Emax = min(Emax,gen%sumEgen%max)
     
     It seems main%gen_weight should be 1 for the deuteron, (without radiative effects), so Phase Space must ONLY be calculated WITHOUT RADIATIVE EFFECTS ! ! !
     
                                                                          |---------- dE'---------|
     2) Normfac / n_acc =  (luminosity / n_tried) * domega_e * domega_p * (gen%e%E%max-gen%e%E%min)
     
     where n_acc -> nentries (in this code), and n_tried is total number of MC events generated (thrown into a volume)
     
   
     Combining definirions from  1) and 2):   

                  |---------------------------------------------------------------------------------- main%weight ----------------------------------------------------------------------------|
                  |------- main%sig -------|  |------------------------------------------- main%gen_weight -------------------------------------------------|                                    |------------------------------'Normfac / n_acc'----------------------------|
     FullWeight = (sigma_theory * Jacobian) * ( (min(gen%sumEgen%max,gen%sumEgen%max) - max(gen%sumEgen%min,gen%sumEgen%min)) / (gen%e%E%max - gen%e%E%min) )  * jacobian_corr * coulomb_corr  * ( (luminosity / n_tried) * domega_e * domega_p * ((gen%e%E%max-gen%e%E%min)) )
     
     
     PhaseSpace = FullWeight / (sigma_theory*Jacobian) =  Genweight * Jacobian_corr * coulomb_corr  * ( (luminosity / n_tried) * domega_e * domega_p * (dE') )
     **NOTE: Phase space MUST only be calculated in SIMC without raditative effects: Genweight = 1, Jacobian_corr ~ 1, and coulomb_corr = 1
     
     For estimates of yields (i.e., beam time to be allocated), one only needs to weight the events by the FullWeight, assuming 1 hr beam time, and then the yields can be scaled by the time. 
     Also, remember that the actual statistical error on the yield, say in Pm bins,  is ~ 1/sqrt(N_counts) per Pm kinematic bin,  where N_counts is the number of weighted counts

     ----------------------------------------------------------
  */
  
  //---------------------
  //  LOOP OVER ENTRIES
  //---------------------
  
  for (Long64_t i=0; i < nentries; i++) {

    //Get the ith entry from the SNT TTree
    tree->GetEntry(i);

    //convert to GeV
    Ei = Ei/1000.;  // initial e- energy
    kf = kf/1000.;  // final e- momentum
    ki = sqrt(Ei*Ei - me*me);   // initial e- momentum

    Pf = Pf/1000.;  //final proton momentum (GeV/c)

    //Calculate electron final momentum 3-vector
    SetCentralAngles(the_central, phe_central);
    TransportToLab(kf, e_xptar, e_yptar, kf_vec);
    
    fP0.SetXYZM( 0.0, 0.0, ki, me );          // e- beam 4-momentum (assumed to be along +Z in LAB)
    fP1.SetXYZM( kf_vec.X(),  kf_vec.Y(),  kf_vec.Z(), me );             // e- scattered 4-momentum
    fMp.SetXYZM( 0.0, 0.0, 0.0, MP );       // proton mass (for x_bj calculation)
    fA.SetXYZM( 0.0, 0.0, 0.0, tgt_mass );  // initial target 4-momentum (target at rest assumed)    
    fQ         = fP0 - fP1;    // four-momentum transfer 
    fA1        = fA + fQ;      // final system 4-momentum (detected hadron + recoil)
    fMp1       = fMp + fQ;
    Q2        = -fQ.M2();  // 4-momentum transferred squared
    q         = fQ.P();    // 3-momentum trasnfer |q|
    nu        = fQ.E();    // energy transfer, nu
    W2        = fMp1.M2();  // INVARIANT MASS squared
    if (W2>0) { W = TMath::Sqrt(W2); } // INVARIANT MASS GeV
    the       = fP0.Angle( fP1.Vect() ) / dtr ;   // e- scattering angle [deg]
    th_q      = fQ.Theta() / dtr;  // in-plane angle of q-vector relative to +z (along beam) [deg]
    ph_q      = fQ.Phi() / dtr;    // out-of-plane angle of q-vector relative to +z (along beam)[deg]
    X         = Q2 / (2.*MP*nu);   //x-bjorken

   
    // calculate final proton momentum vector
    // VERY IMPORTANT: IN HCANA, central angle for HMS is set to negative by default,
    // (therefore, make sure to use the same condition in SIMC)
    SetCentralAngles(-thp_central, php_central); 
    TransportToLab(Pf, h_xptar, h_yptar, Pf_vec);

    // four-momentum of the detected (X) particle (for A(e,e'X), our X = p (proton))
    fX.SetXYZM(Pf_vec.X(), Pf_vec.Y(), Pf_vec.Z(), MP);   

    // four-momentum of the undetected, recoil system, B
    fB = fA1 - fX;
    
    // opening Angle of X with scattered primary (e-) particle
    fXangle = fX.Angle( fP1.Vect()) / dtr; //[deg]
    thp =  fXangle - the; // proton scattering angle [deg]
	
    // missing momentum components in LAB frame [GeV]
    Pmx_lab = fB.X();
    Pmy_lab = fB.Y(); 
    Pmz_lab = fB.Z();

    // calculate missing momentum magnitude
    Pm =  sqrt(Pmx_lab*Pmx_lab + Pmy_lab*Pmy_lab + Pmz_lab*Pmz_lab);
    
    //--------Rotate the recoil system from +z to +q-------
    rot_to_q.SetZAxis( fQ.Vect(), fP1.Vect() ).Invert(); // rotate 

    xq = fX.Vect();
    bq = fB.Vect();
    
    xq *= rot_to_q;
    bq *= rot_to_q;

    //Calculate Angles of q relative to x(detected proton) and b(recoil neutron)
    // sense of roation for phi is: +/- 180 deg
    th_pq = xq.Theta() / dtr;   //"thpq" [ deg ]                                       
    ph_pq = xq.Phi()   / dtr;     //"out-of-plane angle", "phi_pq"   [deg]                                                                 
    th_rq = bq.Theta() / dtr;   // theta_rq  [deg]                                                                                                    
    ph_rq = bq.Phi()   / dtr;     //"out-of-plane angle", phi_rq [deg]

    // convert [0, 360] to [-180, 180] deg
    if(ph_pq >= 180.){
      ph_pq = ph_pq - 360.;
    }
	
    p_miss_q = -bq;  // missing momentum vector in q-frame

    //Missing Momentum Components in the q-frame [GeV]
    Pmz_q = p_miss_q.Z();   //parallel component to +z (+z is along q)
    Pmx_q = p_miss_q.X();   //in-plane perpendicular component to +z
    Pmy_q = p_miss_q.Y();   //out-of-plane component (Oop)

    MM = fB.M(); // INVARIANT MASS of recoil system (aka missing mass) GeV
    MM2 = MM*MM;

    // Kinetic energies of detected (X) and recoil (B) in GeV
    Tx = fX.E() - MP;  // need to figure out why does this NOT work
    Tr = fB.E() - MM;
    
    // Standard nuclear physics definition of "missing energy":
    // binding energy of X in the target (= removal energy of X).
    // NB: If X is knocked out of a lower shell, the recoil system carries
    // a significant excitation energy. This excitation is included in Emiss
    // here, as it should, since it results from the binding of X.
    Em_nuc = nu - Tx - Tr;
    
    // particle physics missing energy: Target Mass + Beam Energy - scatt ele energy - hadron total energy
    Em = nu + fA.M() - fX.E();
    
    // In production reactions, the "missing energy" is defined 
    // as the total energy of the undetected recoil system.
    // This is the "missing mass", Mrecoil, plus any kinetic energy. (this is different than E_nuc)
    Erecoil = fB.E();  // final recoil energy
    Ep = fX.E();       // final proton energy
      
    // define difference between HMS/SHMS reconstructed z
    ztar_diff =  htar_z - etar_z;

    //SIMC Collimator (definition based on HCANA collimator)
    htarx_corr = tar_x - h_xptar*htar_z*cos(thp_central*dtr);
    etarx_corr = tar_x - e_xptar*etar_z*cos(the_central*dtr);  
    
    
    //Define Collimator (same as in HCANA)
    hXColl = htarx_corr + h_xptar*168.;   //in cm
    hYColl = h_ytar + h_yptar*168.;
    eXColl = etarx_corr + e_xptar*253.;
    eYColl = e_ytar + e_yptar*253.-(0.019+40.*.01*0.052)*e_delta+(0.00019+40*.01*.00052)*e_delta*e_delta; //correct for HB horizontal bend	  
    
    
    
    //=====================================================================================
    
    //---------Calculate Necessary Vertex Quantities for Average Kinematics----------------
    
    //Convert from MeV to GeV
    Ein_v = Ein_v / 1000.;
    Ef_v = Ef_v / 1000.;
    kf_v = sqrt(Ef_v*Ef_v - me*me);    //final electron momentum at vertex
    ki_v = sqrt(Ein_v*Ein_v - me*me);   //initial electron momentum at vertex

    Pf_v = Pf_v / 1000.; // final proton momentum 

    //Calculate electron final momentum 3-vector
    SetCentralAngles(the_central, phe_central);
    TransportToLab(kf_v, e_xptar_v, e_yptar_v, kf_vec_v);
    
    //Calculate 4-Vectors at the vertex    
    fP0_v.SetXYZM(0.0, 0.0, ki_v, me); 
    fP1_v.SetXYZM( kf_vec_v.X(), kf_vec_v.Y(), kf_vec_v.Z(), me);
    fA_v.SetXYZM(0.0, 0.0, 0.0, tgt_mass );
    fMp_v.SetXYZM( 0.0, 0.0, 0.0, MP );       // proton mass (for x_bj calculation)
    fQ_v = fP0_v - fP1_v;
    fA1_v = fA_v + fQ_v;  
    fMp1_v       = fMp_v + fQ_v;

    Q2_v        = -fQ_v.M2();  // 4-momentum transferred squared
    q_v         = fQ_v.P();    // 3-momentum trasnfer |q|
    nu_v        = fQ_v.E();    // energy transfer, nu
    W2_v        = fMp1_v.M2();  // INVARIANT MASS squared
    if (W2_v>0) { W_v = TMath::Sqrt(W2_v); } // INVARIANT MASS GeV
    the_v       = fP0_v.Angle( fP1_v.Vect() ) / dtr ;   // e- scattering angle [deg]
    thq_v      = fQ_v.Theta() / dtr;  // in-plane angle of q-vector relative to +z (along beam) [deg]
    phq_v      = fQ_v.Phi() / dtr;    // out-of-plane angle of q-vector relative to +z (along beam)[deg]
    X_v       = Q2_v / (2.*MP*nu_v);   //x-bjorken

     // calculate final proton momentum vector 
    // VERY IMPORTANT: IN HCANA, central angle for HMS is set to negative by default,
    // (therefore, make sure to use the same condition in SIMC)
    SetCentralAngles(-thp_central, php_central);
    TransportToLab(Pf_v, h_xptar_v, h_yptar_v, Pf_vec_v);

    // four-momentum of the detected (X) particle (for A(e,e'X), our X = p (proton))
    fX_v.SetXYZM(Pf_vec_v.X(), Pf_vec_v.Y(), Pf_vec_v.Z(), MP);

    // four-momentum of the undetected, recoil system, B
    fB_v = fA1_v - fX_v;                 //4-MOMENTUM OF UNDETECTED PARTICLE 

    // opening Angle of X with scattered primary (e-) particle
    fXangle_v = fX_v.Angle( fP1_v.Vect()) / dtr; //[deg]
    thp_v =  fXangle_v - the_v; // proton scattering angle [deg]

    // missing momentum components in LAB frame [GeV]
    Pmx_lab_v = fB_v.X();
    Pmy_lab_v = fB_v.Y(); 
    Pmz_lab_v = fB_v.Z();

    // calculate missing momentum magnitude
    Pm_v =  sqrt(Pmx_lab_v*Pmx_lab_v + Pmy_lab_v*Pmy_lab_v + Pmz_lab_v*Pmz_lab_v);
    
    //--------Rotate the recoil system from +z to +q-------
    rot_to_q_v.SetZAxis( fQ_v.Vect(), fP1_v.Vect() ).Invert();

    xq_v = fX_v.Vect();
    bq_v = fB_v.Vect();
    
    xq_v *= rot_to_q_v;
    bq_v *= rot_to_q_v;

    //Calculate Angles of q relative to x(detected proton) and b(recoil neutron)
    // sense of roation for phi is: +/- 180 deg
    th_pq_v = xq_v.Theta()  / dtr;   //"thpq"                                       
    ph_pq_v = xq_v.Phi()    / dtr;     //"out-of-plane angle", "phi_pq"                                                                    
    th_rq_v = bq_v.Theta()  / dtr;   // theta_rq                                                                                                     
    ph_rq_v = bq_v.Phi()    / dtr;     //"out-of-plane angle", phi_rq

    // convert [0, 360] to [-180, 180] deg
    if(ph_pq_v >= 180.){
      ph_pq_v = ph_pq_v - 360.;
    }

    p_miss_q_v = -bq_v;  // missing momentum vector in q-frame

    //Missing Momentum Components in the q-frame [GeV]
    Pmz_q_v = p_miss_q_v.Z();   //parallel component to +z (+z is along q)
    Pmx_q_v = p_miss_q_v.X();   //in-plane perpendicular component to +z
    Pmy_q_v = p_miss_q_v.Y();   //out-of-plane component (Oop)

    MM_v = fB_v.M(); // INVARIANT MASS of recoil system (aka missing mass) GeV
    MM2_v = MM_v*MM_v;

    // Kinetic energies of detected (X) and recoil (B) in GeV
    Tx_v = fX_v.E() - MP;  // need to figure out why does this NOT work
    Tr_v = fB_v.E() - MM;


    Em_nuc_v = nu_v - Tx_v - Tr_v;
    
    // particle physics missing energy: Target Mass + Beam Energy - scatt ele energy - hadron total energy
    Em_v = nu_v + fA_v.M() - fX_v.E();
    
    // In production reactions, the "missing energy" is defined 
    // as the total energy of the undetected recoil system.
    // This is the "missing mass", Mrecoil, plus any kinetic energy. (this is different than E_nuc)
    Erecoil_v = fB.E();  // final recoil energy
    Ep_v = fX.E();       // final proton energy
      
    
    /* --- SNIPPET from event.f in SIMC with coordinate definitions ---
   ! CALCULATE ANGLE PHI BETWEEN SCATTERING PLANE AND REACTION PLANE.
   ! Therefore, define a new system with the z axis parallel to q, and
   ! the x axis inside the q-z-plane: z' = q, y' = q X z, x' = y' X q
   ! this gives phi the way it is usually defined, i.e. phi=0 for in-plane
   ! particles closer to the downstream beamline than q.
   ! phi=90 is above the horizontal plane when q points to the right, and
   ! below the horizontal plane when q points to the left.
   ! Also take into account the different definitions of x, y and z in
   ! replay and SIMC:
   ! As seen looking downstream:        replay	SIMC	(old_simc)
   !				x	right	down	(left)  ---> this must not be correct for hallc replay: as lab coord:  +X -> beam left, +Y -> beam up, +z -> downstream
   !				y	down	left	(up)
   !				z	all have z pointing downstream
   !
   ! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc
    */

    //---------Light Cone Variables (at vertex) (C.Y. March 05, 2021)---------
    
    /*
      NOTE: In analogy to binning the exp. cross section in (Pm, th_rq) 2D histo bins,
      one can also bin the exp. cross sections in light-cone (alpha, pt) 2D histo bins, which can
      be later used to extact the light-cone momentum distributions. See https://arxiv.org/abs/1501.05377
      alpha -> light-cone momentum fraction of the struck nucleon, pt -> transverse momentum component of 
      the struck nucleon which is equal to the -pt_n (transverse momentum component of the recoil system 
      (neutron for D(e,e'p)n case)
    */

    //recoil ("missing neutron") momentum component along q-vector
    PmPar_v = Pm_v * cos(th_rq_v * dtr);
    
    //recoil momentum component transverse to q-vector
    PmPerp_v = sqrt(pow(Pm_v,2) - pow(PmPar_v,2));
    
    //light-cone momentum fraction of the recoil neutron
    alpha_n_v = (Erecoil_v - PmPar_v) / (MD/2.);
    
    //momentum fraction of struck nucleon (normalized such that: alpha + alpha_n = 2)
    alpha_v = 2. - alpha_n_v;
    
    
    
    //--------------------------------------------------------------------------------
    


    
    
    //------Define ANALYSIS CUTS-------

    //---Collimator CUTS---
    
    if(hmsCollCut_flag)  { hmsColl_Cut =  hms_Coll_gCut->IsInside(hYColl, hXColl);}
    else{hmsColl_Cut=1;}
    
    if(shmsCollCut_flag) { shmsColl_Cut =  shms_Coll_gCut->IsInside(eYColl, eXColl);}
    else{shmsColl_Cut=1;}
	  
    //----Kinematics Cuts----

    //Q2
    if(Q2_cut_flag){c_Q2 = Q2>=c_Q2_min && Q2<=c_Q2_max;}
    else{c_Q2=1;}
    
    //Nuclear Missing Energy, Em
    if(Em_cut_flag){c_Em = Em_nuc>=c_Em_min && Em_nuc<=c_Em_max;}
    else{c_Em=1;}

 
    //----Acceptance Cuts----
    if(hdelta_cut_flag){c_hdelta = h_delta>=c_hdelta_min && h_delta<=c_hdelta_max;} 
    else{c_hdelta=1;}
    
    if(edelta_cut_flag){c_edelta = e_delta>=c_edelta_min && e_delta<=c_edelta_max;} 
    else{c_edelta=1;} 
    
    if(ztarDiff_cut_flag){c_ztarDiff = ztar_diff>=c_ztarDiff_min && ztar_diff<=c_ztarDiff_max;} 
    else{c_ztarDiff=1;} 

    //Combine All CUTS
    c_accpCuts = c_hdelta && c_edelta && c_ztarDiff && hmsColl_Cut && shmsColl_Cut;
    c_kinCuts = c_Q2 && c_Em;
    c_allCuts =  c_accpCuts && c_kinCuts;

    //Full Weight
    FullWeight = (Normfac * charge_factor * eff_factor * Weight ) / nentries;
    FullWeight_forRates = (Normfac * charge_factor * Weight ) / nentries;

    PhaseSpace =  Normfac * charge_factor * eff_factor * Jacobian_corr  / nentries;    //Phase Space with jacobian corr. factor


    //fill histogram (only meaningful if no cuts, and eff_factor set to 1 for daq rates)
    H_W_noCut->Fill(W, FullWeight_forRates); // for coin daq rates estimates
    H_Pm_noCut->Fill(Pm, FullWeight_forRates);

    
    if(c_allCuts) {

      // ------ This is for calculation of avergaed kinematics ---
      // will just use the vertex  kinematics -------
      H_Pm_vs_thrq_v    ->Fill(th_rq_v, Pm_v, FullWeight);	 
      H_Ein_2Davg       ->Fill(th_rq_v, Pm_v, Ein_v*FullWeight);
      H_kf_2Davg        ->Fill(th_rq_v, Pm_v, kf_v*FullWeight);
      H_the_2Davg       ->Fill(th_rq_v, Pm_v, (the_v)*FullWeight);
      H_thp_2Davg       ->Fill(th_rq_v, Pm_v, (thp_v)*FullWeight);
      H_Pf_2Davg        ->Fill(th_rq_v, Pm_v, Pf_v*FullWeight);
      H_Pm_2Davg        ->Fill(th_rq_v, Pm_v, Pm_v*FullWeight);
      H_thrq_2Davg      ->Fill(th_rq_v, Pm_v, (th_rq_v)*FullWeight);

      H_q_2Davg         ->Fill(th_rq_v, Pm_v, q_v*FullWeight);
      H_theta_q_2Davg   ->Fill(th_rq_v, Pm_v, (thq_v)*FullWeight);
      H_Q2_2Davg        ->Fill(th_rq_v, Pm_v, Q2_v*FullWeight);
      H_nu_2Davg        ->Fill(th_rq_v, Pm_v, nu_v*FullWeight);
      H_xbj_2Davg       ->Fill(th_rq_v, Pm_v, X_v*FullWeight);
      H_theta_pq_2Davg  ->Fill(th_rq_v, Pm_v, (th_pq_v)*FullWeight);
      H_phi_pq_2Davg    ->Fill(th_rq_v, Pm_v, (ph_pq_v)*FullWeight);
      H_cphi_pq_2Davg   ->Fill(th_rq_v, Pm_v, cos(ph_pq_v*dtr)*FullWeight);
      H_sphi_pq_2Davg   ->Fill(th_rq_v, Pm_v, sin(ph_pq_v*dtr)*FullWeight);

      //once they are filled, then divide H_[]_2Davg / H_Pm_vs_thrq_v
	
      //--------------------------------------------------------------
      
      // This is for the 2D cross section Pm vs. thrq binned in thrq 
      H_Pm_vs_thrq->Fill(th_rq, Pm, FullWeight);
      H_Pm_vs_thrq_ps->Fill(th_rq, Pm, PhaseSpace);
      
      //Primary (electron) Kinematics
      H_kf  ->Fill(kf,   FullWeight);
      H_the ->Fill(the,  FullWeight);
      H_Q2  ->Fill(Q2,   FullWeight);
      H_xbj ->Fill(X,    FullWeight);
      H_nu  ->Fill(nu,   FullWeight);
      H_q   ->Fill(q,    FullWeight);
      H_thq ->Fill(th_q, FullWeight);
      H_W   ->Fill(W,    FullWeight);

 
      //Secondary (Hadron) Kinematics
      H_Pf    ->Fill(Pf, FullWeight);
      H_thp   ->Fill(thp, FullWeight);
      H_Em    ->Fill(Em, FullWeight);
      H_Em_nuc->Fill(Em_nuc, FullWeight);
      H_Pm    ->Fill(Pm, FullWeight);

      H_thpq->Fill(th_pq, FullWeight);
	    
      H_thrq->Fill(th_rq, FullWeight);
      H_phi_pq->Fill(ph_pq, FullWeight);
      H_cphi_pq->Fill(cos(ph_pq * dtr), FullWeight);

      H_MM->Fill(MM, FullWeight);
      H_MM2->Fill(MM2, FullWeight);


      // fill vertex quantities (for checks)
      H_Em_nuc_v     ->Fill(Em_nuc_v, FullWeight);    
      H_Pm_v     ->Fill(Pm_v, FullWeight);    
      H_Q2_v     ->Fill(Q2_v, FullWeight);   
      H_the_v    ->Fill(the_v, FullWeight);   
      H_thpq_v   ->Fill(th_pq_v, FullWeight);  
      H_thrq_v   ->Fill(th_rq_v, FullWeight); 
      H_phi_pq_v ->Fill(ph_pq_v, FullWeight);
      H_cphi_pq_v ->Fill(cos(ph_pq_v*dtr), FullWeight);

      
      //Target Reconstruction (Hall Coord. System)
      H_htar_x->Fill(tar_x, FullWeight);
      H_htar_y->Fill(htar_y, FullWeight);
      H_htar_z->Fill(htar_z, FullWeight);
      H_etar_x->Fill(tar_x, FullWeight);
      H_etar_y->Fill(etar_y, FullWeight);
      H_etar_z->Fill(etar_z, FullWeight);
      H_ztar_diff->Fill(ztar_diff, FullWeight);
      
      //Hadron arm Reconstructed Quantities ( xtar, ytar, xptar, yptar, delta) 
      H_hytar->Fill(h_ytar, FullWeight);
      H_hxptar->Fill(h_xptar, FullWeight);
      H_hyptar->Fill(h_yptar, FullWeight);
      H_hdelta->Fill(h_delta, FullWeight);
      
      //Hadron arm Focal Plane Quantities
      H_hxfp->Fill(h_xfp, FullWeight);
      H_hyfp->Fill(h_yfp, FullWeight);
      H_hxpfp->Fill(h_xpfp, FullWeight);
      H_hypfp->Fill(h_ypfp, FullWeight);
      
      //Electron Arm Reconstructed Quantities ( xtar, ytar, xptar, yptar, delta)
      H_eytar->Fill(e_ytar, FullWeight);
      H_exptar->Fill(e_xptar, FullWeight);
      H_eyptar->Fill(e_yptar, FullWeight);
      H_edelta->Fill(e_delta, FullWeight);
      
      //Electron Arm Focal Plane Quantities
      H_exfp->Fill(e_xfp, FullWeight);
      H_eyfp->Fill(e_yfp, FullWeight);
      H_expfp->Fill(e_xpfp, FullWeight);
      H_eypfp->Fill(e_ypfp, FullWeight);

      //---- Fill 2D Histos ----
      
      // Xfp vs Yfp
      H_hxfp_vs_hyfp->Fill(h_yfp, h_xfp, FullWeight);
      H_exfp_vs_eyfp->Fill(e_yfp, e_xfp, FullWeight);

      // 2D Collimator Histos
      H_hXColl_vs_hYColl->Fill(hYColl, hXColl, FullWeight);
      H_eXColl_vs_eYColl->Fill(eYColl, eXColl, FullWeight);
	    
      // 2D HMS v. SHMS Acceptance Correlations
      H_hxptar_vs_exptar->Fill(e_xptar, h_xptar, FullWeight); 
      H_hyptar_vs_eyptar->Fill(e_yptar, h_yptar, FullWeight); 
      H_hdelta_vs_edelta->Fill(e_delta, h_delta, FullWeight);

						 
    }

    //No Self Cut Histos
    if(c_accpCuts &&  c_Q2) { H_Em_nuc_nsc->Fill(Em_nuc, FullWeight); }
    if(c_accpCuts &&  c_Em) { H_Q2_nsc->Fill(Q2, FullWeight); }

    if(c_kinCuts && c_hdelta && c_edelta && hmsColl_Cut && shmsColl_Cut) { H_ztar_diff_nsc->Fill(ztar_diff, FullWeight); }

    if(c_kinCuts && c_edelta && c_ztarDiff && hmsColl_Cut && shmsColl_Cut) {H_hdelta_nsc->Fill(h_delta, FullWeight); }
    if(c_kinCuts && c_hdelta && c_ztarDiff && hmsColl_Cut && shmsColl_Cut) {H_edelta_nsc->Fill(e_delta, FullWeight); }

    if(c_kinCuts && c_hdelta && c_edelta && c_ztarDiff && shmsColl_Cut) { H_hXColl_vs_hYColl_nsc->Fill(hYColl, hXColl, FullWeight); }
    if(c_kinCuts && c_hdelta && c_edelta && c_ztarDiff && hmsColl_Cut) { H_eXColl_vs_eYColl_nsc->Fill(eYColl, eXColl, FullWeight); }

    cout << "SIMC Events Completed: " << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";
    
  } // end entry loop

  //Finish Calculating the 2D Average Kinematics (Divide by the sum of the weight)
  H_Ein_2Davg        ->Divide(H_Pm_vs_thrq_v);
  H_kf_2Davg         ->Divide(H_Pm_vs_thrq_v);
  H_the_2Davg        ->Divide(H_Pm_vs_thrq_v);
  H_thp_2Davg        ->Divide(H_Pm_vs_thrq_v);
  H_Pf_2Davg         ->Divide(H_Pm_vs_thrq_v);
  H_Pm_2Davg         ->Divide(H_Pm_vs_thrq_v);
  H_thrq_2Davg       ->Divide(H_Pm_vs_thrq_v);

  H_q_2Davg          ->Divide(H_Pm_vs_thrq_v);
  H_theta_q_2Davg    ->Divide(H_Pm_vs_thrq_v);
  H_Q2_2Davg         ->Divide(H_Pm_vs_thrq_v);
  H_nu_2Davg         ->Divide(H_Pm_vs_thrq_v);
  H_xbj_2Davg        ->Divide(H_Pm_vs_thrq_v);
  H_theta_pq_2Davg   ->Divide(H_Pm_vs_thrq_v);
  H_phi_pq_2Davg     ->Divide(H_Pm_vs_thrq_v);
  H_cphi_pq_2Davg    ->Divide(H_Pm_vs_thrq_v);
  H_sphi_pq_2Davg    ->Divide(H_Pm_vs_thrq_v);

  
  //Calculate Xsec
  H_Pm_vs_thrq_xsec->Divide(H_Pm_vs_thrq, H_Pm_vs_thrq_ps);
  
  //-------------------------
  // WRITE HISTOS TO FILE
  //-------------------------
  
  //Create Output ROOTfile
  outROOT = new TFile(simc_OutputFileName, "RECREATE");
  
  //Make directories to store histograms based on category
  outROOT->mkdir("kin_plots");
  outROOT->mkdir("accp_plots");
  
  //Write Kinematics histos to kin_plots directory
  outROOT->cd("kin_plots");
  kin_HList->Write();
  
  //Write Acceptance histos to accp_plots directory
  outROOT->cd("accp_plots");
  accp_HList->Write();

  // this is needed so thay it may be used by calc_avg_kin.py
  outROOT->cd("");
  kin_HList->Write();
  
  //Close File
  outROOT->Close();
  
  //** IMPORTANT** Consideration of statistical uncertainty based on counts
  /*
    The uncertainty calculated per kinematic bin in SIMC is representative of the
    number of events simulated, and is therefore not OK to use it as an estimate
    of the experimental stat. uncertainty. In other words, the greater the number
    of Monte-Carlo events, the smaller the error bar per kinematic bin. 

    To get a realistic estimate of the statistical uncertainty per kinematic bin,
    assuming the bin (say in missing momentum, Pm) has been properly weighted (See FullWeight above)
    one has to use the standard error for counting N independent variables: absolute error of bin = sqrt(N),
    and relative error = sqrt(N) / N = 1 / sqrt(N)
    
  */
  
  

  //-----------------------------------------------------------------
  //
  // ------ Write output file for storing rates, counts, etc. -------
  //
  //-----------------------------------------------------------------

  //FileStreams objects to READ/WRITE to a .txt file
  ofstream out_file;
  ifstream in_file;
  

  // variable integrations for determining total counts, (e,e'p) rates, and daq rats
  float Pmcnts = H_Pm->Integral();
  float Pmcnts_noCut = H_Pm_noCut->Integral();     
  float rates = H_Pm->Integral() / (time * 3600);
  float daq_rates = H_Pm_noCut->Integral() / (time * 3600);
  

  //  deuteron FSI studies 
  if( analysis_flag == "d2fsi") {
  
  cout << " ----------------------------------------- " << endl; 
  cout << "    d(e,e'p) FSI Studies Rate Estimates    " << endl;
  cout << " ----------------------------------------- " << endl;  
  cout << "" << endl;  
  cout << Form("Pm Setting: %d", pm_set) << endl;
  cout << Form("thrq Setting: %d", thrq_set) << endl;    
  cout << Form("Model: Laget %s", model.Data()) << endl;
  cout << Form("Ib [uA]     = %.3f ", Ib) << endl;
  cout << Form("time [hr]   = %.3f ", time) << endl;
  cout << Form("charge [mC] = %.3f ", charge_factor)<< endl;  
  cout << Form("Pm counts  = %.3f", Pmcnts) << endl;  
  cout << Form("d(e,e'p) Rates [Hz] = %.3E ", rates) << endl;  
  cout << Form("DAQ Rates [Hz] = %.3f", daq_rates) << endl;
  cout << " ----------------------------- " << endl; 

  in_file.open(output_file.Data());

  if(!in_file.fail()){
    
    cout << Form(" %s already exists, will append to it . . . ", output_file.Data() ) << endl;
    
    //Open Report FIle in append mode
    out_file.open(output_file, ios::out | ios::app);
    out_file << Form("%i,     %d,    %s,    %.1f,       %.3E,        %.3E,        %.3f,      %.3f,     %.3f ", pm_set, thrq_set, model.Data(), Pmcnts, rates, daq_rates, Ib, time, charge_factor ) << endl;

    
  }
  
  // create output file if it does not exist
  if(in_file.fail()){
    
    out_file.open(output_file);
    //set headers
    out_file << "# SIMC rates estimates (deut FSI studies)" << endl;
    out_file << "# " << endl;
    out_file << "# header definitions " << endl;
    out_file << "# pm_set: central missing momentum setting [MeV] " << endl;
    out_file << "# thrq_set: central recoil angle setting [deg] " << endl;
    out_file << "# model: Laget pwia or fsi (paris NN potential)" << endl;      
    out_file << "# pm_counts: integrated missing momentum counts (yield) with all cuts applied \n# (not binned in any particular kinematics)" << endl;
    out_file << "# deep_rates: deuteron break-up rates [Hz] " << endl;
    out_file << "# daq_rates: data acquisition rates (no analysis cuts) [Hz] " << endl;
    out_file << "# current: beam current [uA] " << endl;
    out_file << "# time: beam-on-target time [hrs] " << endl;
    out_file << "# charge: beam charge [mC] " << endl;
    out_file << "#" << endl;
    
    out_file <<"pm_set,thrq_set,model,pm_counts,deep_rates,daq_rates,current,time,charge" << endl;
    out_file << Form("%i,     %d,    %s,    %.1f,       %.3E,        %.3E,        %.3f,      %.3f,     %.3f ", pm_set, thrq_set, model.Data(), Pmcnts, rates, daq_rates, Ib, time, charge_factor ) << endl;
    
  }
  
  }
  
  //  polarized deuteron  studies 
  if( analysis_flag == "d2pol") {
    
    cout << " --------------------------------------- " << endl; 
    cout << "    d(e,e'p) Polarized Rate Estimates    " << endl;
    cout << " --------------------------------------- " << endl;  
    cout << "" << endl;  
    cout << Form("Pm Setting: %d", pm_set) << endl;
    cout << Form("Q2 Setting: %.1f", Q2_set) << endl;    
    cout << Form("Model: Laget %s", model.Data()) << endl;
    cout << Form("Ib [nA]     = %.3f ", Ib*1000) << endl;
    cout << Form("time [hr]   = %.3f ", time) << endl;
    cout << Form("charge [mC] = %.3f ", charge_factor)<< endl;  
    cout << Form("Pm counts  = %.3f", Pmcnts) << endl;  
    cout << Form("d(e,e'p) Rates [Hz] = %.3E ", rates) << endl;  
    cout << Form("DAQ Rates [Hz] = %.3f", daq_rates) << endl;
    cout << " ----------------------------- " << endl; 
    
    
    in_file.open(output_file.Data());

    if(!in_file.fail()){
      
      cout << Form(" %s already exists, will append to it . . . ", output_file.Data() ) << endl;
      
      //Open Report FIle in append mode
      out_file.open(output_file, ios::out | ios::app);
      out_file << Form("%i,     %.2f,    %s,     %.1f,       %.3E,        %.3E,        %.3f,      %.3f,     %.3f ", pm_set, Q2_set, model.Data(), Pmcnts, rates, daq_rates, Ib*1000, time, charge_factor ) << endl;

      
    }

    // create output file if it does not exist
    if(in_file.fail()){

      out_file.open(output_file);
      //set headers
      out_file << "# SIMC rates estimates (polarized deut studies)" << endl;
      out_file << "# " << endl;
      out_file << "# header definitions " << endl;
      out_file << "# pm_set: central missing momentum setting [MeV] " << endl;
      out_file << "# Q2_set: central 4-momentum transfer setting [GeV^2] " << endl;
      out_file << "# model: Laget pwia or fsi (paris NN potential)" << endl;      
      out_file << "# pm_counts: integrated missing momentum counts (yield) with all cuts applied \n# (not binned in any particular kinematics)" << endl;
      out_file << "# deep_rates: deuteron (e,e'p) rates [Hz] " << endl;
      out_file << "# daq_rates: data acquisition rates (no analysis cuts) [Hz] " << endl;
      out_file << "# current: beam current [nA] " << endl;
      out_file << "# time: beam-on-target time [hrs] " << endl;
      out_file << "# charge: beam charge [mC] " << endl;
      out_file << "#" << endl;
      
      out_file <<"pm_set,Q2_set,model,pm_counts,deep_rates,daq_rates,current,time,charge" << endl;
      out_file << Form("%i,     %.2f,    %s,     %.1f,       %.3E,        %.3E,        %.3f,      %.3f,     %.3f ", pm_set, Q2_set, model.Data(), Pmcnts, rates, daq_rates, Ib*1000, time, charge_factor ) << endl;

    }

  }
  

  
    // --------------------------------------------------------
    // Write Histogram to numerical data file for plotting
    // --------------------------------------------------------
    
    TString cmd = "";



     if( analysis_flag == "d2fsi") {
       
       cmd = Form("mkdir -p %s", output_hist_data.Data() );
       gSystem->Exec(cmd); // create  histo dir. if it doesn't exist (it will automatically check if it exists, otherwise, will create it_
     }

     

     if( analysis_flag == "d2pol") {
       cmd = Form("mkdir -p %s", output_hist_data.Data() );
       gSystem->Exec(cmd); // create  histo dir. if it doesn't exist (it will automatically check if it exists, otherwise, will create it_
     }

    
    TString class_name;
    TH1F *h_i = 0;
    TH2F *h2_i = 0;

    //-----------------------------------------------------
    //Lopp over kin_HList of histogram objects 
    //-----------------------------------------------------
    string xlabel;
    string ylabel;
    string title;
    for(int i=0; i<kin_HList->GetEntries(); i++) {
    
      //Get the class name for each element on the list (either "TH1F" or TH2F")
      class_name = kin_HList->At(i)->ClassName();
      //Read ith histograms in the list from current run
      if(class_name=="TH1F") {
	//Get histogram from the list
	h_i = (TH1F *)kin_HList->At(i);
	
	title  =  h_i->GetTitle();
	xlabel =  h_i->GetXaxis()->GetTitle();
	ylabel =  h_i->GetYaxis()->GetTitle();
	
	try{
	  title.replace(title.find("#"),1,"$\\");
	  title.replace(title.find("}"),1,"}$");
	  
	  xlabel.replace(xlabel.find("#"),1,"$\\");
	  xlabel.replace(xlabel.find("}"),1,"}$");
	  
	  ylabel.replace(ylabel.find("#"),1,"$\\");
	  ylabel.replace(ylabel.find("}"),1,"}$");

	}
	catch (std::out_of_range){

	}


	if( analysis_flag == "d2fsi") {
	  extract_1d_hist(h_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2fsi_pm%d_thrq%d.txt", output_hist_data.Data(), h_i->GetName(), pm_set, thrq_set));
	}
	
	if( analysis_flag == "d2pol") {
	  extract_1d_hist(h_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2pol_pm%d_Q2_%.1f.txt", output_hist_data.Data(), h_i->GetName(), pm_set, Q2_set));
	}
	
	   
      }
      
      if(class_name=="TH2F") {

	//Get histogram from the list
	h2_i = (TH2F *)kin_HList->At(i);
	
	title =  h2_i->GetTitle();
	xlabel =  h2_i->GetXaxis()->GetTitle();
	ylabel =  h2_i->GetYaxis()->GetTitle();


	try{
	  title.replace(title.find("#"),1,"$\\");
	  title.replace(title.find("}"),1,"}$");

	  
	  xlabel.replace(xlabel.find("#"),1,"$\\");
	  xlabel.replace(xlabel.find("}"),1,"}$");
	  
	  ylabel.replace(ylabel.find("#"),1,"$\\");
	  ylabel.replace(ylabel.find("}"),1,"}$");

	}
	catch (std::out_of_range){
	 
	}


	if( analysis_flag == "d2fsi") {
	  extract_2d_hist(h2_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2fsi_pm%d_thrq%d.txt", output_hist_data.Data(), h2_i->GetName(), pm_set, thrq_set));
	}
	
	if( analysis_flag == "d2pol") {
	  extract_2d_hist(h2_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2pol_pm%d_Q2_%.1f.txt", output_hist_data.Data(), h2_i->GetName(), pm_set, Q2_set));
	}
	
      }
      
    }//end loop over kin_HList

    //-----------------------------------------------------
    //Lopp over accp_HList of histogram objects 
    //-----------------------------------------------------
    
    for(int i=0; i<accp_HList->GetEntries(); i++) {
      
      //Get the class name for each element on the list (either "TH1F" or TH2F")
      class_name = accp_HList->At(i)->ClassName();
      //Read ith histograms in the list from current run
      if(class_name=="TH1F") {
	//Get histogram from the list
	h_i = (TH1F *)accp_HList->At(i);

	title =  h_i->GetTitle();
	xlabel =  h_i->GetXaxis()->GetTitle();
	ylabel =  h_i->GetYaxis()->GetTitle();

	try{
	  title.replace(title.find("#"),1,"$\\");
	  title.replace(title.find("}"),1,"}$");

	  xlabel.replace(xlabel.find("#"),1,"$\\");
	  xlabel.replace(xlabel.find("}"),1,"}$");
	  
	  ylabel.replace(ylabel.find("#"),1,"$\\");
	  ylabel.replace(ylabel.find("}"),1,"}$");

	  
	}
	catch (std::out_of_range){  }

	if( analysis_flag == "d2fsi") {
	  extract_1d_hist(h_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2fsi_pm%d_thrq%d.txt", output_hist_data.Data(), h_i->GetName(), pm_set, thrq_set));
	}

	if( analysis_flag == "d2pol") {
	  extract_1d_hist(h_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2pol_pm%d_Q2_%.1f.txt", output_hist_data.Data(), h_i->GetName(), pm_set, Q2_set));
	}
		
      }
      
      if(class_name=="TH2F") {
	//Get histogram from the list
	h2_i = (TH2F *)accp_HList->At(i);

	title =  h2_i->GetTitle();
	xlabel =  h2_i->GetXaxis()->GetTitle();
	ylabel =  h2_i->GetYaxis()->GetTitle();

	try{
	  title.replace(title.find("#"),1,"$\\");
	  title.replace(title.find("}"),1,"}$");

	  xlabel.replace(xlabel.find("#"),1,"$\\");
	  xlabel.replace(xlabel.find("}"),1,"}$");
	  
	  ylabel.replace(ylabel.find("#"),1,"$\\");
	  ylabel.replace(ylabel.find("}"),1,"}$");
	
	  
	}
	catch (std::out_of_range){  }


	if( analysis_flag == "d2fsi") {
	  extract_2d_hist(h2_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2fsi_pm%d_thrq_%d.txt", output_hist_data.Data(), h2_i->GetName(), pm_set, thrq_set));
	}
	
	if( analysis_flag == "d2pol") {
	  extract_2d_hist(h2_i, xlabel.c_str(), ylabel.c_str(), title.c_str(), Form("%s/%s_yield_d2pol_pm%d_Q2_%.1f.txt", output_hist_data.Data(), h2_i->GetName(), pm_set, Q2_set));
	}
	
      }
      
    }//end loop over accp_HList
    

}



