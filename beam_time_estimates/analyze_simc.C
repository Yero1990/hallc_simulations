#include "utils/parse_utils.h"
#include "utils/hist_utils.h"

void analyze_simc(Bool_t heep_check=true, int pm_set=0, TString model="", Bool_t rad_flag=true){

  
  /* 
     User Input:
     pm_set     : central missing momentum setting (580, 700, 800, etc,.)
     model      : "pwia" or "fsi"
     rad_flag   :  1 (true - simulate radiative effects), 0 (false - simulate without radiative effects)
     Ib         : beam current (uA)
     time       : beam-on-target time (hours)
  
     If doing heep check (rate estimates of H(e,e'p) elastics, then set heep_check=true, and leave pm_set and model their default values
  */
  
  TString h_arm_name = "HMS";
  TString e_arm_name = "SHMS";

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

  TString rad;
  if(rad_flag) { rad = "rad";}
  else {rad = "norad";}

  if (heep_check){

    //---Read In File Names with cuts and histogram binning information
    input_CutFileName  = "inp/set_basic_heep_cuts.inp";
    input_HBinFileName = "inp/set_basic_heep_histos.inp";

   
    TString ext0 = "kin0";
    TString ext1 = "kin1";
    TString ext2 = "kin2";
    TString ext = ext2;

    
    //Define File Name Patterns
    simc_infile = Form("infiles/cafe_heep_scan_%s_%s_new.data", ext.Data(), rad.Data());
    simc_InputFileName = Form("worksim/cafe_heep_scan_%s_%s_new.root", ext.Data(), rad.Data());
    simc_OutputFileName = Form("cafe_heep_scan_%s_%s_output.root", ext.Data(), rad.Data());
     

    /*
    simc_infile          = "infiles/d2_heep_3288.data";
    simc_InputFileName   = "worksim/d2_heep_3288.root";
    simc_OutputFileName  = "d2_heep_3288_output.root";
    */
  }

  else{

    //---Read In File Names with cuts and histogram binning information
    input_CutFileName  = Form("inp/set_basic_cuts_pm%d.inp", pm_set);
    input_HBinFileName = Form("inp/set_basic_histos_pm%d.inp", pm_set);
    
    //Define File Name Patterns
    simc_infile = Form("infiles/d2_pm%d_laget%s_%s_mod.data",  pm_set, model.Data(), rad.Data());
    simc_InputFileName = Form("worksim/d2_pm%d_laget%s_%s_mod.root", pm_set, model.Data(), rad.Data());
    simc_OutputFileName = Form("d2_pm%d_laget%s_%s_mod_output.root",  pm_set, model.Data(), rad.Data());

  }
  
  //---------------------------------------------------------------------------------------------------------

  //----------------------------
  // READ CENTRAL KIN. SETTINGS
  //----------------------------

  //beam energy (GeV)
  Double_t beam_e = (stod(split(split(FindString("Ebeam", simc_infile.Data())[0], '!')[0], '=')[1]))/1000.; 

  //e- arm central momentum setting (GeV/c)
  Double_t e_Pcen = (stod(split(split(FindString("spec%e%P", simc_infile.Data())[0], '!')[0], '=')[1]))/1000.; 
  //proton arm central momentum setting (GeV/c)
  Double_t h_Pcen = (stod(split(split(FindString("spec%p%P", simc_infile.Data())[0], '!')[0], '=')[1]))/1000.; 

  //e- arm angle (deg)
  Double_t e_angle = (stod(split(split(FindString("spec%e%theta", simc_infile.Data())[0], '!')[0], '=')[1])); 
  //p arm angle (deg)
  Double_t h_angle = (stod(split(split(FindString("spec%p%theta", simc_infile.Data())[0], '!')[0], '=')[1])); 


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
  Double_t shms_hsize = 17.;  //cm
  Double_t shms_vsize = 25.;

  //Scaling factor to scale collimator cuts from original size cut
  Double_t hms_scale=1.;   //Default
  Double_t shms_scale=1.;

  //==Collimator==
  hmsCollCut_flag = stoi(split(FindString("hmsCollCut_flag", input_CutFileName.Data())[0], '=')[1]);
  shmsCollCut_flag = stoi(split(FindString("shmsCollCut_flag", input_CutFileName.Data())[0], '=')[1]);
  
  //==Read Collimator Cut Scale Factor==
  hms_scale  =  stod(split(FindString("hms_scale", input_CutFileName.Data())[0], '=')[1]);
  shms_scale =  stod(split(FindString("shms_scale", input_CutFileName.Data())[0], '=')[1]);
  
  
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
  	   
  Double_t thx_nbins   = stod(split(FindString("thx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  //proton(hadron) angle
  Double_t thx_xmin    = stod(split(FindString("thx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t thx_xmax    = stod(split(FindString("thx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
           
  Double_t Em_nbins   = stod(split(FindString("Em_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t Em_xmin    = stod(split(FindString("Em_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t Em_xmax    = stod(split(FindString("Em_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
           	     		  	       
  Double_t Pm_nbins   = stod(split(FindString("Pm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t Pm_xmin    = stod(split(FindString("Pm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t Pm_xmax    = stod(split(FindString("Pm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);               				                 	       				  	       
  	   	     		  	       
  Double_t MM_nbins   = stod(split(FindString("MM_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t MM_xmin    = stod(split(FindString("MM_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t MM_xmax    = stod(split(FindString("MM_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
           
  Double_t MM2_nbins  = stod(split(FindString("MM2_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  Double_t MM2_xmin   = stod(split(FindString("MM2_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t MM2_xmax   = stod(split(FindString("MM2_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
           
  Double_t thxq_nbins  = stod(split(FindString("thxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  Double_t thxq_xmin   = stod(split(FindString("thxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Double_t thxq_xmax   = stod(split(FindString("thxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       	      	  	       
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
  TH1F *H_kf      = new TH1F("H_kf", "Final e^{-} Momentum", kf_nbins, kf_xmin, kf_xmax);
  H_kf->Sumw2();
  H_kf->SetDefaultSumw2(); 
  TH1F *H_the     = new TH1F("H_the", "Electron Scattering Angle, #theta_{e}", the_nbins, the_xmin, the_xmax);
  TH1F *H_Q2      = new TH1F("H_Q2","4-Momentum Transfer, Q^{2}", Q2_nbins, Q2_xmin, Q2_xmax);
  TH1F *H_Q2_nsc  = new TH1F("H_Q2_nsc","4-Momentum Transfer, Q^{2}", Q2_nbins, Q2_xmin, Q2_xmax);  //nsc stands for no self-cut (i.e, all cuts except on itself)

  TH1F *H_xbj     = new TH1F("H_xbj", "x-Bjorken", X_nbins, X_xmin, X_xmax);  
  TH1F *H_nu      = new TH1F("H_nu","Energy Transfer, #nu", nu_nbins, nu_xmin, nu_xmax); 
  TH1F *H_q       = new TH1F("H_q", "3-Momentum Transfer, |#vec{q}|", q_nbins, q_xmin, q_xmax);
  TH1F *H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}", thq_nbins, thq_xmin, thq_xmax); 
  TH1F *H_W       = new TH1F("H_W", "Invariant Mass, W", W_nbins, W_xmin, W_xmax);  
  TH1F *H_W_noCut       = new TH1F("H_W_noCut", "Invariant Mass, W (no cuts, realistic rates)", W_nbins, W_xmin, W_xmax);  

  //Secondary (Hadron) Kinematics (recoil and missing are used interchageably) ()
  TH1F *H_Pf      = new TH1F("H_Pf", "Final Hadron Momentum (detected), p_{f}", Pf_nbins, Pf_xmin, Pf_xmax);
  TH1F *H_thx     = new TH1F("H_thx", "Hadron Scattering Angle (detected), #theta_{x}", thx_nbins, thx_xmin, thx_xmax);
  TH1F *H_Em      = new TH1F("H_Em","Nuclear Missing Energy", Em_nbins, Em_xmin, Em_xmax); 
  TH1F *H_Em_nsc      = new TH1F("H_Em_nsc","Nuclear Missing Energy", Em_nbins, Em_xmin, Em_xmax); 

  TH1F *H_Pm      = new TH1F("H_Pm","Missing Momentum, P_{miss}", Pm_nbins, Pm_xmin, Pm_xmax); 
  TH1F *H_MM      = new TH1F("H_MM","Missing Mass, M_{miss}", MM_nbins, MM_xmin, MM_xmax);        
  TH1F *H_MM2     = new TH1F("H_MM2","Missing Mass Squared, M^{2}_{miss}", MM2_nbins, MM2_xmin, MM2_xmax); 
  TH1F *H_thxq    = new TH1F("H_thxq", "In-Plane Angle, #theta_{xq}", thxq_nbins, thxq_xmin, thxq_xmax);
  TH1F *H_thrq    = new TH1F("H_thrq", "In-Plane Angle, #theta_{rq}", thrq_nbins, thrq_xmin, thrq_xmax);

  //2D Pm vs. thrq (for cross section calculation)
  TH2F *H_Pm_vs_thrq  = new TH2F("H_Pm_vs_thrq", "Pm vs. #theta_{rq} (yield)", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
  TH2F *H_Pm_vs_thrq_ps  = new TH2F("H_Pm_vs_thrq_ps", "Pm vs. #theta_{rq} (phase space)", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
  TH2F *H_Pm_vs_thrq_xsec  = new TH2F("H_Pm_vs_thrq_xsec", "Pm vs. #theta_{rq} (xsec)", thrq_nbins, thrq_xmin, thrq_xmax, Pm_nbins, Pm_xmin, Pm_xmax);
    
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
  kin_HList->Add( H_thx      );
  kin_HList->Add( H_Em       );
  kin_HList->Add( H_Em_nsc   );
  kin_HList->Add( H_Pm       );
  kin_HList->Add( H_MM       );
  kin_HList->Add( H_MM2      );
  kin_HList->Add( H_thxq     );
  kin_HList->Add( H_thrq     );

  kin_HList->Add( H_Pm_vs_thrq );
  kin_HList->Add( H_Pm_vs_thrq_ps );
  //kin_HList->Add( H_Pm_vs_thrq_xsec );

  //----------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
  //----------------------------------------------------------------------


  //Electron Arm Focal Plane Quantities
  TH1F *H_exfp = new TH1F("H_exfp", Form("%s X_{fp}; X_{fp} [cm]; Counts / mC", e_arm_name.Data()), exfp_nbins, exfp_xmin, exfp_xmax);
  TH1F *H_eyfp = new TH1F("H_eyfp", Form("%s Y_{fp}; Y_{fp} [cm]; Counts / mC", e_arm_name.Data()), eyfp_nbins, eyfp_xmin, eyfp_xmax);
  TH1F *H_expfp = new TH1F("H_expfp", Form("%s X'_{fp}; X'_{fp} [rad]; Counts / mC", e_arm_name.Data()), expfp_nbins, expfp_xmin, expfp_xmax);
  TH1F *H_eypfp = new TH1F("H_eypfp", Form("%s Y'_{fp}; Y'_{fp} [rad]; Counts / mC", e_arm_name.Data()), eypfp_nbins, eypfp_xmin, eypfp_xmax);
  
  //Electron Arm Reconstructed Quantities 
  TH1F *H_eytar = new TH1F("H_eytar", Form("%s Y_{tar}; Y_{tar} [cm]; Counts / mC", e_arm_name.Data()), eytar_nbins, eytar_xmin, eytar_xmax);
  TH1F *H_exptar = new TH1F("H_exptar", Form("%s X'_{tar}; X'_{tar} [rad]; Counts / mC", e_arm_name.Data()), exptar_nbins, exptar_xmin, exptar_xmax);
  TH1F *H_eyptar = new TH1F("H_eyptar", Form("%s Y'_{tar}; Y'_{tar} [rad]; Counts / mC", e_arm_name.Data()), eyptar_nbins, eyptar_xmin, eyptar_xmax);
  TH1F *H_edelta = new TH1F("H_edelta", Form("%s Momentum Acceptance, #delta; #delta [%%]; Counts / mC", e_arm_name.Data()), edelta_nbins, edelta_xmin, edelta_xmax);
  TH1F *H_edelta_nsc = new TH1F("H_edelta_nsc", Form("%s Momentum Acceptance, #delta; #delta [%%]; Counts / mC", e_arm_name.Data()), edelta_nbins, edelta_xmin, edelta_xmax);
  
  //Hadron arm Focal Plane Quantities
  TH1F *H_hxfp = new TH1F("H_hxfp", Form("%s  X_{fp}; X_{fp} [cm]; Counts / mC", h_arm_name.Data()), hxfp_nbins, hxfp_xmin, hxfp_xmax);
  TH1F *H_hyfp = new TH1F("H_hyfp", Form("%s  Y_{fp}; Y_{fp} [cm]; Counts / mC", h_arm_name.Data()), hyfp_nbins, hyfp_xmin, hyfp_xmax);
  TH1F *H_hxpfp = new TH1F("H_hxpfp", Form("%s  X'_{fp}; X'_{fp} [rad]; Counts / mC", h_arm_name.Data()), hxpfp_nbins, hxpfp_xmin, hxpfp_xmax );
  TH1F *H_hypfp = new TH1F("H_hypfp", Form("%s  Y'_{fp}; Y'_{fp} [rad]; Counts / mC", h_arm_name.Data()), hypfp_nbins, hypfp_xmin, hypfp_xmax);

  //Hadron arm Reconstructed Quantities 
  TH1F *H_hytar = new TH1F("H_hytar", Form("%s  Y_{tar}; Y_{tar} [cm]; Counts / mC", h_arm_name.Data()), hytar_nbins, hytar_xmin, hytar_xmax);
  TH1F *H_hxptar = new TH1F("H_hxptar", Form("%s  X'_{tar}; X'_{tar} [rad]; Counts / mC", h_arm_name.Data()), hxptar_nbins, hxptar_xmin, hxptar_xmax);
  TH1F *H_hyptar = new TH1F("H_hyptar", Form("%s  Y'_{tar}; Y'_{tar} [rad]; Counts / mC", h_arm_name.Data()), hyptar_nbins, hyptar_xmin, hyptar_xmax );
  TH1F *H_hdelta = new TH1F("H_hdelta", Form("%s  Momentum Acceptance, #delta; #delta [%%]; Counts / mC", h_arm_name.Data()), hdelta_nbins, hdelta_xmin, hdelta_xmax);
  TH1F *H_hdelta_nsc = new TH1F("H_hdelta_nsc", Form("%s  Momentum Acceptance, #delta; #delta [%%]; Counts / mC", h_arm_name.Data()), hdelta_nbins, hdelta_xmin, hdelta_xmax);
  

  //Target Reconstruction (Hall Coord. System) 
  TH1F *H_htar_x = new TH1F("H_htar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts / mC", h_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  TH1F *H_htar_y = new TH1F("H_htar_y", Form("%s y_Target (Lab); y-Target [cm]; Counts / mC", h_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  TH1F *H_htar_z = new TH1F("H_htar_z", Form("%s z_Target (Lab); z-Target [cm]; Counts / mC", h_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);
  TH1F *H_etar_x = new TH1F("H_etar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts / mC", e_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  TH1F *H_etar_y = new TH1F("H_etar_y", Form("%s y-Target (Lab); y-Target [cm]; Counts / mC", e_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  TH1F *H_etar_z = new TH1F("H_etar_z", Form("%s z-Target (Lab); z-Target [cm]; Counts / mC", e_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);

  //difference in reaction vertex z (user-defined)
  TH1F *H_ztar_diff = new TH1F("H_ztar_diff", "Ztar Difference; z-Target Difference [cm]; Counts / mC", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);
  TH1F *H_ztar_diff_nsc = new TH1F("H_ztar_diff_nsc", "Ztar Difference; z-Target Difference [cm]; Counts / mC", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);


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
  TH2F *H_hdelta_vs_edelta_nsc = new TH2F("H_hdelta_vs_edelta_nsc", "HMS vs. SHMS, #delta",   edelta_nbins, edelta_xmin, edelta_xmax, hdelta_nbins, hdelta_xmin, hdelta_xmax);
  
  

  
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
  accp_HList->Add( H_hdelta_vs_edelta_nsc);
  
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

  //--- Define Leaf Variable Names ----

  //Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)
  Double_t ki;                    //initial electron momentum  
  Double_t kf;                    //final electron momentum
  Double_t theta_e;               //Central electron arm angle relative to +z (hall coord. system)
  Double_t Q2;                   //Four-momentum trasfer
  Double_t X;                    //B-jorken X  scaling variable
  Double_t nu;                   //Energy Transfer
  Double_t q;                  //magnitude of the 3-vector q
  Double_t th_q;                 //angle between q and +z (hall coord. system)
  Double_t W;                    // invariant mass
  
  //Secondary Kinematics (USED BY DATA AND SIMC)
  Double_t Pf;                     //final proton momentum
  Double_t Ep;                      //final proton energy (needs to be calculated)
  Double_t theta_p;               //to be calculated separately (in data)
  Double_t Em;                     //Standard Missing Energy for H(e,e'p)
  Double_t Pm;                     //Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
  Double_t MM;                   //Missing Mass (neutron Mass)
  Double_t MM2;                   //Missing Mass Squared
  Double_t th_xq;                  //detected particle in-plane angle w.r.to q-vector
  Double_t th_rq;                  //recoil particle in-plane angle w.r.to q-vector

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

  
  //----- Set Branch Address ------

  //Primary Kinematics (electron kinematics)
  //ki needs to be calculated (initial e- momentum)
  tree->SetBranchAddress("e_pf",    &kf);
  tree->SetBranchAddress("theta_e", &theta_e);
  tree->SetBranchAddress("Q2", &Q2);  
  //Xbj needs to be calculated in the event loop
  tree->SetBranchAddress("nu", &nu);
  tree->SetBranchAddress("q", &q);
  //th_q needs to be calculated in the event loop
  tree->SetBranchAddress("W", &W);

  //Secondary Kinematics (hadron kinematics)
  tree->SetBranchAddress("h_pf",    &Pf);
  tree->SetBranchAddress("theta_p", &theta_p);
  tree->SetBranchAddress("Em", &Em);
  tree->SetBranchAddress("Pm", &Pm);
  //Missing Mass (MM) and MM2 will be defined in entry loop
  tree->SetBranchAddress("theta_pq", &th_xq); 
  tree->SetBranchAddress("theta_rq", &th_rq);  
 
  
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
  //  DEFINE FULL WEIGHT VARIABLES
  //-------------------------------

  // STEP1: Determine the charge factor:
  // definition: total charge deposited on target over a time period
  // SIMC input files are set to 'events / 1mC'
   
  // Charge factor is the total integrated charge assuming a beam current and run time
  Double_t Ib = 60;       //beam current in (uA) microAmps (micro-Coulombs / sec),   1 mC = 1000 uC
  Double_t time = 1.0;     //estimated time (in hours) a run takes (start - end) of run
  Double_t charge_factor = Ib * time * 3600. / 1000.;

  //target boiling slopes for Hydrofen and Deuterium (during commissioning)
  //Double_t LH2_slope = 0.00063396; 
  //Double_t LD2_slope = 0.00080029;  //NORM. yield loss / uA

  //more realistic target boiling slopes based on improved measurements later on
  Double_t LH2_slope = 0.0002; 
  Double_t LD2_slope = 0.00025;
  
  // STEP2: Estimate Efficiencies (use efficiencies from commissioning experiment)
  // coin. rates were ~ 2.5 Hz in commissioning,
  //Double_t e_trk      = 0.964;
  //Double_t h_trk      = 0.988;
  //Double_t daq_lt     = 0.98;   //(it was 0.926 during commissioning due to large logic windows ~100 ns in HMS, but now is smaller)

  //for heep checks
  Double_t e_trk      = 0.99;
  Double_t h_trk      = 0.99;
  Double_t daq_lt     = 0.99;
  
  Double_t tgt_boil;
  if(heep_check){
    tgt_boil   = 1. - LH2_slope * Ib;
  }
  else{
    tgt_boil   = 1. - LD2_slope * Ib;
  }

  Double_t proton_abs = 1.0; // 0.9534;  let assume no proton absorption thru material (since for heep singles, only electron thru SHMS, and does NOT get absorbed) 

  Double_t eff_factor;

  //assuming heep check singles
  if(heep_check){
    eff_factor = e_trk * daq_lt * tgt_boil;
  }
  
  else{
    eff_factor = e_trk * h_trk * daq_lt * tgt_boil * proton_abs;
  }


  Double_t FullWeight;
  Double_t FullWeight_forRates; //this is the full weight for true rate estimates (so no inefficiencies accounted for)
  Double_t PhaseSpace;

  /* ----------------------------------------------------------
     
     NOTE on FullWeight and Phase Space definitions : 
     
     FullWeight = Weight * Normfac / n_acc   (assuming 1 mC of charge and detector efficiencies are 100 % (in fraction is 1.0 ) )
     
     1)  Weight = (sigma_theory*Jacobian) * main%gen_weight * jacobian_corr * coulomb_corr (see event.f) for 'LAGET_DEUT'
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

    //-----Define Additional Kinematic Variables--------
    
    X = Q2 / (2.*MP*nu);
    th_q = th_xq + theta_p;
    ztar_diff =  htar_z - etar_z;

    Pf = Pf/1000.;  //final proton momentum (GeV/c)
    Ep = sqrt(MP*MP + Pf*Pf);
    
    //Missing Mass
    if(heep_check){
      MM2 = Em*Em - Pm*Pm;
      MM = sqrt(MM2);
    }
    else{
      MM = sqrt( pow(nu+MD-Ep,2) - Pm*Pm );  //recoil mass (neutron missing mass)
      MM2 = MM * MM;
    }

    
    //SIMC Collimator (definition based on HCANA collimator)
    htarx_corr = tar_x - h_xptar*htar_z*cos(h_angle*dtr);
    etarx_corr = tar_x - e_xptar*etar_z*cos(e_angle*dtr);  
    
    
    //Define Collimator (same as in HCANA)
    hXColl = htarx_corr + h_xptar*168.;   //in cm
    hYColl = h_ytar + h_yptar*168.;
    eXColl = etarx_corr + e_xptar*253.;
    eYColl = e_ytar + e_yptar*253.-(0.019+40.*.01*0.052)*e_delta+(0.00019+40*.01*.00052)*e_delta*e_delta; //correct for HB horizontal bend	  

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
    
    //Missing Energy, Em
    if(Em_cut_flag){c_Em = Em>=c_Em_min && Em<=c_Em_max;}
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

    /*
    cout << "FullWeight = " << FullWeight << endl;
    cout << "Normfac = " << Normfac << endl;
    cout << "charge_factor = " << charge_factor << endl;
    cout << "eff_factor = " << eff_factor << endl;
    cout << "Weight = " << Weight << endl;
    cout << "nentries = " << nentries << endl;
    cout << "Jacobian_corr = " << Jacobian_corr << endl;
    */

    //fill histogram (no cuts)
    H_W_noCut->Fill(W, FullWeight_forRates);
     
    if(c_allCuts) {


      // This is for the 2D cross section Pm vs. thnq binned in thnq 
      H_Pm_vs_thrq->Fill(th_rq/dtr, Pm, FullWeight);
      H_Pm_vs_thrq_ps->Fill(th_rq/dtr, Pm, PhaseSpace);
      
      //Primary (electron) Kinematics
      H_kf->Fill(kf/1000., FullWeight);
      H_the->Fill(theta_e/dtr, FullWeight);
      H_Q2->Fill(Q2, FullWeight);
      H_xbj->Fill(X, FullWeight);
      H_nu->Fill(nu, FullWeight);
      H_q->Fill(q, FullWeight);
      H_thq->Fill(th_q/dtr, FullWeight);
      H_W->Fill(W, FullWeight);

 
      //Secondary (Hadron) Kinematics
      H_Pf->Fill(Pf, FullWeight);
      H_thx->Fill(theta_p/dtr, FullWeight);
      H_Em->Fill(Em, FullWeight);
      H_Pm->Fill(Pm, FullWeight);     
      H_thxq->Fill(th_xq/dtr, FullWeight);
      H_thrq->Fill(th_rq/dtr, FullWeight);
      H_MM->Fill(MM, FullWeight);
      H_MM2->Fill(MM2, FullWeight);
   
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
    if(c_accpCuts &&  c_Q2) { H_Em_nsc->Fill(Em, FullWeight); }
    if(c_accpCuts &&  c_Em) { H_Q2_nsc->Fill(Q2, FullWeight); }

    if(c_kinCuts && c_hdelta && c_edelta && hmsColl_Cut && shmsColl_Cut) { H_ztar_diff_nsc->Fill(ztar_diff, FullWeight); }

    if(c_kinCuts && c_edelta && c_ztarDiff && hmsColl_Cut && shmsColl_Cut) {H_hdelta_nsc->Fill(h_delta, FullWeight); }
    if(c_kinCuts && c_hdelta && c_ztarDiff && hmsColl_Cut && shmsColl_Cut) {H_edelta_nsc->Fill(e_delta, FullWeight); }

    if(c_kinCuts && c_hdelta && c_edelta && c_ztarDiff && shmsColl_Cut) { H_hXColl_vs_hYColl_nsc->Fill(hYColl, hXColl, FullWeight); }
    if(c_kinCuts && c_hdelta && c_edelta && c_ztarDiff && hmsColl_Cut) { H_eXColl_vs_eYColl_nsc->Fill(eYColl, eXColl, FullWeight); }

    
    cout << "SIMC Events Completed: " << std::setprecision(2) << double(i) / nentries * 100. << "  % " << std::flush << "\r";
    
  } // end entry loop


  //Calculate Xsec
  //H_Pm_vs_thrq_xsec->Divide(H_Pm_vs_thrq, H_Pm_vs_thrq_ps);
  
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

  
  //------------------------------------------
  // Extract The Yield binned in Pm vs th_rq
  //------------------------------------------
  // extract_2d_hist(H_Pm_vs_thrq, "#theta_{rq} [deg]", "Missing Momentum, P_{m} [GeV/c]", Form("yield_pm%d_model%s_%s_%.1fuA_%.1fhr.txt",  pm_set, model.Data(), rad.Data(), Ib, time));

  //--------
  // Extrack numerical data for histogram plotting
  //--------

  /*
  //extract cuts-related histos
  extract_1d_hist(H_Q2_nsc, "4-Momentum Transfers, Q2 [GeV^2]", Form("yield_Q2_pm%d.txt", pm_set));
  extract_1d_hist(H_Em_nsc, "Missing Energy, Em [GeV]", Form("yield_Em_pm%d.txt", pm_set));
  extract_1d_hist(H_edelta_nsc, "SHMS Delta [%%]", Form("yield_edelta_pm%d.txt", pm_set));
  extract_1d_hist(H_hdelta_nsc, "HMS Delta [%%]", Form("yield_hdelta_pm%d.txt", pm_set));
  extract_1d_hist(H_ztar_diff_nsc, "Ztar diff [cm]", Form("yield_ztardiff_pm%d.txt", pm_set));

  extract_2d_hist(H_hXColl_vs_hYColl_nsc, "HMS Y Collimator [cm]", "HMS X Collimator [cm]", Form("yield_hColl_pm%d.txt",  pm_set));
  extract_2d_hist(H_eXColl_vs_eYColl_nsc, "SHMS Y Collimator [cm]", "SHMS X Collimator [cm]", Form("yield_eColl_pm%d.txt",  pm_set));  
  extract_2d_hist(H_hdelta_vs_edelta_nsc, "SHMS Delta [%%]", "HMS Delta [%%]", Form("yield_delta2d_pm%d.txt",  pm_set));

  //extarct focal plane histos
  extract_2d_hist(H_hxfp_vs_hyfp, "HMS Y Focal Plane [cm]", "HMS X Focal Plane [cm]", Form("yield_hfp_pm%d.txt", pm_set));
  extract_2d_hist(H_exfp_vs_eyfp, "SHMS Y Focal Plane [cm]", "SHMS X Focal Plane [cm]", Form("yield_efp_pm%d.txt", pm_set));

  //extract electron kinematics
  extract_1d_hist(H_kf, "Final Electron Momentum, kf, [GeV]", Form("yield_kf_pm%d.txt", pm_set));
  extract_1d_hist(H_the, "Final Electron Angle, th_e, [deg]", Form("yield_the_pm%d.txt", pm_set));
  extract_1d_hist(H_nu, "Energy Transfer, nu [GeV]", Form("yield_nu_pm%d.txt", pm_set));
  extract_1d_hist(H_q, "3-momentum transfer, |q| [GeV/c]", Form("yield_q_pm%d.txt", pm_set));
  extract_1d_hist(H_thq, "Recoil Angle q relative to +z (lab), th_q [deg]", Form("yield_thq_pm%d.txt", pm_set));
  extract_1d_hist(H_xbj, "x-Bjorken", Form("yield_xbj_pm%d.txt", pm_set));
  */
  
  //extract hadron kinematics
  //extract_1d_hist(H_MM, "Missing Mass, MM [GeV]", Form("yield_MM_pm%d.txt", pm_set));  
  //extract_1d_hist(H_Pm, "Missing Momentum, Pm [GeV/c]", Form("yield_Pm_pm%d_noCUTS.txt", pm_set));
  /*
  extract_1d_hist(H_Pf, "Final Proton Momentum, Pf, [GeV]", Form("yield_Pf_pm%d.txt", pm_set));
  extract_1d_hist(H_thx, "Proton Scattering Angle, th_p, [deg]", Form("yield_thp_pm%d.txt", pm_set));
  extract_1d_hist(H_thxq, "Proton Angle wrto q-vector, thxq, [deg]", Form("yield_thxq_pm%d.txt", pm_set));
  extract_1d_hist(H_thrq, "Neutron Angle wrto q-vector, thrq, [deg]", Form("yield_thrq_pm%d.txt", pm_set));
  */

}
