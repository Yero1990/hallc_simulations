#include "parse_utils.h"

void analyze_simc(int pm_set, TString model="fsi", TString rad = "norad"){


  TString h_arm_name = "HMS";
  TString e_arm_name = "SHMS";
  
  //-------------------
  // READ FILENAMES
  //-------------------

  //Input parameter controls filenames
  TString main_controls_fname;
  TString input_CutFileName;
  TString input_HBinFileName;
  TString input_SetFileName;

  //Declare TFile Pointers (reading/writing ROOTfiles)
  TFile *inROOT;
  TFile *outROOT;

  //Input ROOTfile Name (to be read)
  TString simc_InputFileName;
  TString simc_OutputFileName;
  
  TString temp; //temporary string placeholder


  
  //---Read In File Names with cuts and histogram binning information

  main_controls_fname = "main_controls.inp";

  input_SetFileName = trim(split(FindString("fname_file", main_controls_fname.Data())[0], '=')[1]);
  input_CutFileName  = trim(split(FindString("cuts_file", main_controls_fname.Data())[0], '=')[1]);
  input_HBinFileName = trim(split(FindString("hist_file", main_controls_fname.Data())[0], '=')[1]);

  //Define Input/Output (.root) File Name Patterns
  temp = trim(split(FindString("input_ROOTfilePattern", input_SetFileName.Data())[0], '=')[1]);
  simc_InputFileName = Form(temp.Data(),  pm_set, model.Data(), rad.Data());

  temp = trim(split(FindString("output_ROOTfilePattern", input_SetFileName.Data())[0], '=')[1]);
  simc_OutputFileName = Form(temp.Data(),  pm_set, model.Data(), rad.Data());

  //---------------------------------------------------------------------------------------------------------
  
  //--------------------
  // READ ANALYSIS CUTS
  //-------------------

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
  TH1F *H_the     = new TH1F("H_the", "Electron Scattering Angle, #theta_{e}", the_nbins, the_xmin, the_xmax);
  TH1F *H_Q2      = new TH1F("H_Q2","4-Momentum Transfer, Q^{2}", Q2_nbins, Q2_xmin, Q2_xmax); 
  TH1F *H_xbj     = new TH1F("H_xbj", "x-Bjorken", X_nbins, X_xmin, X_xmax);  
  TH1F *H_nu      = new TH1F("H_nu","Energy Transfer, #nu", nu_nbins, nu_xmin, nu_xmax); 
  TH1F *H_q       = new TH1F("H_q", "3-Momentum Transfer, |#vec{q}|", q_nbins, q_xmin, q_xmax);
  TH1F *H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}", thq_nbins, thq_xmin, thq_xmax); 

  //Secondary (Hadron) Kinematics (recoil and missing are used interchageably) ()
  TH1F *H_Pf      = new TH1F("H_Pf", "Final Hadron Momentum (detected), p_{f}", Pf_nbins, Pf_xmin, Pf_xmax);
  TH1F *H_thx     = new TH1F("H_thx", "Hadron Scattering Angle (detected), #theta_{x}", thx_nbins, thx_xmin, thx_xmax);
  TH1F *H_Em_nuc  = new TH1F("H_Em_nuc","Nuclear Missing Energy", Em_nbins, Em_xmin, Em_xmax); 
  TH1F *H_Pm      = new TH1F("H_Pm","Missing Momentum, P_{miss}", Pm_nbins, Pm_xmin, Pm_xmax); 
  TH1F *H_MM      = new TH1F("H_MM","Missing Mass, M_{miss}", MM_nbins, MM_xmin, MM_xmax);        
  TH1F *H_MM2     = new TH1F("H_MM2","Missing Mass Squared, M^{2}_{miss}", MM2_nbins, MM2_xmin, MM2_xmax); 
  TH1F *H_thxq    = new TH1F("H_thxq", "In-Plane Angle, #theta_{xq}", thxq_nbins, thxq_xmin, thxq_xmax);
  TH1F *H_thrq    = new TH1F("H_thrq", "In-Plane Angle, #theta_{rq}", thrq_nbins, thrq_xmin, thrq_xmax);

    
  //Add Kin Histos to TList

  //Add Primary Kin Histos
  kin_HList->Add( H_kf     );
  kin_HList->Add( H_the    );
  kin_HList->Add( H_Q2     );
  kin_HList->Add( H_xbj    );
  kin_HList->Add( H_nu     );
  kin_HList->Add( H_q      );
  kin_HList->Add( H_thq    );

  //Add Secondary Kin Histos
  kin_HList->Add( H_Pf       );
  kin_HList->Add( H_thx      );
  kin_HList->Add( H_Em_nuc   );
  kin_HList->Add( H_Pm       );
  kin_HList->Add( H_MM       );
  kin_HList->Add( H_MM2      );
  kin_HList->Add( H_thxq     );
  kin_HList->Add( H_thrq     );


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
  

  //Target Reconstruction (Hall Coord. System) 
  TH1F *H_htar_x = new TH1F("H_htar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts / mC", h_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  TH1F *H_htar_y = new TH1F("H_htar_y", Form("%s y_Target (Lab); y-Target [cm]; Counts / mC", h_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  TH1F *H_htar_z = new TH1F("H_htar_z", Form("%s z_Target (Lab); z-Target [cm]; Counts / mC", h_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);
  TH1F *H_etar_x = new TH1F("H_etar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts / mC", e_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  TH1F *H_etar_y = new TH1F("H_etar_y", Form("%s y-Target (Lab); y-Target [cm]; Counts / mC", e_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  TH1F *H_etar_z = new TH1F("H_etar_z", Form("%s z-Target (Lab); z-Target [cm]; Counts / mC", e_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);

  //difference in reaction vertex z (user-defined)
  TH1F *H_ztar_diff = new TH1F("H_ztar_diff", "Ztar Difference; z-Target Difference [cm]; Counts / mC", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);


  //2D Collimator Histos
  TH2F *H_hXColl_vs_hYColl = new TH2F("H_hXColl_vs_hYColl", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", h_arm_name.Data(), h_arm_name.Data(), h_arm_name.Data()), hYColl_nbins, hYColl_xmin, hYColl_xmax,  hXColl_nbins, hXColl_xmin, hXColl_xmax);
  TH2F *H_eXColl_vs_eYColl = new TH2F("H_eXColl_vs_eYColl", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", e_arm_name.Data(), e_arm_name.Data(), e_arm_name.Data()), eYColl_nbins, eYColl_xmin, eYColl_xmax, eXColl_nbins, eXColl_xmin, eXColl_xmax); 
  
  //2D Hour Glass Histos
  TH2F *H_hxfp_vs_hyfp  = new TH2F("H_hxfp_vs_hyfp", Form("%s  X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", h_arm_name.Data()),  hyfp_nbins, hyfp_xmin, hyfp_xmax, hxfp_nbins, hxfp_xmin, hxfp_xmax);
  TH2F *H_exfp_vs_eyfp  = new TH2F("H_exfp_vs_eyfp", Form("%s  X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", e_arm_name.Data()),  eyfp_nbins, eyfp_xmin, eyfp_xmax, exfp_nbins, exfp_xmin, exfp_xmax);

  
  //Add ACCP Histos to TList
  accp_HList->Add( H_exfp       );
  accp_HList->Add( H_eyfp       );
  accp_HList->Add( H_expfp      );
  accp_HList->Add( H_eypfp      );

  accp_HList->Add( H_eytar       );
  accp_HList->Add( H_exptar      );
  accp_HList->Add( H_eyptar      );
  accp_HList->Add( H_edelta      );
  
  accp_HList->Add( H_hxfp       );
  accp_HList->Add( H_hyfp       );
  accp_HList->Add( H_hxpfp      );
  accp_HList->Add( H_hypfp      );
  
  accp_HList->Add( H_hytar       );
  accp_HList->Add( H_hxptar      );
  accp_HList->Add( H_hyptar      );
  accp_HList->Add( H_hdelta      );

  accp_HList->Add( H_htar_x       );
  accp_HList->Add( H_htar_y       );
  accp_HList->Add( H_htar_z       );
  accp_HList->Add( H_etar_x       );
  accp_HList->Add( H_etar_y       );
  accp_HList->Add( H_etar_z       );
  accp_HList->Add( H_ztar_diff    );

  accp_HList->Add( H_hXColl_vs_hYColl  );
  accp_HList->Add( H_eXColl_vs_eYColl  );
  
  accp_HList->Add( H_hxfp_vs_hyfp  );
  accp_HList->Add( H_exfp_vs_eyfp  );


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
  //Electron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t e_xfp;
  Double_t e_xpfp;
  Double_t e_yfp;
  Double_t e_ypfp;
  
  Double_t e_ytar;
  Double_t e_yptar;
  Double_t e_xptar;
  Double_t e_delta;
  Double_t kf;                        //final electron momentum
  Double_t ki;                        //initial electron momentum

  //Hadron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t h_xfp;
  Double_t h_xpfp;
  Double_t h_yfp;
  Double_t h_ypfp;
  
  Double_t h_ytar;
  Double_t h_yptar;
  Double_t h_xptar;
  Double_t h_delta;
  Double_t Pf;                 //final proton momentum
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System (USED BY DATA AND SIMC)
  Double_t tar_x; //For SIMC ONLY (It is the same for HMS/SHMS)

  Double_t  htar_x;
  Double_t  htar_y;
  Double_t  htar_z;
  
  Double_t  etar_x;
  Double_t  etar_y;
  Double_t  etar_z;

  Double_t ztar_diff;

  //Collimators
  Double_t hXColl, hYColl, eXColl, eYColl;

  //Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)
  Double_t theta_e;              //Central electron arm angle relative to +z (hall coord. system)
  Double_t W;                    //Invariant Mass W (should be proton mass in H(e,e'p))
  Double_t W2;                    //Invariant mass squared
  Double_t Q2;                   //Four-momentum trasfer
  Double_t X;                    //B-jorken X  scaling variable
  Double_t nu;                   //Energy Transfer
  Double_t q;                  //magnitude of the 3-vector q
  Double_t th_q;                 //angle between q and +z (hall coord. system)

  //Secondary Kinematics (USED BY DATA AND SIMC)
  Double_t Ep;                     //proton energy
  Double_t Em;                    //Standard Missing Energy for H(e,e'p)
  Double_t Em_nuc;                //Nuclear definition of Missing Energy (Used for D(e,e'p): B.E. of deuteron)
  Double_t Pm;                    //Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
  Double_t Pmx_lab;               //X-Component of Missing Momentum (in Lab(or Hall) frame. +X: beam left, +Y: up, +Z: downstream beam) 
  Double_t Pmy_lab;
  Double_t Pmz_lab;
  Double_t Pmx_q;                 //X-Component of Missing Momentum (in frame where +z_lab is rotated to +z_q. Pmz_q is along +z(parallel to q))
  Double_t Pmy_q;
  Double_t Pmz_q;
  Double_t Kp;                    //Kinetic Energy of detected particle (proton)
  Double_t Kn;                    //Kinetic Energy of recoil system (neutron)
  Double_t M_recoil;              //Missing Mass (neutron Mass)
  Double_t MM2;                   //Missing Mass Squared
  Double_t E_recoil;              //Recoil Energy of the system (neutron total energy)
  Double_t En;                    //Same as above
  Double_t th_pq;                  //Polar angle of detected particle with q   ----> th_pq
  Double_t th_nq;                  //Polar angle of recoil system with q (rad)  ---> th_nq (neutreon-q angle. IMPORTANT in D(e,e'p))
  Double_t ph_pq;                  //Azimuth angle of detected particle with q    ----> phi_pq angle between proton and q-vector
  Double_t ph_nq;                  //Azimuth of recoil system with scattering plane (rad) ----> phi_nq angle between neutron and q-vector
  Double_t xangle;                //Angle of detected particle with scattered electron (Used to determine hadron angle)
  Double_t theta_p;               //to be calculated separately (in data)

  //Light-Cone Momentum Variables
  Double_t PmPar;     //parallel component of recoil momentum relative to q-vector
  Double_t PmPerp;    //transverse component of recoil momentum relative to q-vector
  Double_t alpha_n;   //light-cone momentum fraction of the recoil neutron
  Double_t alpha;     //momentum fraction of struck nucleon (normalized such that: alpha + alpha_n = 2)
  
  //SIMC Specific TTree Variable Names
  Double_t Normfac;
  
  //Thrown quantities (Used to determine spec. resolution)
  Double_t h_deltai;
  Double_t h_yptari;
  Double_t h_xptari;
  Double_t h_ytari;
  
  Double_t e_deltai;
  Double_t e_yptari;
  Double_t e_xptari;
  Double_t e_ytari;
  
  Double_t epsilon;
  Double_t corrsing;
  Double_t fry;
  Double_t radphot;
  Double_t sigcc;
  Double_t Weight;               //This Weight has the cross section in it
  Double_t Jacobian;
  Double_t Genweight;
  Double_t SF_weight;
  Double_t Jacobian_corr;
  Double_t sig;
  Double_t sig_recon;
  Double_t sigcc_recon;
  Double_t coul_corr;
  Double_t Ein;                  //single beam energy value (SIMC Uses this energy. If not corr. for energy loss, it should be same as in input file)
  Double_t theta_rq;
  Double_t SF_weight_recon;
  Double_t h_Thf;

  Double_t prob_abs;  // Probability of absorption of particle in the HMS Collimator
                      //(Must be multiplies by the weight. If particle interation is
                      //NOT simulated, it is set to 1.)

  
  //SIMC Collimator
  Double_t htarx_corr;
  Double_t etarx_corr;

  
  //----- Set Branch Address ------
  //Electron Arm Focal Plane / Reconstructed Quantities 
  tree->SetBranchAddress("e_xfp",  &e_xfp);
  tree->SetBranchAddress("e_xpfp", &e_xpfp);
  tree->SetBranchAddress("e_yfp",  &e_yfp);
  tree->SetBranchAddress("e_ypfp", &e_ypfp);
  
  tree->SetBranchAddress("e_ytar",  &e_ytar);
  tree->SetBranchAddress("e_yptar", &e_yptar);
  tree->SetBranchAddress("e_xptar", &e_xptar);
  tree->SetBranchAddress("e_delta", &e_delta);
  tree->SetBranchAddress("e_pf",    &kf);
  
  //Hadron Arm Focal Plane / Reconstructed Quantities 
  tree->SetBranchAddress("h_xfp",  &h_xfp);
  tree->SetBranchAddress("h_xpfp", &h_xpfp);
  tree->SetBranchAddress("h_yfp",  &h_yfp);
  tree->SetBranchAddress("h_ypfp", &h_ypfp);
        
  tree->SetBranchAddress("h_ytar",  &h_ytar);
  tree->SetBranchAddress("h_yptar", &h_yptar);
  tree->SetBranchAddress("h_xptar", &h_xptar);
  tree->SetBranchAddress("h_delta", &h_delta);
  tree->SetBranchAddress("h_pf",    &Pf);
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System
  tree->SetBranchAddress("tar_x", &tar_x);
  tree->SetBranchAddress("h_yv",  &htar_y);
  tree->SetBranchAddress("h_zv",  &htar_z);
  tree->SetBranchAddress("e_yv",  &etar_y);
  tree->SetBranchAddress("e_zv",  &etar_z);
    
  //Primary Kinematics (electron kinematics)
  tree->SetBranchAddress("theta_e", &theta_e);
  tree->SetBranchAddress("W", &W);
  tree->SetBranchAddress("Q2", &Q2);  
  //Xbj needs to be calculated in the event loop
  tree->SetBranchAddress("nu", &nu);
  tree->SetBranchAddress("q", &q);
  //th_q needs to be calculated in the event loop

  //Secondary Kinematics (hadron kinematics)
  tree->SetBranchAddress("Em", &Em);
  tree->SetBranchAddress("Pm", &Pm);
  
  tree->SetBranchAddress("theta_pq", &th_pq); 
  tree->SetBranchAddress("phi_pq", &ph_pq);  
  tree->SetBranchAddress("theta_p", &theta_p);
  
  //SIMC-SPECIFIC LEAF VARIABLES (Not all may be used here)
  tree->SetBranchAddress("Normfac",  &Normfac);
  tree->SetBranchAddress("h_deltai", &h_deltai);
  tree->SetBranchAddress("h_yptari", &h_yptari);
  tree->SetBranchAddress("h_xptari", &h_xptari);
  tree->SetBranchAddress("h_ytari",  &h_ytari);
  
  tree->SetBranchAddress("e_deltai", &e_deltai);
  tree->SetBranchAddress("e_yptari", &e_yptari);
  tree->SetBranchAddress("e_xptari", &e_xptari);
  tree->SetBranchAddress("e_ytari",  &e_ytari);
  
  tree->SetBranchAddress("epsilon",  &epsilon);
  tree->SetBranchAddress("corrsing", &corrsing);
  tree->SetBranchAddress("fry",      &fry);
  tree->SetBranchAddress("radphot",  &radphot);
  tree->SetBranchAddress("sigcc",    &sigcc);
  tree->SetBranchAddress("Weight",   &Weight);
  tree->SetBranchAddress("Jacobian", &Jacobian);
  tree->SetBranchAddress("Genweight",&Genweight);
  tree->SetBranchAddress("SF_weight", &SF_weight);
  tree->SetBranchAddress("Jacobian_corr", &Jacobian_corr);
  tree->SetBranchAddress("sig", &sig);
  tree->SetBranchAddress("sig_recon", &sig_recon);
  tree->SetBranchAddress("sigcc_recon", &sigcc_recon);
  tree->SetBranchAddress("coul_corr", &coul_corr);
  tree->SetBranchAddress("Ein", &Ein);
  tree->SetBranchAddress("theta_rq", &theta_rq);
  tree->SetBranchAddress("SF_weight_recon", &SF_weight_recon);
  tree->SetBranchAddress("h_Thf", &h_Thf);
  tree->SetBranchAddress("Ein_v", &Ein_v);
  tree->SetBranchAddress("Q2_v", &Q2_v);
  tree->SetBranchAddress("nu_v", &nu_v);
  tree->SetBranchAddress("q_lab_v", &q_lab_v);
  tree->SetBranchAddress("pm_v", &Pm_v);
  tree->SetBranchAddress("pm_par_v", &Pm_par_v);
  tree->SetBranchAddress("pf_v", &Pf_v);
  tree->SetBranchAddress("Ep_v", &Ep_v);
  tree->SetBranchAddress("Ef_v", &Ef_v);
  tree->SetBranchAddress("e_xptar_v", &e_xptar_v);
  tree->SetBranchAddress("e_yptar_v", &e_yptar_v);
  tree->SetBranchAddress("h_xptar_v", &h_xptar_v);
  tree->SetBranchAddress("h_yptar_v", &h_yptar_v);
  tree->SetBranchAddress("probabs", &prob_abs);
  
  
}
