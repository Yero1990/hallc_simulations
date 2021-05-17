//
// program to convert a simc ascii n-tuple to a root tree
//
// reads the data from a default file name
// read data ntuple from stdin

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "./utils/parse_utils.h"

using namespace std;

// create a tree from the simc ascii file

int fmake_tree(string ikin_filename)
//int fmake_tree()
{

  cout << "Running fmake_tree.C . . ." << endl;


  //define the kin. file name (without the .extension)
  string kfname = split(ikin_filename, '.')[0];


  //Directly read normfac from outfiles/kfname.hist
  Double_t Normfac = stod(split(FindString("normfac", Form("outfiles/%s.hist", kfname.c_str()))[0], '=')[1]) ;
 
   
  // simc output is in worksim
  ifstream f_simc("worksim/"+ikin_filename);
  if ( !f_simc.is_open() ){
    cout << "worksim/"+ikin_filename+" is not open !" << endl;
    return(-1);
  }

  // files are open work on them;

  TString NTtitle; // NT title
  TString NTname ; // NT name
  string line;
  Int_t Ntags;  // number of tags

  Int_t i;
  

  // Tree variables and branch names
  Double_t val[100];
  //Double_t Normfac;
  
  TString tag_name[100]; // tag names
  TString var_branch[100]; // branch name string
  TString var_name[100]; 
  TString var_type("/D"); // all variables of the same type

  // read file
  // header information
  // this factor MUST be added by hand or by script
  // f_normfac >> Normfac;
  // this is the default output
  f_simc >> NTtitle;
  f_simc >> NTname;
  f_simc >> Ntags;

  cout << NTtitle << endl;
  cout << NTname << endl;
  cout << Ntags << endl;

  // read tag names
  for (i = 0; i<Ntags; i++) {
    f_simc >> tag_name[i];
    cout << tag_name[i] << endl;
  }


  //replace '.data' with '.root' in file name
  size_t pos = ikin_filename.find("data");
  ikin_filename.replace(pos, std::string("data").length(), "root");   // 5 = length( $name )
  
  string fpath="worksim/"+ikin_filename;

  // create output file in worksim
  TFile *f = new TFile(fpath.c_str(), "RECREATE", NTtitle);

  // Create the tree
  TTree *stree = new TTree(NTname,NTtitle);
  
  // the tree variable 
  // 1st brang contains the norm factor
  
  var_branch[0] = "Normfac"; // get rid of blanks
  var_name[0] = "Normfac"; // get rid of blanks
  var_name[0].Append( var_type );
  stree->Branch(var_branch[0].Data(), &Normfac, var_name[0]);
  cout << "adding branch : " << var_branch[0].Data() << endl;

  Int_t offset = 1;
  for (i=0; i<Ntags;i++){
    var_branch[i+offset] = tag_name[i].Strip(); // get rid of blanks
    var_name[i+offset] = tag_name[i].Strip(); // get rid of blanks
    var_name[i+offset].Append( var_type );

    // add the branches to the tree

    stree->Branch(var_branch[i+offset].Data(), &val[i], var_name[i+offset]);

    cout << "adding branch : " << var_branch[i+offset].Data() << endl;
  }

  
  // read the values and fill the tree

  Long64_t n = 0;

  //  while(! f_simc.eof() ){         // loop until eof
  while( getline(f_simc, line) ){         // loop until eof
    if (n == 0) {
      cout << "First Line: " << line << endl;
      n++;
      continue;
    }
    n++;

    if (n%1000 == 0) cout << "Event : " << n << endl;

    istringstream iss(line);

    for (i=0; i<Ntags; i++) {  // read the data for each event
     iss >> val[i];          
     // cout << "event : " << n << " tag : " << i << " value = " << val[i] << endl;
    }
    // fill the tree's branches
    stree->Fill();
  }
  cout << "Read " << n << " events" << endl;

  f->Write(); // write output file
  
  return(0);

}
