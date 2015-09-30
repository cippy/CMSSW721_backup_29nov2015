#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <math.h>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>
#include "EmanTreeAnalysis.h"
#include "AdishTreeAnalysis.h"
#include "whichApplication.h"

using namespace myAnalyzerTEman;
using namespace myAnalyzerTAdish;

using namespace std;


int main(int argc, char* argv[]) {


//================ Creating chain 
 
#if Application == 1 

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;

  char inputFileName[200];
  char buffer1[200];
  char buffer2[200];
  char rootFileToChain[500];
  char rootFriendFileToChain[500];

  //vector< vector<Double_t> > matrix;
  vector< Double_t > yieldsRow;
  vector< Double_t > efficiencyRow;
  Int_t nSample = 0;
  std::vector<std::string> sampleName;
  sampleName.push_back("DYJetsToLL");
  sampleName.push_back("Top");
  sampleName.push_back("QCD");
  sampleName.push_back("GJets");
  sampleName.push_back("WJetsToLNu");
  sampleName.push_back("ZJetsToNuNu");
  std::vector<std::string> selectionDefinition;
  selectionDefinition.push_back("entry point");
  selectionDefinition.push_back("2 lept SF/OS");
  selectionDefinition.push_back("muons");
  selectionDefinition.push_back("tight Tag");
  selectionDefinition.push_back("mll");
  selectionDefinition.push_back("njets");
  selectionDefinition.push_back("jet1pt > 110");
  selectionDefinition.push_back("jetjetdphi");
  selectionDefinition.push_back("electron veto");
  selectionDefinition.push_back("photon veto");

  strcpy(inputFileName,argv[1]);

  ifstream *inputFile = new ifstream(inputFileName);
  while( !(inputFile->eof()) ){
    inputFile->getline(buffer1,200);
    inputFile->getline(buffer2,200);
    if ((!strstr(buffer1,"#") && !(strspn(buffer1," ") == strlen(buffer1))) && 
          (!strstr(buffer2,"#") && !(strspn(buffer2," ") == strlen(buffer2))))
      {
        sscanf(buffer1,"%s",rootFileToChain);
	sscanf(buffer2,"%s",rootFriendFileToChain);
        
	std::cout << "Creating chain ..." << std::endl;
	TChain* chain = new TChain("tree");
	chain->Add(TString(rootFileToChain));

        std::cout << "Adding friend to chain ..." << std::endl;
	TChain* chFriend = new TChain("mjvars/t");
	chain->AddFriend("mjvars/t",TString(rootFriendFileToChain));
	
	if(!chain) {
	  std::cout << "Error: chain not created. End of programme" << std::endl;
	  exit(EXIT_FAILURE);
	}
	nSample++;
	std::cout << "analysing sample n " << nSample << " : " << sampleName[nSample-1]  << std::endl; 
	std::cout<<chain->GetEntries()<<std::endl;      
	//================ Run Analysis
	//zmumujetsAna tree( chain );
	zmumujetsControlSample tree( chain , sampleName[nSample-1].c_str());
	tree.loop(yieldsRow, efficiencyRow); 
	//matrix.push_back(yieldsRow);
	//matrix.push_back(efficiencyRow);
	delete chain;
	delete chFriend;

      }
    	
  }

  if ( yieldsRow.size() != efficiencyRow.size() ) {
    std::cout << "Warning:  different number of steps for yields and efficiencies." << std::endl; 
    std::cout << "Will use the bigger one." << std::endl; 
  }
  Int_t selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/nSample) : (efficiencyRow.size()/nSample);
  FILE* fp;
  const char* finalFileName = "zmumujetsCSandBkg_yieldsEff.dat";
  if ( (fp=fopen(finalFileName,"w")) == NULL) {
    cout<<"Error: '"<<finalFileName<<"' not opened"<<endl;
  } else {
    cout<<"creating file '"<<finalFileName<<"' ..."<<endl;
    fprintf(fp,"#    step         ");
    for(Int_t i = 0; i < nSample; i++) {
      fprintf(fp,"%-16s ",sampleName[i].c_str());
    }
    fprintf(fp,"\n");
    for (Int_t i = 0; i < selectionSize; i++) {
      fprintf(fp,"%-16s",selectionDefinition[i].c_str());
      for(Int_t j = 0; j < nSample; j++) {
	if (yieldsRow.at( i + j * selectionSize) < 1000) fprintf(fp,"%7.2lf ",yieldsRow.at( i + j * selectionSize));	   
	else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));
	fprintf(fp,"%5.1lf%%   ",(100 * efficiencyRow.at( i + j * selectionSize)));
      } 
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  inputFile->close();
  delete inputFile;

#endif

#if Application == 2

  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zmumujetsAna tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif

#if Application == 3

  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zeejetsAna tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif

#if Application == 4 

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;

  char inputFileName[200];
  char buffer1[200];
  char buffer2[200];
  char rootFileToChain[500];
  char rootFriendFileToChain[500];

  //vector< vector<Double_t> > matrix;
  vector< Double_t > yieldsRow;
  vector< Double_t > efficiencyRow;
  Int_t nSample = 0;
  std::vector<std::string> sampleName;
  sampleName.push_back("DYJetsToLL");
  sampleName.push_back("Top");
  sampleName.push_back("QCD");
  //sampleName.push_back("GJets");
  sampleName.push_back("WJetsToLNu");
  sampleName.push_back("ZJetsToNuNu");
  std::vector<std::string> selectionDefinition;
  selectionDefinition.push_back("entry point");
  selectionDefinition.push_back("2 lept SF/OS");
  selectionDefinition.push_back("electrons");
  selectionDefinition.push_back("tight Tag");
  selectionDefinition.push_back("mll");
  selectionDefinition.push_back("elemet > 200");
  selectionDefinition.push_back("njets");
  selectionDefinition.push_back("jet1pt > 110");
  selectionDefinition.push_back("jetjetdphi");
  selectionDefinition.push_back("muon veto");
  selectionDefinition.push_back("photon veto");

  strcpy(inputFileName,argv[1]);

  ifstream *inputFile = new ifstream(inputFileName);
  while( !(inputFile->eof()) ){
    inputFile->getline(buffer1,200);
    inputFile->getline(buffer2,200);
    if ((!strstr(buffer1,"#") && !(strspn(buffer1," ") == strlen(buffer1))) && 
          (!strstr(buffer2,"#") && !(strspn(buffer2," ") == strlen(buffer2))))
      {
        sscanf(buffer1,"%s",rootFileToChain);
	sscanf(buffer2,"%s",rootFriendFileToChain);
        
	std::cout << "Creating chain ..." << std::endl;
	TChain* chain = new TChain("tree");
	chain->Add(TString(rootFileToChain));

        std::cout << "Adding friend to chain ..." << std::endl;
	TChain* chFriend = new TChain("mjvars/t");
	chain->AddFriend("mjvars/t",TString(rootFriendFileToChain));
	
	if(!chain) {
	  std::cout << "Error: chain not created. End of programme" << std::endl;
	  exit(EXIT_FAILURE);
	}
	nSample++;
	std::cout << "analysing sample n " << nSample << " : " << sampleName[nSample-1]  << std::endl; 
	std::cout<<chain->GetEntries()<<std::endl;      
	//================ Run Analysis
	//zeejetsAna tree( chain );
	zeejetsControlSample tree( chain , sampleName[nSample-1].c_str());
	tree.loop(yieldsRow, efficiencyRow); 
	//matrix.push_back(yieldsRow);
	//matrix.push_back(efficiencyRow);
	delete chain;
	delete chFriend;

      }
    	
  }

  if ( yieldsRow.size() != efficiencyRow.size() ) {
    std::cout << "Warning:  different number of steps for yields and efficiencies." << std::endl; 
    std::cout << "Will use the bigger one." << std::endl; 
  }
  Int_t selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/nSample) : (efficiencyRow.size()/nSample);
  FILE* fp;
  const char* finalFileName = "zeejetsCSandBkg_yieldsEff.dat";
  if ( (fp=fopen(finalFileName,"w")) == NULL) {
    cout<<"Error: '"<<finalFileName<<"' not opened"<<endl;
  } else {
    cout<<"creating file '"<<finalFileName<<"' ..."<<endl;
    fprintf(fp,"#    step         ");
    for(Int_t i = 0; i < nSample; i++) {
      fprintf(fp,"%-16s ",sampleName[i].c_str());
    }
    fprintf(fp,"\n");
    for (Int_t i = 0; i < selectionSize; i++) {
      fprintf(fp,"%-16s",selectionDefinition[i].c_str());
      for(Int_t j = 0; j < nSample; j++) {
	if (yieldsRow.at( i + j * selectionSize) < 1000) fprintf(fp,"%7.2lf ",yieldsRow.at( i + j * selectionSize));	   
	else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));
	fprintf(fp,"%5.1lf%%   ",(100 * efficiencyRow.at( i + j * selectionSize)));
      } 
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  inputFile->close();
  delete inputFile;

#endif

#if Application == 5

  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zmumujetsAna_LepSk tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif

#if Application == 6

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/evVarFriend_ZJetsToNuNu.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  znunujetsAna tree( chain );
  tree.loop(configFileName);
  
  delete chain;
  delete chFriend;

#endif

#if Application == 7

#if defined MUON

  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

#elif defined ELECTRON

  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

#endif

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zlljetsAna tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif


#if Application == 8

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  ifstream inputFile(configFileName);
  Double_t muonOrElectronOrNeutrino_PDGID;

   if (inputFile.is_open()) {

     Double_t value;
     string parameterName;
     while (inputFile >> parameterName >> value) {

       if (parameterName == "LEP_PDG_ID") {

	 muonOrElectronOrNeutrino_PDGID = (Int_t) value;
	 std::cout << "lepton_pdgID = " << value <<std::endl;

	 if (fabs(muonOrElectronOrNeutrino_PDGID) == 13) {
	   std::cout << "Analysis of Z->mumu" << std::endl;
	 } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 11) {
	   std::cout << "Analysis of Z->ee" << std::endl;
	 } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 12 || fabs(muonOrElectronOrNeutrino_PDGID) == 14 || fabs(muonOrElectronOrNeutrino_PDGID) == 16) {
	   std::cout << "Analysis of Z->nunu (any flavour)" << std::endl;
	 }

       }

     }

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   TChain* chain = NULL;
   TChain* chFriend = NULL;

   if (fabs(muonOrElectronOrNeutrino_PDGID) == 13) {

     std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
     std::cout << "Creating chain ..." << std::endl;
     chain = new TChain("tree");
     chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/tree.root");
     chFriend = new TChain("mjvars/t");
     chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

   } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 11) {

     std::cout<<"Using Emanuele's trees with skim on 2 leptons (biased Axe)"<<std::endl;
     std::cout << "Creating chain ..." << std::endl;
     chain = new TChain("tree");
     chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
     chFriend = new TChain("mjvars/t");
     chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

   } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 12 || fabs(muonOrElectronOrNeutrino_PDGID) == 14 || fabs(muonOrElectronOrNeutrino_PDGID) == 16) {

     std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
     std::cout << "Creating chain ..." << std::endl;
     chain = new TChain("tree");
     chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/tree.root");
     chFriend = new TChain("mjvars/t");
     chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/evVarFriend_ZJetsToNuNu.root");

   }

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis
  if (fabs(muonOrElectronOrNeutrino_PDGID) == 13 || fabs(muonOrElectronOrNeutrino_PDGID) == 11) {

    zlljetsAna_new tree( chain );
    tree.loop(configFileName);

  } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 12 || fabs(muonOrElectronOrNeutrino_PDGID) == 14 || fabs(muonOrElectronOrNeutrino_PDGID) == 16) {

    znunujetsAna tree( chain );
    tree.loop(configFileName);

  }

  delete chain;
  delete chFriend;

#endif

#if Application == 9

  

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  Double_t lepton_PDGID;
  Int_t isdata_flag;
  Int_t tau_veto_flag;
  string treePath;
  string friendTreePath;
  string option = "";

  Int_t unweighted_event_flag = 0;  // this flag tells the user if the MC uses unit weight (using w = 1 is basically for debugging purposes)
  Int_t adishTree_flag = 0;  // tells user I'm using Adish's tree

  // following are for CS analysis
  Int_t controlSample_flag = 0;
  Int_t signalRegion_flag = 0;
  string fileWithSamplesPath = "";
  string filename_base = "";

  if (argc > 2 ) {

    for (Int_t i = 2; i < argc; i++) {   // look at all possible options passed

      string thisArgument(argv[i]);

      if (thisArgument  == "-nw" ) {
	
	unweighted_event_flag = 1;    //-nw option stands for "no weight"
	cout << "Option " << thisArgument << " passed: using no event weight (w =1) if allowed" << std::endl;

      }

      if (thisArgument  == "-at" ) {   // "at" means Adish's tree
	     
	adishTree_flag = 1;
	cout << "Option " << thisArgument << " passed: using Adish's trees" << std::endl;

      }

    }

  }

  ifstream inputFile(configFileName);

  if (inputFile.is_open()) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    while (inputFile >> parameterType ) {  // read only first object  here

      if (parameterType == "NUMBER") {

	inputFile >> parameterName >> value;  

	if (parameterName == "LEP_PDG_ID") {

	  lepton_PDGID = (Int_t) value;
	  std::cout << "lepton_pdgID = " << value <<std::endl;

	  if (fabs(lepton_PDGID) == 13) {
	    std::cout << "Analysis of Z->mumu" << std::endl;
	  } else if (fabs(lepton_PDGID) == 11) {
	    std::cout << "Analysis of Z->ee" << std::endl;
	  } else if (fabs(lepton_PDGID) == 12 || fabs(lepton_PDGID) == 14 || fabs(lepton_PDGID) == 16) {
	    std::cout << "Analysis of Z->nunu" << std::endl;
	  }

	} 

	if (parameterName == "ISDATA_FLAG") {

	  isdata_flag = (Int_t) value;

	  if (isdata_flag == 0) std::cout << "Running on MonteCarlo" << std::endl;
	  else std::cout << "Running on data" << std::endl;

	}

	if (parameterName == "TAU_VETO_FLAG") {

	  tau_veto_flag = (Int_t) value;

	  if (tau_veto_flag == 0) std::cout << "Not applying tau veto" << std::endl;
	  else std::cout << "Applying tau veto" << std::endl;

	}

      } else if (parameterType == "STRING") {

	inputFile >> parameterName >> name;

	if (parameterName == "TREE_PATH") {

	  treePath = name;
	  std::cout << setw(20) << "tree : " << treePath <<std::endl;

	}  

	if (parameterName == "FRIEND_TREE_PATH") {

	  friendTreePath = name;
	  std::cout << setw(20) << "friend tree : " << friendTreePath <<std::endl;

	}  

	if (parameterName == "OPTION") {

	  option = name;
	  std::cout << setw(20) << "option : " << option <<std::endl;
  

	} 

	if (parameterName == "FILENAME_BASE") {  // path to file with CS samples and some options

	  filename_base = name;
	  std::cout << setw(20) << "filename_base: " << filename_base <<std::endl;

	} 

	if (parameterName == "PATH_TO_SAMPLES_4CS") {  // path to file with CS samples and some options

	  controlSample_flag = 1;
	  fileWithSamplesPath= name;
	  std::cout << "Performing analysis on control samples." <<std::endl;
	  std::cout << setw(20) << "file pointing to CS samples: " << fileWithSamplesPath<<std::endl;

	} 

	if (parameterName == "PATH_TO_SAMPLES_4SR") {  // path to file with CS samples and some options

	  signalRegion_flag = 1;
	  fileWithSamplesPath= name;
	  std::cout << "Performing analysis on signal region." <<std::endl;
	  std::cout << setw(20) << "file pointing to SR samples: " << fileWithSamplesPath<<std::endl;

	} 

      }

    }

    inputFile.close();
                                                                                                                         
  } else {

    cout << "Error: could not open file " << configFileName << endl;
    exit(EXIT_FAILURE);

  }

  if (controlSample_flag == 1 || signalRegion_flag == 1) {

    // char buffer1[200];
    // char buffer2[200];
    // char rootFileToChain[500];
    // char rootFriendFileToChain[500];

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    //vector< vector<Double_t> > matrix;
    vector< Double_t > yieldsRow;
    vector< Double_t > efficiencyRow;
    Int_t nSample = 0;
    std::vector<std::string> sampleName;
    
    std::vector<std::string> selectionDefinition;

    if (signalRegion_flag == 1) {
    
      selectionDefinition.push_back("entry point");  // include MetNoLep
      selectionDefinition.push_back("jet1pt");
      selectionDefinition.push_back("jetjetdphi");
      selectionDefinition.push_back("njets");
      selectionDefinition.push_back("muon veto");
      selectionDefinition.push_back("electron veto");
      if (tau_veto_flag) selectionDefinition.push_back("tau veto");
      selectionDefinition.push_back("photon veto");

    } else if (controlSample_flag == 1) {

      selectionDefinition.push_back("entry point");        
      selectionDefinition.push_back("preselection");   // include genLep, HLT, MetNoLep
      selectionDefinition.push_back("2lep SF/OS");
      selectionDefinition.push_back("2lep loose");
      if (fabs(lepton_PDGID) == 13) selectionDefinition.push_back("muons");
      else if (fabs(lepton_PDGID) == 11) selectionDefinition.push_back("electrons");
      selectionDefinition.push_back("tight Tag");
      selectionDefinition.push_back("mll");
      selectionDefinition.push_back("jet1pt");
      selectionDefinition.push_back("jetjetdphi");
      selectionDefinition.push_back("njets");
      if (fabs(lepton_PDGID) == 13) selectionDefinition.push_back("electron veto");
      else if (fabs(lepton_PDGID) == 11) selectionDefinition.push_back("muon veto");
      selectionDefinition.push_back("photon veto");
      if (tau_veto_flag) selectionDefinition.push_back("tau veto");
      selectionDefinition.push_back("lep match");

    }

    ifstream sampleFile(fileWithSamplesPath.c_str());
    Int_t fileEndReached_flag = 0;

    if (sampleFile.is_open()) {

      while (fileEndReached_flag == 0) {

	//============================================/
  
	while ( (sampleFile >> parameterType) && (!(parameterType == "#")) ) {  // read only first object  here: if it is '#' it signal that another part is starting

	  if (parameterType == "#STOP") fileEndReached_flag = 1;   //tells me that the file is ended and I don't need to go on reading it.

	  if (parameterType == "NUMBER") {

	    sampleFile >> parameterName >> value;  

	    if (parameterName == "ISDATA_FLAG") {

	      isdata_flag = (Int_t) value;

	      if (isdata_flag == 0) std::cout << "Running on MonteCarlo" << std::endl;
	      else std::cout << "Running on data" << std::endl;

	    }

	  } else if (parameterType == "STRING") {

	    sampleFile >> parameterName >> name;

	    if (parameterName == "SAMPLE_NAME") {

	      sampleName.push_back(name.c_str());
	      std::cout << setw(20) << "sample name : " << name <<std::endl;

	    }  

	    if (parameterName == "TREE_PATH") {

	      treePath = name;
	      std::cout << setw(20) << "tree : " << treePath <<std::endl;

	    }   

	    if (parameterName == "FRIEND_TREE_PATH") {

	      friendTreePath = name;
	      std::cout << setw(20) << "friend tree : " << friendTreePath <<std::endl;

	    }  

	  }

	}

	if ( !fileEndReached_flag) {

	  std::cout << "Creating chain ..." << std::endl;
	  TChain* chain = new TChain("tree");
	  chain->Add(TString(treePath.c_str()));

	  std::cout << "Adding friend to chain ..." << std::endl;
	  TChain* chFriend = new TChain("mjvars/t");
	  chain->AddFriend("mjvars/t",TString(friendTreePath.c_str()));
	
	  if(!chain) {
	    std::cout << "Error: chain not created. End of programme" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  std::cout << "analysing sample n " << (nSample+1) << " : " << sampleName[nSample]  << std::endl; 
	  std::cout<<chain->GetEntries()<<std::endl;      
	  //================ Run Analysis
	  //zmumujetsAna tree( chain );
 
	  if (signalRegion_flag == 1) {

	    monojet_SignalRegion tree( chain , sampleName[nSample].c_str());
	    tree.loop(configFileName, isdata_flag, unweighted_event_flag, yieldsRow, efficiencyRow); 

	  } else if (controlSample_flag == 1) {

	    zlljetsControlSample tree( chain , sampleName[nSample].c_str());
	    tree.loop(configFileName, isdata_flag, unweighted_event_flag, yieldsRow, efficiencyRow); 

	  }

	  nSample++;
	  delete chain;
	  delete chFriend;

	}

	//============================================/
 
      }  //end of while(fileEndReanched_flag == 0)

      sampleFile.close();

    } else {   // end of if (sampleFile.is_open())

      cout << "Error: could not open file " << fileWithSamplesPath<< endl;
      exit(EXIT_FAILURE);

    }

    if ( yieldsRow.size() != efficiencyRow.size() ) {
      std::cout << "Warning:  different number of steps for yields and efficiencies." << std::endl; 
      std::cout << "Will use the bigger one." << std::endl; 
    }

    Int_t selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/nSample) : (efficiencyRow.size()/nSample);

    FILE* fp;
    string finalFileName = filename_base;
    finalFileName += "_yieldsTable.dat";
    if (unweighted_event_flag) finalFileName += "_weq1";   //means with weights equal to 1 (for debugging purposes)

    if ( (fp=fopen(finalFileName.c_str(),"w")) == NULL) {
      cout<<"Error: '"<<finalFileName<<"' not opened"<<endl;
    } else {
      cout<<"creating file '"<<finalFileName<<"' to save table with yields ..."<<endl;
      fprintf(fp,"#    step         ");
      for(Int_t i = 0; i < nSample; i++) {
	fprintf(fp,"%-16s ",sampleName[i].c_str());
      }
      fprintf(fp,"\n");
      for (Int_t i = 0; i < selectionSize; i++) {
	fprintf(fp,"%-16s",selectionDefinition[i].c_str());
	for(Int_t j = 0; j < nSample; j++) {
	  
	  if (yieldsRow.at( i + j * selectionSize) < 0) {

	    string space = "//";
	    fprintf(fp,"%7s ",space.c_str());
	    fprintf(fp,"%5s    ",space.c_str());

	  } else { 

	    if (yieldsRow.at( i + j * selectionSize) < 10) fprintf(fp,"%7.1lf ",yieldsRow.at( i + j * selectionSize));  //	j * selectionSize refers to number for a sample, i refers to the selection step   
	    else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));
	    fprintf(fp,"%5.1lf%%   ",(100 * efficiencyRow.at( i + j * selectionSize)));

	  }
	} 
	fprintf(fp,"\n");
      }
      fclose(fp);
    }

    return 0;

  }  //end of  "if (controlSample_flag == 1 || signalRegion_flag == 1)"


  //if (!isdata_flag && unweighted_event_flag) std::cout << "Using unweighted events (w = 1) " << std::endl;
  //std::cout<<"Using Emanuele's trees with no skim"<<std::endl;
  std::cout<< std::endl;
  if (adishTree_flag) std::cout<<"Using Adish's trees "<<std::endl;
  else std::cout<<"Using Emanuele's trees "<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain;
  if (adishTree_flag) chain = new TChain("tree/tree");
  else chain = new TChain("tree");
  //chain->Add("/cmshome/ciprianim/edimarcoTree/tree_noSkim/DYJetsToLL_M50/tree.root");
  chain->Add(treePath.c_str());
  TChain* chFriend = new TChain("mjvars/t");
  //chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_noSkim/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");
  if (!adishTree_flag) chain->AddFriend("mjvars/t",friendTreePath.c_str());

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  if (adishTree_flag) {

    zlljets_resoResp_forAdish tree( chain );
    tree.loop(configFileName, isdata_flag, unweighted_event_flag);

  } else {

    if ( !(std::strcmp("axel",option.c_str()))) {

      zlljets_Axe_noSkim_light tree( chain );
      tree.loop(configFileName);

    } else if ( !(std::strcmp("rerel",option.c_str()))) {

      zlljets_resoResp tree( chain );
      tree.loop(configFileName, isdata_flag, unweighted_event_flag);

    } else if (!(std::strcmp("ana",option.c_str())) ) {

      zlljetsAna_new tree( chain );
      tree.loop(configFileName, isdata_flag);

    } else if (!(std::strcmp("znunuMC",option.c_str())) ) {

      znunujetsAna tree( chain );
      tree.loop(configFileName);

    } else {

      std::cout << "Option '" << option << "' not found. End of programme" << std::endl;
      return 0;

    }

  }
  
  delete chain;
  delete chFriend;

#endif


  return 0;
}

