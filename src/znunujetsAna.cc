#define znunujetsAna_cxx
#include "EmanTreeAnalysis.h"
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TVirtualFitter.h>
//my headers
#include "functionsForAnalysis.h"
#include "myClasses.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef znunujetsAna_cxx

znunujetsAna::znunujetsAna(TTree *tree) : edimarcoTree_v2(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void znunujetsAna::loop(const char* configFileName)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection for muon veto
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean1",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("jetclean2",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJetClean30",1);    // # of jet with pt > 30 & eta < 2.5 and cleaning for against muons misidentified as PFjets   
   fChain->SetBranchStatus("JetClean_pt",1);  
   fChain->SetBranchStatus("JetClean_eta",1); 

   // fChain->SetBranchStatus("met_pt",1);
   // fChain->SetBranchStatus("met_phi",1);
   fChain->SetBranchStatus("metNoMu_pt",1);   // likely this will coincide with the pt of the Z(nunu)
   
   char ROOT_FNAME[100];
   char TXT_FNAME[100];
   char TEX_FNAME[100]; 
   
   Double_t LUMI;
   Int_t NJETS;
   Double_t J1PT;
   Double_t J1ETA;
   Double_t J2PT;
   Double_t J2ETA;
   Double_t J1J2DPHI;
   Int_t LEP_PDG_ID;
   Int_t TAU_VETO_FLAG;
   Double_t METNOLEP_START;
   string FILENAME_BASE;

   ifstream inputFile(configFileName);

   if (inputFile.is_open()) {

     Double_t value;
     string name;
     string parameterName;
     string parameterType;

     mySpaces(cout,2);
     cout << "Printing content of " << configFileName << " file" << endl;
     mySpaces(cout,1);

     while (inputFile >> parameterType ) {

       if (parameterType == "NUMBER") {

	 inputFile >> parameterName >> value;
	 cout << setw(20) << parameterName << setw(7) << value << endl;

	 if (parameterName == "LUMI") LUMI = value;
	 else if (parameterName == "NJETS") NJETS = value;
	 else if (parameterName == "J1PT") J1PT = value;
	 else if (parameterName == "J1ETA") J1ETA = value;
	 else if (parameterName == "J2PT") J2PT = value;
	 else if (parameterName == "J2ETA") J2ETA = value;
	 else if (parameterName == "J1J2DPHI") J1J2DPHI = value;
	 else if (parameterName == "LEP_PDG_ID") LEP_PDG_ID = value;
	 else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = value;
	 else if (parameterName == "METNOLEP_START") METNOLEP_START = value;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 cout << right << setw(20) << parameterName << "  " << left << name << endl;
	 if (parameterName == "FILENAME_BASE") FILENAME_BASE = name; 

       }

     }
     
     mySpaces(cout,2);

     strcpy(ROOT_FNAME,(FILENAME_BASE + ".root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + ".txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + ".tex").c_str());

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   selection jet1C("jet1C",Form("jet1pt > %4.0lf",(Double_t)J1PT),Form("nJetClean30 >= 1 && JetClean1_pt > %4.0lf && abs(JetClean1_eta) < %1.1lf && jetclean1 > 0.5",(Double_t)J1PT,J1ETA));
   selection jjdphiC("jjdphiC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection njetsC("njets","nJetClean30 <= 2");
   selection eLooseVetoC("eLooseVetoC","electrons veto");
   selection muLooseVetoC("muLooseVeto","muons veto");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC("tauLooseVetoC","taus veto");
   selection metNoMuStartC("metNoMuStartC",Form("metNoMu > %2.0lf",METNOLEP_START));
   
   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   mask znunujetsMonojetSelection("monojet selection on Z->nunu sample with selection flow as Emanuele's");

   if (METNOLEP_START) znunujetsMonojetSelection.append(metNoMuStartC.get2ToId());
   znunujetsMonojetSelection.append(jet1C.get2ToId());
   znunujetsMonojetSelection.append(jjdphiC.get2ToId());
   znunujetsMonojetSelection.append(njetsC.get2ToId());
   znunujetsMonojetSelection.append(eLooseVetoC.get2ToId());
   znunujetsMonojetSelection.append(muLooseVetoC.get2ToId());
   if (TAU_VETO_FLAG) znunujetsMonojetSelection.append(tauLooseVetoC.get2ToId());
   znunujetsMonojetSelection.append(gammaLooseVetoC.get2ToId());

   cout << "Opening file " <<ROOT_FNAME<< endl;

   TFile *rootFile = new TFile(ROOT_FNAME,"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }

   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin","yields of Z#nu#nu control sample in bins of met;#slash{E}_{T};# of events",nMetBins,metBinEdges);

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"znunujetsAna::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     nTotalWeightedEvents += newwgt;  // counting events with weights

     eventMask +=jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT && fabs(JetClean_eta[0] < J1ETA && jetclean1 > 0.5));
     eventMask += jjdphiC.addToMask( nJetClean30 == 1 || (nJetClean30 >= NJETS && fabs(dphijj) < J1J2DPHI && jetclean2 > 0.5));
     eventMask += njetsC.addToMask(nJetClean30 <= NJETS);           
     eventMask += eLooseVetoC.addToMask(nEle10V == 0);
     eventMask += muLooseVetoC.addToMask(nMu10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += metNoMuStartC.addToMask(metNoMu_pt > METNOLEP_START);
     
     znunujetsMonojetSelection.countEvents(eventMask, newwgt);
     

     // filling histogram with yields at the end of the selection in bins of met
     if ( ((eventMask & znunujetsMonojetSelection.globalMask.back()) == znunujetsMonojetSelection.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBin->Fill(metNoMu_pt,newwgt);     
     }
     
   }

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelection);
   cout<<endl;   

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBin, metBinEdges, nMetBins);

   cout<<"creating file '"<<TXT_FNAME<<"' ..."<<endl;
   ofstream myfile(TXT_FNAME,ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }
            
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelection);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBin, metBinEdges, nMetBins);

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   FILE *fp;
   fp = fopen(TEX_FNAME,"w");

   if ( fp == NULL)  cout<<"Error: '"<<TEX_FNAME<<"' not opened"<<endl;
   else {

     cout<<"creating file '"<<TEX_FNAME<<"' ..."<<endl;
     myAddDefaultPackages(fp,TEX_FNAME);
     fprintf(fp,"\\begin{document}\n");
     fprintf(fp,"\n");
     string commentInTable;       
     //makeTableTex(fp, LUMI, nTotalWeightedEvents, &mu_Acc_Eff, commentInTable);      
     commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelection,commentInTable);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file


}
