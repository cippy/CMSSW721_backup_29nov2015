#define monojet_SignalRegion_cxx
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
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
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
#include <TVector3.h>
#include <TVirtualFitter.h>
//my headers
#include "functionsForAnalysis.h"
#include "myClasses.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef monojet_SignalRegion_cxx

monojet_SignalRegion::monojet_SignalRegion(TTree *tree, const char* inputSuffix) : edimarcoTree_v2(tree) {
  //cout <<"check in constructor "<<endl;
  suffix = inputSuffix;  // it is the sample name (e.g. QCD, ZJetsToNuNu ecc...)
  Init(tree);

}

#endif

void monojet_SignalRegion::loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag, vector< Double_t > &yRow, vector< Double_t > &eRow, vector< Double_t > &uncRow)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor
   //fChain->SetBranchStatus("LHEorigWeight",1); // contains negative values: the weight in the event is weight*LHEorigWeight

   fChain->SetBranchStatus("genWeight",1); 

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   //fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
   //fChain->SetBranchStatus("nTau18V",1);
   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean1",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("jetclean2",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJetClean30",1);    // # of jet with pt > 30 & eta < 2.5 and cleaning for against muons misidentified as PFjets   
   fChain->SetBranchStatus("JetClean_pt",1);  
   fChain->SetBranchStatus("JetClean_eta",1);  
   // fChain->SetBranchStatus("nJet",1);         // # of jets with pt > 25 && |eta| < 2.5
   // fChain->SetBranchStatus("nJet30",1);         // # of jets with pt > 30 && |eta| < 2.4
   // fChain->SetBranchStatus("nJet30a",1);       // # of jets with pt > 30 && |eta| < 4.7 
   // fChain->SetBranchStatus("Jet_pt",1);  
   // fChain->SetBranchStatus("Jet_eta",1);  
 
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_phi",1);   
   fChain->SetBranchStatus("LepGood_mass",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   // fChain->SetBranchStatus("ngenLep",1);         // not needed, using GenPart to study generator level quantities
   // fChain->SetBranchStatus("genLep_pdgId",1);
   // fChain->SetBranchStatus("genLep_pt",1);
   // fChain->SetBranchStatus("genLep_eta",1);
   // fChain->SetBranchStatus("genLep_phi",1);
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   if (!ISDATA_FLAG) {
     // fChain->SetBranchStatus("nGenPart",1);
     // fChain->SetBranchStatus("GenPart_pdgId",1);
     // fChain->SetBranchStatus("GenPart_motherId",1);
     // fChain->SetBranchStatus("GenPart_pt",1);
     // fChain->SetBranchStatus("GenPart_eta",1);
     // fChain->SetBranchStatus("GenPart_phi",1);
     // fChain->SetBranchStatus("GenPart_mass",1);
     // fChain->SetBranchStatus("GenPart_motherIndex",1);

     fChain->SetBranchStatus("vtxW",1);   // weight to have better agreement between data and MC (will not be always used)
     fChain->SetBranchStatus("xsec",1);   // weight to have better agreement between data and MC
   }

   //fChain->SetBranchStatus("met_pt",1);
   //fChain->SetBranchStatus("met_eta",1);
   //fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   //fChain->SetBranchStatus("metNoMu_eta",1);
   fChain->SetBranchStatus("metNoMu_phi",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices 

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
   Int_t TAU_VETO_FLAG;
   // Int_t HLT_FLAG;                  // not needed: monojet default trigger paths should be 100% efficient wrt offline selection
   // Double_t HLT_LEP1PT;    
   // Double_t HLT_LEP2PT;
   // Double_t HLT_LEP1ETA;
   // Double_t HLT_LEP2ETA;
   Double_t METNOLEP_START;
   string FILENAME_BASE;
   string DIRECTORY_TO_SAVE_FILES;
   string DIRECTORY_NAME;

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
	 cout << right << setw(20) << parameterName << "  " << left << value << endl;

	 if (parameterName == "LUMI") LUMI = value;
	 else if (parameterName == "NJETS") NJETS = value;
	 else if (parameterName == "J1PT") J1PT = value;
	 else if (parameterName == "J1ETA") J1ETA = value;
	 else if (parameterName == "J2PT") J2PT = value;
	 else if (parameterName == "J2ETA") J2ETA = value;
	 else if (parameterName == "J1J2DPHI") J1J2DPHI = value;
	 else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = value;
	 else if (parameterName == "METNOLEP_START") METNOLEP_START = value;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 cout << right << setw(20) << parameterName << "  " << left << name << endl;

	 if (parameterName == "FILENAME_BASE") {

	   FILENAME_BASE = name; 
	   if ( !ISDATA_FLAG && unweighted_event_flag) FILENAME_BASE += "_weq1";  // if using unit weight, add _weq1 to filename (weq1 means weight = 1)

	 }

	 if (parameterName == "DIRECTORY_PATH") {  // name of directory where files are saved

	  DIRECTORY_TO_SAVE_FILES = name;
	  //std::cout << "Files will be saved in '" << name << "' ." <<std::endl;

	} 

	if (parameterName == "DIRECTORY_NAME") {  // name of directory where files are saved

	  DIRECTORY_NAME = name;
	  //std::cout << "Files will be saved in directory named '" << name << "' ." <<std::endl;

	} 

       }

     }
     
     mySpaces(cout,2);

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   string outputFolder =  DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME + "/";

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   // vector<Double_t> metCut;
   // metCut.push_back(250);
   // metCut.push_back(300);
   // metCut.push_back(350);
   // metCut.push_back(400);
   // metCut.push_back(500);

   selection metNoMuC("metNoMuC",Form("metNoMu > %4.0lf",METNOLEP_START),"first cut on met");
   selection jet1C("jet1C",Form("jet1pt > %4.0lf",(Double_t)J1PT),Form("nJetClean >= 1 && JetClean1_pt > %4.0lf && abs(JetClean1_eta) < %1.1lf && jetclean1 > 0.5",(Double_t)J1PT,J1ETA));
   selection jjdphiC("jjdphiC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection njetsC("njets","nJetClean30 <= 2");
   selection muonLooseVetoC("muonLooseVetoC","muonss veto");; 
   selection electronLooseVetoC("electronLooseVetoC","electrons veto");; 
   selection tauLooseVetoC;
   if (TAU_VETO_FLAG) tauLooseVetoC.set("tauLooseVetoC","tau veto");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");


   Double_t SUMWEIGHTS;   // ==============  To be initialized with proper value to compute event weight in MC ========================
   vector<Double_t> sumWeightVector;
   vector<Int_t> eventsInSubsamples;

   Double_t nTotalWeightedEvents = 0.0;     

   Int_t using_phys14_sample_flag = 0;
   if (FILENAME_BASE.find("phys14") != std::string::npos) {
     using_phys14_sample_flag = 1;    
     cout << "Using phys14 samples" << endl;
   }

   Int_t using_spring15_sample_flag = 0;
   if (FILENAME_BASE.find("spring15") != std::string::npos) {
     using_spring15_sample_flag = 1;    
     cout << "Using spring15 samples" << endl;
   }

   Int_t using_spring15_25ns_sample_flag = 0;
   if (FILENAME_BASE.find("spring15_25ns") != std::string::npos) {
     using_spring15_25ns_sample_flag = 1;    
     cout << "Using spring15_25ns samples" << endl;
     mySumWeight_filler_spring15_25ns(suffix, sumWeightVector);  // this function fills the vector with the proper values of sumWeight depending on the sample
     myEventsInSubsamples_filler_spring15_25ns(suffix, eventsInSubsamples); 
   }

   if ( !ISDATA_FLAG && unweighted_event_flag) cout << "Warning: no weight applied to events (w = 1)" << endl;  // if MC with unit weight, make user know

   if (ISDATA_FLAG) {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_DATA.tex").c_str());
   } else {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_" + suffix + ".root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_" + suffix + ".txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_" + suffix + ".tex").c_str());
   }

   selection::checkMaskLength();
   selection::printActiveSelections(cout); 

   UInt_t maskJetsSelection = njetsC.get2ToId() + jet1C.get2ToId() + jjdphiC.get2ToId();

   UInt_t maskMonoJetSelection = maskJetsSelection + electronLooseVetoC.get2ToId() + muonLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   if ( TAU_VETO_FLAG ) maskMonoJetSelection += tauLooseVetoC.get2ToId();

   mask monojet_SignalRegion("monojet signal selection");

   if (METNOLEP_START) monojet_SignalRegion.append(metNoMuC.get2ToId());
   
   monojet_SignalRegion.append(jet1C.get2ToId());
   monojet_SignalRegion.append(jjdphiC.get2ToId());
   monojet_SignalRegion.append(njetsC.get2ToId());
   monojet_SignalRegion.append(muonLooseVetoC.get2ToId());
   monojet_SignalRegion.append(electronLooseVetoC.get2ToId());
   if (TAU_VETO_FLAG) monojet_SignalRegion.append(tauLooseVetoC.get2ToId());
   monojet_SignalRegion.append(gammaLooseVetoC.get2ToId());

   cout << "Opening file " <<ROOT_FNAME<< endl;

   TFile *rootFile = new TFile((outputFolder + ROOT_FNAME).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<(outputFolder + ROOT_FNAME).c_str()<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }
 

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   //Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};       

   TH1D *HYieldsMetBin = new TH1D("HYieldsMetBin","monojet signal region's yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdges);
   
   TH1D *HvtxDistribution = new TH1D("HvtxDistribution","",40,-0.5,39.5);   
   TH1D *HnjetsDistribution = new TH1D("HnjetsDistribution","njets using nJetClean30",10,-0.5,9.5);   
   TH1D *Hj1j2dphiDistribution = new TH1D("Hj1j2dphiDistribution","",30,0.0,3.0);
   TH1D *Hjet1etaDistribution = new TH1D("Hjet1etaDistribution","",60,-3.0,3.0);
   TH1D *Hjet2etaDistribution = new TH1D("Hjet2etaDistribution","",60,-3.0,3.0);
   TH1D *HmetNoLepDistribution;
   TH1D *Hjet1ptDistribution;  
   TH1D *Hjet2ptDistribution;

   //if (using_phys14_sample_flag) {
   HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",100,0.0,1000.0);
   Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",100,0.0,1000.0); 
   Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",100,0.0,1000.0);
   // } else if (using_spring15_sample_flag) {
   //   HmetNoLepDistribution = new TH1D("HmetNoLepDistribution","",60,METNOLEP_START,METNOLEP_START+600);
   //   Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",60,J1PT,J1PT+600); 
   //   Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",60,J2PT,J2PT+600);
   // }

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   // deciding  what is the event weight
   Double_t newwgt;
   Double_t eventCounter = 0;  // support variable: at the beginning it is set to the number of entries of the first subsample (e.g. HT100to200 or whatever): when the number of 
                                                 // event analyzed reaches this value, it's increased by the number of entries in the following subsample and so on. Basically, it's needed to keep track of
                                                 // the specific subsample that is being analyzed (so that the proper value of sumWeight is used)
   Int_t htbin = 0;  // 

   if (ISDATA_FLAG || unweighted_event_flag) newwgt = 1.0;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"monojet_SignalRegion::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 

     if (jentry%500000 == 0) {
       cout << "entry: " << jentry << endl;
     }

     if(!ISDATA_FLAG && !unweighted_event_flag) {

       // the following if statement is used to set the proper value of sumWeight, which changes depending on the HTbin or on the specific subsample.
       // at first eventCounter = 0 so when the loop on events begins the condition is fulfilled:
       // eventCounter is increased by the # of entries in the first subsample, SUMWEIGHTS is set and the htbin index is increased by 1.
       // the if condition will be fulfilled again when the next subsample starts being analyzed (note that jentry starts from 0, not from 1)

       if (using_spring15_25ns_sample_flag) {

	 if (jentry == eventCounter) {
	   eventCounter += eventsInSubsamples[htbin];
	   SUMWEIGHTS = sumWeightVector[htbin];
	   htbin++;
	   cout << endl;
	   cout << "entry = " << jentry << ":   " ;
	   cout << "htbin = " << htbin << "  --->  ";   // it will print 1, 2, 3 ... but as an index it would be 0, 1, 2 ...
	   cout << "sumWeight = " << SUMWEIGHTS << endl;
	   cout << endl;
	 }

	 newwgt = 1000 * LUMI * xsec * genWeight / SUMWEIGHTS; 

       } else if (using_spring15_sample_flag == 1 && using_spring15_25ns_sample_flag == 0) newwgt = 1000 * LUMI * vtxW * xsec * genWeight / SUMWEIGHTS;    
       // 1000 is because LUMI is in fb^-1 and xsec is in pb
       // old wrong one:     newwgt = LUMI * vtxW * weight * LHEorigWeight; 
       else if (using_phys14_sample_flag) newwgt = LUMI * weight;   // for older trees (backward compatibility)
       else newwgt = LUMI * weight;   // for older trees (backward compatibility)

     }

     nTotalWeightedEvents += newwgt;  // counting events with weights

     // beginning of eventMask building
     
     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT && fabs(JetClean_eta[0] < J1ETA && jetclean1 > 0.5));  //could skip cut on eta 
     eventMask += jjdphiC.addToMask( nJetClean30 == 1 || (nJetClean30 >= NJETS && fabs(dphijj) < J1J2DPHI && jetclean2 > 0.5));
     eventMask += njetsC.addToMask(nJetClean30 <= NJETS);
     eventMask += muonLooseVetoC.addToMask(nMu10V == 0);
     eventMask += electronLooseVetoC.addToMask(nEle10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += metNoMuC.addToMask(metNoMu_pt > METNOLEP_START);

     // for (Int_t i = 0; i <  metCut.size(); i++) {
     //   eventMask += metNoLepC[i].addToMask(metNoMu_pt > metCut[i]);
     // }
     
     // end of eventMask building

     monojet_SignalRegion.countEvents(eventMask,newwgt);

     if ( ((eventMask & monojet_SignalRegion.globalMask.back()) == monojet_SignalRegion.globalMask.back()) ) {
       
       // this histogram holds the final yields in bins of MET
	 HYieldsMetBin->Fill(metNoMu_pt,newwgt);
	 
	 HmetNoLepDistribution->Fill(metNoMu_pt,newwgt);
	 HvtxDistribution->Fill(nVert,newwgt);
	 HnjetsDistribution->Fill(nJetClean30,newwgt);
	 Hjet1etaDistribution->Fill(JetClean_eta[0],newwgt);
	 Hjet1ptDistribution->Fill(JetClean_pt[0],newwgt);
	 if (nJetClean30 == 2) {
	   Hj1j2dphiDistribution->Fill(dphijj,newwgt);
	   Hjet2etaDistribution->Fill(JetClean_eta[1],newwgt);
	   Hjet2ptDistribution->Fill(JetClean_pt[1],newwgt);
	 }

     }
       
   }                        // end of loop on entries

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &monojet_SignalRegion);

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HYieldsMetBin, metBinEdges, nMetBins);
 
   cout<<"creating file '"<<TXT_FNAME<<"' in folder "<< outputFolder <<" ..."<<endl;
   ofstream myfile((outputFolder + TXT_FNAME).c_str(),ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }

   //opening inputFile named configFileName again to save content in myfile named TXT_FNAME

   inputFile.open(configFileName);

   if (inputFile.is_open()) {
     
     mySpaces(myfile,2);
     cout << "Saving content of " << configFileName << " file in "<< TXT_FNAME << endl;
     myfile << "Content of " << configFileName << endl;
     mySpaces(myfile,1);

     Double_t value;
     string name;
     string parameterName;
     string parameterType;

     while (inputFile >> parameterType ) {

       if (parameterType == "NUMBER") {

	 inputFile >> parameterName >> value;
	 myfile << setw(20) << parameterName << setw(7) << value << endl;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 myfile << right << setw(20) << parameterName << "  " << left << name << endl;

       }

     }
     
     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << " to save content in "<< TXT_FNAME << endl;
     exit(EXIT_FAILURE);

   }

   mySpaces(myfile,2);
   if (!ISDATA_FLAG && unweighted_event_flag) myfile << "======   Using unweighted events (w = 1)   ======" << endl;
   mySpaces(myfile,3);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &monojet_SignalRegion);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HYieldsMetBin, metBinEdges, nMetBins);

   myfile.close();

   // filling with yields and efficiency: I will use efficiency with respect to total and not to previous step, but I could make this choice in the config file

   // entry point
   yRow.push_back(nTotalWeightedEvents);
   eRow.push_back(1.0000);
   uncRow.push_back(sqrt(nTotalWeightedEvents));
   
   vector<Int_t> selStep;   //array to store index of step to form selection flow (might want to consider two or more steps together and not separated)

   selStep.push_back(monojet_SignalRegion.whichStepHas(metNoMuC.get2ToId()));
   selStep.push_back(monojet_SignalRegion.whichStepHas(jet1C.get2ToId()));
   selStep.push_back(monojet_SignalRegion.whichStepHas(jjdphiC.get2ToId()));
   selStep.push_back(monojet_SignalRegion.whichStepHas(njetsC.get2ToId()));
   selStep.push_back(monojet_SignalRegion.whichStepHas(muonLooseVetoC.get2ToId()));
   selStep.push_back(monojet_SignalRegion.whichStepHas(electronLooseVetoC.get2ToId()));
   if (TAU_VETO_FLAG) selStep.push_back(monojet_SignalRegion.whichStepHas(tauLooseVetoC.get2ToId()));
   selStep.push_back(monojet_SignalRegion.whichStepHas(gammaLooseVetoC.get2ToId()));

   for(Int_t i = 0; i < selStep.size(); i++) {
  
     yRow.push_back(monojet_SignalRegion.nEvents[selStep[i]]);
     uncRow.push_back(sqrt(yRow.back()));
     if (i == 0) eRow.push_back(monojet_SignalRegion.nEvents[selStep[i]]/nTotalWeightedEvents);
     else if( (i != 0) && (monojet_SignalRegion.nEvents[selStep[i]-1] == 0) ) eRow.push_back(1.0000);
     else eRow.push_back(monojet_SignalRegion.nEvents[selStep[i]]/monojet_SignalRegion.nEvents[selStep[i]-1]);

   }

   // filling last bin with overflow
   myAddOverflowInLastBin(HmetNoLepDistribution);
   myAddOverflowInLastBin(Hjet1ptDistribution);
   myAddOverflowInLastBin(Hjet2ptDistribution);

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   FILE *fp;
   fp = fopen((outputFolder + TEX_FNAME).c_str(),"w");

   if ( fp == NULL)  cout<<"Error: '"<<TEX_FNAME<<"' not opened"<<endl;
   else {

     cout<<"creating file '"<<TEX_FNAME<<" in folder " << outputFolder << "' ..."<<endl;
     myAddDefaultPackages(fp,TEX_FNAME);
     fprintf(fp,"\\begin{document}\n");
     fprintf(fp,"\n");
     string commentInTable;       
     commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &monojet_SignalRegion,commentInTable);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file

   mySpaces(cout,2);

}




