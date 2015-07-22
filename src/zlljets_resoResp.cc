#define zlljets_resoResp_cxx
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

#ifdef zlljets_resoResp_cxx

zlljets_resoResp::zlljets_resoResp(TTree *tree) : edimarcoTree_v2(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void zlljets_resoResp::loop(const char* configFileName, const Int_t ISDATA_FLAG)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

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
     fChain->SetBranchStatus("nGenPart",1);
     fChain->SetBranchStatus("GenPart_pdgId",1);
     fChain->SetBranchStatus("GenPart_motherId",1);
     fChain->SetBranchStatus("GenPart_pt",1);
     fChain->SetBranchStatus("GenPart_eta",1);
     fChain->SetBranchStatus("GenPart_phi",1);
     fChain->SetBranchStatus("GenPart_mass",1);
     fChain->SetBranchStatus("GenPart_motherIndex",1);
   }

   fChain->SetBranchStatus("met_pt",1);
   //fChain->SetBranchStatus("met_eta",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   //fChain->SetBranchStatus("metNoMu_eta",1);
   fChain->SetBranchStatus("metNoMu_phi",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices
   fChain->SetBranchStatus("vtxW",1);   // weight to have better agreement between data and MC 

   char ROOT_FNAME[100];
   char TXT_FNAME[100];
   char TEX_FNAME[100];
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"

   Double_t LUMI;
   Int_t NJETS;
   Double_t J1PT;
   Double_t J1ETA;
   Double_t J2PT;
   Double_t J2ETA;
   Double_t J1J2DPHI;
   Double_t RATIO_BR_ZINV_ZLL;
   Double_t UNC_RATIO_BR_ZINV_ZLL;
   Int_t LEP_PDG_ID;
   Int_t NEG_LEP_PDG_ID2;
   Double_t LEP1PT;
   Double_t LEP2PT;
   Double_t LEP1ETA;
   Double_t LEP2ETA;
   Double_t DILEPMASS_LOW;
   Double_t DILEPMASS_UP;
   Double_t LEP_ISO_04;
   Double_t GENLEP1PT;
   Double_t GENLEP2PT;
   Double_t GENLEP1ETA;
   Double_t GENLEP2ETA;
   Double_t GEN_ZMASS_LOW;
   Double_t GEN_ZMASS_UP;
   Int_t TAU_VETO_FLAG;
   Int_t HLT_FLAG;
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP2PT;
   Double_t HLT_LEP1ETA;
   Double_t HLT_LEP2ETA;
   Int_t NVTXS;                           // # of points for study of u_par and u_perp vs # of reconstructed vertices nvtx
   Int_t FIRST_NVTX;                    // starting number of vertices for met study   
   Double_t METNOLEP_START;
   Int_t JETS_SELECTION_RESORESP_FLAG;
   Int_t PHOTON_VETO_RESORESP_FLAG;
   Int_t LEPTON_VETO_RESORESP_FLAG;
   string FILENAME_BASE;

   ifstream inputFile(configFileName);

   if (inputFile.is_open()) {

     Double_t value;
     string name;
     string parameterName;
     string parameterType;
     vector<Double_t> parameterValue;

     mySpaces(cout,2);
     cout << "Printing content of " << configFileName << " file" << endl;
     mySpaces(cout,1);

     while (inputFile >> parameterType ) {

       if (parameterType == "NUMBER") {

	 inputFile >> parameterName >> value;
	 parameterValue.push_back(value);
	 cout << setw(20) << parameterName << setw(7) << value << endl;

       } else if (parameterType == "STRING") {
	 
	 inputFile >> parameterName >> name;
	 cout << right << setw(20) << parameterName << "  " << left << name << endl;
	 if (parameterName == "FILENAME_BASE") FILENAME_BASE = name; 

       }

     }
     
     // following variables are initialized with values in the file configFileName
     LUMI = parameterValue[0];
     NJETS = (Int_t) parameterValue[1];
     J1PT = parameterValue[2];
     J1ETA = parameterValue[3];
     J2PT = parameterValue[4];
     J2ETA = parameterValue[5];
     J1J2DPHI = parameterValue[6];
     RATIO_BR_ZINV_ZLL = parameterValue[7];
     UNC_RATIO_BR_ZINV_ZLL = parameterValue[8];
     LEP_PDG_ID = (Int_t) parameterValue[9];
     NEG_LEP_PDG_ID2 = (Int_t) parameterValue[10];
     LEP1PT = parameterValue[11];
     LEP2PT = parameterValue[12];
     LEP1ETA = parameterValue[13];
     LEP2ETA = parameterValue[14];
     DILEPMASS_LOW = parameterValue[15];
     DILEPMASS_UP = parameterValue[16];
     LEP_ISO_04 = parameterValue[17];
     GENLEP1PT = parameterValue[18];
     GENLEP2PT = parameterValue[19];
     GENLEP1ETA = parameterValue[20];
     GENLEP2ETA = parameterValue[21];
     GEN_ZMASS_LOW = parameterValue[22];
     GEN_ZMASS_UP = parameterValue[23];
     TAU_VETO_FLAG = (Int_t) parameterValue[24];
     HLT_FLAG = (Int_t) parameterValue[25];
     HLT_LEP1PT = parameterValue[26];
     HLT_LEP2PT = parameterValue[27];
     HLT_LEP1ETA = parameterValue[28];
     HLT_LEP2ETA = parameterValue[29];
     NVTXS = parameterValue[30];
     FIRST_NVTX = parameterValue[31];
     METNOLEP_START = parameterValue[32];
     JETS_SELECTION_RESORESP_FLAG = (Int_t) parameterValue[33];
     PHOTON_VETO_RESORESP_FLAG = (Int_t) parameterValue[34];
     LEPTON_VETO_RESORESP_FLAG = (Int_t) parameterValue[35];
     mySpaces(cout,2);

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   // vector<Double_t> metCut;
   // metCut.push_back(250);
   // metCut.push_back(300);
   // metCut.push_back(350);
   // metCut.push_back(400);
   // metCut.push_back(500);

   //selection njetsC("njetsC",Form("njets <= %i",NJETS),"pt > 30; |eta| < 4.7");   // using nJet30a
   //selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   //selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.5");       // using nJetClen30
   //selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   //selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   //selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),Form("only if njets = %i",NJETS));
  
   // selections for monojet selection (it also includes veto on muons or electrons depending on the sample
   selection jet1C("jet1C",Form("jet1pt > %4.0lf",(Double_t)J1PT),Form("nJetClean >= 1 && JetClean1_pt > %4.0lf && abs(JetClean1_eta) < %1.1lf && jetclean1 > 0.5",(Double_t)J1PT,J1ETA));
   selection jjdphiC("jjdphiC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection njetsC("njets","nJetClean30 <= 2");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC("tauLooseVetoC","tau veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection invMassC("invMassC",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
   // following selections are set differently in the next "if" statements depending on the lepton flavour 
   //selection metNoLepC[metCut.size()];
   selection lepLooseVetoC;
   selection twoLeptonsC;;
   selection twoLepLooseC;;
   selection lep1tightIdIso04C;;
   selection twoLepTightC;
   selection lep1ptC;
   selection lep2ptC;
   selection lep1etaC;  
   selection lep2etaC;
   selection genLepC;  
   selection metNoLepStartC;
   selection HLTlepC;
   // the following are only for electrons
   selection lep2tightIdIso04C;

   // TVector3 metNoLepTV3;   // metNoLep 3D vector, TV3 is to make clear that it's a 3D vector and not a TLorentzVector
   // TVector3 ele; 
   // ele is any electron to compute MetNoEle, for muons it's not needed because it's already in the tree
   //TVector3 met, eleVectorSum;

   // the reason to use TVector3 instead of TVector2 (which would be faster) is that I need to compute dphi between metNoLep and the Z vector (TLorentVector).
   // To do that I use the TVector3::DeltaPhi(const TVector3&) method passing TLorentzVector::Vect() which returns the 3D vector from a Lorentz Vector.
   // there isn't a method giving just the transverse 2D vector (although I could do it myself) from a TLorentVector (but it does exist a TVector3::XYvector()). 
   // Maybe I will change because TVector3 is much slower than TVector2

   TVector2 metNoLepTV, ele;

   TLorentzVector l1gen, l2gen, Zgen;     // gen level Z and l1,l2  (Z->(l1 l2)
   TLorentzVector l1reco, l2reco, Zreco;

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Int_t firstIndexGen = 0;
   Int_t secondIndexGen = 1;
   Int_t recoLepFound_flag = 0;
   Int_t genLepFound_flag = 0;
   Int_t Z_index = 0; 

   Double_t nTotalWeightedEvents = 0.0;     
   Int_t HLT_passed_flag = 1; // some computations (for e) require a trigger preselection, while other don't. The former will be done if the flag is set to 1
   // it's set to 1 because if the trigger selection is not applied every event must be considered to be a "good" event having passed all preselections
   // Actually in this code the trigger is necessary, but I keep it like this nonetheless.

   // following 2 variable are used for acceptance and efficiency selection, define below in the loop: if selection is passed they are set to 1, otherwise they are set to 0
   // Int_t acceptanceSelectionDef = 0;
   // Int_t efficiencySelectionDef = 0;

   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches
   Float_t *ptr_nLep10V = NULL;   

   Float_t *ptr_metNoLepPt = NULL;       // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)
   //Float_t *ptr_metNoLepEta = NULL; 
   Float_t *ptr_metNoLepPhi = NULL;  

   Float_t nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such
   Float_t nLep10V = 0.0;
   Double_t metNoLepPt = 0.0;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e
   //Double_t metNoLepEta = 0.0;
   Double_t metNoLepPhi = 0.0;   // same story as above

   if (ISDATA_FLAG) {
     strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + "_DATA.tex").c_str());
   } else {
     strcpy(ROOT_FNAME,(FILENAME_BASE + ".root").c_str());
     strcpy(TXT_FNAME,(FILENAME_BASE + ".txt").c_str());
     strcpy(TEX_FNAME,(FILENAME_BASE + ".tex").c_str());
   }

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
     
     strcpy(FLAVOUR,"mu");
     strcpy(LL_FLAVOUR,"mumu");
     strcpy(CONTROL_SAMPLE,"Z-->mumu");
         
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_metNoLepPt = &metNoMu_pt;               // for muons  get this variable from the tree 
     //ptr_metNoLepEta = &metNoMu_eta;               // for muons  get this variable from the tree 
     ptr_metNoLepPhi = &metNoMu_phi;         // for muons  get this variable from the tree

     // for (Int_t i = 0; i < metCut.size(); i++) {
     //   metNoLepC[i].set(Form("metNoMuC[%i]",i),Form("metNoMu > %3.0lf",metCut.at(i)));
     // }

     lepLooseVetoC.set("eLooseVetoC","electrons veto");
     twoLeptonsC.set("twomuonsC","muons");
     twoLepLooseC.set("twomuLooseC","2 loose muons");
     lep1tightIdIso04C.set("mu1tightIdIso04C","leading muon tight","tight ID + relIso04 (as Emanuele)");
     twoLepTightC.set("twomuTightC","2 tight muons");
     lep1ptC.set("mu1ptC",Form("mu1pt > %3.0lf",LEP1PT),"leading muon pt");
     lep2ptC.set("mu2ptC",Form("mu2pt > %3.0lf",LEP2PT),"trailing muon pt");
     lep1etaC.set("mu1etaC",Form("|mu1eta| < %1.1lf",LEP1ETA),"leading muon eta");  
     lep2etaC.set("mu2etaC",Form("|mu2eta| < %1.1lf",LEP2ETA),"trailing muon eta");
      if (!ISDATA_FLAG) genLepC.set("genMuonsC","muons generated");     
     metNoLepStartC.set("metNoMu200C",Form("metNoMu > %2.0lf",METNOLEP_START));
     HLTlepC.set("HLTmuonC","HLT for muons");
     lep2tightIdIso04C.set("mu2tightIdIso04C","trailing muon tight","tight ID + relIso04 (as Emanuele)");


   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     strcpy(FLAVOUR,"ele");
     strcpy(LL_FLAVOUR,"ee");
     strcpy(CONTROL_SAMPLE,"Z-->ee");

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons   

     // for (Int_t i = 0; i < metCut.size(); i++) {
     //   metNoLepC[i].set(Form("metNoEleC[%i]",i),Form("metNoEle > %3.0lf",metCut.at(i)));
     // }
     
     lepLooseVetoC.set("muLooseVetoC","muons veto");
     twoLeptonsC.set("twoelectronsC","electrons");
     twoLepLooseC.set("twoeleLooseC","2 loose electrons");
     lep1tightIdIso04C.set("ele1tightIdIso04C","leading electron tight","tight ID + relIso04 (as Emanuele)");
     twoLepTightC.set("twoeleTightC","2 tight electrons");
     lep1ptC.set("ele1ptC",Form("ele1pt > %3.0lf",LEP1PT),"leading electron pt");
     lep2ptC.set("ele2ptC",Form("ele2pt > %3.0lf",LEP2PT),"trailing electron pt");
     lep1etaC.set("ele1etaC",Form("|ele1eta| < %1.1lf",LEP1ETA),"leading electron eta");  
     lep2etaC.set("ele2etaC",Form("|ele2eta| < %1.1lf",LEP2ETA),"trailing electron eta");
     if (!ISDATA_FLAG) genLepC.set("genElectronsC","electrons generated");     
     metNoLepStartC.set("metNoEle200C",Form("metNoEle > %2.0lf",METNOLEP_START));
     HLTlepC.set("HLTelectronC","HLT for electrons");
     lep2tightIdIso04C.set("ele2tightIdIso04C","trailing electron tight","tight ID + relIso04 (as Emanuele)");

   }

   //selection genTausC("genTausC","taus generated");                       
   //selection acceptanceC("acceptanceC","acceptance cuts");
   //selection efficiencyC("efficiencyC","efficiency cuts");

   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   UInt_t maskJetsSelection = njetsC.get2ToId() + jet1C.get2ToId() + jjdphiC.get2ToId();

   UInt_t maskMonoJetSelection = maskJetsSelection + lepLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   if ( TAU_VETO_FLAG ) maskMonoJetSelection += tauLooseVetoC.get2ToId();

   UInt_t maskTightTag;   // holds cuts for lepton tight selection, which is different between muons and electrons

   //mask zlljetsControlSample(Form("%s control sample with selection flow as Emanuele's",CONTROL_SAMPLE));

   mask zlljetsControlSampleGenLep(Form("%s control sample (%s gen if MC) with selection flow as Emanuele's",CONTROL_SAMPLE,FLAVOUR));
   if (!ISDATA_FLAG) zlljetsControlSampleGenLep.append(genLepC.get2ToId());
   zlljetsControlSampleGenLep.append(HLTlepC.get2ToId());

   // mask tautaubkgInZll(Form("tau tau background in %s control sample",CONTROL_SAMPLE));
   // tautaubkgInZll.append(genTausC.get2ToId());

   mask resoAndResponse("selection for resolution and response");
   resoAndResponse.append(HLTlepC.get2ToId());
   resoAndResponse.append(oppChargeLeptonsC.get2ToId() + twoLepLooseC.get2ToId());
   

   // if (HLT_FLAG) {

   //   //zlljetsControlSample.append(HLTlepC.get2ToId());
     
   //   //tautaubkgInZll.append(HLTlepC.get2ToId());

   // }

   if (METNOLEP_START) {
       
     zlljetsControlSampleGenLep.append(metNoLepStartC.get2ToId());
     //resoAndResponse.append(metNoLepStartC.get2ToId());

   }

   if (fabs(LEP_PDG_ID) == 13) {  

     maskTightTag = lep1tightIdIso04C.get2ToId() + lep2tightIdIso04C.get2ToId() + lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();  ;  // for now tight requirements on pt and eta are already included in the loose condition because they coincide (not true for electrons)
   
     // zlljetsControlSample.append(metNoLepStartC.get2ToId());
     // zlljetsControlSample.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     // zlljetsControlSample.append(twoLeptonsC.get2ToId());
     // zlljetsControlSample.append(maskTightTag);
     // zlljetsControlSample.append(invMassC.get2ToId());
         
     //zlljetsControlSampleGenLep.append(metNoLepStartC.get2ToId());
     zlljetsControlSampleGenLep.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(twoLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(maskTightTag); 
     zlljetsControlSampleGenLep.append(invMassC.get2ToId());
   
     // tautaubkgInZll.append(metNoLepStartC.get2ToId());
     // tautaubkgInZll.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     // tautaubkgInZll.append(twoLeptonsC.get2ToId());
     // tautaubkgInZll.append(maskTightTag);
     // tautaubkgInZll.append(invMassC.get2ToId());

     resoAndResponse.append(maskTightTag);
     resoAndResponse.append(invMassC.get2ToId());

   } else if (fabs(LEP_PDG_ID) == 11) {  

     maskTightTag = lep1tightIdIso04C.get2ToId() + lep2tightIdIso04C.get2ToId() + lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();

     // zlljetsControlSample.append(metNoLepStartC.get2ToId());
     // zlljetsControlSample.append(oppChargeLeptonsC.get2ToId()); // skip loose requirement because I wil ask the tight one for both
     // zlljetsControlSample.append(twoLeptonsC.get2ToId());
     // zlljetsControlSample.append(maskTightTag);
     // zlljetsControlSample.append(invMassC.get2ToId());
       
     //zlljetsControlSampleGenLep.append(metNoLepStartC.get2ToId());
     zlljetsControlSampleGenLep.append(oppChargeLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(twoLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(maskTightTag);
     zlljetsControlSampleGenLep.append(invMassC.get2ToId());
   
     // tautaubkgInZll.append(metNoLepStartC.get2ToId());
     // tautaubkgInZll.append(oppChargeLeptonsC.get2ToId());
     // tautaubkgInZll.append(twoLeptonsC.get2ToId());
     // tautaubkgInZll.append(maskTightTag);
     // tautaubkgInZll.append(invMassC.get2ToId());
     
     resoAndResponse.append(maskTightTag);
     resoAndResponse.append(invMassC.get2ToId());

   }
 
   // zlljetsControlSample.append(jet1C.get2ToId());
   // zlljetsControlSample.append(jjdphiC.get2ToId());
   // zlljetsControlSample.append(njetsC.get2ToId());
   // zlljetsControlSample.append(lepLooseVetoC.get2ToId());
   // zlljetsControlSample.append(gammaLooseVetoC.get2ToId());
   
   zlljetsControlSampleGenLep.append(jet1C.get2ToId());
   zlljetsControlSampleGenLep.append(jjdphiC.get2ToId());
   zlljetsControlSampleGenLep.append(njetsC.get2ToId());
   zlljetsControlSampleGenLep.append(lepLooseVetoC.get2ToId());
   zlljetsControlSampleGenLep.append(gammaLooseVetoC.get2ToId());

   // tautaubkgInZll.append(jet1C.get2ToId());
   // tautaubkgInZll.append(jjdphiC.get2ToId());
   // tautaubkgInZll.append(njetsC.get2ToId());
   // tautaubkgInZll.append(lepLooseVetoC.get2ToId());
   // tautaubkgInZll.append(gammaLooseVetoC.get2ToId());

   if (JETS_SELECTION_RESORESP_FLAG) resoAndResponse.append(maskJetsSelection);
   if (PHOTON_VETO_RESORESP_FLAG) resoAndResponse.append(gammaLooseVetoC.get2ToId());
   if (LEPTON_VETO_RESORESP_FLAG) resoAndResponse.append(lepLooseVetoC.get2ToId());

   if (TAU_VETO_FLAG) {
     // zlljetsControlSample.append(tauLooseVetoC.get2ToId());
     zlljetsControlSampleGenLep.append(tauLooseVetoC.get2ToId());
     resoAndResponse.append(tauLooseVetoC.get2ToId());
   }

   cout << "Opening file " <<ROOT_FNAME<< endl;

   TFile *rootFile = new TFile(ROOT_FNAME,"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }
 

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   //Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};       

   //TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin",Form("yields of %s control sample in bins of met;#slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinGenLep = new TH1D("HzlljetsYieldsMetBinGenLep",Form("yields of %s control sample (%s gen if MC) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE,CONTROL_SAMPLE),nMetBins,metBinEdges);
   //TH1D *HzlljetsYieldsMetBinGenTau = new TH1D("HzlljetsYieldsMetBinGenTau",Form("yields of %s control sample (Z->#tau#tau gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   
   TH1D *HZtoLLRecoPt = new TH1D("HZtoLLRecoPt","",101,0.,1010);
   TH1D *HZtoLLGenPt ;
   TH1D *HZtoLLPt_RecoGenRatio;                    // this is the histogram with reco/gen
   TH1D *HZtoLLPt_RecoGenRatio_pdf;             // histogram of reco/gen distribution function
   TH1D *HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf;     

   if (!ISDATA_FLAG) {
     HZtoLLGenPt = new TH1D("HZtoLLGenPt","",101,0.,1010);
     HZtoLLPt_RecoGenRatio = new TH1D("HZtoLLPt_RecoGenRatio","",101,0.,1010.);
     HZtoLLPt_RecoGenRatio_pdf = new TH1D("HZtoLLPt_RecoGenRatio_pdf","",100,0.5,1.5);
     HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf = new TH1D("HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf","",100,0.5,1.5);
   }

   // TH1D* Hacc = new TH1D("Hacc","",nMetBins,metBinEdges);
   // TH1D* Heff = new TH1D("Heff","",nMetBins,metBinEdges);
   // TH1D* Hacceff = new TH1D("Hacceff","",nMetBins,metBinEdges);
   // TH1D *HzvvEstimate = new TH1D("HzvvEstimate",Form("yields of Z->#nu#nu estimated as N(%s) * BR_ratio / (A*#varepsilon)",CONTROL_SAMPLE),nMetBins,metBinEdges);

   //TH1D *HinvMass[nMetBins];
   TH1D *HzlljetsInvMassMetBinGenLep[nMetBins];
   //TH1D *HzlljetsInvMassMetBinGenTau[nMetBins];
   TH1D *HZtoLLRecoPt_MetBin[nMetBins];
   TH1D *HZtoLLGenPt_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_pdf_MetBin[nMetBins];

   Float_t invMassBinWidth = 2.0;  // invariant mass histogram's bin width in GeV
   Int_t NinvMassBins = (DILEPMASS_UP - DILEPMASS_LOW) / invMassBinWidth;

   for (Int_t i = 0; i < nMetBins; i++) {

     //HinvMass[i] = new TH1D(Form("HinvMass_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
     HzlljetsInvMassMetBinGenLep[i] = new TH1D(Form("HzlljetsInvMassMetBinGenLep_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
     //HzlljetsInvMassMetBinGenTau[i] = new TH1D(Form("HzlljetsInvMassMetBinGenTau_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);
     HZtoLLRecoPt_MetBin[i] = new TH1D(Form("HZtoLLRecoPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);

     if (!ISDATA_FLAG) {
       HZtoLLGenPt_MetBin[i] = new TH1D(Form("HZtoLLGenPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
       HZtoLLPt_RecoGenRatio_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
       HZtoLLPt_RecoGenRatio_pdf_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_pdf_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",100,0.5,1.5);
     }

   } 

   TH1D *H_uPerp_VS_Nvtx[NVTXS];
   TH1D *H_uPar_VS_Nvtx[NVTXS]; 
   TH1D *H_uPar_VS_Nvtx_lowZpT[NVTXS];
 
   for (Int_t i = 0; i < NVTXS; i++) {

     H_uPerp_VS_Nvtx[i] = new TH1D(Form("H_uPerp_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"",80,-200,200);  // 5 GeV bins 
     H_uPar_VS_Nvtx[i] = new TH1D(Form("H_uPar_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"ZpT in [250,500]GeV",80,-200,200);  // 5 GeV bins
     H_uPar_VS_Nvtx_lowZpT[i] = new TH1D(Form("H_uPar_VS_Nvtx_lowZpT_nvtx%i",FIRST_NVTX+i),"ZpT in [50,250]GeV",80,-200,200);  // 5 GeV bins

   }

   // zpt bin edges for respose studies (which are independent on the CS selection, they are done on the whole sample)
   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 390., 410., 430., 450., 470., 500., 530., 560, 600., 640., 700., 800.};
   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   //Double_t ZptBinEdges[] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   Double_t ZptBinEdges[] = {20., 40., 60., 80., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};

   Int_t nBinsForResponse = sizeof(ZptBinEdges)/sizeof(Double_t) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements

   TH1D *H_uPerp_VS_ZpT[nBinsForResponse];  
   TH1D *H_uPar_VS_ZpT[nBinsForResponse]; 
   TH1D *H_uPar_ZpT_ratio[nBinsForResponse];  // for the response curve
   TH1D *HZptBinned[nBinsForResponse];
   //the following histograms will give the distribution of met|| / wzpt. The mean value will be used to create the response curve, that is (<met|| / wzpt>) vs wzpt
   // for each point, wzpt will be taken as the average wzpt in the range considered
 
   for (Int_t i = 0; i < nBinsForResponse; i++) {   

     //HZptBinned[i] are histograms with 5 bins in the range given by ZptBinEdges[i] and ZptBinEdges[i+1]
     // the mean wzpt in each bin will be computed as the histogram's mean
     HZptBinned[i] = new TH1D(Form("HZptBinned_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",5,ZptBinEdges[i],ZptBinEdges[i+1]); 
     H_uPar_ZpT_ratio[i] = new TH1D(Form("H_uPar_ZpT_ratio_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",50,0.0,2.0); 
     H_uPerp_VS_ZpT[i] = new TH1D(Form("H_uPerp_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 
     H_uPar_VS_ZpT[i] = new TH1D(Form("H_uPar_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 

   }

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   TH1D *HZptBinEdges = new TH1D("HZptBinEdges","bin edges for ZpT distributions",nBinsForResponse+1,0.0,nBinsForResponse+1);
   for (Int_t i = 0; i <= nBinsForResponse; i++) {
     HZptBinEdges->SetBinContent(i+1,ZptBinEdges[i]);
   }

   TH1D *Hnvtx = new TH1D("Hnvtx","# of vertices for studies of variables as a function of nvtx",NVTXS,0.0,NVTXS);
   for (Int_t i = 0; i < NVTXS; i++) {              // watch out: differently from above, i < NVTXS, not <=, because if NVTXS = 3 I need 3 points, not 4
     Hnvtx->SetBinContent(i+1,FIRST_NVTX+i);
   }

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljets_resoResp::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     if (jentry%500000 == 0) cout << jentry << endl;

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     nTotalWeightedEvents += newwgt;  // counting events with weights

     nLepLoose = *ptr_nLepLoose;          
     nLep10V = *ptr_nLep10V;

     Double_t ZgenMass; 
     Double_t ZtoLLGenPt;    // could do Double_t ZtoLLGenPt = GenPart_pt[Z_index];

     if (!ISDATA_FLAG) {

       genLepFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, 23, firstIndexGen, secondIndexGen, Z_index, GenPart_motherIndex); 
       if (!genLepFound_flag) continue;

       eventMask += genLepC.addToMask( genLepFound_flag );
       l1gen.SetPtEtaPhiM(GenPart_pt[firstIndexGen],GenPart_eta[firstIndexGen],GenPart_phi[firstIndexGen],GenPart_mass[firstIndexGen]);
       l2gen.SetPtEtaPhiM(GenPart_pt[secondIndexGen],GenPart_eta[secondIndexGen],GenPart_phi[secondIndexGen],GenPart_mass[secondIndexGen]);
       Zgen = l1gen + l2gen;
       Z_index = GenPart_motherIndex[firstIndexGen];   //could do Z_index = myGetPartIndex(23, nGenPart, GenPart_pdgId);  
       ZgenMass = Zgen.Mag(); 
       ZtoLLGenPt = Zgen.Pt();    // could do Double_t ZtoLLGenPt = GenPart_pt[Z_index];

     }

     recoLepFound_flag = myGetPairIndexInArray(LEP_PDG_ID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);  
     if (!recoLepFound_flag) continue;     

     l1reco.SetPtEtaPhiM(LepGood_pt[firstIndex],LepGood_eta[firstIndex],LepGood_phi[firstIndex],LepGood_mass[firstIndex]);
     l2reco.SetPtEtaPhiM(LepGood_pt[secondIndex],LepGood_eta[secondIndex],LepGood_phi[secondIndex],LepGood_mass[secondIndex]);
     Zreco = l1reco + l2reco;
     Double_t ZtoLLRecoPt = Zreco.Pt();

     if (fabs(LEP_PDG_ID) == 13) { 

       if ( HLT_FLAG ) {

	 // use the dimuon trigger, not the metNoLep trigger
       	 if ( recoLepFound_flag && (fabs(LepGood_eta[firstIndex]) < HLT_LEP1ETA) && (fabs(LepGood_eta[secondIndex]) < HLT_LEP2ETA) && 
       	      (LepGood_pt[firstIndex] > HLT_LEP1PT) && (LepGood_pt[secondIndex] > HLT_LEP2PT) ) HLT_passed_flag = 1; 	 
       	 else continue;

       }  // end of   if ( HLT_FLAG )

       // if ( genLepFound_flag && (GenPart_pt[firstIndexGen] > GENLEP1PT) && (GenPart_pt[secondIndexGen] > GENLEP2PT) && ( fabs(GenPart_eta[firstIndexGen]) < GENLEP1ETA) && ( fabs(GenPart_eta[secondIndexGen]) < GENLEP2ETA) && (ZgenMass > GEN_ZMASS_LOW) && (ZgenMass < GEN_ZMASS_UP) )  acceptanceSelectionDef = 1;
       // else acceptanceSelectionDef = 0;

       // if (recoLepFound_flag && (nLepLoose == 2) && (LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < LEP_ISO_04)) efficiencySelectionDef = 1;
       // else efficiencySelectionDef = 0;

       metNoLepPt = *ptr_metNoLepPt;       
       //metNoLepEta = *ptr_metNoLepEta; 
       metNoLepPhi = *ptr_metNoLepPhi; 
       //metNoLepTV3.SetPtEtaPhi(metNoLepPt,metNoLepEta,metNoLepPhi);   // will use this 3D vector below
       metNoLepTV.SetMagPhi(metNoLepPt,metNoLepPhi);

     } else if (fabs(LEP_PDG_ID) == 11) { 

       if ( HLT_FLAG ) {

       	 if ( recoLepFound_flag && (LepGood_tightId[firstIndex] == 1) && (LepGood_tightId[secondIndex] == 1) && 
       	      (fabs(LepGood_eta[firstIndex]) < HLT_LEP1ETA) && (fabs(LepGood_eta[secondIndex]) < HLT_LEP2ETA) && 
       	      (LepGood_pt[firstIndex] > HLT_LEP1PT) && (LepGood_pt[secondIndex] > HLT_LEP2PT) ) HLT_passed_flag = 1; 	 
	 else continue;

       }  // end of   if ( HLT_FLAG )

       // metNoLepTV3.SetPtEtaPhi(met_pt,met_eta,met_phi);
       // // summing just electrons from Z if found
       // ele.SetPtEtaPhi(LepGood_pt[firstIndex],LepGood_eta[firstIndex],LepGood_phi[firstIndex]);
       // metNoLepTV3 += ele;
       // ele.SetPtEtaPhi(LepGood_pt[secondIndex],LepGood_eta[secondIndex],LepGood_phi[secondIndex]);
       // metNoLepTV3 += ele;

       metNoLepTV.SetMagPhi(met_pt,met_phi);
       // summing just electrons from Z if found
       ele.SetMagPhi(LepGood_pt[firstIndex],LepGood_phi[firstIndex]);
       metNoLepTV += ele;
       ele.SetMagPhi(LepGood_pt[secondIndex],LepGood_phi[secondIndex]);
       metNoLepTV += ele;

       metNoLepPt = metNoLepTV.Mod();

       //metNoLepTV3 = met + eleVectorSum;  // metNoLep vector created summing real met vector and vector sum of all electrons (actually only those from Z)
       // for muons I would sum all the muons because in any case the analysis would require exactly 2 muons (for control sample) and no muons for the signal
       //metNoLepPt = metNoLepTV3.Pt();  // for electrons we define components by hand, for muons we used the variable in the tree to form the vector
       //metNoLepEta = metNoLepTV3.Eta();  // not needed, I'll just use the TVector directly
       //metNoLepPhi = metNoLepTV3.Phi();  // not needed, I'll just use the TVector directly

       // if ( genLepFound_flag && (GenPart_pt[firstIndexGen] > GENLEP1PT) && (GenPart_pt[secondIndexGen] > GENLEP2PT) &&
       // 	    ( fabs(GenPart_eta[firstIndexGen]) < GENLEP1ETA) && ( fabs(GenPart_eta[secondIndexGen]) < GENLEP2ETA) &&
       // 	    (ZgenMass > GEN_ZMASS_LOW) && (ZgenMass < GEN_ZMASS_UP) ) acceptanceSelectionDef = 1;
       // else acceptanceSelectionDef = 0;

       // if ( recoLepFound_flag && (LepGood_tightId[firstIndex] == 1) && (LepGood_tightId[secondIndex] == 1) &&
       // 	    (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ) ) efficiencySelectionDef = 1;
       // else efficiencySelectionDef = 0;

     }

     // beginning of eventMask building

     
     // genLepC added to mask above if ISDATA_FLAG == false (in order not to repeat here the check 

     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT && fabs(JetClean_eta[0] < J1ETA && jetclean1 > 0.5));  //could skip cut on eta
     eventMask += jjdphiC.addToMask( nJetClean30 == 1 || (nJetClean30 >= NJETS && fabs(dphijj) < J1J2DPHI && jetclean2 > 0.5));
     eventMask += njetsC.addToMask(nJetClean30 <= NJETS);
     eventMask += lepLooseVetoC.addToMask(nLep10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += metNoLepStartC.addToMask(metNoLepPt > METNOLEP_START);

     // for (Int_t i = 0; i <  metCut.size(); i++) {
     //   eventMask += metNoLepC[i].addToMask(metNoLepPt > metCut[i]);
     // }
     
     // the following make sense only if recoLepFound_flag == 1 (i.e. flag is true), which means that fabs(LepGood_pdgId[firstIndex/secondIndex]) == LEP_PDG_ID) is 
     // true
     // also, 2 OS/SF leptons are present
     if (recoLepFound_flag) {

       eventMask += HLTlepC.addToMask(HLT_passed_flag);     
       eventMask += oppChargeLeptonsC.addToMask( 1);
       eventMask += twoLeptonsC.addToMask(1);
       eventMask += twoLepLooseC.addToMask(nLepLoose == 2);
       eventMask += lep1ptC.addToMask((LepGood_pt[firstIndex] > LEP1PT)); 
       eventMask += lep1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < LEP1ETA) );
       eventMask += lep2ptC.addToMask((LepGood_pt[secondIndex] > LEP2PT) );
       eventMask += lep2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < LEP2ETA) );
       eventMask += invMassC.addToMask((mZ1 > DILEPMASS_LOW) && (mZ1 < DILEPMASS_UP));     
       eventMask += lep1tightIdIso04C.addToMask((LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) );
       eventMask += lep2tightIdIso04C.addToMask((LepGood_tightId[secondIndex] == 1) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ));
      
     }

     // end of eventMask building


       // this part, for electrons, is done only if trigger is passed
       // nTotalWeightedEvents counts the number of events passing trigger selection, if any

       //zlljetsControlSample.countEvents(eventMask,newwgt);
       zlljetsControlSampleGenLep.countEvents(eventMask,newwgt);
       //tautaubkgInZll.countEvents(eventMask, newwgt);
       resoAndResponse.countEvents(eventMask, newwgt);

       // filling histogram with yields and invariant mass at the end of the selection in bins of met
       // if ( ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) ) {  
       // 	 // this histogram holds the final yields in bins of MET
       // 	 HzlljetsYieldsMetBin->Fill(metNoLepPt,newwgt);    
       // }

       if ( ((eventMask & zlljetsControlSampleGenLep.globalMask.back()) == zlljetsControlSampleGenLep.globalMask.back()) ) {

	 // this histogram holds the final yields in bins of MET
	 HzlljetsYieldsMetBinGenLep->Fill(metNoLepPt,newwgt);

	 if (!ISDATA_FLAG) {

	   //enter this part if 2 OS/SF leptons were found among gen and reco particles. Now checking compatibilities between pairs
	   // e.g. l1gen = e+, l2gen = e- ; l1reco = e+, l2reco = e- (but the charge order might not coincide)
	   // now we require a DeltaR cut between them to assess that lreco comes from lgen
	   // since 2 OS/SF were found to get inside here, if !(l1gen->l1reco && l2gen->l2reco) then for sure l1gen->l2reco && l2gen->l1reco
	   
	   Double_t DeltaR_lreco_lgen_pair1 = 0.0;
	   Double_t DeltaR_lreco_lgen_pair2 = 0.0;

	   if(LepGood_pdgId[firstIndex] == GenPart_pdgId[firstIndexGen] && LepGood_pdgId[secondIndex] == GenPart_pdgId[secondIndexGen]) {

	     DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l1gen);
	     DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l2gen);

	   } else {
	    
	     DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l2gen);
	     DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l1gen);

	   }

	   if (DeltaR_lreco_lgen_pair1 < 0.1 && DeltaR_lreco_lgen_pair2 < 0.1) {
	     
	     HZtoLLRecoPt->Fill(ZtoLLRecoPt,newwgt);
	     HZtoLLGenPt->Fill(ZtoLLGenPt,newwgt);
	     if (ZtoLLGenPt != 0) {

	       HZtoLLPt_RecoGenRatio_pdf->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);
	       if (ZtoLLRecoPt > 600) HZtoLLPt_RecoGenRatio_pdf_ZpT600ToInf->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);

	     }

	   }

	 } else HZtoLLRecoPt->Fill(ZtoLLRecoPt,newwgt);	 // if running on data just do this

       }

       // if ( ((eventMask & tautaubkgInZll.globalMask.back()) == tautaubkgInZll.globalMask.back()) ) {  
       // 	 // this histogram holds the final yields in bins of MET
       // 	 HzlljetsYieldsMetBinGenTau->Fill(metNoLepPt,newwgt);  
       // }

       if ( ((eventMask & resoAndResponse.globalMask.back()) == resoAndResponse.globalMask.back()) ) {  
	   
	 // following is done if two OS leptons are found (otherwise there would be no Z)
	 // moreover, if we have electron, trigger selection must be passed 
	 //metNoLepPt cut is also included: it's part of the trigger for muons and we also apply it to electrons for consistency
		
	 //metNoLepTV3.SetPtEtaPhi((Double_t)metNoLepPt,(Double_t)metNoLepEta,(Double_t)metNoLepPhi);  // already initialized above
	 //Double_t dphiMetNoLepZ = metNoLepTV3.DeltaPhi(Zreco.Vect());
	 TVector3 Zreco3D = Zreco.Vect();
	 Double_t dphiMetNoLepZ = metNoLepTV.DeltaPhi(Zreco3D.XYvector());

	 Double_t u_par = metNoLepPt * TMath::Cos(dphiMetNoLepZ);  // actually u_par is minus this quantity, but then I do u_par-ZpT instead of u_par+ZpT
	 Double_t u_perp = metNoLepPt * TMath::Sin(dphiMetNoLepZ);

	 if (ZtoLLRecoPt > ZptBinEdges[0]) {  

	   Int_t nvtxBin = nVert - FIRST_NVTX;
	   Int_t lastnvtx = NVTXS + FIRST_NVTX;

	   Double_t uparMinusZrecoPt = u_par - ZtoLLRecoPt;

	   if ((nvtxBin >= 0) && (nVert < lastnvtx)) {

	     H_uPerp_VS_Nvtx[nvtxBin]->Fill(u_perp,newwgt);
	     
	     if (ZtoLLRecoPt < 250 ) {

	       H_uPar_VS_Nvtx_lowZpT[nvtxBin]->Fill(uparMinusZrecoPt,newwgt);
	     
	     } else if (ZtoLLRecoPt < 500) {                       // (met||-wzpt) distribution's width depends on Zpt, thus I use this range

	       H_uPar_VS_Nvtx[nvtxBin]->Fill(uparMinusZrecoPt,newwgt);
	 
	     }       
    
	   }  // end of   if ((nvtxBin >= 0) && (nVert < lastnvtx))

	   /**************************************************/
	   // computing met responses
	   /**************************************************/

	   // first of all I make sure that wzpt is in the appropriate range
	   if ( ZtoLLRecoPt < ZptBinEdges[nBinsForResponse] ) {

	     Int_t respBin = myGetBin(ZtoLLRecoPt,ZptBinEdges,nBinsForResponse);
	     //cout<<"bin = "<<bin<<endl;
	     HZptBinned[respBin]->Fill(ZtoLLRecoPt,newwgt);        
	     H_uPar_ZpT_ratio[respBin]->Fill(u_par/ZtoLLRecoPt,newwgt);     //the mean value of this histogram is the response
	     H_uPerp_VS_ZpT[respBin]->Fill(u_perp,newwgt);
	     H_uPar_VS_ZpT[respBin]->Fill(uparMinusZrecoPt,newwgt);
 
	   }

	 }            // end of if (ZtoLLRecoPt > ZptBinEdges[0])

       }
	

       // now entering analysis in bins of met

       if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) {

	 Int_t bin = myGetBin(metNoLepPt,metBinEdges,nMetBins);
       
	 // if ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) {
	 //   // this histogram holds the invariant mass distribution (one for each met bin)
	 //   HinvMass[bin]->Fill(mZ1,newwgt);   
	 // }

	 if ( ((eventMask & zlljetsControlSampleGenLep.globalMask.back()) == zlljetsControlSampleGenLep.globalMask.back()) ) { 
 
	   HzlljetsInvMassMetBinGenLep[bin]->Fill(mZ1,newwgt); 

	   if (!ISDATA_FLAG) {

	     //enter this part if 2 OS/SF leptons were found among gen and reco particles. Now checking compatibilities between pairs
	     // e.g. l1gen = e+, l2gen = e- ; l1reco = e+, l2reco = e- (but the charge order might not coincide)
	     // now we require a DeltaR cut between them to assess that lreco comes from lgen
	     // since 2 OS/SF were found to get inside here, if !(l1gen->l1reco && l2gen->l2reco) then for sure l1gen->l2reco && l2gen->l1reco

	     Double_t DeltaR_lreco_lgen_pair1 = 0.0;
	     Double_t DeltaR_lreco_lgen_pair2 = 0.0;

	     if(LepGood_pdgId[firstIndex] == GenPart_pdgId[firstIndexGen] && LepGood_pdgId[secondIndex] == GenPart_pdgId[secondIndexGen]) {

	       DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l1gen);
	       DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l2gen);

	     } else {
	    
	       DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l2gen);
	       DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l1gen);

	     }

	     if (DeltaR_lreco_lgen_pair1 < 0.1 && DeltaR_lreco_lgen_pair2 < 0.1) {
	     
	       HZtoLLRecoPt_MetBin[bin]->Fill(ZtoLLRecoPt,newwgt);
	       HZtoLLGenPt_MetBin[bin]->Fill(ZtoLLGenPt,newwgt);
	       if (ZtoLLGenPt != 0) HZtoLLPt_RecoGenRatio_pdf_MetBin[bin]->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);

	     }

	   }

	 } else HZtoLLRecoPt_MetBin[bin]->Fill(ZtoLLRecoPt,newwgt);  // if running on data just do this
	 
       }                      // end of    if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) 
       
   }                        // end of loop on entries

   
   /************************************/
   //                    MET|| & MET_|_ VS NVTX & ZpT
   /************************************/

   //resolution vs nvtx

   Double_t xValues[NVTXS];
   Double_t yValues[NVTXS];
   Double_t yValuesErr[NVTXS];

   for (Int_t i = 0; i < NVTXS; i++) {
     xValues[i] = i + FIRST_NVTX;
     yValues[i] = H_uPar_VS_Nvtx_lowZpT[i]->GetRMS();
     yValuesErr[i] = H_uPar_VS_Nvtx_lowZpT[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx_lowZpT = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetTitle(Form("resolution || from histogram's RMS, ZpT in [%2.0lf,250] GeV",ZptBinEdges[0]));
   GresolutionMetNoLepParZvsNvtx_lowZpT->Draw("AP");
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx_lowZpT->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx_lowZpT->SetName("gr_resolution_uPar_vs_Nvtx_lowZpT");
   GresolutionMetNoLepParZvsNvtx_lowZpT->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uPar_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uPar_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx->SetTitle("resolution || from histogram's RMS, ZpT in [250,500] GeV");
   GresolutionMetNoLepParZvsNvtx->Draw("AP");
   GresolutionMetNoLepParZvsNvtx->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx->SetName("gr_resolution_uPar_vs_Nvtx");
   GresolutionMetNoLepParZvsNvtx->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uPerp_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uPerp_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepOrtZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepOrtZvsNvtx->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsNvtx->Draw("AP");
   GresolutionMetNoLepOrtZvsNvtx->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepOrtZvsNvtx->SetName("gr_resolution_uPerp_vs_Nvtx");
   GresolutionMetNoLepOrtZvsNvtx->Write();

   // response curve

   Double_t response[nBinsForResponse];
   Double_t responseErr[nBinsForResponse];
   Double_t meanZpt[nBinsForResponse];
   Double_t meanZptErr[nBinsForResponse];

   for (Int_t i = 0; i < nBinsForResponse; i++) {
     meanZpt[i] = HZptBinned[i]->GetMean();
     meanZptErr[i] = HZptBinned[i]->GetMeanError();
     response[i] = H_uPar_ZpT_ratio[i]->GetMean();
     responseErr[i] = H_uPar_ZpT_ratio[i]->GetMeanError();
     //cout<<i<<" meanZpt = "<<meanZpt[i]<<" +/- "<<meanZptErr[i]<<"    response = "<<response[i]<<" +/- "<<responseErr[i]<<endl;
   }

   TGraphErrors *GresponseCurve = new TGraphErrors(nBinsForResponse,meanZpt,response,meanZptErr,responseErr);
   GresponseCurve->SetTitle("response curve");
   GresponseCurve->Draw("AP");
   GresponseCurve->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve->GetYaxis()->SetTitle(" < u_{||} / ZpT >");
   GresponseCurve->GetYaxis()->SetRangeUser(0.6, 1.1);
   GresponseCurve->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve->SetName("gr_responseCurve");
   GresponseCurve->Write();

   // resolution vs ZpT

   Double_t resoMetNoLepParZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepParZvsZptErr[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZptErr[nBinsForResponse];

   for (Int_t i = 0; i < nBinsForResponse; i++) {
     resoMetNoLepParZvsZpt[i] = H_uPar_VS_ZpT[i]->GetRMS();
     resoMetNoLepParZvsZptErr[i] = H_uPar_VS_ZpT[i]->GetRMSError();
     resoMetNoLepOrtZvsZpt[i] = H_uPerp_VS_ZpT[i]->GetRMS();
     resoMetNoLepOrtZvsZptErr[i] = H_uPerp_VS_ZpT[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepParZvsZpt,meanZptErr,resoMetNoLepParZvsZptErr);
   GresolutionMetNoLepParZvsZpt->SetTitle("resolution || from histogram's RMS");
   GresolutionMetNoLepParZvsZpt->Draw("AP");
   GresolutionMetNoLepParZvsZpt->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepParZvsZpt->SetName("gr_resolution_uPar_vs_ZpT");
   GresolutionMetNoLepParZvsZpt->Write();

   TGraphErrors *GresolutionMetNoLepOrtZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepOrtZvsZpt,meanZptErr,resoMetNoLepOrtZvsZptErr);
   GresolutionMetNoLepOrtZvsZpt->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsZpt->Draw("AP");
   GresolutionMetNoLepOrtZvsZpt->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepOrtZvsZpt->SetName("gr_resolution_uPerp_vs_ZpT");
   GresolutionMetNoLepOrtZvsZpt->Write();

   // end of TGraphs

   mySpaces(cout,2);
//selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zlljetsControlSampleGenLep);
//selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &tautaubkgInZll);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &resoAndResponse);
   // for (Int_t i = 0; i < nMetBins; i++) {
   //   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEventsNoHLT, lep_acc_eff[i] );
   // }

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBinGenLep, metBinEdges, nMetBins);
 
   cout<<"creating file '"<<TXT_FNAME<<"' ..."<<endl;
   ofstream myfile(TXT_FNAME,ios::out);

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

   mySpaces(myfile,3);
   //selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zlljetsControlSampleGenLep);
   //selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &tautaubkgInZll);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &resoAndResponse);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBinGenLep, metBinEdges, nMetBins);

   // TH1D *HevtPassMonoJetSel = new TH1D("HevtPassMonoJetSel","events passing monojet selection for A#times#epsilon",nMetBins,metBinEdges);
   // TH1D *HevtPassAccSel = new TH1D("HevtPassAccSel","events passing acceptance selection for A#times#epsilon",nMetBins,metBinEdges);
   // TH1D *HevtPassEffSel = new TH1D("HevtPassEffSel","events passing efficiency selection for A#times#epsilon",nMetBins,metBinEdges);

   // // using [0] element to find step, all elements are equivalent for this purpose
   // Int_t stepMonojetSelection_In_lepAccEff = lep_acc_eff[0]->whichStepHas(maskMonoJetSelection);
   // Int_t stepAcceptance_In_lepAccEff = lep_acc_eff[0]->whichStepHas(acceptanceC.get2ToId());
   // Int_t stepEfficiency_In_lepAccEff = lep_acc_eff[0]->whichStepHas(efficiencyC.get2ToId());
   // // cout<<"step: MJ     acc     eff"<<endl;
   // // cout<<stepMonojetSelection_In_lepAccEff<<stepAcceptance_In_lepAccEff<<stepEfficiency_In_lepAccEff<<endl;
   // Double_t acc, eff, accStatErr, effStatErr, acceff, acceffStatErr;

   // mySpaces(cout,2);
   // mySpaces(myfile,2);
   // cout << "Printing acceptance and efficiency." << endl;
   // cout << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;
   // myfile << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;

   // for (Int_t i = 0; i < nMetBins; i++) {
  
   //   HevtPassMonoJetSel->SetBinContent(i+1,lep_acc_eff[i]->getEvents(stepMonojetSelection_In_lepAccEff));
   //   HevtPassMonoJetSel->SetBinError(  i+1,lep_acc_eff[i]->getEventsErr(stepMonojetSelection_In_lepAccEff));
   //   HevtPassAccSel->SetBinContent(    i+1,lep_acc_eff[i]->getEvents(stepAcceptance_In_lepAccEff));
   //   HevtPassAccSel->SetBinError(      i+1,lep_acc_eff[i]->getEventsErr(stepAcceptance_In_lepAccEff));
   //   HevtPassEffSel->SetBinContent(    i+1,lep_acc_eff[i]->getEvents(stepEfficiency_In_lepAccEff));
   //   HevtPassEffSel->SetBinError(      i+1,lep_acc_eff[i]->getEventsErr(stepEfficiency_In_lepAccEff));

   //   acc = lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
   //   eff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff];
   //   accStatErr = sqrt(acc * (1 - acc) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
   //   effStatErr = sqrt(eff * (1 - eff) / lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]);
   //   Hacc->SetBinContent(i+1,acc);
   //   Hacc->SetBinError(i+1,accStatErr);
   //   Heff->SetBinContent(i+1,eff);
   //   Heff->SetBinError(i+1,effStatErr);

   //   cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;
   //   myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;

   // }

   // mySpaces(cout,2);
   // mySpaces(myfile,2);
   // cout << "Printing acceptance * efficiency" << endl;
   // cout << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;
   // myfile << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;

   // for (Int_t i = 0; i < nMetBins; i++) {
   //   // do not merge with previous loop: I want to print them after the previous loop because I might copy and paste this output to make acc * eff table
     
   //   acceff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
   //   acceffStatErr = sqrt(acceff * (1 - acceff) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
   //   Hacceff->SetBinContent(i+1,acceff);
   //   Hacceff->SetBinError(i+1,acceffStatErr);

   //   cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
   //   myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
  
   // }

   myfile.close();

   // TEfficiency *Acceptance = new TEfficiency(*HevtPassAccSel,*HevtPassMonoJetSel);
   // TEfficiency *Efficiency = new TEfficiency(*HevtPassEffSel,*HevtPassAccSel);
   // TEfficiency *AccTimesEff = new TEfficiency(*HevtPassEffSel,*HevtPassMonoJetSel);

   // Acceptance->SetUseWeightedEvents();
   // Efficiency->SetUseWeightedEvents();
   // AccTimesEff->SetUseWeightedEvents();

   // TGraphAsymmErrors *grAE_Acc = Acceptance->CreateGraph();
   // TGraphAsymmErrors *grAE_Eff = Efficiency->CreateGraph();
   // TGraphAsymmErrors *grAE_AccTimesEff = AccTimesEff->CreateGraph();

   // grAE_Acc->Write("grAE_Acc");
   // grAE_Eff->Write("grAE_Eff");
   // grAE_AccTimesEff->Write("grAE_AccTimesEff");

   // now get Z(inv) estimate as N_Zvv = N_Zll * R / (A*e), R being BR(Zvv)/BR(Zll) where l is either mu or e (R ~ 6)
   // TH1D *H_BR_ratio = new TH1D("H_BR_ratio",Form("BR(Z#nu#nu/BR(%s)",CONTROL_SAMPLE),nMetBins,metBinEdges);

   // for(Int_t i = 0; i <= nMetBins; i++) {
   //   H_BR_ratio->SetBinContent(i,RATIO_BR_ZINV_ZLL);
   //   H_BR_ratio->SetBinError(i,UNC_RATIO_BR_ZINV_ZLL);
   // }

   // HzvvEstimate->Multiply(HzlljetsYieldsMetBinGenLep,H_BR_ratio);
   // HzvvEstimate->Divide(Hacceff);
   // delete H_BR_ratio; //no need to save it

   // I add overflow bin's content in the last bin for all histograms where that is needed
   // for those histogram filled with Divide() method, it's not done as long as it was already done on the histograms given as
   // argument to the Divide() method
   myAddOverflowInLastBin(HZtoLLRecoPt);
   if (!ISDATA_FLAG) {
     myAddOverflowInLastBin(HZtoLLGenPt);  
     HZtoLLPt_RecoGenRatio->Divide(HZtoLLRecoPt,HZtoLLGenPt);
   }

   for (Int_t i = 0; i < nMetBins; i++) {
     myAddOverflowInLastBin(HZtoLLRecoPt_MetBin[i]);
   }

   
   if (!ISDATA_FLAG) {
     for (Int_t i = 0; i < nMetBins; i++) {
       myAddOverflowInLastBin(HZtoLLGenPt_MetBin[i]);
       HZtoLLPt_RecoGenRatio_MetBin[i]->Divide(HZtoLLRecoPt_MetBin[i],HZtoLLGenPt_MetBin[i]);  
     }
   }

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
     commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &zlljetsControlSampleGenLep,commentInTable);
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &resoAndResponse);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file

}



