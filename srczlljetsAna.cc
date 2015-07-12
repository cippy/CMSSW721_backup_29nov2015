#define zlljetsAna_cxx
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
#include "thresholds.h"
//#include "edimarcoTreeFriend.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef zlljetsAna_cxx

zlljetsAna::zlljetsAna(TTree *tree) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void zlljetsAna::loop()
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJet",1);         // # of jets with pt > 25 && |eta| < 2.5
   fChain->SetBranchStatus("nJet30",1);         // # of jets with pt > 30 && |eta| < 2.4
   fChain->SetBranchStatus("nJet30a",1);       // # of jets with pt > 30 && |eta| < 4.7 
   fChain->SetBranchStatus("Jet_pt",1);  
   fChain->SetBranchStatus("Jet_eta",1);  
 
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_phi",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   fChain->SetBranchStatus("ngenLep",1);
   fChain->SetBranchStatus("genLep_pdgId",1);
   fChain->SetBranchStatus("genLep_pt",1);
   fChain->SetBranchStatus("genLep_eta",1);
   //fChain->SetBranchStatus("m2l",1);  // m(ll)  (I can compute it myself, maybe it's better)
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_motherId",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_mass",1);
   fChain->SetBranchStatus("GenPart_motherIndex",1);

   fChain->SetBranchStatus("met_pt",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);

   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   vector<Double_t> metCut;
   metCut.push_back(250);
   metCut.push_back(300);
   metCut.push_back(350);
   metCut.push_back(400);
   metCut.push_back(500);
   selection njetsC("njetsC",Form("njets <= %i",NJETS),"pt > 30; |eta| < 4.7");   // using nJet30a
   selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),Form("only if njets = %i",NJETS));
   selection jjdphiEmanC("jjdphiEmanC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection invMassC("invMassC",Form("mass in [%i,%i]",DILEPMASS_LOW,DILEPMASS_UP));

#ifdef MUON

   Int_t &nLep10V = nMu10V;

   selection metNoMuC[metCut.size()];
   selection &lepmetC[metCut.size()];
   for (Int_t i = 0; i < metCut.size(); i++) {
     mumetC[i].set(Form("metNoMuC[%i]",i),Form("metNoMu > %3.0lf",metCut.at(i)));
     lepmetC[i] = mumetC[i];
   }
   selection eLooseVetoC("eLooseVetoC","electrons veto");
   selectio &lepLooseVetoC = eLooseVetoC;
   selection twomuonsC("twomuonsC","muons");
   selection &twoLeptonsC = twomuonsC;
   selection twomuLooseC("twomuLooseC","2 loose muons");
   selection &twoLepLooseC = twomuLooseC;
   selection mu1tightIdC("mu1tightIdC","leading muon tight","tight ID + relIso04 (as Emanuele)");
   selection &lep1tightIdC = mu1tightIdC;
   selection twomuTightC("twomuTightC","2 tight muons");
   selection &twoLepTightC = twomuTightC;
   selection mu1ptC("mu1ptC",Form("mu1pt > %i",LEP1PT),"leading muon pt");
   selection &lep1ptC = mu1ptC;
   selection mu2ptC("mu2ptC",Form("mu2pt > %i",LEP2PT),"trailing muon pt");
   selection &lep2ptC = mu2ptC;
   selection mu1etaC("mu1etaC",Form("|mu1eta| < %1.1lf",LEP1ETA),"leading muon eta");  
   selection &lep1etaC = mu1etaC;
   selection mu2etaC("mu2etaC",Form("|mu2eta| < %1.1lf",LEP2ETA),"trailing muon eta");
   selection &lep2etaC = mu2etaC;
   selection genMuonsC("genMuonsC","muons generated");     
   selection &genLepC = genMuonsC;

#elseif ELECTRON

   Int_t &nLep10V = nEle10V;

   selection metNoEleC[metCut.size()];
   selection &lepmetC[metCut.size()];
   for (Int_t i = 0; i < metCut.size(); i++) {
     elemetC[i].set(Form("metNoEleC[%i]",i),Form("metNoEle > %3.0lf",metCut.at(i)));
     lepmetC[i] = elemetC[i];
   }
   selection muLooseVetoC("muLooseVetoC","muons veto");
   selection &lepLooseVetoC = muLooseVetoC;
   selection twoelectronsC("twoelectronsC","electrons");
   selection &twoLeptonsC = twoelectronsC;
   selection twoeleLooseC("twoeleLooseC","2 loose electrons");
   selection &twoLepLooseC = twoeleLooseC;
   selection ele1tightIdC("ele1tightIdC","leading electron tight","tight ID + relIso04 (as Emanuele)");
   selection &lep1tightIdC = ele1tightIdC;
   selection twoeleTightC("twoeleTightC","2 tight electrons");
   selection &twoLepTightC = twoeleTightC;
   selection ele1ptC("ele1ptC","ele1pt > 20","leading pt electron");
   selection ele1ptC("ele1ptC",Form("ele1pt > %i",LEP1PT),"leading electron pt");
   selection &lep1ptC = ele1ptC;
   selection ele2ptC("ele2ptC",Form("ele2pt > %i",LEP2PT),"trailing electron pt");
   selection &lep2ptC = ele2ptC;
   selection ele1etaC("ele1etaC",Form("|ele1eta| < %1.1lf",LEP1ETA),"leading electron eta");  
   selection &lep1etaC = ele1etaC;
   selection ele2etaC("ele2etaC",Form("|ele2eta| < %1.1lf",LEP2ETA),"trailing electron eta");
   selection genElectronsC("genElectronsC","electrons generated");     
   selection &genLepC = genElectronsC;

   // the following are only for electrons
   selection MetNoEle200C("MetNoEle200C","MetNoEle > 200");
   selection ele2tightIdC("ele2tightIdC","trailing electron tight","tight ID + relIso04 (as Emanuele)");

#endif
   selection genTausC("genTausC","taus generated");                       
   selection acceptanceC("acceptanceC","acceptance cuts");
   selection efficiencyC("efficiencyC","efficiency cuts");

   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   UInt_t maskMonoJetSelection = njetsEmanC.get2ToId() + jet1ptC.get2ToId() + jjdphiEmanC.get2ToId() +
                                   lepLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   mask lep_Acc_Eff(Form("%s acceptance and efficiency",LL_FLAVOUR)); 
   lep_Acc_Eff.append(genLepC.get2ToId());
   lep_Acc_Eff.append(maskMonoJetSelection);
   lep_Acc_Eff.append(acceptanceC.get2ToId());
   lep_Acc_Eff.append(efficiencyC.get2ToId());

   mask zlljetsControlSample(Form("%s control sample with selection flow as Emanuele's",CONTROL_SAMPLE));
   zlljetsControlSample.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zlljetsControlSample.append(twoLeptonsC.get2ToId());
   zlljetsControlSample.append(lep1tightIdC.get2ToId());
   zlljetsControlSample.append(invMassC.get2ToId());
   zlljetsControlSample.append(njetsEmanC.get2ToId());
   zlljetsControlSample.append(jet1ptC.get2ToId());
   zlljetsControlSample.append(jjdphiEmanC.get2ToId());
   zlljetsControlSample.append(lepLooseVetoC.get2ToId());
   zlljetsControlSample.append(gammaLooseVetoC.get2ToId());

   mask zlljetsControlSampleGenLep(Form("%s control sample (%s gen ) with selection flow as Emanuele's",CONTROL_SAMPLE,FLAVOUR));
   zlljetsControlSampleGenLep.append(genLepC.get2ToId());
   zlljetsControlSampleGenLep.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zlljetsControlSampleGenLep.append(twoLeptonsC.get2ToId());
   zlljetsControlSampleGenLep.append(lep1tightIdC.get2ToId());
   zlljetsControlSampleGenLep.append(invMassC.get2ToId());
   zlljetsControlSampleGenLep.append(njetsEmanC.get2ToId());
   zlljetsControlSampleGenLep.append(jet1ptC.get2ToId());
   zlljetsControlSampleGenLep.append(jjdphiEmanC.get2ToId());
   zlljetsControlSampleGenLep.append(lepLooseVetoC.get2ToId());
   zlljetsControlSampleGenLep.append(gammaLooseVetoC.get2ToId());

   mask tautaubkgInZll(Form("tau tau background in %s control sample",CONTROL_SAMPLE));
   tautaubkgInZll.append(genTausC.get2ToId());
   tautaubkgInZll.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   tautaubkgInZll.append(twoleptonsC.get2ToId());
   tautaubkgInZll.append(lep1tightIdC.get2ToId());
   tautaubkgInZll.append(invMassC.get2ToId());
   tautaubkgInZll.append(njetsEmanC.get2ToId());
   tautaubkgInZll.append(jet1ptC.get2ToId());
   tautaubkgInZll.append(jjdphiEmanC.get2ToId());
   tautaubkgInZll.append(lepLooseVetoC.get2ToId());
   tautaubkgInZll.append(gammaLooseVetoC.get2ToId());

   mask *lep_acc_eff[nMetBins];
   for ( Int_t i = 0; i < nMetBins; i++) {
     lep_acc_eff[i] = new mask;
     lep_acc_eff[i]->setName(Form("%s_acc_eff:  %3.0lf < met < %3.0lf",FLAVOUR,metBinEdges[i], metBinEdges[i+1]));
     lep_acc_eff[i].append(genLepC.get2ToId());
     lep_acc_eff[i].append(maskMonoJetSelection);
     lep_acc_eff[i].append(acceptanceC.get2ToId());
     lep_acc_eff[i].append(efficiencyC.get2ToId());
   }

   string name;

   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};    

   TH1D *HinvMass[nMetBins];
   TH1D *HzlljetsInvMassMetBinGenLep[nMetBins];
   TH1D *HzlljetsInvMassMetBinGenTau[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HinvMass[i] = new TH1D(Form("HinvMass[%i]",i),"",30,60.,120.);
     HzlljetsInvMassMetBinGenLep[i] = new TH1D(Form("HzlljetsInvMassMetBinGenLep[%i]",i),"",30,60.,120.);
     HzlljetsInvMassMetBinGenTau[i] = new TH1D(Form("HzlljetsInvMassMetBinGenTau[%i]",i),"",30,60.,120.);
   } 

   TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin",Form("yields of %s control sample in bins of met;#slash{E}_{T};# of events",
								     CONTROL_SAMPLE), nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinGenLep = new TH1D("HzlljetsYieldsMetBinGenLep",Form("yields of %s control sample (Z->#mu#mu gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinGenTau = new TH1D("HzlljetsYieldsMetBinGenTau",Form("yields of %s control sample (Z->#tau#tau gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);

   TH1D *HZtoLLRecoPt = new TH1D("HZtoLLRecoPt","",101,0.,1010);
   TH1D *HZtoLLGenPt = new TH1D("HZtoLLGenPt","",101,0.,1010);
   // this is the histogram with reco/gen
   TH1D *HZtoLLPt_RecoGenRatio = new TH1D("HZtoLLPt_RecoGenRatio","",101,0.,1010.);
   // histogram of reco/gen distribution function
   TH1D *HZtoLLPt_RecoGenRatio_pdf = new TH1D("HZtoLLPt_RecoGenRatio_pdf","",100,0.,2.);
   TH1D *HZtoLLRecoPt_MetBin[nMetBins];
   TH1D *HZtoLLGenPt_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_pdf_MetBin[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HZtoLLRecoPt_MetBin[i] = new TH1D(Form("HZtoLLRecoPt_MetBin[%i]",i),"",101,0.,1010.);
     HZtoLLGenPt_MetBin[i] = new TH1D(Form("HZtoLLGenPt_MetBin[%i]",i),"",101,0.,1010.);
     HZtoLLPt_RecoGenRatio_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_MetBin[%i]",i),"",101,0.,1010.);
     HZtoLLPt_RecoGenRatio_pdf_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_pdf_MetBin[%i]",i),"",100,0.,2.);
   } 

   TVector2 l1, l2;

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljetsAna::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     // I find the indices corresponding to the 2 leading lepton
     //cout<<"entry : "<<jentry<<endl;
     myGetPairIndexInArray(MU_PDGID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);

     mu1.SetMagPhi(LepGood_pt[firstIndex],LepGood_phi[firstIndex]);
     mu2.SetMagPhi(LepGood_pt[secondIndex],LepGood_phi[secondIndex]);

     nTotalWeightedEvents += newwgt;  // counting events with weights

     eventMask += njetsC.addToMask(nJet30a <= NJETS);
     eventMask += njetsEmanC.addToMask( ((nJet30a == 1 ) || (nJet30a == 2 && Jet_eta[1] < JET2ETA)) && Jet_eta[0] < JET1ETA && jetclean > 0.5);
     eventMask += jjdphiEmanC.addToMask( (nJet30a == 1 && Jet_eta[0] < 2.5) || (nJet30a == 2 && abs(dphijj) < J1J2DPHI));
     eventMask += jet1ptC.addToMask(Jet_pt[0] > J1PT);               
     //in Emanuele's tree we have vectors: [0] is the first jet, [1] is the second and so on (ordered in pt)
     eventMask += jet1etaC.addToMask(fabs(Jet_eta[0]) < J1ETA);
     eventMask += jetCleaningC.addToMask(jetclean > 0.5);     
     // jetclean is 1 if cleaning is passed, 0 otherwise. It's applied to first jet and , if any, to the second
     if (nJet30a == 2) {
       eventMask += jet2etaC.addToMask(fabs(Jet_eta[1]) < J2ETA);
       eventMask += jet1jet2dphiC.addToMask(fabs(dphijj) < J1J2DPHI);
     } else if (nJet30a == 1) {
       eventMask += jet2etaC.get2ToId();
       eventMask += jet1jet2dphiC.get2ToId();
     }
     eventMask += eLooseVetoC.addToMask(nEle10V == 0);
     //eventMask += ntausC.addToMask(ntaus);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     //eventMask += mumet200C.addToMask(mumet);
     for (Int_t i = 0; i <  metCut.size(); i++) {
       eventMask += mumetC[i].addToMask(metNoMu_pt > metCut[i]);
     }
     eventMask += oppChargeMuonsC.addToMask( (LepGood_pdgId[0] * LepGood_pdgId[1]) == -169);  // |pdgID| = 13 for muons
     eventMask += oppChargeLeptonsC.addToMask( (LepGood_pdgId[0] + LepGood_pdgId[1]) == 0);
     eventMask += twomuonsC.addToMask((fabs(LepGood_pdgId[0]) == 13) && (fabs(LepGood_pdgId[1]) == 13));
     eventMask += twomuLooseC.addToMask(nMu10V == 2);
     eventMask += muLooseVetoC.addToMask(nMu10V == 0);
     eventMask += mu1tightC.addToMask((LepGood_tightId[0] == 1) && (fabs(LepGood_pdgId[0]) == 13) && (LepGood_relIso04[0] < 0.12 ) && (fabs(LepGood_eta[0]) < 2.4) && (LepGood_pt[0] > 20));
     eventMask += mu1tightIdC.addToMask((LepGood_tightId[0] == 1) && (LepGood_relIso04[0] < 0.12 ) && (fabs(LepGood_pdgId[0]) == 13));
     eventMask += mu1ptC.addToMask((LepGood_pt[0] > 20) && (fabs(LepGood_pdgId[0]) == 13)); 
     eventMask += mu1etaC.addToMask( (fabs(LepGood_eta[0]) < 2.4) && (fabs(LepGood_pdgId[0]) == 13) );
     eventMask += twomuTightC.addToMask(nMu20T == 2);
     eventMask += mu2ptC.addToMask((LepGood_pt[1] > 10) && (fabs(LepGood_pdgId[1]) == 13));
     eventMask += mu2etaC.addToMask((fabs(LepGood_eta[1]) < 2.4) && (fabs(LepGood_pdgId[1]) == 13));
     eventMask += mumuInvMassC.addToMask((mZ1 > 60) && (mZ1 < 120));
     // eventMask += genMuonsC.addToMask( ((genLep_pdgId[0] * genLep_pdgId[1]) == -169) && (genLep_motherId[0] == 23) && (genLep_motherId[1] == 23) );
     eventMask += genMuonsC.addToMask(myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 13, 23) );
     eventMask += genTausC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23) );  

     mu_Acc_Eff.countEvents(eventMask,newwgt);
     zmumujetsControlSample.countEvents(eventMask,newwgt);
     zmumujetsControlSampleGenMu.countEvents(eventMask,newwgt);
     tautaubkgInZmumu.countEvents(eventMask, newwgt);

     Double_t ZtoMuMuRecoPt = (mu1 + mu2).Mod();
     // since when I use the following index I will ask for 2 mu from Z, it's enough to use directly the genZ instead of 2 genMu
     // thus I look for a Z whose pdgId is 23
     Int_t Z_index = myGetPartIndex(23,nGenPart,GenPart_pdgId);
     Double_t ZtoMuMuGenPt = GenPart_pt[Z_index];

     // filling histogram with yields and invariant mass at the end of the selection in bins of met
     if ( ((eventMask & zmumujetsControlSample.globalMask.back()) == zmumujetsControlSample.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBin->Fill(metNoMu_pt,newwgt);     
     }
     if ( ((eventMask & zmumujetsControlSampleGenMu.globalMask.back()) == zmumujetsControlSampleGenMu.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBinGenMu->Fill(metNoMu_pt,newwgt);
       HZtoMuMuRecoPt->Fill(ZtoMuMuRecoPt,newwgt);
       HZtoMuMuGenPt->Fill(ZtoMuMuGenPt,newwgt);
       if (ZtoMuMuGenPt != 0) HZtoMuMuPt_RecoGenRatio_pdf->Fill(ZtoMuMuRecoPt/ZtoMuMuGenPt,newwgt);
     }
     if ( ((eventMask & tautaubkgInZmumu.globalMask.back()) == tautaubkgInZmumu.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBinGenTau->Fill(metNoMu_pt,newwgt);  
     }

     if ((metNoMu_pt > metBinEdges[0]) && (metNoMu_pt < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoMu_pt,metBinEdges,nMetBins);
       mu_acc_eff[bin]->countEvents(eventMask,newwgt);
       //find step where cut on invariant mass is added the first time
       //Int_t index = zmumujetsControlSample.whichStepHas(mumuInvMassC.get2ToId()); 
       //if (((eventMask & zmumujetsControlSample.globalMask[index]) == zmumujetsControlSample.globalMask[index]) && (index < zmumujetsControlSample.getMaskSize())) {
       if ((eventMask & zmumujetsControlSample.globalMask.back()) == zmumujetsControlSample.globalMask.back()) {
	 // this histogram holds the invariant mass distribution (one for each met bin)
   	 HmumuInvMass[bin]->Fill(mZ1,newwgt);   
       }
       if ( ((eventMask & zmumujetsControlSampleGenMu.globalMask.back()) == zmumujetsControlSampleGenMu.globalMask.back()) ) {  
	 HzmumujetsInvMassMetBinGenMu[bin]->Fill(mZ1,newwgt); 
	 HZtoMuMuRecoPt_MetBin[bin]->Fill(ZtoMuMuRecoPt,newwgt);
	 HZtoMuMuGenPt_MetBin[bin]->Fill(ZtoMuMuGenPt,newwgt);
	 if (ZtoMuMuGenPt != 0) HZtoMuMuPt_RecoGenRatio_pdf_MetBin[bin]->Fill(ZtoMuMuRecoPt/ZtoMuMuGenPt,newwgt);
       }
       if ( ((eventMask & tautaubkgInZmumu.globalMask.back()) == tautaubkgInZmumu.globalMask.back()) ) {  
       	 HzmumujetsInvMassMetBinGenTau[bin]->Fill(mZ1,newwgt);   
       }
     
     } 

   }

   cout<<endl;   
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &mu_Acc_Eff);
   for (Int_t i = 0; i < nMetBins; i++ ) {
     selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, mu_acc_eff[i] );
   }
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zmumujetsControlSampleGenMu);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &tautaubkgInZmumu);
   cout<<endl;   

   for (Int_t i = 0; i < nMetBins; i++) {
     cout<<" mumet in ["<<metBinEdges[i]<<" , "<<metBinEdges[i+1]<<"] :     HmumuInvMass["<<i<<"]->GetSumOfWeights() = ";
     cout<<HmumuInvMass[i]->GetSumOfWeights()<<endl;
   }

   cout << "Printing acceptance and efficiency." << endl;
   TH1D* Hacc = new TH1D("Hacc","",nMetBins,metBinEdges);
   TH1D* Heff = new TH1D("Heff","",nMetBins,metBinEdges);
   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I will copy and paste this output to make acc & eff table
     Double_t acc = mu_acc_eff[i]->nEvents[2]/mu_acc_eff[i]->nEvents[1];
     Double_t eff = mu_acc_eff[i]->nEvents[3]/mu_acc_eff[i]->nEvents[2];
     Double_t accStatErr = sqrt(acc * (1 - acc) / mu_acc_eff[i]->nEvents[1]);
     Double_t effStatErr = sqrt(eff * (1 - eff) / mu_acc_eff[i]->nEvents[2]);
     Hacc->SetBinContent(i+1,acc);
     Hacc->SetBinError(i+1,accStatErr);
     Heff->SetBinContent(i+1,eff);
     Heff->SetBinError(i+1,effStatErr);
     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     cout<<acc<<" "<<accStatErr<<" ";
     cout<<eff<<" "<<effStatErr<<endl;
   }
   Hacc->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Hacc->GetYaxis()->SetTitle("Acceptance");
   Hacc->GetYaxis()->CenterTitle();
   Heff->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Heff->GetYaxis()->SetTitle("efficiency");
   Heff->GetYaxis()->CenterTitle();

   TH1D* Hacceff = new TH1D("Hacceff","",nMetBins,metBinEdges);
   cout << "Printing acceptance * efficiency" << endl;
   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I will copy and paste this output to make acc * eff table
     Double_t acceff = mu_acc_eff[i]->nEvents[3]/mu_acc_eff[i]->nEvents[1];
     Double_t acceffStatErr = sqrt(acceff * (1 - acceff) / mu_acc_eff[i]->nEvents[1]);
     Hacceff->SetBinContent(i+1,acceff);
     Hacceff->SetBinError(i+1,acceffStatErr);
     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     cout<<acceff<<" "<<acceffStatErr<<endl;;
   }
   Hacceff->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Hacceff->GetYaxis()->SetTitle("A #times #epsilon");
   Hacceff->GetYaxis()->CenterTitle();

   name = "hist_"+ ll_flavour +"_AccEff_metBin.root";
   Hacceff->SaveAs(name.c_str());

   name = "hist_" + ll_flavour +"_accEff_all_metBin.root";
   TFile *accAndEff_file = new TFile(name.c_str(),"RECREATE");
   if (!accAndEff_file->IsOpen()) {
     cout<<"Error: file \""<<name<<"\" was not opened."<<endl;
   } else {
     for (Int_t i = 0; i < nMetBins; i++) {
       Hacc->Write();
       Heff->Write();
       Hacceff->Write();
     }
     accAndEff_file->Close();
   }
   delete accAndEff_file;

   char answer = '\0';
   name = "z"+ ll_flavour +"jetsAnaYields.txt";

   answer = myAskToSaveFile(name.c_str());

   if (answer == 'y') {
     
     ofstream myfile(name.c_str(),ios::app);

     if ( !myfile.is_open() ) {
       cout<<"Error: unable to open file "<<name.c_str()<<" !"<<endl;

     } else {
       string command = "date>>" + name;
       //when writing on file, the date is printed as well unless an error occurs
       if ( system(command.c_str()) != 0) cout<<"Error during \"system("<<command<<")\" call"<<endl;  
       myfile<<endl;
       myfile<<endl;
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &mu_Acc_Eff);      
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zmumujetsControlSampleGenMu);
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &tautaubkgInZmumu);
       myfile<<endl;
       myfile<<endl;
       myfile.close();

     }
 
   }

   name = "z"+ ll_flavour + "jetsAnaYields.tex";
   
   answer = myAskToSaveFile(name.c_str());

   if (answer == 'y') {

     //creating a .tex file to build tables with data
     FILE *fp;	 
     //cout<<"name given to file: "<<texfname<<endl;
     if ( (fp=fopen(name.c_str(),"w")) == NULL) {
       cout<<"Error: '"<<name.c_str()<<"' not opened"<<endl;
     } else {
       cout<<"creating file '"<<name.c_str()<<"' ..."<<endl;
       myAddDefaultPackages(fp,name.c_str());
       fprintf(fp,"\\begin{document}\n");
       fprintf(fp,"\n");
       string commentInTable;       
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &mu_Acc_Eff, commentInTable);      
       commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zmumujetsControlSample,commentInTable);
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zmumujetsControlSampleGenMu,commentInTable);
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &tautaubkgInZmumu,commentInTable);
       fprintf(fp,"\\end{document}\n");      
       fclose(fp);
     }

   }

   name = "histZ" + ll_flavour + "jetsAnaInvMass.root"
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile = new TFile(name.c_str(),"RECREATE");

   if (!histFile->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HmumuInvMass[i]->Write();
     }

     histFile->Close();
     
   }

   delete histFile;

   name = "histZ"+ ll_flavour +"jetsInvMassCS_ZJetsTo" + LL_flavour + ".root"
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile2 = new TFile(name.c_str(),"RECREATE");

   if (!histFile2->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HzmumujetsInvMassMetBinGenMu[i]->Write();
     }

     histFile2->Close();
     
   }

   delete histFile2;


   name = "histZ"+ ll_flavour +"jetsInvMassCS_ZJetsToTauTau.root";
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile3 = new TFile(name.c_str(),"RECREATE");

   if (!histFile3->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HzmumujetsInvMassMetBinGenTau[i]->Write();
     }

     histFile3->Close();
     
   }

   delete histFile3;				       
				  
   name = "histZ"+ ll_flavour +"jetsAnaYieldsMetBin.root";
   cout<<"Saving histogram \""<<HzmumujetsYieldsMetBin->GetName()<<"\" in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *YieldsFile = new TFile(name.c_str(),"RECREATE");

   if (!YieldsFile->IsOpen()) cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;
   else HzmumujetsYieldsMetBin->Write();

   cout<<"MET < "<<metBinEdges[0]<<" : yield = ";
   cout<<HzmumujetsYieldsMetBin->GetBinContent(0)<<" +/- "<<HzmumujetsYieldsMetBin->GetBinError(0)<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     cout<<"MET in ["<<metBinEdges[i]<<","<<metBinEdges[i+1]<<"] : yield = ";
     cout<<HzmumujetsYieldsMetBin->GetBinContent(i+1)<<" +/- "<<HzmumujetsYieldsMetBin->GetBinError(i+1)<<endl;
   }
   cout<<"MET > "<<metBinEdges[nMetBins]<<" : yield = ";
   cout<<HzmumujetsYieldsMetBin->GetBinContent(nMetBins + 1)<<" +/- "<<HzmumujetsYieldsMetBin->GetBinError(nMetBins + 1)<<endl;

   YieldsFile->Close();
     
   delete YieldsFile;

   name = "histZ" +ll_flavour+" jetsAnaYieldsMetBinGen"+ L_flavour +".root"
   HzmumujetsYieldsMetBinGenMu->SaveAs(name.c_str());
   cout<<"HzmumujetsYieldsMetBinGenMu: yield per bin"<<endl;
   cout<<"<"<<metBinEdges[0]<<" ";
   cout<<HzmumujetsYieldsMetBinGenMu->GetBinContent(0)<<" "<<HzmumujetsYieldsMetBinGenMu->GetBinError(0)<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     cout<<metBinEdges[i]<<"-"<<metBinEdges[i+1]<<" ";
     cout<<HzmumujetsYieldsMetBinGenMu->GetBinContent(i+1)<<" "<<HzmumujetsYieldsMetBinGenMu->GetBinError(i+1)<<endl;
   }
   cout<<">"<<metBinEdges[nMetBins]<<" ";
   cout<<HzmumujetsYieldsMetBinGenMu->GetBinContent(nMetBins + 1)<<" "<<HzmumujetsYieldsMetBinGenMu->GetBinError(nMetBins + 1)<<endl;

   name = "histZ" +ll_flavour+" jetsAnaYieldsMetBinGenTau.root"
     HzmumujetsYieldsMetBinGenTau->SaveAs(name.c_str());

   // I add overflow bin's content in the last bin for all histograms where that is needed
   // for those histogram filled with Divide() method, it's not done as long as it was already done on the histograms given as
   // argument to the Divide() method
   myAddOverflowInLastBin(HZtoMuMuRecoPt);
   myAddOverflowInLastBin(HZtoMuMuGenPt);

   // saving results about PtZReco/PtZGen
   name = "histZ" + ll_flavour + "RecoPt.root";
   HZtoMuMuRecoPt->SaveAs(name.c_str());
   name = "histZ" + ll_flavour + "GenPt.root";
   HZtoMuMuGenPt->SaveAs(name.c_str());
   HZtoMuMuPt_RecoGenRatio->Divide(HZtoMuMuRecoPt,HZtoMuMuGenPt);
   name = "histZ" + ll_flavour + "Pt_RecoGenRatio.root";
   HZtoMuMuPt_RecoGenRatio->SaveAs(name.c_str());
   name = "histZ" + ll_flavour + "Pt_RecoGenRatio_pdf.root";
   HZtoMuMuPt_RecoGenRatio_pdf->SaveAs(name.c_str());

   // save binned histograms in another file
   name = "histZ"+ ll_flavour +"Pt_RecoGenRatio_MetBin.root";
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile4 = new TFile(name.c_str(),"RECREATE");

   if (!histFile4->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       myAddOverflowInLastBin(HZtoMuMuRecoPt_MetBin[i]);
       myAddOverflowInLastBin(HZtoMuMuGenPt_MetBin[i]);
       HZtoMuMuPt_RecoGenRatio_MetBin[i]->Divide(HZtoMuMuRecoPt_MetBin[i],HZtoMuMuGenPt_MetBin[i]);
       HZtoMuMuPt_RecoGenRatio_MetBin[i]->Write();
     }

     histFile4->Close();
     
   }

   delete histFile4;

   name = "histZ" + ll_flavour + "Pt_RecoGenRatio_pdf_MetBin.root";
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile5 = new TFile(name.c_str(),"RECREATE");

   if (!histFile5->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoMuMuPt_RecoGenRatio_pdf_MetBin[i]->Write();
     }

     histFile5->Close();
     
   }

   delete histFile5;

   name  = "histZ" + ll_flavour + "RecoPt_MetBin.root";
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile6 = new TFile(name.c_str(),"RECREATE");

   if (!histFile6->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoMuMuRecoPt_MetBin[i]->Write();
     }

     histFile6->Close();
     
   }

   delete histFile6;

   name = "histZ" + ll_flavour + "GenPt_MetBin.root";
   cout<<"Saving histograms in file \""<<name.c_str()<<"\" ..."<<endl;
   TFile *histFile7 = new TFile(name.c_str(),"RECREATE");

   if (!histFile7->IsOpen()) {

     cout<<"Error: file \""<<name.c_str()<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoMuMuGenPt_MetBin[i]->Write();
     }

     histFile7->Close();
     
   }

   delete histFile7;

}
