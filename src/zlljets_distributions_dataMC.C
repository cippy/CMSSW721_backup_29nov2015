//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

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
#include <iomanip> 

#include <Rtypes.h> // to use kColor

using namespace std;

void zlljets_distributions_dataMC(const char* suffix = "_GLTcutZmass80to100_mumu", const Int_t lepFlavour = 0, const Int_t MCuncBand_flag = 0) {

  // in this macro lepFlavour = 0 means muons, while 1 means electrons

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/distributions/dataMC_comparison/";
  string plotFileExtension = ".pdf";

  vector<string> MC_TexLabel;
  MC_TexLabel.push_back("Z(#tau#tau)+jets MC");
  MC_TexLabel.push_back("W(l#nu)+jets MC");
  MC_TexLabel.push_back("TT+Jets MC");
  if (lepFlavour == 0) MC_TexLabel.push_back("Z(#mu#mu)+jets MC");
  else if (lepFlavour == 1) MC_TexLabel.push_back("Z(ee)+jets MC");

  string data_TexLabel = "data"; 

  string DataFileName; 
  if (lepFlavour == 0) DataFileName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";
  else if (lepFlavour == 1) DataFileName = "zeejets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";

  vector<TH1D*> hMCinvMass;
  vector<TH1D*> hMCmetNoLep;
  vector<TH1D*> hMCjet1Pt;
  vector<TH1D*> hMCuParMinusZpt;
  vector<TH1D*> hMCuPerp;
  vector<TH1D*> hMCnVertices;
  vector<TH1D*> hMCzPtSpectrum;
  vector<TH1D*> hMCnjets;

  TH1D* hDatainvMass;
  TH1D* hDatametNoLep;
  TH1D* hDatajet1Pt;
  TH1D* hDatauParMinusZpt;
  TH1D* hDatauPerp;
  TH1D* hDatanVertices;
  TH1D* hDatazPtSpectrum;
  TH1D* hDatanjets;

  Int_t invMassLogy_flag = 0;
  Int_t metLogy_flag = 1;
  Int_t j1ptLogy_flag = 1;
  Int_t uParLogy_flag = 1;
  Int_t uPerpLogy_flag = 1;
  Int_t nvtxLogy_flag = 0;
  Int_t zptLogy_flag = 1;
  Int_t njetsLogy_flag = 1;

  vector<string> MCfileName;
  if (lepFlavour == 0) {
    MCfileName.push_back("ztautaujets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("wjetslnu_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
  } else if (lepFlavour == 1) {
    MCfileName.push_back("ztautaujets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("wjetslnu_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("ttjets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("zeejets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
  }

  Int_t nFiles = (Int_t)MCfileName.size();
  Int_t histColor[] = {kYellow,kRed,kBlue,kCyan};

  for(Int_t i = 0; i < nFiles; i++) {
    cout<<"fileName : "<<MCfileName[i]<<endl;

    TFile* f = TFile::Open(MCfileName[i].c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<MCfileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    
    hMCinvMass.push_back((TH1D*)f->Get("HinvMass_NoJetCuts"));
    hMCmetNoLep.push_back((TH1D*)f->Get("HmetNoLepDistribution_NoJetCuts"));
    hMCjet1Pt.push_back((TH1D*)f->Get("Hjet1ptDistribution_NoJetCuts"));
    hMCuParMinusZpt.push_back((TH1D*)f->Get("H_uParMinusZpT_Distribution"));
    hMCuPerp.push_back((TH1D*)f->Get("H_uPerp_Distribution"));
    hMCnVertices.push_back((TH1D*)f->Get("HvtxDistribution_NoJetCuts"));
    hMCzPtSpectrum.push_back((TH1D*)f->Get("HzptDistribution_NoJetCuts"));
    hMCnjets.push_back((TH1D*)f->Get("HnjetsDistribution_NoJetCuts"));
   
  }

  cout<<"fileName : "<< DataFileName<<endl;
  TFile* f = TFile::Open(DataFileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<DataFileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  
  hDatainvMass = (TH1D*)f->Get("HinvMass_NoJetCuts");
  hDatametNoLep = (TH1D*)f->Get("HmetNoLepDistribution_NoJetCuts");
  hDatajet1Pt = (TH1D*)f->Get("Hjet1ptDistribution_NoJetCuts");
  hDatauParMinusZpt = (TH1D*)f->Get("H_uParMinusZpT_Distribution");
  hDatauPerp = (TH1D*)f->Get("H_uPerp_Distribution");
  hDatanVertices = (TH1D*)f->Get("HvtxDistribution_NoJetCuts");
  hDatazPtSpectrum = (TH1D*)f->Get("HzptDistribution_NoJetCuts");
  hDatanjets = (TH1D*)f->Get("HnjetsDistribution_NoJetCuts");

  THStack* hstack_invMass = new THStack("hstack_invMass",""); 
  THStack* hstack_uPerp = new THStack("hstack_uPerp",""); 
  THStack* hstack_uParMinusZpt = new THStack("hstack_uParMinusZpt",""); 
  THStack* hstack_nVertices = new THStack("hstack_nVertices","");
  THStack* hstack_jet1Pt = new THStack("hstack_jet1Pt","");
  THStack* hstack_zPtSpectrum = new THStack("hstack_zPtSpectrum","");
  THStack* hstack_metNoLep = new THStack("hstack_metNoLep","");
  THStack* hstack_njets = new THStack("hstack_njets","");

  for (Int_t j = 0; j < nFiles; j++) {

    // set bin error on MC to sqrt(bincontent) for DYJetsToLL and to 50% to others

    if ( j == (nFiles-1)) {

      for (Int_t i = 1; i <= hMCinvMass[j]->GetNbinsX(); i++) {
	hMCinvMass[j]->SetBinError(i,sqrt(hMCinvMass[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCuPerp[j]->GetNbinsX(); i++) {
	hMCuPerp[j]->SetBinError(i,sqrt(hMCuPerp[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCuParMinusZpt[j]->GetNbinsX(); i++) {
	hMCuParMinusZpt[j]->SetBinError(i,sqrt(hMCuParMinusZpt[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCnVertices[j]->GetNbinsX(); i++) {
	hMCnVertices[j]->SetBinError(i,sqrt(hMCnVertices[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCjet1Pt[j]->GetNbinsX(); i++) {
	hMCjet1Pt[j]->SetBinError(i,sqrt(hMCjet1Pt[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCzPtSpectrum[j]->GetNbinsX(); i++) {
	hMCzPtSpectrum[j]->SetBinError(i,sqrt(hMCzPtSpectrum[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {
	hMCmetNoLep[j]->SetBinError(i,sqrt(hMCmetNoLep[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCnjets[j]->GetNbinsX(); i++) {
	hMCnjets[j]->SetBinError(i,sqrt(hMCnjets[j]->GetBinContent(i)));
      }

    } else {

      for (Int_t i = 1; i <= hMCinvMass[j]->GetNbinsX(); i++) {
	hMCinvMass[j]->SetBinError(i,(hMCinvMass[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCuPerp[j]->GetNbinsX(); i++) {
	hMCuPerp[j]->SetBinError(i,(hMCuPerp[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCuParMinusZpt[j]->GetNbinsX(); i++) {
	hMCuParMinusZpt[j]->SetBinError(i,(hMCuParMinusZpt[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCnVertices[j]->GetNbinsX(); i++) {
	hMCnVertices[j]->SetBinError(i,(hMCnVertices[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCjet1Pt[j]->GetNbinsX(); i++) {
	hMCjet1Pt[j]->SetBinError(i,(hMCjet1Pt[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCzPtSpectrum[j]->GetNbinsX(); i++) {
	hMCzPtSpectrum[j]->SetBinError(i,(hMCzPtSpectrum[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {
	hMCmetNoLep[j]->SetBinError(i,(hMCmetNoLep[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCnjets[j]->GetNbinsX(); i++) {
	hMCnjets[j]->SetBinError(i,(hMCnjets[j]->GetBinContent(i)/2.));
      }

    }

    hMCinvMass[j]->SetFillColor(histColor[j]);
    hstack_invMass->Add(hMCinvMass[j]);
    hMCuPerp[j]->SetFillColor(histColor[j]);
    hstack_uPerp->Add(hMCuPerp[j]);
    hMCuParMinusZpt[j]->SetFillColor(histColor[j]);
    hstack_uParMinusZpt->Add(hMCuParMinusZpt[j]);
    hMCnVertices[j]->SetFillColor(histColor[j]);
    hstack_nVertices->Add(hMCnVertices[j]);
    hMCjet1Pt[j]->SetFillColor(histColor[j]);
    hstack_jet1Pt->Add(hMCjet1Pt[j]);
    hMCzPtSpectrum[j]->SetFillColor(histColor[j]);
    hstack_zPtSpectrum->Add(hMCzPtSpectrum[j]);
    hMCmetNoLep[j]->SetFillColor(histColor[j]);
    hstack_metNoLep->Add(hMCmetNoLep[j]);
    hMCnjets[j]->SetFillColor(histColor[j]);
    hstack_njets->Add(hMCnjets[j]);

  }

  TH1D * ratio = NULL; // will use it for the ratio plots

  TPad *subpadMass_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadMass_2 = NULL; 

  TCanvas *cInvMass = new TCanvas("invMass","",700,700);
  TLegend *legInvMass = new TLegend(0.7,0.7,0.89,0.89);

  subpadMass_1 = new TPad("padMass_1","",0.0,0.28,1.0,1.0);
  if (invMassLogy_flag) subpadMass_1->SetLogy();
  //subpadMass_1->SetBottomMargin(0);
  subpadMass_2 = new TPad("padMass_2","",0.0,0.0,1.0,0.27);
  subpadMass_2->SetGridy();
  //subpadMass_2->SetTopMargin(0);
  subpadMass_2->SetBottomMargin(0.3);
  subpadMass_1->Draw();
  subpadMass_2->Draw();

  subpadMass_1->cd();
  hstack_invMass->Draw("HIST");
  //hstack_invMass->SetMinimum(0.3);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_invMass->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatainvMass->SetMarkerStyle(8); // large dot
  hDatainvMass->Draw("EX0 SAME"); //X0 doesn't draw x error
  if (lepFlavour == 0) hstack_invMass->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  else if (lepFlavour == 1) hstack_invMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hstack_invMass->GetXaxis()->SetTitleSize(0.04);
  hstack_invMass->GetYaxis()->SetTitle("events / 1 GeV");
  hstack_invMass->GetYaxis()->SetTitleSize(0.04);
  hstack_invMass->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legInvMass->AddEntry(hMCinvMass[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legInvMass->AddEntry(hDatainvMass,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legInvMass->Draw(); 
  legInvMass->SetMargin(0.3); 
  legInvMass->SetBorderSize(0);
  subpadMass_2->cd();
  ratio = new TH1D(*hDatainvMass);
  ratio->Divide((TH1D*)hstack_invMass->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  if (lepFlavour == 0) ratio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  else if (lepFlavour == 1)  ratio->GetXaxis()->SetTitle("m_{ee} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (invMassLogy_flag) cInvMass->SaveAs( (plotDirectoryPath + cInvMass->GetName() + suffix + plotFileExtension).c_str() );
  else cInvMass->SaveAs( (plotDirectoryPath + cInvMass->GetName() + suffix + "noLogY" + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadUpar_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadUpar_2 = NULL; 

  TCanvas *cuParMinusZpt = new TCanvas("uParMinusZpt","",700,700);
  TLegend *leguParMinusZpt = new TLegend(0.7,0.7,0.89,0.89);

  subpadUpar_1 = new TPad("padUpar_1","",0.0,0.28,1.0,1.0);
  if (uParLogy_flag) subpadUpar_1->SetLogy();
  //subpadUpar_1->SetBottomMargin(0);
  subpadUpar_2 = new TPad("padUpar_2","",0.0,0.0,1.0,0.27);
  subpadUpar_2->SetGridy();
  //subpadUpar_2->SetTopMargin(0);
  subpadUpar_2->SetBottomMargin(0.3);
  subpadUpar_1->Draw();
  subpadUpar_2->Draw();
  
  subpadUpar_1->cd();
  hstack_uParMinusZpt->Draw("HIST");
  hstack_uParMinusZpt->SetMinimum(0.03);
  hstack_uParMinusZpt->SetMaximum(3000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_uParMinusZpt->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatauParMinusZpt->SetMarkerStyle(8); // large dot
  hDatauParMinusZpt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hDatauParMinusZpt->GetYaxis()->SetRangeUser(0.03,3050.0);
  hstack_uParMinusZpt->GetXaxis()->SetTitle("u_{||} - Z_{pT} [GeV]");
  hstack_uParMinusZpt->GetXaxis()->SetTitleSize(0.04);
  hstack_uParMinusZpt->GetYaxis()->SetTitle("events / 8 GeV");
  hstack_uParMinusZpt->GetYaxis()->SetTitleSize(0.04);
  hstack_uParMinusZpt->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    leguParMinusZpt->AddEntry(hMCuParMinusZpt[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  leguParMinusZpt->AddEntry(hDatauParMinusZpt,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  leguParMinusZpt->Draw(); 
  leguParMinusZpt->SetMargin(0.3); 
  leguParMinusZpt->SetBorderSize(0);
  subpadUpar_2->cd();
  ratio = new TH1D(*hDatauParMinusZpt);
  ratio->Divide((TH1D*)hstack_uParMinusZpt->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("u_{||} - Z_{pT} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (uParLogy_flag) cuParMinusZpt->SaveAs( (plotDirectoryPath + cuParMinusZpt->GetName() + suffix + plotFileExtension).c_str() );
  else cuParMinusZpt->SaveAs( (plotDirectoryPath + cuParMinusZpt->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadUperp_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadUperp_2 = NULL; 

  TCanvas *cuPerp = new TCanvas("uPerp","",700,700);
  TLegend *leguPerp = new TLegend(0.7,0.7,0.89,0.89);

  subpadUperp_1 = new TPad("padUperp_1","",0.0,0.28,1.0,1.0);
  if (uPerpLogy_flag) subpadUperp_1->SetLogy();
  //subpadUperp_1->SetBottomMargin(0);
  subpadUperp_2 = new TPad("padUperp_2","",0.0,0.0,1.0,0.27);
  subpadUperp_2->SetGridy();
  //subpadUperp_2->SetTopMargin(0);
  subpadUperp_2->SetBottomMargin(0.3);
  subpadUperp_1->Draw();
  subpadUperp_2->Draw();

  subpadUperp_1->cd();
  hstack_uPerp->Draw("HIST");
  hstack_uPerp->SetMinimum(0.03);
  hstack_uPerp->SetMaximum(3000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_uPerp->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatauPerp->SetMarkerStyle(8); // large dot
  hDatauPerp->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_uPerp->GetXaxis()->SetTitle("u_{#perp}  [GeV]");
  hstack_uPerp->GetXaxis()->SetTitleSize(0.04);
  hstack_uPerp->GetYaxis()->SetTitle("events / 8 GeV");
  hstack_uPerp->GetYaxis()->SetTitleSize(0.04);
  hstack_uPerp->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    leguPerp->AddEntry(hMCuPerp[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  leguPerp->AddEntry(hDatauPerp,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  leguPerp->Draw(); 
  leguPerp->SetMargin(0.3); 
  leguPerp->SetBorderSize(0);
  subpadUperp_2->cd();
  ratio = new TH1D(*hDatauPerp);
  ratio->Divide((TH1D*)hstack_uPerp->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("u_{#perp}  [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (uPerpLogy_flag) cuPerp->SaveAs( (plotDirectoryPath + cuPerp->GetName() + suffix + plotFileExtension).c_str() );
  else cuPerp->SaveAs( (plotDirectoryPath + cuPerp->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadNvtx_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadNvtx_2 = NULL; 

  TCanvas *cnVertices = new TCanvas("nVertices","",700,700);
  TLegend *legnVertices = new TLegend(0.7,0.7,0.89,0.89);
 
  subpadNvtx_1 = new TPad("padNvtx_1","",0.0,0.28,1.0,1.0);
  if (nvtxLogy_flag) subpadNvtx_1->SetLogy();
  //subpadNvtx_1->SetBottomMargin(0);
  subpadNvtx_2 = new TPad("padNvtx_2","",0.0,0.0,1.0,0.27);
  subpadNvtx_2->SetGridy();
  //subpadNvtx_2->SetTopMargin(0);
  subpadNvtx_2->SetBottomMargin(0.3);
  subpadNvtx_1->Draw();
  subpadNvtx_2->Draw();

  subpadNvtx_1->cd();
  hstack_nVertices->Draw("HIST");
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_nVertices->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatanVertices->SetMarkerStyle(8); // large dot
  hDatanVertices->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_nVertices->GetXaxis()->SetTitle("N_{vtx} [GeV]");
  hstack_nVertices->GetXaxis()->SetTitleSize(0.04);
  hstack_nVertices->GetYaxis()->SetTitle("events");
  hstack_nVertices->GetYaxis()->SetTitleSize(0.04);
  hstack_nVertices->GetYaxis()->CenterTitle();
  //hstack_nVertices->SetMinimum(0.3);
  for (Int_t j = 0; j < nFiles; j++) {
    legnVertices->AddEntry(hMCnVertices[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legnVertices->AddEntry(hDatanVertices,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legnVertices->Draw(); 
  legnVertices->SetMargin(0.3); 
  legnVertices->SetBorderSize(0);
  subpadNvtx_2->cd();
  ratio = new TH1D(*hDatanVertices);
  ratio->Divide((TH1D*)hstack_nVertices->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("N_{vtx} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (nvtxLogy_flag) cnVertices->SaveAs( (plotDirectoryPath + cnVertices->GetName() + suffix + plotFileExtension).c_str() );
  else cnVertices->SaveAs( (plotDirectoryPath + cnVertices->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadJ1pt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ1pt_2 = NULL; 

  TCanvas *cjet1Pt = new TCanvas("jet1Pt","",700,700);
  TLegend *legjet1Pt = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ1pt_1 = new TPad("padJ1pt_1","",0.0,0.28,1.0,1.0);
  if (j1ptLogy_flag) subpadJ1pt_1->SetLogy();
  //subpadJ1pt_1->SetBottomMargin(0);
  subpadJ1pt_2 = new TPad("padJ1pt_2","",0.0,0.0,1.0,0.27);
  subpadJ1pt_2->SetGridy();
  //subpadJ1pt_2->SetTopMargin(0);
  subpadJ1pt_2->SetBottomMargin(0.3);
  subpadJ1pt_1->Draw();
  subpadJ1pt_2->Draw();

  subpadJ1pt_1->cd();
  hstack_jet1Pt->Draw("HIST");
  hstack_jet1Pt->SetMinimum(0.3);
  hstack_jet1Pt->SetMaximum(4000.0);
  // hstack_jet1Pt->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_jet1Pt->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatajet1Pt->SetMarkerStyle(8); // large dot
  hDatajet1Pt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_jet1Pt->GetXaxis()->SetTitle("leading jet pT [GeV]");
  hstack_jet1Pt->GetXaxis()->SetTitleSize(0.04);
  //hstack_jet1Pt->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_jet1Pt->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_jet1Pt->GetYaxis()->SetTitleSize(0.04);
  hstack_jet1Pt->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legjet1Pt->AddEntry(hMCjet1Pt[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legjet1Pt->AddEntry(hDatajet1Pt,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legjet1Pt->Draw(); 
  legjet1Pt->SetMargin(0.3); 
  legjet1Pt->SetBorderSize(0);
  subpadJ1pt_2->cd();
  ratio = new TH1D(*hDatajet1Pt);
  ratio->Divide((TH1D*)hstack_jet1Pt->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("leading jet pT [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (j1ptLogy_flag) cjet1Pt->SaveAs( (plotDirectoryPath + cjet1Pt->GetName() + suffix + plotFileExtension).c_str() );
  else cjet1Pt->SaveAs( (plotDirectoryPath + cjet1Pt->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadZpt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadZpt_2 = NULL; 

  TCanvas *czPtSpectrum = new TCanvas("zPtSpectrum","",700,700);
  TLegend *legzPtSpectrum = new TLegend(0.7,0.7,0.89,0.89);

  subpadZpt_1 = new TPad("padZpt_1","",0.0,0.28,1.0,1.0);
  if (zptLogy_flag) subpadZpt_1->SetLogy();
  //subpadZpt_1->SetBottomMargin(0);
  subpadZpt_2 = new TPad("padZpt_2","",0.0,0.0,1.0,0.27);
  subpadZpt_2->SetGridy();
  //subpadZpt_2->SetTopMargin(0);
  subpadZpt_2->SetBottomMargin(0.3);
  subpadZpt_1->Draw();
  subpadZpt_2->Draw();

  subpadZpt_1->cd();
  hstack_zPtSpectrum->Draw("HIST");
  hstack_zPtSpectrum->SetMinimum(0.3);
  hstack_zPtSpectrum->SetMaximum(4000.0);
  //if (MCuncBand_flag) hstack_zPtSpectrum->Draw("E SAME");
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_zPtSpectrum->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatazPtSpectrum->SetMarkerStyle(8); // large dot
  hDatazPtSpectrum->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_zPtSpectrum->GetXaxis()->SetTitle("Z_{pT} [GeV]");
  hstack_zPtSpectrum->GetXaxis()->SetTitleSize(0.04);
  //hstack_zPtSpectrum->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_zPtSpectrum->GetYaxis()->SetTitle("events / 5 GeV");
  hstack_zPtSpectrum->GetYaxis()->SetTitleSize(0.04);
  hstack_zPtSpectrum->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legzPtSpectrum->AddEntry(hMCzPtSpectrum[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legzPtSpectrum->AddEntry(hDatazPtSpectrum,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legzPtSpectrum->Draw(); 
  legzPtSpectrum->SetMargin(0.3); 
  legzPtSpectrum->SetBorderSize(0);
  //czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );
  subpadZpt_2->cd();
  ratio = new TH1D(*hDatazPtSpectrum);
  ratio->Divide((TH1D*)hstack_zPtSpectrum->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
   if (zptLogy_flag) czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );
   else czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadMet_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadMet_2 = NULL; 

  TCanvas *cmetNoLep = new TCanvas("metNoLep","",700,700);
  TLegend *legmetNoLep = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadMet_1 = new TPad("padMet_1","",0.0,0.28,1.0,1.0);
  if (metLogy_flag) subpadMet_1->SetLogy();
  //subpadMet_1->SetBottomMargin(0);
  subpadMet_2 = new TPad("padMet_2","",0.0,0.0,1.0,0.27);
  subpadMet_2->SetGridy();
  //subpadMet_2->SetTopMargin(0);
  subpadMet_2->SetBottomMargin(0.3);
  subpadMet_1->Draw();
  subpadMet_2->Draw();

  subpadMet_1->cd();
  hstack_metNoLep->Draw("HIST");
  hstack_metNoLep->SetMinimum(0.3);
  hstack_metNoLep->SetMaximum(4000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_metNoLep->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatametNoLep->SetMarkerStyle(8); // large dot
  hDatametNoLep->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_metNoLep->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  hstack_metNoLep->GetXaxis()->SetTitleSize(0.04);
  hstack_metNoLep->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_metNoLep->GetYaxis()->SetTitleSize(0.04);
  hstack_metNoLep->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legmetNoLep->AddEntry(hMCmetNoLep[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legmetNoLep->AddEntry(hDatametNoLep,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legmetNoLep->Draw(); 
  legmetNoLep->SetMargin(0.3); 
  legmetNoLep->SetBorderSize(0);
  subpadMet_2->cd();
  ratio = new TH1D(*hDatametNoLep);
  ratio->Divide((TH1D*)hstack_metNoLep->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (metLogy_flag) cmetNoLep->SaveAs( (plotDirectoryPath + cmetNoLep->GetName() + suffix + plotFileExtension).c_str() );
  else cmetNoLep->SaveAs( (plotDirectoryPath + cmetNoLep->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadNjets_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadNjets_2 = NULL; 

  TCanvas *cnjets = new TCanvas("njets","",700,700);
  TLegend *legnjets = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadNjets_1 = new TPad("padNjets_1","",0.0,0.28,1.0,1.0);
  if (njetsLogy_flag) subpadNjets_1->SetLogy();
  //subpadNjets_1->SetBottomMargin(0);
  subpadNjets_2 = new TPad("padNjets_2","",0.0,0.0,1.0,0.27);
  subpadNjets_2->SetGridy();
  //subpadNjets_2->SetTopMargin(0);
  subpadNjets_2->SetBottomMargin(0.3);
  subpadNjets_1->Draw();
  subpadNjets_2->Draw();

  subpadNjets_1->cd();
  hstack_njets->Draw("HIST");
  hstack_njets->SetMinimum(0.3);
  hstack_njets->SetMaximum(15000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_njets->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatanjets->SetMarkerStyle(8); // large dot
  hDatanjets->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_njets->GetXaxis()->SetTitle("N_{jets}");
  hstack_njets->GetXaxis()->SetTitleSize(0.04);
  hstack_njets->GetXaxis()->SetRange(1,5);
  hstack_njets->GetYaxis()->SetTitle("events ");
  hstack_njets->GetYaxis()->SetTitleSize(0.04);
  hstack_njets->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legnjets->AddEntry(hMCnjets[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legnjets->AddEntry(hDatanjets,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legnjets->Draw(); 
  legnjets->SetMargin(0.3); 
  legnjets->SetBorderSize(0);
  subpadNjets_2->cd();
  ratio = new TH1D(*hDatanjets);
  ratio->Divide((TH1D*)hstack_njets->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("N_{jets}");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetXaxis()->SetRange(1,5);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (njetsLogy_flag) cnjets->SaveAs( (plotDirectoryPath + cnjets->GetName() + suffix + plotFileExtension).c_str() );
  else cnjets->SaveAs( (plotDirectoryPath + cnjets->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;


}

void zlljets_distributions_dataMC_monojetSelection(const char* suffix = "_MonoJetCuts_J1pt30_Zmass80to100_mumu", const Int_t lepFlavour = 0, const Int_t MCuncBand_flag = 0, const Double_t XaxisMaxPt = -1) {

  // in this macro lepFlavour = 0 means muons, while 1 means electrons
  // if XaxisMaxPt = -1 (default), use default value, otherwise use the value chosen by user (to be used when the default range is very different from that of the data

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/distributions/dataMC_comparison/";
  string plotFileExtension = ".pdf";

  vector<string> MC_TexLabel;
  MC_TexLabel.push_back("Z(#tau#tau)+jets MC");
  MC_TexLabel.push_back("W(l#nu)+jets MC");
  MC_TexLabel.push_back("TT+Jets MC");
  if (lepFlavour == 0) MC_TexLabel.push_back("Z(#mu#mu)+jets MC");
  else if (lepFlavour == 1) MC_TexLabel.push_back("Z(ee)+jets MC");

  string data_TexLabel = "data"; 

  string DataFileName; 
  if (lepFlavour == 0) DataFileName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";
  else if (lepFlavour == 1) DataFileName = "zeejets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";

  vector<TH1D*> hMCinvMass;
  vector<TH1D*> hMCmetNoLep;
  vector<TH1D*> hMCjet1Pt;
  vector<TH1D*> hMCjet2Pt;
  vector<TH1D*> hMCjet1eta;
  vector<TH1D*> hMCjet2eta;
  vector<TH1D*> hMCj1j2dphi;
  vector<TH1D*> hMCnVertices;
  vector<TH1D*> hMCzPtSpectrum;
  vector<TH1D*> hMCnjets;

  TH1D* hDatainvMass;
  TH1D* hDatametNoLep;
  TH1D* hDatajet1Pt;
  TH1D* hDatajet2Pt;
  TH1D* hDatajet1eta;
  TH1D* hDatajet2eta;
  TH1D* hDataj1j2dphi;
  TH1D* hDatanVertices;
  TH1D* hDatazPtSpectrum;
  TH1D* hDatanjets;

  Int_t invMassLogy_flag = 0;
  Int_t metLogy_flag = 1;
  Int_t j1ptLogy_flag = 1;
  Int_t j2ptLogy_flag = 1;
  Int_t j1etaLogy_flag = 1;
  Int_t j2etaLogy_flag = 1;
  Int_t j1j2dphiLogy_flag = 1;
  Int_t nvtxLogy_flag = 0;
  Int_t zptLogy_flag = 1;
  Int_t njetsLogy_flag = 1;

  vector<string> MCfileName;
  if (lepFlavour == 0) {
    MCfileName.push_back("ztautaujets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("wjetslnu_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
  } else if (lepFlavour == 1) {
    MCfileName.push_back("ztautaujets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("wjetslnu_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("ttjets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
    MCfileName.push_back("zeejets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
  }

  Int_t nFiles = (Int_t)MCfileName.size();
  Int_t histColor[] = {kYellow,kRed,kBlue,kCyan};

  for(Int_t i = 0; i < nFiles; i++) {
    cout<<"fileName : "<<MCfileName[i]<<endl;

    TFile* f = TFile::Open(MCfileName[i].c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<MCfileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    
    hMCinvMass.push_back((TH1D*)f->Get("HinvMass"));
    hMCmetNoLep.push_back((TH1D*)f->Get("HmetNoLepDistribution"));
    hMCnVertices.push_back((TH1D*)f->Get("HvtxDistribution"));
    hMCzPtSpectrum.push_back((TH1D*)f->Get("HzptDistribution"));
    hMCnjets.push_back((TH1D*)f->Get("HnjetsDistribution"));
    hMCjet1Pt.push_back((TH1D*)f->Get("Hjet1ptDistribution"));
    hMCjet2Pt.push_back((TH1D*)f->Get("Hjet2ptDistribution"));
    hMCjet1eta.push_back((TH1D*)f->Get("Hjet1etaDistribution"));
    hMCjet2eta.push_back((TH1D*)f->Get("Hjet2etaDistribution"));
    hMCj1j2dphi.push_back((TH1D*)f->Get("Hj1j2dphiDistribution"));

  }

  cout<<"fileName : "<< DataFileName<<endl;
  TFile* f = TFile::Open(DataFileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<DataFileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  
  hDatainvMass = (TH1D*)f->Get("HinvMass");
  hDatametNoLep = (TH1D*)f->Get("HmetNoLepDistribution");
  hDatajet1Pt = (TH1D*)f->Get("Hjet1ptDistribution");
  hDatajet2Pt = (TH1D*)f->Get("Hjet2ptDistribution");
  hDatajet1eta = (TH1D*)f->Get("Hjet1etaDistribution");
  hDatajet2eta = (TH1D*)f->Get("Hjet2etaDistribution");
  hDataj1j2dphi = (TH1D*)f->Get("Hj1j2dphiDistribution");
  hDatanVertices = (TH1D*)f->Get("HvtxDistribution");
  hDatazPtSpectrum = (TH1D*)f->Get("HzptDistribution");
  hDatanjets = (TH1D*)f->Get("HnjetsDistribution");

  THStack* hstack_invMass = new THStack("hstack_invMass",""); 
  THStack* hstack_metNoLep = new THStack("hstack_metNoLep","");
  THStack* hstack_jet1Pt = new THStack("hstack_jet1Pt","");
  THStack* hstack_jet2Pt = new THStack("hstack_jet2Pt","");
  THStack* hstack_jet1eta = new THStack("hstack_jet1eta",""); 
  THStack* hstack_jet2eta = new THStack("hstack_jet2eta",""); 
  THStack* hstack_j1j2dphi = new THStack("hstack_j1j2dphi","");
  THStack* hstack_nVertices = new THStack("hstack_nVertices","");
  THStack* hstack_zPtSpectrum = new THStack("hstack_zPtSpectrum",""); 
  THStack* hstack_njets = new THStack("hstack_njets","");

  for (Int_t j = 0; j < nFiles; j++) {

    // set bin error on MC to sqrt(bincontent) for DYJetsToLL and to 50% to others

    if ( j == (nFiles-1)) {

      for (Int_t i = 1; i <= hMCinvMass[j]->GetNbinsX(); i++) {
	hMCinvMass[j]->SetBinError(i,sqrt(hMCinvMass[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {
	hMCmetNoLep[j]->SetBinError(i,sqrt(hMCmetNoLep[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCjet1Pt[j]->GetNbinsX(); i++) {
	hMCjet1Pt[j]->SetBinError(i,sqrt(hMCjet1Pt[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCjet2Pt[j]->GetNbinsX(); i++) {
	hMCjet2Pt[j]->SetBinError(i,sqrt(hMCjet2Pt[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCjet1eta[j]->GetNbinsX(); i++) {
	hMCjet1eta[j]->SetBinError(i,sqrt(hMCjet1eta[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCjet2eta[j]->GetNbinsX(); i++) {
	hMCjet2eta[j]->SetBinError(i,sqrt(hMCjet2eta[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCj1j2dphi[j]->GetNbinsX(); i++) {
	hMCj1j2dphi[j]->SetBinError(i,sqrt(hMCj1j2dphi[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCnVertices[j]->GetNbinsX(); i++) {
	hMCnVertices[j]->SetBinError(i,sqrt(hMCnVertices[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCzPtSpectrum[j]->GetNbinsX(); i++) {
	hMCzPtSpectrum[j]->SetBinError(i,sqrt(hMCzPtSpectrum[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCnjets[j]->GetNbinsX(); i++) {
	hMCnjets[j]->SetBinError(i,sqrt(hMCnjets[j]->GetBinContent(i)));
      }

    } else {

      for (Int_t i = 1; i <= hMCinvMass[j]->GetNbinsX(); i++) {
	hMCinvMass[j]->SetBinError(i,(hMCinvMass[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {
	hMCmetNoLep[j]->SetBinError(i,(hMCmetNoLep[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCjet1Pt[j]->GetNbinsX(); i++) {
	hMCjet1Pt[j]->SetBinError(i,(hMCjet1Pt[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCjet2Pt[j]->GetNbinsX(); i++) {
	hMCjet2Pt[j]->SetBinError(i,(hMCjet2Pt[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCjet1eta[j]->GetNbinsX(); i++) {
	hMCjet1eta[j]->SetBinError(i,(hMCjet1eta[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCjet2eta[j]->GetNbinsX(); i++) {
	hMCjet2eta[j]->SetBinError(i,(hMCjet2eta[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCj1j2dphi[j]->GetNbinsX(); i++) {
	hMCj1j2dphi[j]->SetBinError(i,(hMCj1j2dphi[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCnVertices[j]->GetNbinsX(); i++) {
	hMCnVertices[j]->SetBinError(i,(hMCnVertices[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCzPtSpectrum[j]->GetNbinsX(); i++) {
	hMCzPtSpectrum[j]->SetBinError(i,(hMCzPtSpectrum[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCnjets[j]->GetNbinsX(); i++) {
	hMCnjets[j]->SetBinError(i,(hMCnjets[j]->GetBinContent(i)/2.));
      }

    }

    hMCinvMass[j]->SetFillColor(histColor[j]);
    hstack_invMass->Add(hMCinvMass[j]);
    hMCmetNoLep[j]->SetFillColor(histColor[j]);
    hstack_metNoLep->Add(hMCmetNoLep[j]);
    hMCjet1Pt[j]->SetFillColor(histColor[j]);
    hstack_jet1Pt->Add(hMCjet1Pt[j]);
    hMCjet2Pt[j]->SetFillColor(histColor[j]);
    hstack_jet2Pt->Add(hMCjet2Pt[j]);
    hMCjet1eta[j]->SetFillColor(histColor[j]);
    hstack_jet1eta->Add(hMCjet1eta[j]);
    hMCjet2eta[j]->SetFillColor(histColor[j]);
    hstack_jet2eta->Add(hMCjet2eta[j]);
    hMCj1j2dphi[j]->SetFillColor(histColor[j]);
    hstack_j1j2dphi->Add(hMCj1j2dphi[j]);
    hMCnVertices[j]->SetFillColor(histColor[j]);
    hstack_nVertices->Add(hMCnVertices[j]);    
    hMCzPtSpectrum[j]->SetFillColor(histColor[j]);
    hstack_zPtSpectrum->Add(hMCzPtSpectrum[j]);  
    hMCnjets[j]->SetFillColor(histColor[j]);
    hstack_njets->Add(hMCnjets[j]);

  }

  TH1D * ratio = NULL; // will use it for the ratio plots

  TPad *subpadMass_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadMass_2 = NULL; 

  TCanvas *cInvMass = new TCanvas("invMass","",700,700);
  TLegend *legInvMass = new TLegend(0.7,0.7,0.89,0.89);

  subpadMass_1 = new TPad("padMass_1","",0.0,0.28,1.0,1.0);
  if (invMassLogy_flag) subpadMass_1->SetLogy();
  //subpadMass_1->SetBottomMargin(0);
  subpadMass_2 = new TPad("padMass_2","",0.0,0.0,1.0,0.27);
  subpadMass_2->SetGridy();
  //subpadMass_2->SetTopMargin(0);
  subpadMass_2->SetBottomMargin(0.3);
  subpadMass_1->Draw();
  subpadMass_2->Draw();

  subpadMass_1->cd();
  hstack_invMass->Draw("HIST");
  //hstack_invMass->SetMinimum(0.3);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_invMass->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatainvMass->SetMarkerStyle(8); // large dot
  hDatainvMass->Draw("EX0 SAME"); //X0 doesn't draw x error
  if (lepFlavour == 0) hstack_invMass->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  else if (lepFlavour == 1) hstack_invMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hstack_invMass->GetXaxis()->SetTitleSize(0.04);
  hstack_invMass->GetYaxis()->SetTitle("events / 1 GeV");
  hstack_invMass->GetYaxis()->SetTitleSize(0.04);
  hstack_invMass->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legInvMass->AddEntry(hMCinvMass[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legInvMass->AddEntry(hDatainvMass,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legInvMass->Draw(); 
  legInvMass->SetMargin(0.3); 
  legInvMass->SetBorderSize(0);
  subpadMass_2->cd();
  ratio = new TH1D(*hDatainvMass);
  ratio->Divide((TH1D*)hstack_invMass->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  if (lepFlavour == 0) ratio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  else if (lepFlavour == 1)  ratio->GetXaxis()->SetTitle("m_{ee} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (invMassLogy_flag) cInvMass->SaveAs( (plotDirectoryPath + cInvMass->GetName() + suffix + plotFileExtension).c_str() );
  else cInvMass->SaveAs( (plotDirectoryPath + cInvMass->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadMet_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadMet_2 = NULL; 

  TCanvas *cmetNoLep = new TCanvas("metNoLep","",700,700);
  TLegend *legmetNoLep = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadMet_1 = new TPad("padMet_1","",0.0,0.28,1.0,1.0);
  if (metLogy_flag) subpadMet_1->SetLogy();
  //subpadMet_1->SetBottomMargin(0);
  subpadMet_2 = new TPad("padMet_2","",0.0,0.0,1.0,0.27);
  subpadMet_2->SetGridy();
  //subpadMet_2->SetTopMargin(0);
  subpadMet_2->SetBottomMargin(0.3);
  subpadMet_1->Draw();
  subpadMet_2->Draw();

  subpadMet_1->cd();
  hstack_metNoLep->Draw("HIST");
  hstack_metNoLep->SetMinimum(0.3);
  hstack_metNoLep->SetMaximum(1000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_metNoLep->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatametNoLep->SetMarkerStyle(8); // large dot
  hDatametNoLep->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_metNoLep->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  hstack_metNoLep->GetXaxis()->SetTitleSize(0.04);
  if (XaxisMaxPt > 0) hstack_metNoLep->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  hstack_metNoLep->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_metNoLep->GetYaxis()->SetTitleSize(0.04);
  hstack_metNoLep->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legmetNoLep->AddEntry(hMCmetNoLep[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legmetNoLep->AddEntry(hDatametNoLep,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legmetNoLep->Draw(); 
  legmetNoLep->SetMargin(0.3); 
  legmetNoLep->SetBorderSize(0);
  subpadMet_2->cd();
  ratio = new TH1D(*hDatametNoLep);
  ratio->Divide((TH1D*)hstack_metNoLep->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  if (XaxisMaxPt > 0) ratio->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (metLogy_flag) cmetNoLep->SaveAs( (plotDirectoryPath + cmetNoLep->GetName() + suffix + plotFileExtension).c_str() );
  else cmetNoLep->SaveAs( (plotDirectoryPath + cmetNoLep->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;

  
  
  TPad *subpadJ1pt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ1pt_2 = NULL; 

  TCanvas *cjet1Pt = new TCanvas("jet1Pt","",700,700);
  TLegend *legjet1Pt = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ1pt_1 = new TPad("padJ1pt_1","",0.0,0.28,1.0,1.0);
  if (j1ptLogy_flag) subpadJ1pt_1->SetLogy();
  //subpadJ1pt_1->SetBottomMargin(0);
  subpadJ1pt_2 = new TPad("padJ1pt_2","",0.0,0.0,1.0,0.27);
  subpadJ1pt_2->SetGridy();
  //subpadJ1pt_2->SetTopMargin(0);
  subpadJ1pt_2->SetBottomMargin(0.3);
  subpadJ1pt_1->Draw();
  subpadJ1pt_2->Draw();

  subpadJ1pt_1->cd();
  hstack_jet1Pt->Draw("HIST");
  hstack_jet1Pt->SetMinimum(0.3);
  hstack_jet1Pt->SetMaximum(1000.0);
  // hstack_jet1Pt->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_jet1Pt->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatajet1Pt->SetMarkerStyle(8); // large dot
  hDatajet1Pt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_jet1Pt->GetXaxis()->SetTitle("leading jet pT [GeV]");
  hstack_jet1Pt->GetXaxis()->SetTitleSize(0.04);
  if (XaxisMaxPt > 0) hstack_jet1Pt->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  //hstack_jet1Pt->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_jet1Pt->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_jet1Pt->GetYaxis()->SetTitleSize(0.04);
  hstack_jet1Pt->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legjet1Pt->AddEntry(hMCjet1Pt[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legjet1Pt->AddEntry(hDatajet1Pt,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legjet1Pt->Draw(); 
  legjet1Pt->SetMargin(0.3); 
  legjet1Pt->SetBorderSize(0);
  subpadJ1pt_2->cd();
  ratio = new TH1D(*hDatajet1Pt);
  ratio->Divide((TH1D*)hstack_jet1Pt->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("leading jet pT [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  if (XaxisMaxPt > 0) ratio->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (j1ptLogy_flag) cjet1Pt->SaveAs( (plotDirectoryPath + cjet1Pt->GetName() + suffix + plotFileExtension).c_str() );
  else cjet1Pt->SaveAs( (plotDirectoryPath + cjet1Pt->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;

 

  TPad *subpadJ2pt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ2pt_2 = NULL; 

  TCanvas *cjet2Pt = new TCanvas("jet2Pt","",700,700);
  TLegend *legjet2Pt = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ2pt_1 = new TPad("padJ2pt_1","",0.0,0.28,1.0,1.0);
  if (j2ptLogy_flag) subpadJ2pt_1->SetLogy();
  //subpadJ2pt_1->SetBottomMargin(0);
  subpadJ2pt_2 = new TPad("padJ2pt_2","",0.0,0.0,1.0,0.27);
  subpadJ2pt_2->SetGridy();
  //subpadJ2pt_2->SetTopMargin(0);
  subpadJ2pt_2->SetBottomMargin(0.3);
  subpadJ2pt_1->Draw();
  subpadJ2pt_2->Draw();

  subpadJ2pt_1->cd();
  hstack_jet2Pt->Draw("HIST");
  hstack_jet2Pt->SetMinimum(0.3);
  hstack_jet2Pt->SetMaximum(1000.0);
  // hstack_jet2Pt->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_jet2Pt->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatajet2Pt->SetMarkerStyle(8); // large dot
  hDatajet2Pt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_jet2Pt->GetXaxis()->SetTitle("second jet pT [GeV]");
  hstack_jet2Pt->GetXaxis()->SetTitleSize(0.04);
  if (XaxisMaxPt > 0) hstack_jet2Pt->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  //hstack_jet2Pt->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_jet2Pt->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_jet2Pt->GetYaxis()->SetTitleSize(0.04);
  hstack_jet2Pt->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legjet2Pt->AddEntry(hMCjet2Pt[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legjet2Pt->AddEntry(hDatajet2Pt,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legjet2Pt->Draw(); 
  legjet2Pt->SetMargin(0.3); 
  legjet2Pt->SetBorderSize(0);
  subpadJ2pt_2->cd();
  ratio = new TH1D(*hDatajet2Pt);
  ratio->Divide((TH1D*)hstack_jet2Pt->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("second jet pT [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  if (XaxisMaxPt > 0) ratio->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (j2ptLogy_flag) cjet2Pt->SaveAs( (plotDirectoryPath + cjet2Pt->GetName() + suffix + plotFileExtension).c_str() );
  else cjet2Pt->SaveAs( (plotDirectoryPath + cjet2Pt->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadJ1eta_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ1eta_2 = NULL; 

  TCanvas *cjet1eta = new TCanvas("jet1eta","",700,700);
  TLegend *legjet1eta = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ1eta_1 = new TPad("padJ1eta_1","",0.0,0.28,1.0,1.0);
  if (j1etaLogy_flag) subpadJ1eta_1->SetLogy();
  //subpadJ1eta_1->SetBottomMargin(0);
  subpadJ1eta_2 = new TPad("padJ1eta_2","",0.0,0.0,1.0,0.27);
  subpadJ1eta_2->SetGridy();
  //subpadJ1eta_2->SetTopMargin(0);
  subpadJ1eta_2->SetBottomMargin(0.3);
  subpadJ1eta_1->Draw();
  subpadJ1eta_2->Draw();

  subpadJ1eta_1->cd();
  hstack_jet1eta->Draw("HIST");
  hstack_jet1eta->SetMinimum(0.3);
  hstack_jet1eta->SetMaximum(1000.0);
  // hstack_jet1eta->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_jet1eta->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatajet1eta->SetMarkerStyle(8); // large dot
  hDatajet1eta->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_jet1eta->GetXaxis()->SetTitle("leading jet #eta");
  hstack_jet1eta->GetXaxis()->SetTitleSize(0.04);
  //hstack_jet1eta->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_jet1eta->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_jet1eta->GetYaxis()->SetTitleSize(0.04);
  hstack_jet1eta->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legjet1eta->AddEntry(hMCjet1eta[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legjet1eta->AddEntry(hDatajet1eta,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legjet1eta->Draw(); 
  legjet1eta->SetMargin(0.3); 
  legjet1eta->SetBorderSize(0);
  subpadJ1eta_2->cd();
  ratio = new TH1D(*hDatajet1eta);
  ratio->Divide((TH1D*)hstack_jet1eta->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("leading jet #eta");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (j1etaLogy_flag) cjet1eta->SaveAs( (plotDirectoryPath + cjet1eta->GetName() + suffix + plotFileExtension).c_str() );
  else cjet1eta->SaveAs( (plotDirectoryPath + cjet1eta->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadJ2eta_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ2eta_2 = NULL; 

  TCanvas *cjet2eta = new TCanvas("jet2eta","",700,700);
  TLegend *legjet2eta = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ2eta_1 = new TPad("padJ2eta_1","",0.0,0.28,1.0,1.0);
  if (j2etaLogy_flag) subpadJ2eta_1->SetLogy();
  //subpadJ2eta_1->SetBottomMargin(0);
  subpadJ2eta_2 = new TPad("padJ2eta_2","",0.0,0.0,1.0,0.27);
  subpadJ2eta_2->SetGridy();
  //subpadJ2eta_2->SetTopMargin(0);
  subpadJ2eta_2->SetBottomMargin(0.3);
  subpadJ2eta_1->Draw();
  subpadJ2eta_2->Draw();

  subpadJ2eta_1->cd();
  hstack_jet2eta->Draw("HIST");
  hstack_jet2eta->SetMinimum(0.3);
  hstack_jet2eta->SetMaximum(1000.0);
  // hstack_jet2eta->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_jet2eta->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatajet2eta->SetMarkerStyle(8); // large dot
  hDatajet2eta->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_jet2eta->GetXaxis()->SetTitle("second jet #eta");
  hstack_jet2eta->GetXaxis()->SetTitleSize(0.04);
  //hstack_jet2eta->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_jet2eta->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_jet2eta->GetYaxis()->SetTitleSize(0.04);
  hstack_jet2eta->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legjet2eta->AddEntry(hMCjet2eta[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legjet2eta->AddEntry(hDatajet2eta,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legjet2eta->Draw(); 
  legjet2eta->SetMargin(0.3); 
  legjet2eta->SetBorderSize(0);
  subpadJ2eta_2->cd();
  ratio = new TH1D(*hDatajet2eta);
  ratio->Divide((TH1D*)hstack_jet2eta->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("second jet #eta");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (j2etaLogy_flag) cjet2eta->SaveAs( (plotDirectoryPath + cjet2eta->GetName() + suffix + plotFileExtension).c_str() );
  else cjet2eta->SaveAs( (plotDirectoryPath + cjet2eta->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;

  

  TPad *subpadJ1j2dphi_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ1j2dphi_2 = NULL; 

  TCanvas *cj1j2dphi = new TCanvas("j1j2dphi","",700,700);
  TLegend *legj1j2dphi = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ1j2dphi_1 = new TPad("padJ1j2dphi_1","",0.0,0.28,1.0,1.0);
  if (j1j2dphiLogy_flag) subpadJ1j2dphi_1->SetLogy();
  //subpadJ1j2dphi_1->SetBottomMargin(0);
  subpadJ1j2dphi_2 = new TPad("padJ1j2dphi_2","",0.0,0.0,1.0,0.27);
  subpadJ1j2dphi_2->SetGridy();
  //subpadJ1j2dphi_2->SetTopMargin(0);
  subpadJ1j2dphi_2->SetBottomMargin(0.3);
  subpadJ1j2dphi_1->Draw();
  subpadJ1j2dphi_2->Draw();

  subpadJ1j2dphi_1->cd();
  hstack_j1j2dphi->Draw("HIST");
  hstack_j1j2dphi->SetMinimum(0.03);
  hstack_j1j2dphi->SetMaximum(1000.0);
  // hstack_j1j2dphi->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_j1j2dphi->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDataj1j2dphi->SetMarkerStyle(8); // large dot
  hDataj1j2dphi->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_j1j2dphi->GetXaxis()->SetTitle("#Delta#phi (jet1,jet2)");
  hstack_j1j2dphi->GetXaxis()->SetTitleSize(0.04);
  //hstack_j1j2dphi->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_j1j2dphi->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_j1j2dphi->GetYaxis()->SetTitleSize(0.04);
  hstack_j1j2dphi->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legj1j2dphi->AddEntry(hMCj1j2dphi[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legj1j2dphi->AddEntry(hDataj1j2dphi,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legj1j2dphi->Draw(); 
  legj1j2dphi->SetMargin(0.3); 
  legj1j2dphi->SetBorderSize(0);
  subpadJ1j2dphi_2->cd();
  ratio = new TH1D(*hDataj1j2dphi);
  ratio->Divide((TH1D*)hstack_j1j2dphi->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("#Delta#phi (jet1,jet2)");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (j1j2dphiLogy_flag) cj1j2dphi->SaveAs( (plotDirectoryPath + cj1j2dphi->GetName() + suffix + plotFileExtension).c_str() );
  else cj1j2dphi->SaveAs( (plotDirectoryPath + cj1j2dphi->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadNvtx_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadNvtx_2 = NULL; 

  TCanvas *cnVertices = new TCanvas("nVertices","",700,700);
  TLegend *legnVertices = new TLegend(0.7,0.7,0.89,0.89);
 
  subpadNvtx_1 = new TPad("padNvtx_1","",0.0,0.28,1.0,1.0);
  if (nvtxLogy_flag) subpadNvtx_1->SetLogy();
  //subpadNvtx_1->SetBottomMargin(0);
  subpadNvtx_2 = new TPad("padNvtx_2","",0.0,0.0,1.0,0.27);
  subpadNvtx_2->SetGridy();
  //subpadNvtx_2->SetTopMargin(0);
  subpadNvtx_2->SetBottomMargin(0.3);
  subpadNvtx_1->Draw();
  subpadNvtx_2->Draw();

  subpadNvtx_1->cd();
  hstack_nVertices->Draw("HIST");
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_nVertices->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatanVertices->SetMarkerStyle(8); // large dot
  hDatanVertices->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_nVertices->GetXaxis()->SetTitle("N_{vtx} [GeV]");
  hstack_nVertices->GetXaxis()->SetTitleSize(0.04);
  hstack_nVertices->GetYaxis()->SetTitle("events");
  hstack_nVertices->GetYaxis()->SetTitleSize(0.04);
  hstack_nVertices->GetYaxis()->CenterTitle();
  //hstack_nVertices->SetMinimum(0.3);
  for (Int_t j = 0; j < nFiles; j++) {
    legnVertices->AddEntry(hMCnVertices[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legnVertices->AddEntry(hDatanVertices,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legnVertices->Draw(); 
  legnVertices->SetMargin(0.3); 
  legnVertices->SetBorderSize(0);
  subpadNvtx_2->cd();
  ratio = new TH1D(*hDatanVertices);
  ratio->Divide((TH1D*)hstack_nVertices->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("N_{vtx} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (nvtxLogy_flag) cnVertices->SaveAs( (plotDirectoryPath + cnVertices->GetName() + suffix + plotFileExtension).c_str() );
  else cnVertices->SaveAs( (plotDirectoryPath + cnVertices->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadZpt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadZpt_2 = NULL; 

  TCanvas *czPtSpectrum = new TCanvas("zPtSpectrum","",700,700);
  TLegend *legzPtSpectrum = new TLegend(0.7,0.7,0.89,0.89);

  subpadZpt_1 = new TPad("padZpt_1","",0.0,0.28,1.0,1.0);
  if (zptLogy_flag) subpadZpt_1->SetLogy();
  //subpadZpt_1->SetBottomMargin(0);
  subpadZpt_2 = new TPad("padZpt_2","",0.0,0.0,1.0,0.27);
  subpadZpt_2->SetGridy();
  //subpadZpt_2->SetTopMargin(0);
  subpadZpt_2->SetBottomMargin(0.3);
  subpadZpt_1->Draw();
  subpadZpt_2->Draw();

  subpadZpt_1->cd();
  hstack_zPtSpectrum->Draw("HIST");
  hstack_zPtSpectrum->SetMinimum(0.3);
  hstack_zPtSpectrum->SetMaximum(1000.0);
  //if (MCuncBand_flag) hstack_zPtSpectrum->Draw("E SAME");
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_zPtSpectrum->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatazPtSpectrum->SetMarkerStyle(8); // large dot
  hDatazPtSpectrum->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_zPtSpectrum->GetXaxis()->SetTitle("Z_{pT} [GeV]");
  hstack_zPtSpectrum->GetXaxis()->SetTitleSize(0.04);
  if (XaxisMaxPt > 0) hstack_zPtSpectrum->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  //hstack_zPtSpectrum->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_zPtSpectrum->GetYaxis()->SetTitle("events / 5 GeV");
  hstack_zPtSpectrum->GetYaxis()->SetTitleSize(0.04);
  hstack_zPtSpectrum->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legzPtSpectrum->AddEntry(hMCzPtSpectrum[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legzPtSpectrum->AddEntry(hDatazPtSpectrum,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legzPtSpectrum->Draw(); 
  legzPtSpectrum->SetMargin(0.3); 
  legzPtSpectrum->SetBorderSize(0);
  //czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );
  subpadZpt_2->cd();
  ratio = new TH1D(*hDatazPtSpectrum);
  ratio->Divide((TH1D*)hstack_zPtSpectrum->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  if (XaxisMaxPt > 0) ratio->GetXaxis()->SetRangeUser(0,XaxisMaxPt);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (zptLogy_flag) czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );
  else czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadNjets_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadNjets_2 = NULL; 

  TCanvas *cnjets = new TCanvas("njets","",700,700);
  TLegend *legnjets = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadNjets_1 = new TPad("padNjets_1","",0.0,0.28,1.0,1.0);
  if (njetsLogy_flag) subpadNjets_1->SetLogy();
  //subpadNjets_1->SetBottomMargin(0);
  subpadNjets_2 = new TPad("padNjets_2","",0.0,0.0,1.0,0.27);
  subpadNjets_2->SetGridy();
  //subpadNjets_2->SetTopMargin(0);
  subpadNjets_2->SetBottomMargin(0.3);
  subpadNjets_1->Draw();
  subpadNjets_2->Draw();

  subpadNjets_1->cd();
  hstack_njets->Draw("HIST");
  hstack_njets->SetMinimum(0.3);
  hstack_njets->SetMaximum(5000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_njets->GetStack()->Last())->DrawCopy("E2 SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    //stackCopy->Draw("E2 SAME");
  }
  hDatanjets->SetMarkerStyle(8); // large dot
  hDatanjets->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_njets->GetXaxis()->SetTitle("N_{jets}");
  hstack_njets->GetXaxis()->SetTitleSize(0.04);
  hstack_njets->GetXaxis()->SetRange(1,5);
  hstack_njets->GetYaxis()->SetTitle("events ");
  hstack_njets->GetYaxis()->SetTitleSize(0.04);
  hstack_njets->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legnjets->AddEntry(hMCnjets[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legnjets->AddEntry(hDatanjets,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legnjets->Draw(); 
  legnjets->SetMargin(0.3); 
  legnjets->SetBorderSize(0);
  subpadNjets_2->cd();
  ratio = new TH1D(*hDatanjets);
  ratio->Divide((TH1D*)hstack_njets->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("N_{jets}");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetXaxis()->SetRange(1,5);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  if (njetsLogy_flag) cnjets->SaveAs( (plotDirectoryPath + cnjets->GetName() + suffix + plotFileExtension).c_str() );
  else cnjets->SaveAs( (plotDirectoryPath + cnjets->GetName() + suffix + "_noLogY" + plotFileExtension).c_str() );
  delete ratio;

}



void zlljets_distributions_dataMC_new(const char* suffix = "_GLTcutZmass80to100_mumu") {

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/distributions/dataMC_comparison/";
  string plotFileExtension = ".pdf";

  vector<string> MC_TexLabel;
  MC_TexLabel.push_back("TT+Jets MC");
  MC_TexLabel.push_back("Z(#mu#mu)+jets MC");
  string data_TexLabel = "data"; 

  vector<TH1D*> hMC;
  vector<TH1D*> hData;

  vector<string> hName;
  hName.push_back("HmetNoLepDistribution");
  hName.push_back("HzptDistribution");
  hName.push_back("Hjet1ptDistribution");
  hName.push_back("HinvMass");
  hName.push_back("HvtxDistribution");
  hName.push_back("H_uPerp_Distribution");
  hName.push_back("H_uParMinusZpT_Distribution");

  vector<string> cNameSuffix;
  cNameSuffix.push_back("metNoLep");
  cNameSuffix.push_back("zpt");
  cNameSuffix.push_back("jet1pt");
  cNameSuffix.push_back("invMass");
  cNameSuffix.push_back("Nvtx");
  cNameSuffix.push_back("uPerp");
  cNameSuffix.push_back("uParMinusZpT");

  Int_t logyScale_flag[] = {1,1,1,0,0,1,1};  // decide whether corresponding distribution is to be drwan with log scale on y axis

  Int_t nDistr = hName.size();

  vector<string> xAxisLabel;
  xAxisLabel.push_back("#slash{E}_{T} [GeV]");
  xAxisLabel.push_back("ZpT [GeV]");
  xAxisLabel.push_back("leading jet pT [GeV]");
  xAxisLabel.push_back("m_{#mu#mu} [GeV]");
  xAxisLabel.push_back("N_{vtx}");
  xAxisLabel.push_back("u_{#perp} [GeV]");
  xAxisLabel.push_back("u_{||} - Z_{pT} [GeV]");

  vector<string> yAxisLabel;
  yAxisLabel.push_back("events / 10 GeV");
  yAxisLabel.push_back("events / 5 GeV");
  yAxisLabel.push_back("events / 10 GeV");
  yAxisLabel.push_back("events / 1 GeV");
  yAxisLabel.push_back("events");
  yAxisLabel.push_back("events / 8 GeV");
  yAxisLabel.push_back("events / 8 GeV");

  TCanvas *c[nDistr];
  TLegend *leg[nDistr];
  THStack *hstack[nDistr];

  for (Int_t distr = 0; distr < nDistr; distr++) { 

    c[distr] = new TCanvas(Form("c_%s",cNameSuffix[distr].c_str()),"");
    if (logyScale_flag[distr]) c[distr]->SetLogy();
    hstack[distr] = new THStack(Form("hstack_%s",cNameSuffix[distr].c_str()),"");
    leg[distr] = new TLegend(0.7,0.7,0.89,0.89);
			
  }

  string DataFileName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";
  vector<string> MCfileName;
  MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
  MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  

  Int_t nFiles = (Int_t)MCfileName.size();
  Int_t histColor[] = {kRed,kCyan};

  TFile* f[nFiles];
  
  for (Int_t i = 0; i < nFiles; i++) {

    f[i] = TFile::Open(MCfileName[i].c_str(),"READ");
      if (!f[i] || !f[i]->IsOpen()) {
	cout<<"*******************************"<<endl;
	cout<<"Error opening file \""<<MCfileName[i]<<"\".\nApplication will be terminated."<<endl;
	cout<<"*******************************"<<endl;
	exit(EXIT_FAILURE);
      }

  }

  TFile* fData = TFile::Open(DataFileName.c_str(),"READ");
  if (!fData || !fData->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<DataFileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  for (Int_t distr = 0; distr < nDistr; distr++) { 

    cout << "Distribution --> " << cNameSuffix[distr] << endl; 

    for(Int_t i = 0; i < nFiles; i++) {
      
      //cout<<"fileName : "<<MCfileName[i]<<endl 
      
      hMC.push_back((TH1D*)f[i]->Get(hName[distr].c_str()));
      if (!hMC[distr + i]) {
	cout << "Error: coud not get histogram "<< hName[distr].c_str() << " from file " << MCfileName[i] << endl;
	exit(EXIT_FAILURE);
      }    
   
    }
   
    hData.push_back((TH1D*)fData->Get(hName[distr].c_str()));
    
    for (Int_t j = 0; j < nFiles; j++) {
      
      // set bin error on MC to sqrt(bincontent) for DYJetsToLL and to 50% to others
      
      if ( j == (nFiles-1)) {

	for (Int_t i = 1; i <= hMC[distr + j]->GetNbinsX(); i++) {
	  hMC[distr + j]->SetBinError(i,sqrt(hMC[distr + j]->GetBinContent(i)));
	}
	
      } else {

	for (Int_t i = 1; i <= hMC[distr + j]->GetNbinsX(); i++) {
	  hMC[distr + j]->SetBinError(i,(hMC[distr + j]->GetBinContent(i)/2.));
	}	

      }

      hMC[distr + j]->SetFillColor(histColor[j]);
      hstack[distr]->Add(hMC[distr + j]);

    }

    hData[distr]->SetMarkerStyle(8);

    // filling canvases and drawing them

    c[distr]->cd();
    hstack[distr]->Draw("HE");
    hstack[distr]->Draw("E SAME");                        // must draw twice for having error bars in the filled area
    hData[distr]->Draw("EX0 SAME");                    //X0 doesn't draw x error
    hstack[distr]->GetXaxis()->SetTitle(xAxisLabel[distr].c_str());
    hstack[distr]->GetXaxis()->SetTitleSize(0.04);
    hstack[distr]->GetYaxis()->SetTitle(yAxisLabel[distr].c_str());
    hstack[distr]->GetYaxis()->SetTitleSize(0.04);
    hstack[distr]->GetYaxis()->CenterTitle();
    for (Int_t j = 0; j < nFiles; j++) {
      leg[distr]->AddEntry(hMC[distr +j],Form("%s",MC_TexLabel[j].c_str()),"lf");
    }
    leg[distr]->AddEntry(hData[distr],Form("%s",data_TexLabel.c_str()),"p");
    gStyle->SetStatStyle(0);
    leg[distr]->Draw(); 
    leg[distr]->SetMargin(0.3); 
    leg[distr]->SetBorderSize(0);
    c[distr]->SaveAs( (plotDirectoryPath + cNameSuffix[distr] + suffix + plotFileExtension).c_str() );

  }


}
			     

//THE FOLLOWING MUST BE COMPLETED
/*
void zlljets_distributions_onlyMC(const char* suffix = "_GLTcutZmass80to100_onlyMC_mumu", const Int_t lepFlavour = 0, const Int_t MCuncBand_flag = 0, const Int_t spring15_flag = 0) {

  // in this macro lepFlavour = 0 means muons, while 1 means electrons

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/distributions/onlyMC/";
  string plotFileExtension = ".pdf";

  vector<string> MC_TexLabel;
  MC_TexLabel.push_back("Z(#tau#tau)+jets MC");
  MC_TexLabel.push_back("W(l#nu)+jets MC");
  MC_TexLabel.push_back("TT+Jets MC");
  if (lepFlavour == 0) MC_TexLabel.push_back("Z(#mu#mu)+jets MC");
  else if (lepFlavour == 1) MC_TexLabel.push_back("Z(ee)+jets MC");

  string data_TexLabel = "data"; 

  vector<TH1D*> hMCinvMass;
  vector<TH1D*> hMCmetNoLep;
  vector<TH1D*> hMCjet1Pt;
  vector<TH1D*> hMCuParMinusZpt;
  vector<TH1D*> hMCuPerp;
  vector<TH1D*> hMCnVertices;
  vector<TH1D*> hMCzPtSpectrum;
  vector<TH1D*> hMCnjets;

  TH1D* hDatainvMass;
  TH1D* hDatametNoLep;
  TH1D* hDatajet1Pt;
  TH1D* hDatauParMinusZpt;
  TH1D* hDatauPerp;
  TH1D* hDatanVertices;
  TH1D* hDatazPtSpectrum;
  TH1D* hDatanjets;

  vector<string> MCfileName;
  if (spring15_flag) {
    if (lepFlavour == 0) {
      MCfileName.push_back("ztautaujets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("wjetslnu_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
    } else if (lepFlavour == 1) {
      MCfileName.push_back("ztautaujets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("wjetslnu_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("ttjets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("zeejets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
    }
  } else {   //use root file for phys14 
    if (lepFlavour == 0) {
      MCfileName.push_back("ztautaujets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("wjetslnu_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
    } else if (lepFlavour == 1) {
      MCfileName.push_back("ztautaujets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("wjetslnu_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("ttjets_ee_resoResp_spring15_50ns_GLTcutZmass80to100.root");
      MCfileName.push_back("zeejets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  
    }
  }

  Int_t nFiles = (Int_t)MCfileName.size();
  Int_t histColor[] = {kYellow,kRed,kBlue,kCyan};

  for(Int_t i = 0; i < nFiles; i++) {
    cout<<"fileName : "<<MCfileName[i]<<endl;

    TFile* f = TFile::Open(MCfileName[i].c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<MCfileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    
    hMCinvMass.push_back((TH1D*)f->Get("HinvMass_NoJetCuts"));
    hMCmetNoLep.push_back((TH1D*)f->Get("HmetNoLepDistribution_NoJetCuts"));
    hMCjet1Pt.push_back((TH1D*)f->Get("Hjet1ptDistribution_NoJetCuts"));
    hMCuParMinusZpt.push_back((TH1D*)f->Get("H_uParMinusZpT_Distribution"));
    hMCuPerp.push_back((TH1D*)f->Get("H_uPerp_Distribution"));
    hMCnVertices.push_back((TH1D*)f->Get("HvtxDistribution_NoJetCuts"));
    hMCzPtSpectrum.push_back((TH1D*)f->Get("HzptDistribution_NoJetCuts"));
    hMCnjets.push_back((TH1D*)f->Get("HnjetsDistribution_NoJetCuts"));
   
  }

  cout<<"fileName : "<< DataFileName<<endl;
  TFile* f = TFile::Open(DataFileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<DataFileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  
  hDatainvMass = (TH1D*)f->Get("HinvMass_NoJetCuts");
  hDatametNoLep = (TH1D*)f->Get("HmetNoLepDistribution_NoJetCuts");
  hDatajet1Pt = (TH1D*)f->Get("Hjet1ptDistribution_NoJetCuts");
  hDatauParMinusZpt = (TH1D*)f->Get("H_uParMinusZpT_Distribution");
  hDatauPerp = (TH1D*)f->Get("H_uPerp_Distribution");
  hDatanVertices = (TH1D*)f->Get("HvtxDistribution_NoJetCuts");
  hDatazPtSpectrum = (TH1D*)f->Get("HzptDistribution_NoJetCuts");
  hDatanjets = (TH1D*)f->Get("HnjetsDistribution_NoJetCuts");

  THStack* hstack_invMass = new THStack("hstack_invMass",""); 
  THStack* hstack_uPerp = new THStack("hstack_uPerp",""); 
  THStack* hstack_uParMinusZpt = new THStack("hstack_uParMinusZpt",""); 
  THStack* hstack_nVertices = new THStack("hstack_nVertices","");
  THStack* hstack_jet1Pt = new THStack("hstack_jet1Pt","");
  THStack* hstack_zPtSpectrum = new THStack("hstack_zPtSpectrum","");
  THStack* hstack_metNoLep = new THStack("hstack_metNoLep","");
  THStack* hstack_njets = new THStack("hstack_njets","");

  for (Int_t j = 0; j < nFiles; j++) {

    // set bin error on MC to sqrt(bincontent) for DYJetsToLL and to 50% to others

    if ( j == (nFiles-1)) {

      for (Int_t i = 1; i <= hMCinvMass[j]->GetNbinsX(); i++) {
	hMCinvMass[j]->SetBinError(i,sqrt(hMCinvMass[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCuPerp[j]->GetNbinsX(); i++) {
	hMCuPerp[j]->SetBinError(i,sqrt(hMCuPerp[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCuParMinusZpt[j]->GetNbinsX(); i++) {
	hMCuParMinusZpt[j]->SetBinError(i,sqrt(hMCuParMinusZpt[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCnVertices[j]->GetNbinsX(); i++) {
	hMCnVertices[j]->SetBinError(i,sqrt(hMCnVertices[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCjet1Pt[j]->GetNbinsX(); i++) {
	hMCjet1Pt[j]->SetBinError(i,sqrt(hMCjet1Pt[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCzPtSpectrum[j]->GetNbinsX(); i++) {
	hMCzPtSpectrum[j]->SetBinError(i,sqrt(hMCzPtSpectrum[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {
	hMCmetNoLep[j]->SetBinError(i,sqrt(hMCmetNoLep[j]->GetBinContent(i)));
      }
      for (Int_t i = 1; i <= hMCnjets[j]->GetNbinsX(); i++) {
	hMCnjets[j]->SetBinError(i,sqrt(hMCnjets[j]->GetBinContent(i)));
      }

    } else {

      for (Int_t i = 1; i <= hMCinvMass[j]->GetNbinsX(); i++) {
	hMCinvMass[j]->SetBinError(i,(hMCinvMass[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCuPerp[j]->GetNbinsX(); i++) {
	hMCuPerp[j]->SetBinError(i,(hMCuPerp[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCuParMinusZpt[j]->GetNbinsX(); i++) {
	hMCuParMinusZpt[j]->SetBinError(i,(hMCuParMinusZpt[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCnVertices[j]->GetNbinsX(); i++) {
	hMCnVertices[j]->SetBinError(i,(hMCnVertices[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCjet1Pt[j]->GetNbinsX(); i++) {
	hMCjet1Pt[j]->SetBinError(i,(hMCjet1Pt[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCzPtSpectrum[j]->GetNbinsX(); i++) {
	hMCzPtSpectrum[j]->SetBinError(i,(hMCzPtSpectrum[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {
	hMCmetNoLep[j]->SetBinError(i,(hMCmetNoLep[j]->GetBinContent(i)/2.));
      }
      for (Int_t i = 1; i <= hMCnjets[j]->GetNbinsX(); i++) {
	hMCnjets[j]->SetBinError(i,(hMCnjets[j]->GetBinContent(i)/2.));
      }

    }

    hMCinvMass[j]->SetFillColor(histColor[j]);
    hstack_invMass->Add(hMCinvMass[j]);
    hMCuPerp[j]->SetFillColor(histColor[j]);
    hstack_uPerp->Add(hMCuPerp[j]);
    hMCuParMinusZpt[j]->SetFillColor(histColor[j]);
    hstack_uParMinusZpt->Add(hMCuParMinusZpt[j]);
    hMCnVertices[j]->SetFillColor(histColor[j]);
    hstack_nVertices->Add(hMCnVertices[j]);
    hMCjet1Pt[j]->SetFillColor(histColor[j]);
    hstack_jet1Pt->Add(hMCjet1Pt[j]);
    hMCzPtSpectrum[j]->SetFillColor(histColor[j]);
    hstack_zPtSpectrum->Add(hMCzPtSpectrum[j]);
    hMCmetNoLep[j]->SetFillColor(histColor[j]);
    hstack_metNoLep->Add(hMCmetNoLep[j]);
    hMCnjets[j]->SetFillColor(histColor[j]);
    hstack_njets->Add(hMCnjets[j]);

  }

  TH1D * ratio = NULL; // will use it for the ratio plots

  TPad *subpadMass_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadMass_2 = NULL; 

  TCanvas *cInvMass = new TCanvas("invMass","",700,700);
  TLegend *legInvMass = new TLegend(0.7,0.7,0.89,0.89);

  subpadMass_1 = new TPad("padMass_1","",0.0,0.28,1.0,1.0);
  //subpadMass_1->SetLogy();
  //subpadMass_1->SetBottomMargin(0);
  subpadMass_2 = new TPad("padMass_2","",0.0,0.0,1.0,0.27);
  subpadMass_2->SetGridy();
  //subpadMass_2->SetTopMargin(0);
  subpadMass_2->SetBottomMargin(0.3);
  subpadMass_1->Draw();
  subpadMass_2->Draw();

  subpadMass_1->cd();
  hstack_invMass->Draw("HIST");
  //hstack_invMass->SetMinimum(0.3);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_invMass->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatainvMass->SetMarkerStyle(8); // large dot
  hDatainvMass->Draw("EX0 SAME"); //X0 doesn't draw x error
  if (lepFlavour == 0) hstack_invMass->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  else if (lepFlavour == 1) hstack_invMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hstack_invMass->GetXaxis()->SetTitleSize(0.04);
  hstack_invMass->GetYaxis()->SetTitle("events / 1 GeV");
  hstack_invMass->GetYaxis()->SetTitleSize(0.04);
  hstack_invMass->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legInvMass->AddEntry(hMCinvMass[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legInvMass->AddEntry(hDatainvMass,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legInvMass->Draw(); 
  legInvMass->SetMargin(0.3); 
  legInvMass->SetBorderSize(0);
  subpadMass_2->cd();
  ratio = new TH1D(*hDatainvMass);
  ratio->Divide((TH1D*)hstack_invMass->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  if (lepFlavour == 0) ratio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  else if (lepFlavour == 1)  ratio->GetXaxis()->SetTitle("m_{ee} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cInvMass->SaveAs( (plotDirectoryPath + cInvMass->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadUpar_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadUpar_2 = NULL; 

  TCanvas *cuParMinusZpt = new TCanvas("uParMinusZpt","",700,700);
  TLegend *leguParMinusZpt = new TLegend(0.7,0.7,0.89,0.89);

  subpadUpar_1 = new TPad("padUpar_1","",0.0,0.28,1.0,1.0);
  subpadUpar_1->SetLogy();
  //subpadUpar_1->SetBottomMargin(0);
  subpadUpar_2 = new TPad("padUpar_2","",0.0,0.0,1.0,0.27);
  subpadUpar_2->SetGridy();
  //subpadUpar_2->SetTopMargin(0);
  subpadUpar_2->SetBottomMargin(0.3);
  subpadUpar_1->Draw();
  subpadUpar_2->Draw();
  
  subpadUpar_1->cd();
  hstack_uParMinusZpt->Draw("HIST");
  hstack_uParMinusZpt->SetMinimum(0.03);
  hstack_uParMinusZpt->SetMaximum(3000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_uParMinusZpt->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatauParMinusZpt->SetMarkerStyle(8); // large dot
  hDatauParMinusZpt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hDatauParMinusZpt->GetYaxis()->SetRangeUser(0.03,3050.0);
  hstack_uParMinusZpt->GetXaxis()->SetTitle("u_{||} - Z_{pT} [GeV]");
  hstack_uParMinusZpt->GetXaxis()->SetTitleSize(0.04);
  hstack_uParMinusZpt->GetYaxis()->SetTitle("events / 8 GeV");
  hstack_uParMinusZpt->GetYaxis()->SetTitleSize(0.04);
  hstack_uParMinusZpt->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    leguParMinusZpt->AddEntry(hMCuParMinusZpt[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  leguParMinusZpt->AddEntry(hDatauParMinusZpt,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  leguParMinusZpt->Draw(); 
  leguParMinusZpt->SetMargin(0.3); 
  leguParMinusZpt->SetBorderSize(0);
  subpadUpar_2->cd();
  ratio = new TH1D(*hDatauParMinusZpt);
  ratio->Divide((TH1D*)hstack_uParMinusZpt->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("u_{||} - Z_{pT} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cuParMinusZpt->SaveAs( (plotDirectoryPath + cuParMinusZpt->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadUperp_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadUperp_2 = NULL; 

  TCanvas *cuPerp = new TCanvas("uPerp","",700,700);
  TLegend *leguPerp = new TLegend(0.7,0.7,0.89,0.89);

  subpadUperp_1 = new TPad("padUperp_1","",0.0,0.28,1.0,1.0);
  subpadUperp_1->SetLogy();
  //subpadUperp_1->SetBottomMargin(0);
  subpadUperp_2 = new TPad("padUperp_2","",0.0,0.0,1.0,0.27);
  subpadUperp_2->SetGridy();
  //subpadUperp_2->SetTopMargin(0);
  subpadUperp_2->SetBottomMargin(0.3);
  subpadUperp_1->Draw();
  subpadUperp_2->Draw();

  subpadUperp_1->cd();
  hstack_uPerp->Draw("HIST");
  hstack_uPerp->SetMinimum(0.03);
  hstack_uPerp->SetMaximum(3000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_uPerp->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatauPerp->SetMarkerStyle(8); // large dot
  hDatauPerp->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_uPerp->GetXaxis()->SetTitle("u_{#perp}  [GeV]");
  hstack_uPerp->GetXaxis()->SetTitleSize(0.04);
  hstack_uPerp->GetYaxis()->SetTitle("events / 8 GeV");
  hstack_uPerp->GetYaxis()->SetTitleSize(0.04);
  hstack_uPerp->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    leguPerp->AddEntry(hMCuPerp[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  leguPerp->AddEntry(hDatauPerp,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  leguPerp->Draw(); 
  leguPerp->SetMargin(0.3); 
  leguPerp->SetBorderSize(0);
  subpadUperp_2->cd();
  ratio = new TH1D(*hDatauPerp);
  ratio->Divide((TH1D*)hstack_uPerp->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("u_{#perp}  [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cuPerp->SaveAs( (plotDirectoryPath + cuPerp->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadNvtx_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadNvtx_2 = NULL; 

  TCanvas *cnVertices = new TCanvas("nVertices","",700,700);
  TLegend *legnVertices = new TLegend(0.7,0.7,0.89,0.89);
 
  subpadNvtx_1 = new TPad("padNvtx_1","",0.0,0.28,1.0,1.0);
  subpadNvtx_1->SetLogy();
  //subpadNvtx_1->SetBottomMargin(0);
  subpadNvtx_2 = new TPad("padNvtx_2","",0.0,0.0,1.0,0.27);
  subpadNvtx_2->SetGridy();
  //subpadNvtx_2->SetTopMargin(0);
  subpadNvtx_2->SetBottomMargin(0.3);
  subpadNvtx_1->Draw();
  subpadNvtx_2->Draw();

  subpadNvtx_1->cd();
  hstack_nVertices->Draw("HIST");
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_nVertices->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatanVertices->SetMarkerStyle(8); // large dot
  hDatanVertices->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_nVertices->GetXaxis()->SetTitle("N_{vtx} [GeV]");
  hstack_nVertices->GetXaxis()->SetTitleSize(0.04);
  hstack_nVertices->GetYaxis()->SetTitle("events");
  hstack_nVertices->GetYaxis()->SetTitleSize(0.04);
  hstack_nVertices->GetYaxis()->CenterTitle();
  //hstack_nVertices->SetMinimum(0.3);
  for (Int_t j = 0; j < nFiles; j++) {
    legnVertices->AddEntry(hMCnVertices[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legnVertices->AddEntry(hDatanVertices,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legnVertices->Draw(); 
  legnVertices->SetMargin(0.3); 
  legnVertices->SetBorderSize(0);
  subpadNvtx_2->cd();
  ratio = new TH1D(*hDatanVertices);
  ratio->Divide((TH1D*)hstack_nVertices->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("N_{vtx} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cnVertices->SaveAs( (plotDirectoryPath + cnVertices->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadJ1pt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadJ1pt_2 = NULL; 

  TCanvas *cjet1Pt = new TCanvas("jet1Pt","",700,700);
  TLegend *legjet1Pt = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadJ1pt_1 = new TPad("padJ1pt_1","",0.0,0.28,1.0,1.0);
  subpadJ1pt_1->SetLogy();
  //subpadJ1pt_1->SetBottomMargin(0);
  subpadJ1pt_2 = new TPad("padJ1pt_2","",0.0,0.0,1.0,0.27);
  subpadJ1pt_2->SetGridy();
  //subpadJ1pt_2->SetTopMargin(0);
  subpadJ1pt_2->SetBottomMargin(0.3);
  subpadJ1pt_1->Draw();
  subpadJ1pt_2->Draw();

  subpadJ1pt_1->cd();
  hstack_jet1Pt->Draw("HIST");
  hstack_jet1Pt->SetMinimum(0.3);
  hstack_jet1Pt->SetMaximum(4000.0);
  // hstack_jet1Pt->SetMinimum(5000.0);
 if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_jet1Pt->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatajet1Pt->SetMarkerStyle(8); // large dot
  hDatajet1Pt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_jet1Pt->GetXaxis()->SetTitle("leading jet pT [GeV]");
  hstack_jet1Pt->GetXaxis()->SetTitleSize(0.04);
  //hstack_jet1Pt->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_jet1Pt->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_jet1Pt->GetYaxis()->SetTitleSize(0.04);
  hstack_jet1Pt->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legjet1Pt->AddEntry(hMCjet1Pt[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legjet1Pt->AddEntry(hDatajet1Pt,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legjet1Pt->Draw(); 
  legjet1Pt->SetMargin(0.3); 
  legjet1Pt->SetBorderSize(0);
  subpadJ1pt_2->cd();
  ratio = new TH1D(*hDatajet1Pt);
  ratio->Divide((TH1D*)hstack_jet1Pt->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("leading jet pT [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cjet1Pt->SaveAs( (plotDirectoryPath + cjet1Pt->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadZpt_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadZpt_2 = NULL; 

  TCanvas *czPtSpectrum = new TCanvas("zPtSpectrum","",700,700);
  TLegend *legzPtSpectrum = new TLegend(0.7,0.7,0.89,0.89);

  subpadZpt_1 = new TPad("padZpt_1","",0.0,0.28,1.0,1.0);
  subpadZpt_1->SetLogy();
  //subpadZpt_1->SetBottomMargin(0);
  subpadZpt_2 = new TPad("padZpt_2","",0.0,0.0,1.0,0.27);
  subpadZpt_2->SetGridy();
  //subpadZpt_2->SetTopMargin(0);
  subpadZpt_2->SetBottomMargin(0.3);
  subpadZpt_1->Draw();
  subpadZpt_2->Draw();

  subpadZpt_1->cd();
  hstack_zPtSpectrum->Draw("HIST");
  hstack_zPtSpectrum->SetMinimum(0.3);
  hstack_zPtSpectrum->SetMaximum(4000.0);
  //if (MCuncBand_flag) hstack_zPtSpectrum->Draw("E SAME");
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_zPtSpectrum->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatazPtSpectrum->SetMarkerStyle(8); // large dot
  hDatazPtSpectrum->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_zPtSpectrum->GetXaxis()->SetTitle("Z_{pT} [GeV]");
  hstack_zPtSpectrum->GetXaxis()->SetTitleSize(0.04);
  //hstack_zPtSpectrum->GetYaxis()->SetRangeUser(0.3,4000.0);
  gPad->Update();
  hstack_zPtSpectrum->GetYaxis()->SetTitle("events / 5 GeV");
  hstack_zPtSpectrum->GetYaxis()->SetTitleSize(0.04);
  hstack_zPtSpectrum->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legzPtSpectrum->AddEntry(hMCzPtSpectrum[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legzPtSpectrum->AddEntry(hDatazPtSpectrum,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legzPtSpectrum->Draw(); 
  legzPtSpectrum->SetMargin(0.3); 
  legzPtSpectrum->SetBorderSize(0);
  //czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );
  subpadZpt_2->cd();
  ratio = new TH1D(*hDatazPtSpectrum);
  ratio->Divide((TH1D*)hstack_zPtSpectrum->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;



  TPad *subpadMet_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadMet_2 = NULL; 

  TCanvas *cmetNoLep = new TCanvas("metNoLep","",700,700);
  TLegend *legmetNoLep = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadMet_1 = new TPad("padMet_1","",0.0,0.28,1.0,1.0);
  subpadMet_1->SetLogy();
  //subpadMet_1->SetBottomMargin(0);
  subpadMet_2 = new TPad("padMet_2","",0.0,0.0,1.0,0.27);
  subpadMet_2->SetGridy();
  //subpadMet_2->SetTopMargin(0);
  subpadMet_2->SetBottomMargin(0.3);
  subpadMet_1->Draw();
  subpadMet_2->Draw();

  subpadMet_1->cd();
  hstack_metNoLep->Draw("HIST");
  hstack_metNoLep->SetMinimum(0.3);
  hstack_metNoLep->SetMaximum(4000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_metNoLep->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatametNoLep->SetMarkerStyle(8); // large dot
  hDatametNoLep->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_metNoLep->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  hstack_metNoLep->GetXaxis()->SetTitleSize(0.04);
  hstack_metNoLep->GetYaxis()->SetTitle("events / 10 GeV");
  hstack_metNoLep->GetYaxis()->SetTitleSize(0.04);
  hstack_metNoLep->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legmetNoLep->AddEntry(hMCmetNoLep[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legmetNoLep->AddEntry(hDatametNoLep,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legmetNoLep->Draw(); 
  legmetNoLep->SetMargin(0.3); 
  legmetNoLep->SetBorderSize(0);
  subpadMet_2->cd();
  ratio = new TH1D(*hDatametNoLep);
  ratio->Divide((TH1D*)hstack_metNoLep->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cmetNoLep->SaveAs( (plotDirectoryPath + cmetNoLep->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


  TPad *subpadNjets_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpadNjets_2 = NULL; 

  TCanvas *cnjets = new TCanvas("njets","",700,700);
  TLegend *legnjets = new TLegend(0.7,0.7,0.89,0.89);
  
  subpadNjets_1 = new TPad("padNjets_1","",0.0,0.28,1.0,1.0);
  subpadNjets_1->SetLogy();
  //subpadNjets_1->SetBottomMargin(0);
  subpadNjets_2 = new TPad("padNjets_2","",0.0,0.0,1.0,0.27);
  subpadNjets_2->SetGridy();
  //subpadNjets_2->SetTopMargin(0);
  subpadNjets_2->SetBottomMargin(0.3);
  subpadNjets_1->Draw();
  subpadNjets_2->Draw();

  subpadNjets_1->cd();
  hstack_njets->Draw("HIST");
  hstack_njets->SetMinimum(0.3);
  hstack_njets->SetMaximum(15000.0);
  if (MCuncBand_flag) {
    TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_njets->GetStack()->Last())->DrawCopy("SAME"));
    stackCopy->SetFillColor(kGray+2);
    stackCopy->SetFillStyle(3017);
    stackCopy->Draw("E2 SAME");
  }
  hDatanjets->SetMarkerStyle(8); // large dot
  hDatanjets->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_njets->GetXaxis()->SetTitle("N_{jets}");
  hstack_njets->GetXaxis()->SetTitleSize(0.04);
  hstack_njets->GetXaxis()->SetRange(1,5);
  hstack_njets->GetYaxis()->SetTitle("events ");
  hstack_njets->GetYaxis()->SetTitleSize(0.04);
  hstack_njets->GetYaxis()->CenterTitle();
  for (Int_t j = 0; j < nFiles; j++) {
    legnjets->AddEntry(hMCnjets[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  legnjets->AddEntry(hDatanjets,Form("%s",data_TexLabel.c_str()),"p");
  gStyle->SetStatStyle(0);
  legnjets->Draw(); 
  legnjets->SetMargin(0.3); 
  legnjets->SetBorderSize(0);
  subpadNjets_2->cd();
  ratio = new TH1D(*hDatanjets);
  ratio->Divide((TH1D*)hstack_njets->GetStack()->Last());
  ratio->SetStats(0);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitle("N_{jets}");
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.07);
  ratio->GetYaxis()->SetTitle("data/MC");
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->DrawCopy("EX0");
  ratio->SetMarkerStyle(8);  //medium dot
  cnjets->SaveAs( (plotDirectoryPath + cnjets->GetName() + suffix + plotFileExtension).c_str() );
  delete ratio;


}
*/
