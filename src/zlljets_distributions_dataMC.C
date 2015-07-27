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

void zlljets_distributions_dataMC(const char* suffix = "_GLTcutZmass80to100_mumu") {

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/distributions/dataMC_comparison/";
  string plotFileExtension = ".pdf";

  vector<string> MC_TexLabel;
  MC_TexLabel.push_back("Z(#tau#tau)+jets MC");
  MC_TexLabel.push_back("TT+Jets MC");
  MC_TexLabel.push_back("Z(#mu#mu)+jets MC");

  string data_TexLabel = "data"; 

  string DataFileName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";

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
  MCfileName.push_back("ztautaujets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
  MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
  MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  

  Int_t nFiles = (Int_t)MCfileName.size();
  Int_t histColor[] = {kRed,kBlue,kCyan};

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
    hMCjet1Pt.push_back((TH1D*)f->Get("Hjet1ptDistribution"));
    hMCuParMinusZpt.push_back((TH1D*)f->Get("H_uParMinusZpT_Distribution"));
    hMCuPerp.push_back((TH1D*)f->Get("H_uPerp_Distribution"));
    hMCnVertices.push_back((TH1D*)f->Get("HvtxDistribution"));
    hMCzPtSpectrum.push_back((TH1D*)f->Get("HzptDistribution"));
    hMCnjets.push_back((TH1D*)f->Get("HnjetsDistribution"));
   
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
  hDatauParMinusZpt = (TH1D*)f->Get("H_uParMinusZpT_Distribution");
  hDatauPerp = (TH1D*)f->Get("H_uPerp_Distribution");
  hDatanVertices = (TH1D*)f->Get("HvtxDistribution");
  hDatazPtSpectrum = (TH1D*)f->Get("HzptDistribution");
  hDatanjets = (TH1D*)f->Get("HnjetsDistribution");

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

  TCanvas *cInvMass = new TCanvas("invMass","");
  TLegend *legInvMass = new TLegend(0.7,0.7,0.89,0.89);
  //cInvMass->SetLogy();
  hstack_invMass->Draw("HE");
  //hstack_invMass->SetMinimum(0.3);
  hstack_invMass->Draw("E SAME");
  hDatainvMass->SetMarkerStyle(8); // large dot
  hDatainvMass->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_invMass->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
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
  cInvMass->SaveAs( (plotDirectoryPath + cInvMass->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cuParMinusZpt = new TCanvas("uParMinusZpt","");
  TLegend *leguParMinusZpt = new TLegend(0.7,0.7,0.89,0.89);
  cuParMinusZpt->SetLogy();
  hstack_uParMinusZpt->Draw("HE");
  hstack_uParMinusZpt->SetMinimum(0.3);
  hstack_uParMinusZpt->SetMaximum(3000.0);
  hstack_uParMinusZpt->Draw("E SAME");
  hDatauParMinusZpt->SetMarkerStyle(8); // large dot
  hDatauParMinusZpt->Draw("EX0 SAME"); //X0 doesn't draw x error
  hDatauParMinusZpt->GetYaxis()->SetRangeUser(0.3,3050.0);
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
  cuParMinusZpt->SaveAs( (plotDirectoryPath + cuParMinusZpt->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cuPerp = new TCanvas("uPerp","");
  TLegend *leguPerp = new TLegend(0.7,0.7,0.89,0.89);
  cuPerp->SetLogy();
  hstack_uPerp->Draw("HE");
  hstack_uPerp->SetMinimum(0.3);
  hstack_uPerp->SetMaximum(3000.0);
  hstack_uPerp->Draw("E SAME");
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
  cuPerp->SaveAs( (plotDirectoryPath + cuPerp->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cnVertices = new TCanvas("nVertices","");
  TLegend *legnVertices = new TLegend(0.7,0.7,0.89,0.89);
  cnVertices->SetLogy();
  hstack_nVertices->Draw("HE");
  hstack_nVertices->Draw("E SAME");
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
  cnVertices->SaveAs( (plotDirectoryPath + cnVertices->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cjet1Pt = new TCanvas("jet1Pt","");
  TLegend *legjet1Pt = new TLegend(0.7,0.7,0.89,0.89);
  cjet1Pt->SetLogy();
  hstack_jet1Pt->Draw("HE");
  hstack_jet1Pt->SetMinimum(0.3);
  // hstack_jet1Pt->SetMinimum(5000.0);
  hstack_jet1Pt->Draw("E SAME");
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
  cjet1Pt->SaveAs( (plotDirectoryPath + cjet1Pt->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *czPtSpectrum = new TCanvas("zPtSpectrum","");
  TLegend *legzPtSpectrum = new TLegend(0.7,0.7,0.89,0.89);
  czPtSpectrum->SetLogy();
  hstack_zPtSpectrum->Draw("HE");
  hstack_zPtSpectrum->SetMinimum(0.3);
  // hstack_zPtSpectrum->SetMinimum(5000.0);
  hstack_zPtSpectrum->Draw("E SAME");
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
  czPtSpectrum->SaveAs( (plotDirectoryPath + czPtSpectrum->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cmetNoLep = new TCanvas("metNoLep","");
  TLegend *legmetNoLep = new TLegend(0.7,0.7,0.89,0.89);
  cmetNoLep->SetLogy();
  hstack_metNoLep->Draw("HE");
  hstack_metNoLep->SetMinimum(0.3);
  hstack_metNoLep->SetMaximum(4000.0);
  hstack_metNoLep->Draw("E SAME");
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
  cmetNoLep->SaveAs( (plotDirectoryPath + cmetNoLep->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cnjets = new TCanvas("njets","");
  TLegend *legnjets = new TLegend(0.7,0.7,0.89,0.89);
  cnjets->SetLogy();
  hstack_njets->Draw("HE");
  hstack_njets->SetMinimum(0.3);
  hstack_njets->SetMaximum(4000.0);
  hstack_njets->Draw("E SAME");
  hDatanjets->SetMarkerStyle(8); // large dot
  hDatanjets->Draw("EX0 SAME"); //X0 doesn't draw x error
  hstack_njets->GetXaxis()->SetTitle("N_{jets}");
  hstack_njets->GetXaxis()->SetTitleSize(0.04);
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
  cnjets->SaveAs( (plotDirectoryPath + cnjets->GetName() + suffix + plotFileExtension).c_str() );

  


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
			     
