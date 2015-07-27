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

void resoAnaZtoLLjets_dataMC(const Double_t max_zpt_xaxis = 800, const char* fMCName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root", const char* fDataName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root", const char* suffix = "_spring15_50ns_GLTcutZmass80to100_mumu", const Double_t min_zpt_xaxis = 0) {

  // const char* fMCName = "zmumujets_resoResp_noSkim_light.root";
  // const char* fDataName = "zeejets_resoResp_noSkim_light.root";

  // suffix is a string used to name files without overwriting already existing ones

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/resoRespAnalysis/dataMC_comparison/";
  string plotFileExtension = ".pdf";

  string MC_TexLabel = "#mu#mu MC";
  string data_TexLabel = "data"; 

  Int_t default_YaxisRange = 1;

  TFile* fMC = TFile::Open(fMCName,"READ");

  if (!fMC->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fMCName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  TGraphErrors *grMC_resoPar_vs_nvtx = (TGraphErrors*)fMC->Get("gr_resolution_uPar_vs_Nvtx");
  TGraphErrors *grMC_resoPar_vs_nvtx_lowZpT = (TGraphErrors*)fMC->Get("gr_resolution_uPar_vs_Nvtx_lowZpT");
  TGraphErrors *grMC_resoPerp_vs_nvtx = (TGraphErrors*)fMC->Get("gr_resolution_uPerp_vs_Nvtx");
  TGraphErrors *grMC_resoPar_vs_ZpT = (TGraphErrors*)fMC->Get("gr_resolution_uPar_vs_ZpT");
  TGraphErrors *grMC_resoPerp_vs_ZpT = (TGraphErrors*)fMC->Get("gr_resolution_uPerp_vs_ZpT");
  TGraphErrors *grMC_responseCurve = (TGraphErrors*)fMC->Get("gr_responseCurve");
  TGraphErrors *grMC_responseCurve_0jets = (TGraphErrors*)fMC->Get("gr_responseCurve_0jets");

  if (!grMC_resoPar_vs_nvtx || !grMC_resoPar_vs_nvtx_lowZpT || !grMC_resoPerp_vs_nvtx || !grMC_resoPar_vs_ZpT || !grMC_resoPerp_vs_ZpT || !grMC_responseCurve || grMC_responseCurve_0jets) {
    cout <<" Error: one ore more objects not read from file" << fMCName << endl;
    exit(EXIT_FAILURE);
  }

  grMC_resoPar_vs_nvtx->SetMarkerColor(kRed);
  grMC_resoPar_vs_nvtx_lowZpT->SetMarkerColor(kRed);
  grMC_resoPerp_vs_nvtx->SetMarkerColor(kRed);
  grMC_resoPar_vs_ZpT->SetMarkerColor(kRed);
  grMC_resoPerp_vs_ZpT->SetMarkerColor(kRed);
  grMC_responseCurve->SetMarkerColor(kRed);
  grMC_responseCurve_0jets->SetMarkerColor(kRed);
  grMC_resoPar_vs_nvtx->SetMarkerStyle(21);
  grMC_resoPar_vs_nvtx_lowZpT->SetMarkerStyle(20);
  grMC_resoPerp_vs_nvtx->SetMarkerStyle(21);
  grMC_resoPar_vs_ZpT->SetMarkerStyle(21);
  grMC_resoPerp_vs_ZpT->SetMarkerStyle(21);
  grMC_responseCurve->SetMarkerStyle(21);
  grMC_responseCurve_0jets->SetMarkerStyle(21);

  TFile* fData = TFile::Open(fDataName,"READ");

  if (!fData->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fDataName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  TGraphErrors *grData_resoPar_vs_nvtx = (TGraphErrors*)fData->Get("gr_resolution_uPar_vs_Nvtx");
  TGraphErrors *grData_resoPar_vs_nvtx_lowZpT = (TGraphErrors*)fData->Get("gr_resolution_uPar_vs_Nvtx_lowZpT");
  TGraphErrors *grData_resoPerp_vs_nvtx = (TGraphErrors*)fData->Get("gr_resolution_uPerp_vs_Nvtx");
  TGraphErrors *grData_resoPar_vs_ZpT = (TGraphErrors*)fData->Get("gr_resolution_uPar_vs_ZpT");
  TGraphErrors *grData_resoPerp_vs_ZpT = (TGraphErrors*)fData->Get("gr_resolution_uPerp_vs_ZpT");
  TGraphErrors *grData_responseCurve = (TGraphErrors*)fData->Get("gr_responseCurve");
  TGraphErrors *grData_responseCurve_0jets = (TGraphErrors*)fData->Get("gr_responseCurve_0jets");

  if (!grData_resoPar_vs_nvtx || !grData_resoPar_vs_nvtx_lowZpT || !grData_resoPerp_vs_nvtx || !grData_resoPar_vs_ZpT || !grData_resoPerp_vs_ZpT || !grData_responseCurve || grData_responseCurve_0jets) {
    cout <<" Error: one ore more objects not read from file" << fDataName << endl;
    exit(EXIT_FAILURE);
  }

  grData_resoPar_vs_nvtx->SetMarkerColor(kBlue);
  grData_resoPar_vs_nvtx_lowZpT->SetMarkerColor(kBlue);
  grData_resoPerp_vs_nvtx->SetMarkerColor(kBlue);
  grData_resoPar_vs_ZpT->SetMarkerColor(kBlue);
  grData_resoPerp_vs_ZpT->SetMarkerColor(kBlue);
  grData_responseCurve->SetMarkerColor(kBlue);
  grData_responseCurve_0jets->SetMarkerColor(kBlue);
  grData_resoPar_vs_nvtx->SetMarkerStyle(22);
  grData_resoPar_vs_nvtx_lowZpT->SetMarkerStyle(23);
  grData_resoPerp_vs_nvtx->SetMarkerStyle(22);
  grData_resoPar_vs_ZpT->SetMarkerStyle(22);
  grData_resoPerp_vs_ZpT->SetMarkerStyle(22);
  grData_responseCurve->SetMarkerStyle(22);
  grData_responseCurve_0jets->SetMarkerStyle(22);

  TCanvas *c_resoPar_vs_nvtx = new TCanvas("resoPar_vs_nvtx","");
  TLegend *leg_resoPar_vs_nvtx = new TLegend(0.38,0.78,0.49,0.89);
  TMultiGraph *mg_resoPar_vs_nvtx = new TMultiGraph();
  mg_resoPar_vs_nvtx->Add(grMC_resoPar_vs_nvtx,"P");
  mg_resoPar_vs_nvtx->Add(grData_resoPar_vs_nvtx,"P");
  mg_resoPar_vs_nvtx->Draw("A");
  mg_resoPar_vs_nvtx->GetXaxis()->SetTitle("nvtx");
  mg_resoPar_vs_nvtx->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_nvtx->GetYaxis()->SetTitleOffset(1.2);
  mg_resoPar_vs_nvtx->GetYaxis()->SetTitleSize(0.04);
  if (default_YaxisRange) mg_resoPar_vs_nvtx->GetYaxis()->SetRangeUser(10.0,60.0);
  TPaveText *Pave_c_resoPar_vs_nvtx = new TPaveText(0.38,0.67,0.65,0.72,"NB NDC"); //NB abort borders drawing (no shadows), NDC is to use NDC coordinates
  Pave_c_resoPar_vs_nvtx->AddText("250 < ZpT < 500 GeV");  
  Pave_c_resoPar_vs_nvtx->Draw();
  leg_resoPar_vs_nvtx->AddEntry(grMC_resoPar_vs_nvtx,MC_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx->AddEntry(grData_resoPar_vs_nvtx,data_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx->Draw(); 
  leg_resoPar_vs_nvtx->SetMargin(0.3); 
  leg_resoPar_vs_nvtx->SetBorderSize(0);
  c_resoPar_vs_nvtx->SaveAs( (plotDirectoryPath + c_resoPar_vs_nvtx->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPar_vs_nvtx_lowZpT = new TCanvas("resoPar_vs_nvtx_lowZpT","");
  TLegend *leg_resoPar_vs_nvtx_lowZpT = new TLegend(0.38,0.78,0.49,0.89);
  TMultiGraph *mg_resoPar_vs_nvtx_lowZpT = new TMultiGraph();
  mg_resoPar_vs_nvtx_lowZpT->Add(grMC_resoPar_vs_nvtx_lowZpT,"P");
  mg_resoPar_vs_nvtx_lowZpT->Add(grData_resoPar_vs_nvtx_lowZpT,"P");
  mg_resoPar_vs_nvtx_lowZpT->Draw("A");
  mg_resoPar_vs_nvtx_lowZpT->GetXaxis()->SetTitle("nvtx");
  mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetTitleOffset(1.2);
  mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetTitleSize(0.04);
  if (default_YaxisRange) mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetRangeUser(8.0,35.0);
  TPaveText *Pave_c_resoPar_vs_nvtx_lowZpT = new TPaveText(0.50,0.84,0.80,0.89,"NB NDC"); //NB abort borders drawing (no shadows), NDC is to use NDC coordinates
  Pave_c_resoPar_vs_nvtx_lowZpT->SetShadowColor(0);   // shadow in white (probably it's not done at all)
  Pave_c_resoPar_vs_nvtx_lowZpT->AddText("ZpT < 250 GeV");  
  Pave_c_resoPar_vs_nvtx_lowZpT->Draw();
  leg_resoPar_vs_nvtx_lowZpT->AddEntry(grMC_resoPar_vs_nvtx_lowZpT,MC_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx_lowZpT->AddEntry(grData_resoPar_vs_nvtx_lowZpT,data_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx_lowZpT->Draw(); 
  leg_resoPar_vs_nvtx_lowZpT->SetMargin(0.3); 
  leg_resoPar_vs_nvtx_lowZpT->SetBorderSize(0);
  c_resoPar_vs_nvtx_lowZpT->SaveAs( (plotDirectoryPath + c_resoPar_vs_nvtx_lowZpT->GetName() + suffix + plotFileExtension).c_str());

  TGraphErrors *grMC_resoPar_vs_nvtx_clone = (TGraphErrors*)grMC_resoPar_vs_nvtx->Clone("gr_resolution_uPar_vs_Nvtx_clone");
  TGraphErrors *grMC_resoPar_vs_nvtx_lowZpT_clone = (TGraphErrors*)grMC_resoPar_vs_nvtx_lowZpT->Clone("gr_resolution_uPar_vs_Nvtx_lowZpT_clone");
  TGraphErrors *grData_resoPar_vs_nvtx_clone = (TGraphErrors*)grData_resoPar_vs_nvtx->Clone("gr_resolution_uPar_vs_Nvtx_clone");
  TGraphErrors *grData_resoPar_vs_nvtx_lowZpT_clone = (TGraphErrors*)grData_resoPar_vs_nvtx_lowZpT->Clone("gr_resolution_uPar_vs_Nvtx_lowZpT_clone");

  TCanvas *c_resoPar_vs_nvtx_comp = new TCanvas("resoPar_vs_nvtx_comp","");
  TLegend *leg_resoPar_vs_nvtx_comp = new TLegend(0.28,0.67,0.58,0.89);
  TMultiGraph *mg_resoPar_vs_nvtx_comp = new TMultiGraph();
  mg_resoPar_vs_nvtx_comp->Add(grMC_resoPar_vs_nvtx_lowZpT_clone,"P");
  mg_resoPar_vs_nvtx_comp->Add(grData_resoPar_vs_nvtx_lowZpT_clone,"P");
  mg_resoPar_vs_nvtx_comp->Add(grMC_resoPar_vs_nvtx_clone,"P");
  mg_resoPar_vs_nvtx_comp->Add(grData_resoPar_vs_nvtx_clone,"P");
  mg_resoPar_vs_nvtx_comp->Draw("A");
  mg_resoPar_vs_nvtx_comp->GetXaxis()->SetTitle("nvtx");
  mg_resoPar_vs_nvtx_comp->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_nvtx_comp->GetYaxis()->SetTitleOffset(1.2);
  mg_resoPar_vs_nvtx_comp->GetYaxis()->SetTitleSize(0.04);
  if (default_YaxisRange) mg_resoPar_vs_nvtx_comp->GetYaxis()->SetRangeUser(8.0,50.0);
  leg_resoPar_vs_nvtx_comp->AddEntry((TObject*)0,"ZpT < 250 GeV","");
  leg_resoPar_vs_nvtx_comp->AddEntry(grMC_resoPar_vs_nvtx_lowZpT_clone,MC_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx_comp->AddEntry(grData_resoPar_vs_nvtx_lowZpT_clone,data_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx_comp->AddEntry((TObject*)0,"","");
  leg_resoPar_vs_nvtx_comp->AddEntry((TObject*)0,"250 < ZpT < 500 GeV","");
  leg_resoPar_vs_nvtx_comp->AddEntry(grMC_resoPar_vs_nvtx_clone,MC_TexLabel.c_str(),"lp");
  leg_resoPar_vs_nvtx_comp->AddEntry(grData_resoPar_vs_nvtx_clone,data_TexLabel.c_str(),"lp");  
  leg_resoPar_vs_nvtx_comp->Draw(); 
  leg_resoPar_vs_nvtx_comp->SetMargin(0.3); 
  leg_resoPar_vs_nvtx_comp->SetBorderSize(0);
  c_resoPar_vs_nvtx_comp->SaveAs( (plotDirectoryPath + c_resoPar_vs_nvtx_comp->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPerp_vs_nvtx = new TCanvas("resoPerp_vs_nvtx","");
  TMultiGraph *mg_resoPerp_vs_nvtx = new TMultiGraph();
  TLegend *leg_resoPerp_vs_nvtx = new TLegend(0.38,0.78,0.49,0.89);
  mg_resoPerp_vs_nvtx->Add(grMC_resoPerp_vs_nvtx,"AP");
  mg_resoPerp_vs_nvtx->Add(grData_resoPerp_vs_nvtx,"P");
  mg_resoPerp_vs_nvtx->Draw("A");
  mg_resoPerp_vs_nvtx->GetXaxis()->SetTitle("nvtx");
  mg_resoPerp_vs_nvtx->GetYaxis()->SetTitle("#sigma (u_{#perp}) [GeV]");
  mg_resoPerp_vs_nvtx->GetYaxis()->SetTitleOffset(1.2);
  mg_resoPerp_vs_nvtx->GetYaxis()->SetTitleSize(0.04);
  if (default_YaxisRange) mg_resoPerp_vs_nvtx->GetYaxis()->SetRangeUser(8.0,35.0);
  leg_resoPerp_vs_nvtx->AddEntry(grMC_resoPerp_vs_nvtx,MC_TexLabel.c_str(),"lp");
  leg_resoPerp_vs_nvtx->AddEntry(grData_resoPerp_vs_nvtx,data_TexLabel.c_str(),"lp");
  leg_resoPerp_vs_nvtx->Draw(); 
  leg_resoPerp_vs_nvtx->SetMargin(0.3); 
  leg_resoPerp_vs_nvtx->SetBorderSize(0);
  c_resoPerp_vs_nvtx->SaveAs( (plotDirectoryPath + c_resoPerp_vs_nvtx->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPar_vs_ZpT = new TCanvas("resoPar_vs_ZpT","");
  TMultiGraph *mg_resoPar_vs_ZpT = new TMultiGraph();
  TLegend *leg_resoPar_vs_ZpT = new TLegend(0.38,0.78,0.49,0.89);
  mg_resoPar_vs_ZpT->Add(grMC_resoPar_vs_ZpT,"AP");
  mg_resoPar_vs_ZpT->Add(grData_resoPar_vs_ZpT,"P");
  mg_resoPar_vs_ZpT->Draw("A");
  mg_resoPar_vs_ZpT->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_resoPar_vs_ZpT->GetXaxis()->SetRangeUser(min_zpt_xaxis, max_zpt_xaxis);
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitleOffset(1.2);
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitleSize(0.04);
  mg_resoPar_vs_ZpT->GetYaxis()->SetRangeUser(10.0,40.0);
  leg_resoPar_vs_ZpT->AddEntry(grMC_resoPar_vs_ZpT,MC_TexLabel.c_str(),"lp");
  leg_resoPar_vs_ZpT->AddEntry(grData_resoPar_vs_ZpT,data_TexLabel.c_str(),"lp");
  leg_resoPar_vs_ZpT->Draw(); 
  leg_resoPar_vs_ZpT->SetMargin(0.3); 
  leg_resoPar_vs_ZpT->SetBorderSize(0);
  c_resoPar_vs_ZpT->SaveAs( (plotDirectoryPath + c_resoPar_vs_ZpT->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPerp_vs_ZpT = new TCanvas("resoPerp_vs_ZpT","");
  TMultiGraph *mg_resoPerp_vs_ZpT = new TMultiGraph();
  TLegend *leg_resoPerp_vs_ZpT = new TLegend(0.38,0.78,0.49,0.89);
  mg_resoPerp_vs_ZpT->Add(grMC_resoPerp_vs_ZpT,"AP");
  mg_resoPerp_vs_ZpT->Add(grData_resoPerp_vs_ZpT,"P");
  mg_resoPerp_vs_ZpT->Draw("A");
  mg_resoPerp_vs_ZpT->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_resoPerp_vs_ZpT->GetXaxis()->SetRangeUser(min_zpt_xaxis, max_zpt_xaxis);
  mg_resoPerp_vs_ZpT->GetYaxis()->SetTitle("#sigma (u_{#perp}) [GeV]");
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitleOffset(1.2);
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitleSize(0.04);
  if (default_YaxisRange) mg_resoPerp_vs_ZpT->GetYaxis()->SetRangeUser(12.0,35.0);
  leg_resoPerp_vs_ZpT->AddEntry(grMC_resoPerp_vs_ZpT,MC_TexLabel.c_str(),"lp");
  leg_resoPerp_vs_ZpT->AddEntry(grData_resoPerp_vs_ZpT,data_TexLabel.c_str(),"lp");
  leg_resoPerp_vs_ZpT->Draw(); 
  leg_resoPerp_vs_ZpT->SetMargin(0.3); 
  leg_resoPerp_vs_ZpT->SetBorderSize(0);
  c_resoPerp_vs_ZpT->SaveAs( (plotDirectoryPath + c_resoPerp_vs_ZpT->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_responseCurve = new TCanvas("responseCurve","");
  TMultiGraph *mg_responseCurve = new TMultiGraph();
  TLegend *leg_responseCurve = new TLegend(0.38,0.78,0.49,0.89);
  mg_responseCurve->Add(grMC_responseCurve,"AP");
  mg_responseCurve->Add(grData_responseCurve,"P");
  mg_responseCurve->Draw("A");
  mg_responseCurve->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_responseCurve->GetXaxis()->SetRangeUser(min_zpt_xaxis, max_zpt_xaxis);
  mg_responseCurve->GetYaxis()->SetTitle("< u_{||} / ZpT >");
  mg_responseCurve->GetYaxis()->SetTitleOffset(1.2);
  mg_responseCurve->GetYaxis()->SetTitleSize(0.04);
  leg_responseCurve->AddEntry(grMC_responseCurve,MC_TexLabel.c_str(),"lp");
  leg_responseCurve->AddEntry(grData_responseCurve,data_TexLabel.c_str(),"lp");
  leg_responseCurve->Draw(); 
  leg_responseCurve->SetMargin(0.3); 
  leg_responseCurve->SetBorderSize(0);
  c_responseCurve->SaveAs( (plotDirectoryPath + c_responseCurve->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_responseCurve_largeScale = new TCanvas("responseCurve_largeScale","");
  TMultiGraph *mg_responseCurve_largeScale = new TMultiGraph();
  TLegend *leg_responseCurve_largeScale = new TLegend(0.38,0.78,0.49,0.89);
  mg_responseCurve_largeScale->Add(grMC_responseCurve,"AP");
  mg_responseCurve_largeScale->Add(grData_responseCurve,"P");
  mg_responseCurve_largeScale->Draw("A");
  mg_responseCurve_largeScale->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_responseCurve_largeScale->GetXaxis()->SetRangeUser(min_zpt_xaxis, max_zpt_xaxis);
  mg_responseCurve_largeScale->GetYaxis()->SetTitle("< u_{||} / ZpT >");
  mg_responseCurve->GetYaxis()->SetTitleOffset(1.2);
  mg_responseCurve->GetYaxis()->SetTitleSize(0.04);
  mg_responseCurve_largeScale->GetYaxis()->SetRangeUser(0.8,1.1);
  leg_responseCurve_largeScale->AddEntry(grMC_responseCurve,MC_TexLabel.c_str(),"lp");
  leg_responseCurve_largeScale->AddEntry(grData_responseCurve,data_TexLabel.c_str(),"lp");
  leg_responseCurve_largeScale->Draw(); 
  leg_responseCurve_largeScale->SetMargin(0.3); 
  leg_responseCurve_largeScale->SetBorderSize(0);
  c_responseCurve_largeScale->SaveAs( (plotDirectoryPath + c_responseCurve_largeScale->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_responseCurve_0jets = new TCanvas("responseCurve_0jets","");
  TMultiGraph *mg_responseCurve_0jets = new TMultiGraph();
  TLegend *leg_responseCurve_0jets = new TLegend(0.38,0.78,0.49,0.89);
  mg_responseCurve_0jets->Add(grMC_responseCurve_0jets,"AP");
  mg_responseCurve_0jets->Add(grData_responseCurve_0jets,"P");
  mg_responseCurve_0jets->Draw("A");
  mg_responseCurve_0jets->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_responseCurve_0jets->GetXaxis()->SetRangeUser(min_zpt_xaxis, max_zpt_xaxis);
  mg_responseCurve_0jets->GetYaxis()->SetTitle("< u_{||} / ZpT >");
  mg_responseCurve_0jets->GetYaxis()->SetTitleOffset(1.2);
  mg_responseCurve_0jets->GetYaxis()->SetTitleSize(0.04);
  leg_responseCurve_0jets->AddEntry((TObject*)0,"0 jets events","");
  leg_responseCurve_0jets->AddEntry(grMC_responseCurve_0jets,MC_TexLabel.c_str(),"lp");
  leg_responseCurve_0jets->AddEntry(grData_responseCurve_0jets,data_TexLabel.c_str(),"lp");
  leg_responseCurve_0jets->Draw(); 
  leg_responseCurve_0jets->SetMargin(0.3); 
  leg_responseCurve_0jets->SetBorderSize(0);
  c_responseCurve_0jets->SaveAs( (plotDirectoryPath + c_responseCurve_0jets->GetName() + suffix + plotFileExtension).c_str());

  fMC->Close();
  fData->Close();

}
