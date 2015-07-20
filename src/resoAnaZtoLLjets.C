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

void resoAnaZtoLLjets(const char* fmumuName = "zmumujets_resoResp_noSkim_light.root", const char* feeName = "zeejets_resoResp_noSkim_light.root", const char* suffix = "") {

  // const char* fmumuName = "zmumujets_resoResp_noSkim_light.root";
  // const char* feeName = "zeejets_resoResp_noSkim_light.root";

  // suffix is a string used to name files without overwriting already existing ones

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/resoRespAnalysis/";
  string plotFileExtension = ".pdf";

  TFile* fmumu = TFile::Open(fmumuName,"READ");

  if (!fmumu->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fmumuName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  TGraphErrors *grmumu_resoPar_vs_nvtx = (TGraphErrors*)fmumu->Get("gr_resolution_uPar_vs_Nvtx");
  TGraphErrors *grmumu_resoPar_vs_nvtx_lowZpT = (TGraphErrors*)fmumu->Get("gr_resolution_uPar_vs_Nvtx_lowZpT");
  TGraphErrors *grmumu_resoPerp_vs_nvtx = (TGraphErrors*)fmumu->Get("gr_resolution_uPerp_vs_Nvtx");
  TGraphErrors *grmumu_resoPar_vs_ZpT = (TGraphErrors*)fmumu->Get("gr_resolution_uPar_vs_ZpT");
  TGraphErrors *grmumu_resoPerp_vs_ZpT = (TGraphErrors*)fmumu->Get("gr_resolution_uPerp_vs_ZpT");
  TGraphErrors *grmumu_responseCurve = (TGraphErrors*)fmumu->Get("gr_responseCurve");

  if (!grmumu_resoPar_vs_nvtx || !grmumu_resoPar_vs_nvtx_lowZpT || !grmumu_resoPerp_vs_nvtx || !grmumu_resoPar_vs_ZpT || !grmumu_resoPerp_vs_ZpT || !grmumu_responseCurve) {
    cout <<" Error: one ore more objects not read from file" << fmumuName << endl;
    exit(EXIT_FAILURE);
  }

  grmumu_resoPar_vs_nvtx->SetMarkerColor(kRed);
  grmumu_resoPar_vs_nvtx_lowZpT->SetMarkerColor(kRed);
  grmumu_resoPerp_vs_nvtx->SetMarkerColor(kRed);
  grmumu_resoPar_vs_ZpT->SetMarkerColor(kRed);
  grmumu_resoPerp_vs_ZpT->SetMarkerColor(kRed);
  grmumu_responseCurve->SetMarkerColor(kRed);
  grmumu_resoPar_vs_nvtx->SetMarkerStyle(21);
  grmumu_resoPar_vs_nvtx_lowZpT->SetMarkerStyle(20);
  grmumu_resoPerp_vs_nvtx->SetMarkerStyle(21);
  grmumu_resoPar_vs_ZpT->SetMarkerStyle(21);
  grmumu_resoPerp_vs_ZpT->SetMarkerStyle(21);
  grmumu_responseCurve->SetMarkerStyle(21);

  TFile* fee = TFile::Open(feeName,"READ");

  if (!fee->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<feeName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  TGraphErrors *gree_resoPar_vs_nvtx = (TGraphErrors*)fee->Get("gr_resolution_uPar_vs_Nvtx");
  TGraphErrors *gree_resoPar_vs_nvtx_lowZpT = (TGraphErrors*)fee->Get("gr_resolution_uPar_vs_Nvtx_lowZpT");
  TGraphErrors *gree_resoPerp_vs_nvtx = (TGraphErrors*)fee->Get("gr_resolution_uPerp_vs_Nvtx");
  TGraphErrors *gree_resoPar_vs_ZpT = (TGraphErrors*)fee->Get("gr_resolution_uPar_vs_ZpT");
  TGraphErrors *gree_resoPerp_vs_ZpT = (TGraphErrors*)fee->Get("gr_resolution_uPerp_vs_ZpT");
  TGraphErrors *gree_responseCurve = (TGraphErrors*)fee->Get("gr_responseCurve");

  if (!gree_resoPar_vs_nvtx || !gree_resoPar_vs_nvtx_lowZpT || !gree_resoPerp_vs_nvtx || !gree_resoPar_vs_ZpT || !gree_resoPerp_vs_ZpT || !gree_responseCurve) {
    cout <<" Error: one ore more objects not read from file" << feeName << endl;
    exit(EXIT_FAILURE);
  }

  gree_resoPar_vs_nvtx->SetMarkerColor(kBlue);
  gree_resoPar_vs_nvtx_lowZpT->SetMarkerColor(kBlue);
  gree_resoPerp_vs_nvtx->SetMarkerColor(kBlue);
  gree_resoPar_vs_ZpT->SetMarkerColor(kBlue);
  gree_resoPerp_vs_ZpT->SetMarkerColor(kBlue);
  gree_responseCurve->SetMarkerColor(kBlue);
  gree_resoPar_vs_nvtx->SetMarkerStyle(22);
  gree_resoPar_vs_nvtx_lowZpT->SetMarkerStyle(23);
  gree_resoPerp_vs_nvtx->SetMarkerStyle(22);
  gree_resoPar_vs_ZpT->SetMarkerStyle(22);
  gree_resoPerp_vs_ZpT->SetMarkerStyle(22);
  gree_responseCurve->SetMarkerStyle(22);

  TCanvas *c_resoPar_vs_nvtx = new TCanvas("c_resoPar_vs_nvtx","");
  TLegend *leg_resoPar_vs_nvtx = new TLegend(0.38,0.78,0.49,0.89);
  TMultiGraph *mg_resoPar_vs_nvtx = new TMultiGraph();
  mg_resoPar_vs_nvtx->Add(grmumu_resoPar_vs_nvtx,"P");
  mg_resoPar_vs_nvtx->Add(gree_resoPar_vs_nvtx,"P");
  mg_resoPar_vs_nvtx->Draw("A");
  mg_resoPar_vs_nvtx->GetXaxis()->SetTitle("nvtx");
  mg_resoPar_vs_nvtx->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_nvtx->GetYaxis()->SetTitleOffset(1.4);
  mg_resoPar_vs_nvtx->GetYaxis()->SetRangeUser(15.0,60.0);
  TPaveText *Pave_c_resoPar_vs_nvtx = new TPaveText(0.38,0.67,0.65,0.72,"NB NDC"); //NB abort borders drawing (no shadows), NDC is to use NDC coordinates
  Pave_c_resoPar_vs_nvtx->AddText("250 < ZpT < 500 GeV");  
  Pave_c_resoPar_vs_nvtx->Draw();
  leg_resoPar_vs_nvtx->AddEntry(grmumu_resoPar_vs_nvtx,"#mu#mu","lp");
  leg_resoPar_vs_nvtx->AddEntry(gree_resoPar_vs_nvtx,"ee","lp");
  leg_resoPar_vs_nvtx->Draw(); 
  leg_resoPar_vs_nvtx->SetMargin(0.3); 
  leg_resoPar_vs_nvtx->SetBorderSize(0);
  c_resoPar_vs_nvtx->SaveAs( (plotDirectoryPath + c_resoPar_vs_nvtx->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPar_vs_nvtx_lowZpT = new TCanvas("c_resoPar_vs_nvtx_lowZpT","");
  TLegend *leg_resoPar_vs_nvtx_lowZpT = new TLegend(0.38,0.78,0.49,0.89);
  TMultiGraph *mg_resoPar_vs_nvtx_lowZpT = new TMultiGraph();
  mg_resoPar_vs_nvtx_lowZpT->Add(grmumu_resoPar_vs_nvtx_lowZpT,"P");
  mg_resoPar_vs_nvtx_lowZpT->Add(gree_resoPar_vs_nvtx_lowZpT,"P");
  mg_resoPar_vs_nvtx_lowZpT->Draw("A");
  mg_resoPar_vs_nvtx_lowZpT->GetXaxis()->SetTitle("nvtx");
  mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetTitleOffset(1.4);
  mg_resoPar_vs_nvtx_lowZpT->GetYaxis()->SetRangeUser(15.0,35.0);
  TPaveText *Pave_c_resoPar_vs_nvtx_lowZpT = new TPaveText(0.38,0.67,0.65,0.72,"NB NDC"); //NB abort borders drawing (no shadows), NDC is to use NDC coordinates
  Pave_c_resoPar_vs_nvtx_lowZpT->AddText("20 < ZpT < 250 GeV");  
  Pave_c_resoPar_vs_nvtx_lowZpT->Draw();
  leg_resoPar_vs_nvtx_lowZpT->AddEntry(grmumu_resoPar_vs_nvtx_lowZpT,"#mu#mu","lp");
  leg_resoPar_vs_nvtx_lowZpT->AddEntry(gree_resoPar_vs_nvtx_lowZpT,"ee","lp");
  leg_resoPar_vs_nvtx_lowZpT->Draw(); 
  leg_resoPar_vs_nvtx_lowZpT->SetMargin(0.3); 
  leg_resoPar_vs_nvtx_lowZpT->SetBorderSize(0);
  c_resoPar_vs_nvtx_lowZpT->SaveAs( (plotDirectoryPath + c_resoPar_vs_nvtx_lowZpT->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPar_vs_nvtx_comp = new TCanvas("c_resoPar_vs_nvtx_comp","");
  TLegend *leg_resoPar_vs_nvtx_comp = new TLegend(0.28,0.67,0.58,0.89);
  TMultiGraph *mg_resoPar_vs_nvtx_comp = new TMultiGraph();
  mg_resoPar_vs_nvtx_comp->Add(grmumu_resoPar_vs_nvtx_lowZpT,"P");
  mg_resoPar_vs_nvtx_comp->Add(gree_resoPar_vs_nvtx_lowZpT,"P");
  mg_resoPar_vs_nvtx_comp->Add(grmumu_resoPar_vs_nvtx,"P");
  mg_resoPar_vs_nvtx_comp->Add(gree_resoPar_vs_nvtx,"P");
  mg_resoPar_vs_nvtx_comp->Draw("A");
  mg_resoPar_vs_nvtx_comp->GetXaxis()->SetTitle("nvtx");
  mg_resoPar_vs_nvtx_comp->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_nvtx_comp->GetYaxis()->SetTitleOffset(1.4);
  mg_resoPar_vs_nvtx_comp->GetYaxis()->SetRangeUser(15.0,50.0);
  leg_resoPar_vs_nvtx_comp->AddEntry((TObject*)0,"20 < ZpT < 250 GeV","");
  leg_resoPar_vs_nvtx_comp->AddEntry(grmumu_resoPar_vs_nvtx_lowZpT,"#mu#mu","lp");
  leg_resoPar_vs_nvtx_comp->AddEntry(gree_resoPar_vs_nvtx_lowZpT,"ee","lp");
  leg_resoPar_vs_nvtx_comp->AddEntry((TObject*)0,"","");
  leg_resoPar_vs_nvtx_comp->AddEntry((TObject*)0,"250 < ZpT < 500 GeV","");
  leg_resoPar_vs_nvtx_comp->AddEntry(grmumu_resoPar_vs_nvtx,"#mu#mu","lp");
  leg_resoPar_vs_nvtx_comp->AddEntry(gree_resoPar_vs_nvtx,"ee","lp");  
  leg_resoPar_vs_nvtx_comp->Draw(); 
  leg_resoPar_vs_nvtx_comp->SetMargin(0.3); 
  leg_resoPar_vs_nvtx_comp->SetBorderSize(0);
  c_resoPar_vs_nvtx_comp->SaveAs( (plotDirectoryPath + c_resoPar_vs_nvtx_comp->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPerp_vs_nvtx = new TCanvas("c_resoPerp_vs_nvtx","");
  TMultiGraph *mg_resoPerp_vs_nvtx = new TMultiGraph();
  TLegend *leg_resoPerp_vs_nvtx = new TLegend(0.38,0.78,0.49,0.89);
  mg_resoPerp_vs_nvtx->Add(grmumu_resoPerp_vs_nvtx,"AP");
  mg_resoPerp_vs_nvtx->Add(gree_resoPerp_vs_nvtx,"P");
  mg_resoPerp_vs_nvtx->Draw("A");
  mg_resoPerp_vs_nvtx->GetXaxis()->SetTitle("nvtx");
  mg_resoPerp_vs_nvtx->GetYaxis()->SetTitle("#sigma (u_{#perp}) [GeV]");
  mg_resoPerp_vs_nvtx->GetYaxis()->SetTitleOffset(1.4);
  mg_resoPerp_vs_nvtx->GetYaxis()->SetRangeUser(12.0,30.0);
  leg_resoPerp_vs_nvtx->AddEntry(grmumu_resoPerp_vs_nvtx,"#mu#mu","lp");
  leg_resoPerp_vs_nvtx->AddEntry(gree_resoPerp_vs_nvtx,"ee","lp");
  leg_resoPerp_vs_nvtx->Draw(); 
  leg_resoPerp_vs_nvtx->SetMargin(0.3); 
  leg_resoPerp_vs_nvtx->SetBorderSize(0);
  c_resoPerp_vs_nvtx->SaveAs( (plotDirectoryPath + c_resoPerp_vs_nvtx->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPar_vs_ZpT = new TCanvas("c_resoPar_vs_ZpT","");
  TMultiGraph *mg_resoPar_vs_ZpT = new TMultiGraph();
  TLegend *leg_resoPar_vs_ZpT = new TLegend(0.38,0.78,0.49,0.89);
  mg_resoPar_vs_ZpT->Add(grmumu_resoPar_vs_ZpT,"AP");
  mg_resoPar_vs_ZpT->Add(gree_resoPar_vs_ZpT,"P");
  mg_resoPar_vs_ZpT->Draw("A");
  mg_resoPar_vs_ZpT->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
  mg_resoPar_vs_ZpT->GetYaxis()->SetTitleOffset(1.4);
  mg_resoPar_vs_ZpT->GetYaxis()->SetRangeUser(15.0,60.0);
  leg_resoPar_vs_ZpT->AddEntry(grmumu_resoPar_vs_ZpT,"#mu#mu","lp");
  leg_resoPar_vs_ZpT->AddEntry(gree_resoPar_vs_ZpT,"ee","lp");
  leg_resoPar_vs_ZpT->Draw(); 
  leg_resoPar_vs_ZpT->SetMargin(0.3); 
  leg_resoPar_vs_ZpT->SetBorderSize(0);
  c_resoPar_vs_ZpT->SaveAs( (plotDirectoryPath + c_resoPar_vs_ZpT->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_resoPerp_vs_ZpT = new TCanvas("c_resoPerp_vs_ZpT","");
  TMultiGraph *mg_resoPerp_vs_ZpT = new TMultiGraph();
  TLegend *leg_resoPerp_vs_ZpT = new TLegend(0.38,0.78,0.49,0.89);
  mg_resoPerp_vs_ZpT->Add(grmumu_resoPerp_vs_ZpT,"AP");
  mg_resoPerp_vs_ZpT->Add(gree_resoPerp_vs_ZpT,"P");
  mg_resoPerp_vs_ZpT->Draw("A");
  mg_resoPerp_vs_ZpT->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_resoPerp_vs_ZpT->GetYaxis()->SetTitle("#sigma (u_{#perp}) [GeV]");
  mg_resoPerp_vs_ZpT->GetYaxis()->SetTitleOffset(1.4);
  mg_resoPerp_vs_ZpT->GetYaxis()->SetRangeUser(12.0,30.0);
  leg_resoPerp_vs_ZpT->AddEntry(grmumu_resoPerp_vs_ZpT,"#mu#mu","lp");
  leg_resoPerp_vs_ZpT->AddEntry(gree_resoPerp_vs_ZpT,"ee","lp");
  leg_resoPerp_vs_ZpT->Draw(); 
  leg_resoPerp_vs_ZpT->SetMargin(0.3); 
  leg_resoPerp_vs_ZpT->SetBorderSize(0);
  c_resoPerp_vs_ZpT->SaveAs( (plotDirectoryPath + c_resoPerp_vs_ZpT->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_responseCurve = new TCanvas("c_responseCurve","");
  TMultiGraph *mg_responseCurve = new TMultiGraph();
  TLegend *leg_responseCurve = new TLegend(0.38,0.78,0.49,0.89);
  mg_responseCurve->Add(grmumu_responseCurve,"AP");
  mg_responseCurve->Add(gree_responseCurve,"P");
  mg_responseCurve->Draw("A");
  mg_responseCurve->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_responseCurve->GetYaxis()->SetTitle("< u_{||} / ZpT >");
  mg_responseCurve->GetYaxis()->SetTitleOffset(1.4);
  leg_responseCurve->AddEntry(grmumu_responseCurve,"#mu#mu","lp");
  leg_responseCurve->AddEntry(gree_responseCurve,"ee","lp");
  leg_responseCurve->Draw(); 
  leg_responseCurve->SetMargin(0.3); 
  leg_responseCurve->SetBorderSize(0);
  c_responseCurve->SaveAs( (plotDirectoryPath + c_responseCurve->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *c_responseCurve_largeScale = new TCanvas("c_responseCurve_largeScale","");
  TMultiGraph *mg_responseCurve_largeScale = new TMultiGraph();
  TLegend *leg_responseCurve_largeScale = new TLegend(0.38,0.78,0.49,0.89);
  mg_responseCurve_largeScale->Add(grmumu_responseCurve,"AP");
  mg_responseCurve_largeScale->Add(gree_responseCurve,"P");
  mg_responseCurve_largeScale->Draw("A");
  mg_responseCurve_largeScale->GetXaxis()->SetTitle("ZpT [GeV]");
  mg_responseCurve_largeScale->GetYaxis()->SetTitle("< u_{||} / ZpT >");
  mg_responseCurve_largeScale->GetYaxis()->SetTitleOffset(1.4);
  mg_responseCurve_largeScale->GetYaxis()->SetRangeUser(0.5,1.1);
  leg_responseCurve_largeScale->AddEntry(grmumu_responseCurve,"#mu#mu","lp");
  leg_responseCurve_largeScale->AddEntry(gree_responseCurve,"ee","lp");
  leg_responseCurve_largeScale->Draw(); 
  leg_responseCurve_largeScale->SetMargin(0.3); 
  leg_responseCurve_largeScale->SetBorderSize(0);
  c_responseCurve_largeScale->SaveAs( (plotDirectoryPath + c_responseCurve_largeScale->GetName() + suffix + plotFileExtension).c_str());

  fmumu->Close();
  fee->Close();

}
