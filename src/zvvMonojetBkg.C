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

void zvvMonojetBkg(const Int_t spring15_25ns_flag = 1){

  // if spring15_25ns_flag == 1, Axe computed with spring15_25ns sample will be used. For any othe value, phys14 will be used for Axe (the latter was found to be lower)

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)
  gStyle->SetFillColor(10);

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/znunuEstimate/SRoverCR/";  
  string suffix = "_mumu";
  string plotFileExtension = ".pdf";

  if (spring15_25ns_flag == 1)  suffix += "_AxeFrom_spring15_25ns";
  else suffix += "_AxeFrom_phys14";

  // ZvvMC
  string f1name = "monojet_SR_spring15_25ns_ZJetsToNuNu.root";
  TH1D *h1 = NULL;
  string h1name = "HYieldsMetBin";

  //ZmumuMC
  string f2name = "zmumujets_CS_spring15_25ns_DYJetsToLL.root";
  TH1D *h2 = NULL;
  string h2name = "HYieldsMetBin";

  //Axe mumu
  string f3name;
  if (spring15_25ns_flag == 1)  {
    f3name = "zmumujets_Axe_spring15_25ns_metNoMuSkim200.root";
    cout << "using Axe computed with spring15_25ns_metNoMuSkim200 sample" << endl;
  } else {
    f3name = "zmumujets_Axe_phys14_noSkim.root";
    cout << "using Axe computed with phys14 sample (no skim)" << endl;
  }
  TH1D *h3 = NULL;
  string h3name = "HacceffW";

  Double_t Rbr = 5.942;  //ratio between BR(Z->vv) and BR(Z->mumu)
  Double_t RbrError = 0.019; // error on BR ratio

  TFile* f1 = TFile::Open(f1name.c_str(),"READ");

  if (!f1 || !f1->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<f1name<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  h1 = (TH1D*)f1->Get(h1name.c_str());
  if (!h1) {
    cout << "Error: h1 not found. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  TH1D *hzvvMC = (TH1D*) h1->Clone();

  TFile* f2 = TFile::Open(f2name.c_str(),"READ");

  if (!f2 || !f2->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<f2name<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  h2 = (TH1D*)f2->Get(h2name.c_str());
  if (!h2) {
    cout << "Error: h2 not found. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  TH1D *hzmumuMC = (TH1D*) h2->Clone();

  TFile* f3 = TFile::Open(f3name.c_str(),"READ");

  if (!f3 || !f3->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<f3name<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  h3 = (TH1D*)f3->Get(h3name.c_str());
  if (!h3) {
    cout << "Error: h3 not found. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  TH1D *hAxeMC = (TH1D*) h3->Clone();
  // ==========

  // znunu /zmumu from MC
  TH1D* hratio = (TH1D*) hzvvMC->Clone("hratio");
  hratio->Divide(hzmumuMC);

  TH1D* hBRratio = (TH1D*) hAxeMC->Clone("hBRratio");  // I clone this to get the same binning
  
  for (Int_t i = 1; i <= hBRratio->GetNbinsX(); i++) {  // now I fill it with BR(Z->vv)/BR(Z->ll)

    hBRratio->SetBinContent(i, Rbr);
    hBRratio->SetBinError(i, RbrError);

  }

  TH1D* hBRratioOverAxe = (TH1D*) hBRratio->Clone("hBRratioOverAxe");
  hBRratioOverAxe->Divide(hAxeMC);

  // now here we go with the canvas
  TH1D * ratioplot = NULL; // will use it for the ratio plots

  TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpad_2 = NULL; 
  TCanvas *c = new TCanvas("zvvBkg_scaleFactor","MC scale factor",700,700);

  TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
  
  subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  //subpad_1->SetBottomMargin(0);
  subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  subpad_2->SetGridy();
  //subpad_2->SetTopMargin(0);
  subpad_2->SetBottomMargin(0.3);
  subpad_1->Draw();
  subpad_2->Draw();

  subpad_1->cd();
  hratio->SetTitle("");
  hratio->SetStats(0);  // no statistics box
  hratio->Draw("HIST E");
  hratio->SetLineColor(kBlue);
  hBRratioOverAxe->SetLineColor(kRed);
  hBRratioOverAxe->Draw("HE SAME"); 
  //hratio->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  hratio->GetXaxis()->SetTitle("");
  hratio->GetXaxis()->SetTitleSize(0.06);
  hratio->GetYaxis()->SetTitle("Z(#nu#nu)/Z(#mu#mu) scale factor");
  hratio->GetYaxis()->SetTitleSize(0.06);
  hratio->GetYaxis()->SetTitleOffset(0.7);
  hratio->GetYaxis()->CenterTitle();
  hratio->GetYaxis()->SetRangeUser(4.5,11.0);
  leg->AddEntry(hratio,"MC ratio","lf");
  leg->AddEntry(hBRratioOverAxe,"Axe","lf");
  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  subpad_2->cd();
  ratioplot = new TH1D(*hratio);
  ratioplot->Divide(hBRratioOverAxe);
  ratioplot->SetStats(0);
  ratioplot->GetXaxis()->SetLabelSize(0.10);
  ratioplot->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  ratioplot->GetXaxis()->SetTitleSize(0.15);
  ratioplot->GetXaxis()->SetTitleOffset(0.8);
  ratioplot->GetYaxis()->SetLabelSize(0.10);
  ratioplot->GetYaxis()->SetTitle("ratio");
  ratioplot->GetYaxis()->SetTitleSize(0.15);
  ratioplot->GetYaxis()->SetTitleOffset(0.3);
  ratioplot->GetYaxis()->CenterTitle();
  ratioplot->GetYaxis()->SetRangeUser(0.5,1.5);
  ratioplot->GetYaxis()->SetNdivisions(011);
  ratioplot->DrawCopy("HE");       // use DrawCopy so that I can delete ratioplot pointer and use it again in the next canvas
  ratioplot->SetMarkerStyle(8);  //medium dot
  c->SaveAs( (plotDirectoryPath + c->GetName() + suffix + plotFileExtension).c_str() );
  delete ratioplot;


  TH1D *ZvvEstimatedFromZll = (TH1D*) hzmumuMC->Clone("ZvvEstimatedFromZll");
  ZvvEstimatedFromZll->Multiply(hBRratioOverAxe);

  TH1D *hzvvNormalized = NULL;
  TH1D *hzllNormalized = NULL;

  TPad *subpad2_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpad2_2 = NULL; 
  TCanvas *c2 = new TCanvas("metShape_zvvANDzll","met shape",700,700);

  TLegend *leg2 = new TLegend(0.6,0.7,0.89,0.89);
  
  subpad2_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  subpad2_1->SetLogy();
  //subpad2_1->SetBottomMargin(0);
  subpad2_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  subpad2_2->SetGridy();
  //subpad2_2->SetTopMargin(0);
  subpad2_2->SetBottomMargin(0.3);
  subpad2_1->Draw();
  subpad2_2->Draw();

  subpad2_1->cd();
  hzvvNormalized = (TH1D*) hzvvMC->DrawNormalized("HIST E");
  hzvvNormalized->SetLineColor(kBlue);
  hzvvNormalized->SetTitle("");
  hzvvNormalized->SetStats(0);  // no statistics box
  hzllNormalized = (TH1D*) ZvvEstimatedFromZll->DrawNormalized("HE SAME"); 
  hzllNormalized->SetLineColor(kRed);
  //hzvvMC->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  hzvvNormalized->GetXaxis()->SetTitle("");
  hzvvNormalized->GetXaxis()->SetTitleSize(0.06);
  hzvvNormalized->GetYaxis()->SetTitle("a.u.");
  hzvvNormalized->GetYaxis()->SetTitleSize(0.06);
  hzvvNormalized->GetYaxis()->SetTitleOffset(0.7);
  hzvvNormalized->GetYaxis()->CenterTitle();
  leg2->AddEntry(hzvvNormalized,"Z(#nu#nu)+jets MC","lf");
  leg2->AddEntry(hzllNormalized,"Z(#mu#mu)+jets MC","lf");
  gStyle->SetStatStyle(0);
  leg2->Draw(); 
  leg2->SetMargin(0.3); 
  leg2->SetBorderSize(0);
  subpad2_2->cd();
  ratioplot = new TH1D(*hzvvNormalized);
  ratioplot->Divide(hzllNormalized);
  ratioplot->SetStats(0);
  ratioplot->GetXaxis()->SetLabelSize(0.10);
  ratioplot->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  ratioplot->GetXaxis()->SetTitleSize(0.14);
  ratioplot->GetXaxis()->SetTitleOffset(0.8);
  ratioplot->GetYaxis()->SetLabelSize(0.10);
  ratioplot->GetYaxis()->SetTitle("Z(#nu#nu)/Z(#mu#mu)");
  ratioplot->GetYaxis()->SetTitleSize(0.15);
  ratioplot->GetYaxis()->SetTitleOffset(0.3);
  ratioplot->GetYaxis()->CenterTitle();
  ratioplot->GetYaxis()->SetRangeUser(0.5,1.5);
  ratioplot->GetYaxis()->SetNdivisions(011);
  ratioplot->DrawCopy("HE");   // use DrawCopy so that I can delete ratioplot pointer and use it again in the next canvas
  ratioplot->SetMarkerStyle(8);  //medium dot
  c2->SaveAs( (plotDirectoryPath + c2->GetName() + suffix + plotFileExtension).c_str() );


}
//end of macro

