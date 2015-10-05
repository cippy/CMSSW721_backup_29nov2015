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
#include <algorithm>  // to use the "reverse" function to everse the order in the array

#include <Rtypes.h> // to use kColor

using namespace std;

void setSampleName(const Int_t signalRegion_flag, vector<string> &sampleName, vector<string> &MC_TexLabel) {

  if (signalRegion_flag == 1) {

    sampleName.push_back("GJets");
    sampleName.push_back("DYJetsToLL");
    sampleName.push_back("QCD");
    sampleName.push_back("Diboson");
    sampleName.push_back("Top");
    sampleName.push_back("WJetsToLNu");
    sampleName.push_back("ZJetsToNuNu");
       
    MC_TexLabel.push_back("#gamma + jets");
    MC_TexLabel.push_back("Z(ll)+jets");
    MC_TexLabel.push_back("QCD multijets");
    MC_TexLabel.push_back("Diboson");
    MC_TexLabel.push_back("t#bar{t},single t");
    MC_TexLabel.push_back("W(l#nu)+jets");
    MC_TexLabel.push_back("Z(#nu#nu)+jets");
    
  } else {

    sampleName.push_back("ZJetsToNuNu");
    sampleName.push_back("GJets");
    sampleName.push_back("QCD");
    sampleName.push_back("WJetsToLNu");
    sampleName.push_back("Diboson");
    sampleName.push_back("Top");
    sampleName.push_back("DYJetsToLL");

    MC_TexLabel.push_back("Z(#nu#nu)+jets");
    MC_TexLabel.push_back("#gamma + jets");
    MC_TexLabel.push_back("QCD multijets");
    MC_TexLabel.push_back("W(l#nu)+jets");
    MC_TexLabel.push_back("Diboson");
    MC_TexLabel.push_back("t#bar{t},single t");
    MC_TexLabel.push_back("Z(#mu#mu)+jets");

  }

}


void setHistColor(vector<Int_t> &histColor, const Int_t nSamples) {

  Int_t colorList[] = {kCyan, kViolet, kBlue, kRed, kOrange+1, kYellow, kGreen};  // the first color is for the main object. This array may contain more values than nSamples

  for (Int_t i = 0; i < nSamples; i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)

    histColor.push_back(colorList[i]);

  }

  std::reverse(histColor.begin(), histColor.end());   // reverse order in the array

}

void metDistribution(const Int_t signalRegion_flag = 1) {

  // if signalRegion_flag == 1 (default), will do met distribution in the monojet signal region, else it will do the control region

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)
  gStyle->SetFillColor(10);

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/monojet/met_distribution/";
  //string plotDirectoryPath = "./";
  string plotFileExtension = ".pdf";
  string suffix = "_mumu";

  TH1D* hmet = NULL;
  vector<TH1D*> hMCmetNoLep;

  vector<string> sampleName;
  vector<string> MC_TexLabel;
  setSampleName(signalRegion_flag, sampleName, MC_TexLabel);

  Int_t nFiles = (Int_t)sampleName.size();

  vector<Int_t> histColor;
   setHistColor(histColor, nFiles);

   string filenameBase;
   string canvasName;
   if (signalRegion_flag == 1) {
     filenameBase = "monojet_SR_spring15_25ns_";
     canvasName = "metDistribution_monojetSR";
   } else {
     filenameBase = "zmumujets_CS_spring15_25ns_";
     canvasName = "metDistribution_zjetsCS" + suffix;
   }
 
  string filenameExtension = ".root";

  vector<string> MCfileName;
  for (Int_t i = 0; i < nFiles; i++){
    MCfileName.push_back(filenameBase + sampleName[i] + filenameExtension);
  }

  for(Int_t i = 0; i < nFiles; i++) {

    cout<<"fileName : "<<MCfileName[i]<<endl;

    TFile* f = TFile::Open(MCfileName[i].c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<MCfileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }

    //cout << "check 1 " << endl;    

    hmet = (TH1D*)f->Get("HYieldsMetBin");
    if (!hmet) {
      cout << "Error: histogram not found in file ' " << MCfileName[i] << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hMCmetNoLep.push_back( (TH1D*)hmet->Clone() );

  } 

  THStack* hstack_metNoLep = new THStack("hstack_metNoLep","");

  for (Int_t j = 0; j < nFiles; j++) {

    for (Int_t i = 1; i <= hMCmetNoLep[j]->GetNbinsX(); i++) {

      hMCmetNoLep[j]->SetBinError(i,sqrt(hMCmetNoLep[j]->GetBinContent(i)));

    }

    hMCmetNoLep[j]->SetFillColor(histColor[j]);
    hstack_metNoLep->Add(hMCmetNoLep[j]);

  }

  // now here we go with the canvas
  // TH1D * ratioplot = NULL; // will use it for the ratio plots

  // TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas
  // TPad *subpad_2 = NULL; 
  TCanvas *c = new TCanvas(canvasName.c_str(),"met distribution");
  c->SetLogy();

  TLegend *leg = new TLegend(0.6,0.55,0.89,0.89);
  
  // subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  // //subpad_1->SetBottomMargin(0);
  // subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  // subpad_2->SetGridy();
  // //subpad_2->SetTopMargin(0);
  // subpad_2->SetBottomMargin(0.3);
  // subpad_1->Draw();
  // subpad_2->Draw();

  //subpad_1->cd();
  hstack_metNoLep->Draw("HIST");
  //hstack_metNoLep->SetMinimum(0.3);
  //hstack_metNoLep->SetMaximum(4000.0);
  TH1D* stackCopy = (TH1D*)(((TH1D*)hstack_metNoLep->GetStack()->Last())->DrawCopy("E2 SAME"));
  stackCopy->SetFillColor(kBlack);
  stackCopy->SetFillStyle(3017);
  hstack_metNoLep->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  hstack_metNoLep->GetXaxis()->SetTitleSize(0.06);
  hstack_metNoLep->GetXaxis()->SetTitleOffset(0.6);
  hstack_metNoLep->GetYaxis()->SetTitle("events");
  hstack_metNoLep->GetYaxis()->SetTitleSize(0.06);
  hstack_metNoLep->GetYaxis()->SetTitleOffset(0.8);
  hstack_metNoLep->GetYaxis()->CenterTitle();
  for (Int_t j = (nFiles-1); j >= 0; j--) {
    leg->AddEntry(hMCmetNoLep[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }
  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);

  // subpad_2->cd();
  // ratioplot = new TH1D(*hratio);
  // ratioplot->Divide(hBRratioOverAxe);
  // ratioplot->SetStats(0);
  // ratioplot->GetXaxis()->SetLabelSize(0.10);
  // ratioplot->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
  // ratioplot->GetXaxis()->SetTitleSize(0.15);
  // ratioplot->GetXaxis()->SetTitleOffset(0.8);
  // ratioplot->GetYaxis()->SetLabelSize(0.10);
  // ratioplot->GetYaxis()->SetTitle("ratio");
  // ratioplot->GetYaxis()->SetTitleSize(0.15);
  // ratioplot->GetYaxis()->SetTitleOffset(0.3);
  // ratioplot->GetYaxis()->CenterTitle();
  // ratioplot->GetYaxis()->SetRangeUser(0.5,1.5);
  // ratioplot->GetYaxis()->SetNdivisions(011);
  // ratioplot->DrawCopy("HE");
  // ratioplot->SetMarkerStyle(8);  //medium dot
  c->SaveAs( (plotDirectoryPath + c->GetName() + plotFileExtension).c_str() );

}
