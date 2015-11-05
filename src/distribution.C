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

void setSampleName(const Int_t signalRegion0_controlRegion1, vector<string> &sampleName, vector<string> &MC_TexLabel, const Int_t mumu0_ee1 = 0) {

  // for Control region, mumu0_ee1 says if we use muon (0) or electron (1)

  if (signalRegion0_controlRegion1 == 0) {

    sampleName.push_back("GJets");
    sampleName.push_back("DYJetsToLL");
    //sampleName.push_back("QCD");
    sampleName.push_back("Diboson");
    sampleName.push_back("Top");
    sampleName.push_back("WJetsToLNu");
    sampleName.push_back("ZJetsToNuNu");
       
    MC_TexLabel.push_back("#gamma + jets");
    MC_TexLabel.push_back("Z(ll)+jets");
    //MC_TexLabel.push_back("QCD multijets");
    MC_TexLabel.push_back("Diboson");
    MC_TexLabel.push_back("t#bar{t},single t");
    MC_TexLabel.push_back("W(l#nu)+jets");
    MC_TexLabel.push_back("Z(#nu#nu)+jets");
    
  } else {

    if (mumu0_ee1 == 0) {

      sampleName.push_back("ZJetsToNuNu");
      sampleName.push_back("GJets");
      //sampleName.push_back("QCD");
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");
      sampleName.push_back("DYJetsToLL");

      MC_TexLabel.push_back("Z(#nu#nu)+jets");
      MC_TexLabel.push_back("#gamma + jets");
      //MC_TexLabel.push_back("QCD multijets");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");
      MC_TexLabel.push_back("Z(#mu#mu)+jets");

    } else {

      // sampleName.push_back("ZJetsToNuNu");
      // sampleName.push_back("GJets");
      // sampleName.push_back("QCD");
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");
      sampleName.push_back("DYJetsToLL");

      // MC_TexLabel.push_back("Z(#nu#nu)+jets");
      // MC_TexLabel.push_back("#gamma + jets");
      // MC_TexLabel.push_back("QCD multijets");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");
      MC_TexLabel.push_back("Z(ee)+jets");

    }

  }

}


void setSampleName2lepSkim(const Int_t signalRegion0_controlRegion1, vector<string> &sampleName, vector<string> &MC_TexLabel, const Int_t mumu0_ee1 = 0) {

  // for Control region, mumu0_ee1 says if we use muon (0) or electron (1)

  if (signalRegion0_controlRegion1 == 0) {

    sampleName.push_back("GJets");
    sampleName.push_back("DYJetsToLL");
    //sampleName.push_back("QCD");
    sampleName.push_back("Diboson");
    sampleName.push_back("Top");
    sampleName.push_back("WJetsToLNu");
    sampleName.push_back("ZJetsToNuNu");
       
    MC_TexLabel.push_back("#gamma + jets");
    MC_TexLabel.push_back("Z(ll)+jets");
    //MC_TexLabel.push_back("QCD multijets");
    MC_TexLabel.push_back("Diboson");
    MC_TexLabel.push_back("t#bar{t},single t");
    MC_TexLabel.push_back("W(l#nu)+jets");
    MC_TexLabel.push_back("Z(#nu#nu)+jets");
    
  } else {

    if (mumu0_ee1 == 0) {

      // sampleName.push_back("ZJetsToNuNu");
      // sampleName.push_back("GJets");
      // sampleName.push_back("QCD");
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");
      sampleName.push_back("DYJetsToLL");

      // MC_TexLabel.push_back("Z(#nu#nu)+jets");
      // MC_TexLabel.push_back("#gamma + jets");
      // MC_TexLabel.push_back("QCD multijets");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");
      MC_TexLabel.push_back("Z(#mu#mu)+jets");

    } else {

      // sampleName.push_back("ZJetsToNuNu");
      // sampleName.push_back("GJets");
      // sampleName.push_back("QCD");
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");
      sampleName.push_back("DYJetsToLL");

      // MC_TexLabel.push_back("Z(#nu#nu)+jets");
      // MC_TexLabel.push_back("#gamma + jets");
      // MC_TexLabel.push_back("QCD multijets");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");
      MC_TexLabel.push_back("Z(ee)+jets");

    }

  }

}


void setHistColor(vector<Int_t> &histColor, const Int_t nSamples) {

  Int_t colorList[] = {kCyan, kViolet, kBlue, kRed, /*kOrange+1, */kYellow, kGreen};  // the first color is for the main object. This array may contain more values than nSamples

  for (Int_t i = 0; i < nSamples; i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)

    histColor.push_back(colorList[i]);

  }

  std::reverse(histColor.begin(), histColor.end());   // reverse order in the array

}



void setDistribution(const Int_t signalRegion0_controlRegion1, const Int_t mumu0_ee1, const string var, string &hvarName, string &xAxisName) {

  if ( !(strcmp("metBin",var.c_str())) ) { 

    hvarName = "HYieldsMetBin";
    xAxisName = "#slash{E}_{T} [GeV]";

  } else if ( !(strcmp("met",var.c_str())) ) {

    hvarName = "HmetNoLepDistribution"; 
    xAxisName = "#slash{E}_{T} [GeV]";

  } else if ( !(strcmp("nvtx",var.c_str())) ) {

    hvarName = "HvtxDistribution"; 
    xAxisName = "N(vertices)";

  } else if ( !(strcmp("njets",var.c_str())) ) {

    hvarName = "HnjetsDistribution"; 
    xAxisName = "N(jets)";

  } else if ( !(strcmp("j1pt",var.c_str())) ) {

    hvarName = "Hjet1ptDistribution"; 
    xAxisName = "leading jet p_{T} [GeV]";

  } else if ( !(strcmp("j2pt",var.c_str())) ) {

    hvarName = "Hjet2ptDistribution"; 
    xAxisName = "second jet p_{T} [GeV]";

  } else if ( !(strcmp("invMass",var.c_str())) ) {

    hvarName = "HinvMass";
    if (mumu0_ee1 == 0) xAxisName = "m(#mu#mu) [GeV]";
    else if (mumu0_ee1 == 1) xAxisName = "m(ee) [GeV]";

  } else if ( !(strcmp("zpt",var.c_str())) ) {

    hvarName = "HzptDistribution";
    xAxisName = "p_{T}(Z) [GeV]";

  } else {

    cout << " Variable " << var << " not available. End of programme." << endl;
    exit(EXIT_FAILURE);

  }

}

// ================================================




void distribution(const string folderNameWithRootFiles = "", const Int_t signalRegion0_controlRegion1 = 0, const Int_t mumu0_ee1 = 0, const Int_t data0_noData1 = 0, const string var = "met", const Int_t yAxisLog_flag = 0, const Int_t MCpoissonUncertainty_flag = 0, const Double_t xAxisMin = 0, const Double_t xAxisMax = -1, const Double_t yAxisMin = 0, const Double_t yAxisMax = -1, const Int_t binDensity_flag = 0, const Int_t MCnormalizedToData_flag = 0) {

  // if signalRegion0_controlRegion1 == 0 (default), will do met distribution in the monojet signal region, else it will do the control region

  // mumu0_ee1 is for lepton flavour in CS (not used in SR)

  // data0_noData1 is to use or not a data file to compared with MC

  // yAxisLog_flag is to choose whether or not to use log scale in Y axis (default is 0, that is, normal scale)
  // xAxisMin and xAxisMax are the ranges for x Axis. Default values are used if xAxisMin > xAxisMax (otherwise user values are used)

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)
  gStyle->SetFillColor(10);

  string filenameExtension = ".root";
  // string fileDirectoryPath = "spring15_25ns_rootfiles/";
  string fileDirectoryPath = "/cmshome/ciprianim/CMSSW721/output/" + folderNameWithRootFiles + "/";

  string plotDirectoryPath = fileDirectoryPath;
  // string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/monojet/met_distribution/";
  //string plotDirectoryPath = "./distributions/";
  string plotFileExtension = ".pdf";
  string suffix;

  if (mumu0_ee1 == 0) suffix = "_mumu";
  else if (mumu0_ee1 == 1) suffix = "_ee";
  else {
    cout << "Error: mumu0_ee1 must be 0 or 1. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  TH1D* hvar = NULL;   // to get histogram from file
  string hvarName;          // name of histogram to take from file
  string xAxisName;        // name of X axis when plotting distribution. It is a tex string (with ROOT standard), e.g. "#slash{E}_{T} [GeV]" for MET

  // ===== TO BE MODIFIED =====

  // hvarName = "HmetNoLepDistribution";
  // xAxisName = "#slash{E}_{T} [GeV]";

  setDistribution(signalRegion0_controlRegion1, mumu0_ee1, var, hvarName, xAxisName);
    
  // =====================

  vector<TH1D*> hMC;
  TH1D* hdata = NULL;

  vector<string> sampleName;
  vector<string> MC_TexLabel;
  if (data0_noData1 == 1) setSampleName(signalRegion0_controlRegion1, sampleName, MC_TexLabel, mumu0_ee1);
  else setSampleName2lepSkim(signalRegion0_controlRegion1, sampleName, MC_TexLabel, mumu0_ee1);

  string data_TexLabel = "data";

  Int_t nFiles = (Int_t) sampleName.size();

  vector<Int_t> histColor;
   setHistColor(histColor, nFiles);

   string filenameBase;
   string canvasName;

   if (signalRegion0_controlRegion1 == 0) {

     filenameBase = "monojet_SR_spring15_25ns_";
     canvasName = var + "_monojetSR_";

   } else {

     canvasName = var + "_zjetsCS" + suffix;

     if (mumu0_ee1 == 0) filenameBase = "zmumujets_CS_spring15_25ns_";
     else if (mumu0_ee1 == 1) filenameBase = "zeejets_CS_spring15_25ns_";

   }
 

  vector<string> MCfileName;
  for (Int_t i = 0; i < nFiles; i++) {
    MCfileName.push_back(fileDirectoryPath + filenameBase + sampleName[i] + filenameExtension);
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

    hvar = (TH1D*)f->Get(hvarName.c_str());
    if (!hvar) {
      cout << "Error: histogram not found in file ' " << MCfileName[i] << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hMC.push_back( (TH1D*)hvar->Clone() );

  } 

  // ==== FILE NAME WILL HAVE TO BE MODIFIED, NOW IT IS JUST A TEST =====

  string datafileName = fileDirectoryPath;

  if (data0_noData1 == 0) {

    if (signalRegion0_controlRegion1 == 0) {

      datafileName += "monojet_SR_spring15_25ns_DATA.root";

    } else {

      if (mumu0_ee1 == 0) datafileName += "zmumujets_CS_spring15_25ns_DATA.root";
      else if (mumu0_ee1 == 1) datafileName += "zeejets_CS_spring15_25ns_DATA.root";

    }

  }

  // ==== opening data file ======

  if (data0_noData1 == 0) {

    cout<<"fileName : "<<datafileName<<endl;

    TFile* f = TFile::Open(datafileName.c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<datafileName<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }

    hvar = (TH1D*)f->Get(hvarName.c_str());

    if (!hvar) {
      cout << "Error: histogram not found in file ' " << datafileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hdata = (TH1D*)hvar->Clone();

  }

  // ===============================

  THStack* hMCstack = new THStack("hMCstack","");
  Double_t stackNorm = 0.0;

  for (Int_t j = 0; j < nFiles; j++) {

    for (Int_t i = 1; i <= hMC[j]->GetNbinsX(); i++) {

      if (MCpoissonUncertainty_flag == 1) {

	hMC[j]->SetBinError(i,sqrt(hMC[j]->GetBinContent(i)));
	
      }

    }

    hMC[j]->SetFillColor(histColor[j]);
    stackNorm += hMC[j]->Integral();

  }

  // loop again on MC histograms to scale them and then fill the thstack

  for (Int_t j = 0; j < nFiles; j++) {

    if (data0_noData1 == 0 && MCnormalizedToData_flag != 0) {

      Double_t dataNorm = hdata->Integral();

      if (binDensity_flag != 0) hMC[j]->Scale(dataNorm/stackNorm,"width");
      else hMC[j]->Scale(dataNorm/stackNorm);

    } else if (binDensity_flag != 0) hMC[j]->Scale(1.0,"width");  // option width divides by bin width and manages the correct error setting

    hMCstack->Add(hMC[j]);

  }

  if (data0_noData1 == 0) {

    cout << "Events in data: " << hdata->Integral() << endl;
    cout << "Events in MC   : " << ((TH1D*) hMCstack->GetStack()->Last())->Integral() << endl;

  }
  
  if (data0_noData1 == 0 && binDensity_flag != 0) hdata->Scale(1.0,"width");
  //if (data0_noData1 == 0 && MCnormalizedToData_flag != 0)

  // now here we go with the canvas

  TH1D * ratioplot = NULL; // will use it for the ratio plots

  TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpad_2 = NULL; 

  TCanvas *c;
  if (data0_noData1 == 0) c = new TCanvas(canvasName.c_str(), (var + " distribution").c_str(), 700, 700);
  else c = new TCanvas(canvasName.c_str(), (var + " distribution").c_str());
  TLegend *leg = new TLegend(0.7,0.6,0.99,0.94);  

  // if there are data, split canvas to draw the dta/MC ratio plot

  if (data0_noData1 == 0) {

    subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
    if (yAxisLog_flag) subpad_1->SetLogy();
    //subpad_1->SetBottomMargin(0);
    subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
    subpad_2->SetGridy();
    //subpad_2->SetTopMargin(0);
    subpad_2->SetBottomMargin(0.3);
    subpad_1->Draw();
    subpad_2->Draw();

    subpad_1->cd();

  } else if (yAxisLog_flag) c->SetLogy();

  
  hMCstack->Draw("HIST");
  //if (yAxisMin > 0) hMCstack->SetMinimum(yAxisMin);

  if (yAxisMin < yAxisMax) {
    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
    else c->Update();  
    hMCstack->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);
  }

  // if (data0_noData1 == 0) {

  //   if (yAxisMin > 0) hMCstack->GetYaxis()->SetRangeUser(yAxisMin, subpad_1->GetY2());
  //   else hMCstack->GetYaxis()->SetRangeUser(yAxisMin, c->GetY2());

  // }

  //hMCstack->SetMaximum(4000.0);
  TH1D* stackCopy = (TH1D*)(((TH1D*)hMCstack->GetStack()->Last())->DrawCopy("E2 SAME"));
  stackCopy->SetFillColor(kBlack);
  stackCopy->SetFillStyle(3144);

  if (data0_noData1 == 1) {    //  when using data ( == 0) the x axis will not have labels (they will only be below in the ratio plot
    hMCstack->GetXaxis()->SetTitle(xAxisName.c_str());
    hMCstack->GetXaxis()->SetTitleSize(0.06);
    hMCstack->GetXaxis()->SetTitleOffset(0.6);
  }

  if (xAxisMin < xAxisMax) {
    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
    else c->Update();  
    hMCstack->GetXaxis()->SetRangeUser(xAxisMin,xAxisMax);
  }

  if (binDensity_flag == 0) hMCstack->GetYaxis()->SetTitle("events");
  else hMCstack->GetYaxis()->SetTitle("events / GeV");
  hMCstack->GetYaxis()->SetTitleSize(0.06);
  hMCstack->GetYaxis()->SetTitleOffset(0.8);
  hMCstack->GetYaxis()->CenterTitle();
  for (Int_t j = (nFiles-1); j >= 0; j--) {
    leg->AddEntry(hMC[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }

  if (data0_noData1 == 0) {    
    hdata->SetMarkerStyle(8); // large dot
    hdata->Draw("EX0 SAME"); //X0 doesn't draw x error
    leg->AddEntry(hdata,Form("%s",data_TexLabel.c_str()),"p");
  }

  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);

  if (data0_noData1 == 0) { // if using data, draw the ratio plot

    subpad_2->cd();
    ratioplot = new TH1D(*hdata);
    ratioplot->Divide(stackCopy);
    ratioplot->SetStats(0);
    ratioplot->SetTitle("");
    ratioplot->GetXaxis()->SetLabelSize(0.10);
    ratioplot->GetXaxis()->SetTitle(xAxisName.c_str());
    ratioplot->GetXaxis()->SetTitleSize(0.14);
    ratioplot->GetXaxis()->SetTitleOffset(0.8);
    ratioplot->GetYaxis()->SetLabelSize(0.10);
    ratioplot->GetYaxis()->SetTitle("data / MC");
    ratioplot->GetYaxis()->SetTitleSize(0.15);
    ratioplot->GetYaxis()->SetTitleOffset(0.3);
    ratioplot->GetYaxis()->CenterTitle();
    ratioplot->GetYaxis()->SetRangeUser(0.5,1.5);
    ratioplot->GetYaxis()->SetNdivisions(011);
    if (xAxisMin < xAxisMax) {
      //subpad_2->Update();  // to be done after Draw() to access pad parameters such as default axis range  
      ratioplot->GetXaxis()->SetRangeUser(xAxisMin,xAxisMax);
    }
    ratioplot->SetMarkerStyle(8);  //medium dot
    ratioplot->DrawCopy("E");  

  }

  c->SaveAs( (plotDirectoryPath + c->GetName() + plotFileExtension).c_str() );

}
