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
  MC_TexLabel.push_back("TT+Jets MC");
  MC_TexLabel.push_back("Z(#mu#mu)+jets MC");

  string data_TexLabel = "data"; 

  Int_t default_YaxisRange = 1;

  string DataFileName = "zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100_DATA.root";

  vector<TH1D*> hMCinvMass;
  vector<TH1D*> hMCmetNoLep;
  vector<TH1D*> hMCjet1Pt;
  vector<TH1D*> hMCuParMinusZpt;
  vector<TH1D*> hMCuPerp;
  vector<TH1D*> hMCnVertices;
  vector<TH1D*> hMCzPtSpectrum;

  TH1D* hDatainvMass;
  TH1D* hDatametNoLep;
  TH1D* hDatajet1Pt;
  TH1D* hDatauParMinusZpt;
  TH1D* hDatauPerp;
  TH1D* hDatanVertices;
  TH1D* hDatazPtSpectrum;

  vector<string> MCfileName;
  MCfileName.push_back("ttjets_mumu_resoResp_spring15_50ns_GLTcutZmass80to100.root");
  MCfileName.push_back("zmumujets_resoResp_spring15_50ns_GLTcutZmass80to100.root");  

  Int_t nFiles = (Int_t)MCfileName.size();
  Int_t histColor[] = {kRed,kCyan};

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

  

  THStack* hstack_invMass = new THStack("hstack_invMass",""); 
  THStack* hstack_uPerp = new THStack("hstack_uPerp",""); 
  THStack* hstack_uParMinusZpt = new THStack("hstack_uParMinusZpt",""); 
  
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
    }

    hMCinvMass[j]->SetFillColor(histColor[j]);
    hstack_invMass->Add(hMCinvMass[j]);
    hMCuPerp[j]->SetFillColor(histColor[j]);
    hstack_uPerp->Add(hMCuPerp[j]);
    hMCuParMinusZpt[j]->SetFillColor(histColor[j]);
    hstack_uParMinusZpt->Add(hMCuParMinusZpt[j]);

  }

  TCanvas *cInvMass = new TCanvas("invMass","");
  TLegend *legInvMass = new TLegend(0.7,0.7,0.89,0.89);
  //cInvMass->SetLogy();
  hstack_invMass->Draw("HE");
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
  hstack_uParMinusZpt->Draw("E SAME");
  hDatauParMinusZpt->SetMarkerStyle(8); // large dot
  hDatauParMinusZpt->Draw("EX0 SAME"); //X0 doesn't draw x error
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
  hstack_uPerp->Draw("E SAME");
  hDatauPerp->SetMarkerStyle(8); // large dot
  hDatauPerp->Draw("EX0 SAME"); //X0 doesn't draw x error
  hDatauPerp->GetYaxis()->SetRangeUser(0.3,3050.0);
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
			     
