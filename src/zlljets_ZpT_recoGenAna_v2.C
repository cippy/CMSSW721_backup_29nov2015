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

void zlljets_ZpT_recoGenAna_v2(const string fmumuName = "zmumujets_resoResp_noSkim_light.root", const string feeName = "zeejets_resoResp_noSkim_light.root") {

  TH1::SetDefaultSumw2();   

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/ZpTAnalysis/";
  string plotFileExtension = ".pdf";

  // index 0 of array is for muons, index 1 is for electrons, thus recoZpt[0] refers to Z->mumu ((maybe it will be better to change it)

  vector<string> fileName;
  fileName.push_back(fmumuName);
  fileName.push_back(feeName);

  vector<string> sample;
  sample.push_back("Z(#mu#mu)+jets");
  sample.push_back("Z(ee)+jets");

  vector<TH1D*> hZpTreco;
  vector<TH1D*> hZpTgen;
  vector<TH1D*> hZpTrecoGenRatio;
  vector<TH1D*> hZpTrecoGenRatio_pdf;

  Int_t nFiles = (Int_t)fileName.size();
  Int_t histColor[] = {kBlue,kOrange,kRed,kGreen};

  for(Int_t i = 0; i < nFiles; i++) {

    cout<<"fileName : "<<fileName[i]<<endl;

    TFile* f = TFile::Open(fileName[i].c_str(),"READ");
    
    if (!f || !f->IsOpen()) {

      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<fileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);

    }

    hZpTreco.push_back((TH1D*)f->Get("HZtoLLRecoPt"));
    hZpTgen.push_back((TH1D*)f->Get("HZtoLLGenPt"));
    hZpTrecoGenRatio.push_back((TH1D*)f->Get("HZtoLLPt_RecoGenRatio"));
    hZpTrecoGenRatio_pdf.push_back((TH1D*)f->Get("HZtoLLPt_RecoGenRatio_pdf"));

    if (!(hZpTreco[i] && hZpTgen[i] && hZpTrecoGenRatio[i] && hZpTrecoGenRatio_pdf[i])) {

      cout << "Error: could not get histograms from file " << fileName[i] << endl;
      exit(EXIT_FAILURE);

    }
	 
  }

  for (Int_t j = 0; j < nFiles; j++) {

    hZpTreco[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hZpTreco[j]->SetLineColor(histColor[j]);
    // no need to use TH1::Scale() between reco and gen of Z->mumu or Z->ee, because histogram are filled for the same set of events, so integral
    // should be the same in principle
    //hZpTreco[j]->SetFillColorAlpha(histColor[j], 0.80);
    hZpTgen[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hZpTgen[j]->SetLineColor(histColor[j+2]);
    hZpTrecoGenRatio[j]->SetStats(kFALSE);
    hZpTrecoGenRatio[j]->SetLineColor(histColor[j]);
    hZpTrecoGenRatio_pdf[j]->SetStats(kFALSE);
    hZpTrecoGenRatio_pdf[j]->SetLineColor(histColor[j+2]);
    

  }

  TPad *subpad = NULL;  // will use it to access specific subpad in canvas

  TCanvas *cZtomumuPtRecoGen = new TCanvas("cZtomumuPtRecoGen","");
  TLegend *legZtomumuPtRecoGen = new TLegend(0.70,0.7,0.89,0.89); 
  cZtomumuPtRecoGen->Divide(1,2);
  cZtomumuPtRecoGen->cd(1);
  subpad = (TPad*)cZtomumuPtRecoGen->GetPad(1);
  subpad->SetLogy();
  //cout <<"Drawing histograms and ratio"<<endl;
  hZpTreco[0]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  hZpTreco[0]->GetXaxis()->SetTitleSize(0.04);
  hZpTreco[0]->GetYaxis()->SetTitle("# events");
  hZpTreco[0]->GetYaxis()->SetTitleSize(0.045);
  hZpTreco[0]->GetYaxis()->CenterTitle();
  hZpTreco[0]->Draw("HE");
  hZpTgen[0]->Draw("HE SAME");
  legZtomumuPtRecoGen->AddEntry(hZpTreco[0],"Z(#mu#mu) reco","l");
  legZtomumuPtRecoGen->AddEntry(hZpTgen[0],"Z(#mu#mu) gen","l");
  gStyle->SetStatStyle(0);
  legZtomumuPtRecoGen->Draw(); 
  legZtomumuPtRecoGen->SetMargin(0.3); 
  legZtomumuPtRecoGen->SetBorderSize(0);
  cZtomumuPtRecoGen->cd(2);
  hZpTrecoGenRatio[0]->Draw("E");
  hZpTrecoGenRatio[0]->SetMarkerStyle(7);  //medium dot
  cZtomumuPtRecoGen->SaveAs( (plotDirectoryPath + cZtomumuPtRecoGen->GetName() + plotFileExtension).c_str() );

  TCanvas *cZtoeePtRecoGen = new TCanvas("cZtoeePtRecoGen","");
  TLegend *legZtoeePtRecoGen = new TLegend(0.70,0.7,0.89,0.89);
  cZtoeePtRecoGen->Divide(1,2);
  subpad = (TPad*)cZtoeePtRecoGen->GetPad(1);
  cZtoeePtRecoGen->cd(1);
  subpad->SetLogy();
  //cout <<"Drawing histograms and ratio"<<endl;
  hZpTreco[1]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  hZpTreco[1]->GetXaxis()->SetTitleSize(0.04);
  hZpTreco[1]->GetYaxis()->SetTitle("# events");
  hZpTreco[1]->GetYaxis()->SetTitleSize(0.045);
  hZpTreco[1]->GetYaxis()->CenterTitle();
  hZpTreco[1]->Draw("HE");
  hZpTgen[1]->Draw("HE SAME");
  legZtoeePtRecoGen->AddEntry(hZpTreco[1],"Z(ee) reco","l");
  legZtoeePtRecoGen->AddEntry(hZpTgen[1],"Z(ee) gen","l");
  gStyle->SetStatStyle(0);
  legZtoeePtRecoGen->Draw(); 
  legZtoeePtRecoGen->SetMargin(0.3); 
  legZtoeePtRecoGen->SetBorderSize(0);
  cZtoeePtRecoGen->cd(2);
  hZpTrecoGenRatio[1]->Draw("E");
  hZpTrecoGenRatio[1]->SetMarkerStyle(7);  //medium dot
  cZtoeePtRecoGen->SaveAs( (plotDirectoryPath + cZtoeePtRecoGen->GetName() + plotFileExtension).c_str() );

  // now I normalize histograms to same area (here is 1) to compare Z->mumu and Z->ee which have different numbers of events

  for (Int_t j = 0; j < nFiles; j++) {

    hZpTreco[j]->Scale(1./hZpTreco[j]->Integral(0,1 + hZpTreco[j]->GetNbinsX()));  // normalize to unity, use integral including underflow & overflow bin
    hZpTgen[j]->Scale(1./hZpTgen[j]->Integral(0,1 + hZpTgen[j]->GetNbinsX()));

  }

  TCanvas *cZtollPtRecoGen = new TCanvas("cZtollPtRecoGen","");
  TLegend *legZtollPtRecoGen = new TLegend(0.70,0.7,0.89,0.89); 
  cZtollPtRecoGen->SetLogy();
  //cout <<"Drawing histograms and ratio"<<endl;
  hZpTreco[0]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  hZpTreco[0]->GetXaxis()->SetTitleSize(0.04);
  hZpTreco[0]->GetYaxis()->SetTitle("# events");
  hZpTreco[0]->GetYaxis()->SetTitleSize(0.045);
  hZpTreco[0]->GetYaxis()->CenterTitle();
  hZpTreco[0]->Draw("HE");
  hZpTgen[0]->Draw("HE SAME");
  hZpTreco[1]->Draw("HE SAME");
  hZpTgen[1]->Draw("HE SAME");
  legZtollPtRecoGen->AddEntry(hZpTreco[0],"Z(#mu#mu) reco","l");
  legZtollPtRecoGen->AddEntry(hZpTgen[0],"Z(#mu#mu) gen","l");
  legZtollPtRecoGen->AddEntry(hZpTreco[1],"Z(ee) reco","l");
  legZtollPtRecoGen->AddEntry(hZpTgen[1],"Z(ee) gen","l");
  gStyle->SetStatStyle(0);
  legZtollPtRecoGen->Draw(); 
  legZtollPtRecoGen->SetMargin(0.3); 
  legZtollPtRecoGen->SetBorderSize(0);
  cZtollPtRecoGen->SaveAs( (plotDirectoryPath + cZtollPtRecoGen->GetName() + plotFileExtension).c_str() );

  Int_t bins;
  Double_t firstEdge;
  Double_t lastEdge;
  TH1D *hZpTMuEleRecoRatio;
  TCanvas *cZpTMuEleRecoRatio;

  if( hZpTreco[0]->GetNbinsX() == hZpTreco[1]->GetNbinsX() ) {

    bins = hZpTreco[0]->GetNbinsX();
    firstEdge = hZpTreco[0]->GetXaxis()->GetBinLowEdge(1);
    lastEdge = hZpTreco[0]->GetXaxis()->GetBinUpEdge(bins);
    hZpTMuEleRecoRatio = new TH1D("hZpTMuEleRecoRatio","",bins,firstEdge,lastEdge);
    

    // cout << "hZpTreco[0]->GetNbinsX() : "<< hZpTreco[0]->GetNbinsX() << endl;
    // cout << "hZpTreco[1]->GetNbinsX() : "<< hZpTreco[1]->GetNbinsX() << endl;
    // cout << "hZpTMuEleGenRatio->GetNbinsX() : "<< hZpTMuEleGenRatio->GetNbinsX() << endl;

    if (!hZpTMuEleRecoRatio->Divide(hZpTreco[0],hZpTreco[1])) cout << " Error in hZpTMuEleRecoRatio->Divide(hZpTreco[0],hZpTreco[1])" << endl;
    hZpTMuEleRecoRatio->SetStats(kFALSE);
   
    cZpTMuEleRecoRatio = new TCanvas("cZpTMuEleRecoRatio","");
    hZpTMuEleRecoRatio->GetXaxis()->SetTitle("Z_{pT}[GeV]");
    hZpTMuEleRecoRatio->GetXaxis()->SetTitleSize(0.04);
    hZpTMuEleRecoRatio->GetYaxis()->SetTitle("reco pT_{Z(#mu#mu)} / pT_{Z(ee)} ");
    hZpTMuEleRecoRatio->GetYaxis()->SetTitleSize(0.045);
    hZpTMuEleRecoRatio->GetYaxis()->SetRangeUser(0.5,1.5);
    hZpTMuEleRecoRatio->GetYaxis()->CenterTitle();
    hZpTMuEleRecoRatio->Draw("E");
    gStyle->SetStatStyle(0);
    
  }

  TH1D *hZpTMuEleGenRatio;
  TCanvas *cZpTMuEleGenRatio;

  if (hZpTgen[0]->GetNbinsX() == hZpTgen[1]->GetNbinsX()) {

    bins = hZpTgen[0]->GetNbinsX();
    firstEdge = hZpTgen[0]->GetXaxis()->GetBinLowEdge(1);
    lastEdge = hZpTgen[0]->GetXaxis()->GetBinUpEdge(bins);
    hZpTMuEleGenRatio = new TH1D("hZpTMuEleGenRatio","",bins,firstEdge,lastEdge);

    if (!hZpTMuEleGenRatio->Divide(hZpTgen[0],hZpTgen[1])) cout << " Error in hZpTMuEleGenRatio->Divide(hZpTgen[0],hZpTgen[1])" << endl;
    hZpTMuEleGenRatio->SetStats(kFALSE);

    cZpTMuEleGenRatio = new TCanvas("cZpTMuEleGenRatio","");
    hZpTMuEleGenRatio->GetXaxis()->SetTitle("Z_{pT}[GeV]");
    hZpTMuEleGenRatio->GetXaxis()->SetTitleSize(0.04);
    hZpTMuEleGenRatio->GetYaxis()->SetTitle("gen pT_{Z(#mu#mu)} / pT_{Z(ee)} ");
    hZpTMuEleGenRatio->GetYaxis()->SetTitleSize(0.045);
    hZpTMuEleGenRatio->GetYaxis()->SetRangeUser(0.5,1.5);
    hZpTMuEleGenRatio->GetYaxis()->CenterTitle();
    hZpTMuEleGenRatio->Draw("E");
    gStyle->SetStatStyle(0);

  }


}
