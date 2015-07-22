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
#include <TPaveLabel.h>
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

void zlljets_ZpT_recoGenAna(const string fmumuName = "zmumujets_resoResp_noSkim_light.root", const string feeName = "zeejets_resoResp_noSkim_light.root", const char* suffix = "") {

  TH1::SetDefaultSumw2();   

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/ZpTAnalysis/";
  //string plotDirectoryPath = "./tmpPlots/";  used when working in my pc
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

  Double_t zptStart = 400.0;
  Double_t zptEnd = 800.0;
  TH1D *HZtoLLPt_RecoGenRatio_pdf_inRange[2];

  for(Int_t i = 0; i < 2; i++) {
	
    HZtoLLPt_RecoGenRatio_pdf_inRange[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_pdf_%d_Zpt%2.0lfTo%2.0lf",i+1,zptStart,zptEnd),"",100,0.5,1.5);

  }

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

    if ( !(hZpTreco[i] && hZpTgen[i] && hZpTrecoGenRatio[i] && hZpTrecoGenRatio_pdf[i] ) ) {

      cout << "Error: could not get histograms from file " << fileName[i] << endl;
      exit(EXIT_FAILURE);

    }
	 
  }

  Int_t startingBin = hZpTrecoGenRatio[0]->FindBin(zptStart) ;  // will consider ratio pdf for Zpt > 600 GeV  (they should be equal for Z->ee and Z->nunu)
  Int_t endingBin = hZpTrecoGenRatio[0]->FindBin(zptEnd) ;     // will consider ratio pdf for Zpt < 1000 GeV  (they should be equal for Z->ee and Z->nunu)
  cout << "ZpT bins from "<<startingBin<<" to "<<endingBin<< endl;

  for (Int_t j = 0; j < 2; j++) {

    for (Int_t i = startingBin; i < endingBin; i++) {

      HZtoLLPt_RecoGenRatio_pdf_inRange[j]->Fill(hZpTrecoGenRatio[j]->GetBinContent(i));

    }

    HZtoLLPt_RecoGenRatio_pdf_inRange[j]->Scale(1./HZtoLLPt_RecoGenRatio_pdf_inRange[j]->Integral(0,1 + HZtoLLPt_RecoGenRatio_pdf_inRange[j]->GetNbinsX()));

  }   // note that this distribution will be almost empty because it's the distribution of the mean of various ZpT points.

  //cout << " ****************  check ************** " << endl;

  for (Int_t j = 0; j < nFiles; j++) {

    hZpTreco[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hZpTreco[j]->SetLineColor(histColor[j]);
    // no need to use TH1::Scale() between reco and gen of Z->mumu or Z->ee, because histogram are filled for the same set of events, so integral
    // should be the same in principle
    //hZpTreco[j]->SetFillColorAlpha(histColor[j], 0.80);
    hZpTgen[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hZpTgen[j]->SetLineColor(histColor[j+2]);
    hZpTrecoGenRatio[j]->SetStats(kFALSE);
    //hZpTrecoGenRatio[j]->SetLineColor(histColor[j]);
    hZpTrecoGenRatio_pdf[j]->SetStats(kFALSE);
    hZpTrecoGenRatio_pdf[j]->SetLineColor(histColor[j+2]);
    HZtoLLPt_RecoGenRatio_pdf_inRange[j]->SetStats(kFALSE);
    HZtoLLPt_RecoGenRatio_pdf_inRange[j]->SetLineColor(histColor[j*2]);
    

  }

  TPad *subpad1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpad2 = NULL; 

  TCanvas *cZtomumuPtRecoGen = new TCanvas("cZtomumuPtRecoGen","",700,700);
  TLegend *legZtomumuPtRecoGen = new TLegend(0.70,0.7,0.89,0.89); 
  subpad1 = new TPad("pad1","",0.0,0.36,1.0,1.0);
  subpad1->SetLogy();
  subpad1->SetBottomMargin(0);
  subpad2 = new TPad("pad2","",0.0,0.0,1.0,0.36);
  subpad2->SetGrid();
  subpad2->SetTopMargin(0);
  subpad2->SetBottomMargin(0.2);
  subpad1->Draw();
  subpad2->Draw();
  //cout <<"Drawing histograms and ratio"<<endl;
  subpad1->cd();
  // hZpTreco[0]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  // hZpTreco[0]->GetXaxis()->SetTitleSize(0.04);
  //hZpTreco[0]->GetXaxis()->SetLabelSize(0.04);
  hZpTreco[0]->GetYaxis()->SetTitle("# events");
  hZpTreco[0]->GetYaxis()->SetTitleSize(0.04);
  hZpTreco[0]->GetYaxis()->CenterTitle();
  //hZpTreco[0]->GetYaxis()->SetLabelSize(0.04);
  hZpTreco[0]->Draw("HE");
  hZpTgen[0]->Draw("HE SAME");
  legZtomumuPtRecoGen->AddEntry(hZpTreco[0],"Z(#mu#mu) reco","l");
  legZtomumuPtRecoGen->AddEntry(hZpTgen[0],"Z(#mu#mu) gen","l");
  gStyle->SetStatStyle(0);
  legZtomumuPtRecoGen->Draw(); 
  legZtomumuPtRecoGen->SetMargin(0.3); 
  legZtomumuPtRecoGen->SetBorderSize(0);
  subpad2->cd();
  hZpTrecoGenRatio[0]->GetXaxis()->SetLabelSize(0.08);
  hZpTrecoGenRatio[0]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  hZpTrecoGenRatio[0]->GetXaxis()->SetTitleSize(0.07);
  hZpTrecoGenRatio[0]->GetXaxis()->SetTitleOffset(1.2);
  hZpTrecoGenRatio[0]->GetYaxis()->SetLabelSize(0.08);
  hZpTrecoGenRatio[0]->GetYaxis()->SetTitle("reco/gen ZpT");
  hZpTrecoGenRatio[0]->GetYaxis()->SetTitleSize(0.1);
  hZpTrecoGenRatio[0]->GetYaxis()->CenterTitle();
  hZpTrecoGenRatio[0]->Draw("E");
  hZpTrecoGenRatio[0]->SetMarkerStyle(7);  //medium dot
  cZtomumuPtRecoGen->SaveAs( (plotDirectoryPath + cZtomumuPtRecoGen->GetName() + suffix + plotFileExtension).c_str() );

  TPad *subpad1_bis = NULL;
  TPad *subpad2_bis = NULL;

  TCanvas *cZtoeePtRecoGen = new TCanvas("cZtoeePtRecoGen","",700,700);
  TLegend *legZtoeePtRecoGen = new TLegend(0.70,0.7,0.89,0.89); 
  subpad1_bis = new TPad("pad1","",0.0,0.36,1.0,1.0);
  subpad1_bis->SetLogy();
  subpad1_bis->SetBottomMargin(0);
  subpad2_bis = new TPad("pad2","",0.0,0.0,1.0,0.36);
  subpad2_bis->SetGrid();
  subpad2_bis->SetTopMargin(0);
  subpad2_bis->SetBottomMargin(0.2);
  subpad1_bis->Draw();
  subpad2_bis->Draw();
  //cout <<"Drawing histograms and ratio"<<endl;
  subpad1_bis->cd();
  // hZpTreco[1]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  // hZpTreco[1]->GetXaxis()->SetTitleSize(0.04);
  //hZpTreco[1]->GetXaxis()->SetLabelSize(0.04);
  hZpTreco[1]->GetYaxis()->SetTitle("# events");
  hZpTreco[1]->GetYaxis()->SetTitleSize(0.04);
  hZpTreco[1]->GetYaxis()->CenterTitle();
  //hZpTreco[1]->GetYaxis()->SetLabelSize(0.04);
  hZpTreco[1]->Draw("HE");
  hZpTgen[1]->Draw("HE SAME");
  legZtoeePtRecoGen->AddEntry(hZpTreco[0],"Z(ee) reco","l");
  legZtoeePtRecoGen->AddEntry(hZpTgen[0],"Z(ee) gen","l");
  gStyle->SetStatStyle(0);
  legZtoeePtRecoGen->Draw(); 
  legZtoeePtRecoGen->SetMargin(0.3); 
  legZtoeePtRecoGen->SetBorderSize(0);
  subpad2_bis->cd();
  hZpTrecoGenRatio[1]->GetXaxis()->SetLabelSize(0.08);
  hZpTrecoGenRatio[1]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  hZpTrecoGenRatio[1]->GetXaxis()->SetTitleSize(0.07);
  hZpTrecoGenRatio[1]->GetXaxis()->SetTitleOffset(1.2);
  hZpTrecoGenRatio[1]->GetYaxis()->SetLabelSize(0.08);
  hZpTrecoGenRatio[1]->GetYaxis()->SetTitle("reco/gen ZpT");
  hZpTrecoGenRatio[1]->GetYaxis()->SetTitleSize(0.1);
  hZpTrecoGenRatio[1]->GetYaxis()->CenterTitle();
  hZpTrecoGenRatio[1]->Draw("E");
  hZpTrecoGenRatio[1]->SetMarkerStyle(7);  //medium dot
  cZtoeePtRecoGen->SaveAs( (plotDirectoryPath + cZtoeePtRecoGen->GetName() + suffix + plotFileExtension).c_str() );

  // TCanvas *cZtoeePtRecoGen = new TCanvas("cZtoeePtRecoGen","",700,700);
  // TLegend *legZtoeePtRecoGen = new TLegend(0.70,0.7,0.89,0.89);
  // cZtoeePtRecoGen->Divide(1,2,0,0);
  // cZtoeePtRecoGen->cd(1);
  // subpad1 = (TPad*)cZtoeePtRecoGen->GetPad(1);
  // subpad1->SetPad(0.0,0.36,0.98,0.99);
  // subpad1->SetLogy();
  // //cout <<"Drawing histograms and ratio"<<endl;
  // // hZpTreco[1]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  // // hZpTreco[1]->GetXaxis()->SetTitleSize(0.04);
  // hZpTreco[1]->GetYaxis()->SetTitle("# events");
  // hZpTreco[1]->GetYaxis()->SetTitleSize(0.045);
  // hZpTreco[1]->GetYaxis()->CenterTitle();
  // //hZpTreco[1]->GetYaxis()->SetLabelSize(0.1);
  // hZpTreco[1]->Draw("HE");
  // hZpTgen[1]->Draw("HE SAME");
  // legZtoeePtRecoGen->AddEntry(hZpTreco[1],"Z(ee) reco","l");
  // legZtoeePtRecoGen->AddEntry(hZpTgen[1],"Z(ee) gen","l");
  // gStyle->SetStatStyle(0);
  // legZtoeePtRecoGen->Draw(); 
  // legZtoeePtRecoGen->SetMargin(0.3); 
  // legZtoeePtRecoGen->SetBorderSize(0);
  // cZtoeePtRecoGen->cd(2);
  // subpad2 = (TPad*)cZtoeePtRecoGen->GetPad(2);
  // subpad2->SetPad(0.0,0.0,0.98,0.36);
  // subpad2->SetGrid();
  // hZpTrecoGenRatio[1]->GetXaxis()->SetLabelSize(0.08);
  // hZpTrecoGenRatio[1]->GetXaxis()->SetTitle("Z_{pT}[GeV]");
  // hZpTrecoGenRatio[1]->GetXaxis()->SetTitleSize(0.06);
  // hZpTrecoGenRatio[1]->GetXaxis()->SetTitleOffset(1.2);
  // hZpTrecoGenRatio[1]->GetYaxis()->SetLabelSize(0.08);
  // hZpTrecoGenRatio[1]->Draw("E");
  // hZpTrecoGenRatio[1]->SetMarkerStyle(7);  //medium dot
  // cZtoeePtRecoGen->SaveAs( (plotDirectoryPath + cZtoeePtRecoGen->GetName() + suffix + plotFileExtension).c_str() );

  // now I normalize histograms to same area (here is 1) to compare Z->mumu and Z->ee which have different numbers of events

  for (Int_t j = 0; j < nFiles; j++) {

    hZpTreco[j]->Scale(1./hZpTreco[j]->Integral(0,1 + hZpTreco[j]->GetNbinsX()));  // normalize to unity, use integral including underflow & overflow bin
    hZpTgen[j]->Scale(1./hZpTgen[j]->Integral(0,1 + hZpTgen[j]->GetNbinsX()));
    hZpTrecoGenRatio_pdf[j]->Scale(1./hZpTrecoGenRatio_pdf[j]->Integral(0,1 + hZpTrecoGenRatio_pdf[j]->GetNbinsX()));

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
  cZtollPtRecoGen->SaveAs( (plotDirectoryPath + cZtollPtRecoGen->GetName() + suffix + plotFileExtension).c_str() );

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
    cZpTMuEleRecoRatio->SaveAs( (plotDirectoryPath + cZpTMuEleRecoRatio->GetName() + suffix + plotFileExtension).c_str() );
    
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
    cZpTMuEleGenRatio->SaveAs( (plotDirectoryPath + cZpTMuEleGenRatio->GetName() + suffix + plotFileExtension).c_str() );
 
  }


  TCanvas *cZpTRecoGenRatio_pdf = new TCanvas("cZpTRecoGenRatio_pdf","");
  TLegend *legZpTRecoGenRatio_pdf = new TLegend(0.78,0.78,0.89,0.89);
  hZpTrecoGenRatio_pdf[0]->GetXaxis()->SetTitle("reco Z_{pT}/gen Z_{pT}");
  hZpTrecoGenRatio_pdf[0]->GetXaxis()->SetTitleSize(0.04);
  hZpTrecoGenRatio_pdf[0]->GetXaxis()->SetRangeUser(0.7,1.3);
  hZpTrecoGenRatio_pdf[0]->GetYaxis()->SetTitle("a.u. ");
  hZpTrecoGenRatio_pdf[0]->GetYaxis()->SetTitleSize(0.045);
  hZpTrecoGenRatio_pdf[0]->GetYaxis()->SetRangeUser(0.0,0.4);
  hZpTrecoGenRatio_pdf[0]->GetYaxis()->CenterTitle();
  hZpTrecoGenRatio_pdf[0]->Draw("HE");
  hZpTrecoGenRatio_pdf[1]->Draw("HE SAME");
  legZpTRecoGenRatio_pdf->AddEntry(hZpTrecoGenRatio_pdf[0],"Z(#mu#mu)","l");
  legZpTRecoGenRatio_pdf->AddEntry(hZpTrecoGenRatio_pdf[1],"Z(ee)","l");
  gStyle->SetStatStyle(0);
  legZpTRecoGenRatio_pdf->Draw(); 
  legZpTRecoGenRatio_pdf->SetMargin(0.3); 
  legZpTRecoGenRatio_pdf->SetBorderSize(0);
  cZpTRecoGenRatio_pdf->SaveAs( (plotDirectoryPath + cZpTRecoGenRatio_pdf->GetName() + suffix + plotFileExtension).c_str() );

  TCanvas *cZpTRecoGenRatio_pdf_inRange = new TCanvas(Form("ZpTRecoGenRatio_pdf_Zpt%2.0lfTo%2.0lf",zptStart,zptEnd),"");
  TLegend *legZpTRecoGenRatio_pdf_inRange = new TLegend(0.78,0.78,0.89,0.89);
  TPaveLabel *title = new TPaveLabel(0.1,0.94,0.9,0.98,Form("Zpt recoGen ratio: %2.0lf < ZpT < %2.0lf",zptStart,zptEnd));
  title->Draw();
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetXaxis()->SetTitle("reco Z_{pT}/gen Z_{pT}");
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetXaxis()->SetTitleSize(0.04);
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetXaxis()->SetRangeUser(0.7,1.3);
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetYaxis()->SetTitle("a.u. ");
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetYaxis()->SetTitleSize(0.045);
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetYaxis()->SetRangeUser(0.0,0.4);
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->GetYaxis()->CenterTitle();
  HZtoLLPt_RecoGenRatio_pdf_inRange[0]->Draw("HE");
  HZtoLLPt_RecoGenRatio_pdf_inRange[1]->Draw("HE SAME");
  legZpTRecoGenRatio_pdf_inRange->AddEntry(HZtoLLPt_RecoGenRatio_pdf_inRange[0],"Z(#mu#mu)","l");
  legZpTRecoGenRatio_pdf_inRange->AddEntry(HZtoLLPt_RecoGenRatio_pdf_inRange[1],"Z(ee)","l");
  gStyle->SetStatStyle(0);
  legZpTRecoGenRatio_pdf_inRange->Draw(); 
  legZpTRecoGenRatio_pdf_inRange->SetMargin(0.3); 
  legZpTRecoGenRatio_pdf_inRange->SetBorderSize(0);
  cZpTRecoGenRatio_pdf_inRange->SaveAs( (plotDirectoryPath + cZpTRecoGenRatio_pdf_inRange->GetName() + suffix + plotFileExtension).c_str() );

}
