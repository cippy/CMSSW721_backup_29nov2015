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

void met_NoMu_NoEle_comparison() 
{

  TH1::SetDefaultSumw2();   

  vector<TH1D*> hMetShape;
  
  vector<TH1D*> hAxe;                     // only for Z(mumu) and Z(ee)
  vector<TH1D*> hZvvYieldsEstimate;       // only for Z(mumu) and Z(ee)

  TH1D* hZvvYields = NULL;

  vector<string> fileName;
  fileName.push_back("zmumujets_Axe_noSkim.root");
  fileName.push_back("zeejets_Axe_noSkim.root");
  fileName.push_back("znunujetsAna.root");
  vector<string> sample;
  sample.push_back("Z(#mu#mu)+jets");
  sample.push_back("Z(ee)+jets");
  sample.push_back("Z(#nu#nu)+jets");
  Int_t nFiles = (Int_t)fileName.size();
  Int_t histColor[] = {kBlue,kRed,kGreen};

  for(Int_t i = 0; i < nFiles; i++) {

    cout<<"fileName : "<<fileName[i]<<endl;

    TFile* f = TFile::Open(fileName[i].c_str(),"READ");
    
    if (!f || !f->IsOpen()) {

      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<fileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);

    }
    
    if (fileName[i] == "znunujetsAna.root") {

      hMetShape.push_back((TH1D*)f->Get("HzlljetsYieldsMetBin"));
      hZvvYields = (TH1D*)hMetShape[2]->Clone("hZvvYields");
      hZvvYields->SetLineColor(histColor[2]);

    } else {

      hMetShape.push_back((TH1D*)f->Get("HzlljetsYieldsMetBinGenLep"));
      hAxe.push_back((TH1D*)f->Get("Hacceff"));
      hZvvYieldsEstimate.push_back((TH1D*)f->Get("HzvvEstimate"));

    }
	 
  }

  for (Int_t j = 0; j < 2; j++) {
    // correct shape dividing by Axe bin by bin
    if ( !hMetShape[j]->Divide(hAxe[j]) ) cout << "Error in hMetShape[" << j << "]->Divide(hAxe[" <<j << "])" <<endl; 
    hZvvYieldsEstimate[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hZvvYieldsEstimate[j]->SetLineColor(histColor[j]);
  }

  for (Int_t j = 0; j < nFiles; j++) {

    hMetShape[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hMetShape[j]->SetLineColor(histColor[j]);
    hMetShape[j]->Scale(1./hMetShape[j]->Integral(0,1 + hMetShape[j]->GetNbinsX()));  // normalize to unity, use integral including underflow & overflow bin
    //hMetShape[j]->SetFillColorAlpha(histColor[j], 0.80);

  }

  // check if 2 histogram have same number of bins (they should by definition, but better to check)
  if (hMetShape[0]->GetNbinsX() != hMetShape[1]->GetNbinsX()) {
    cout << " Warning: the 2 histograms have different number of bins. The lower number will be used." << endl;
  }

  // in case the number of bins is different, use the lower one (hoping at least that the first low edge is the same for both)
  Int_t lessBinHistoIndex = (hMetShape[0]->GetNbinsX() <= hMetShape[1]->GetNbinsX()) ? 0 : 1; 
  Int_t nMetBins = hMetShape[lessBinHistoIndex]->GetNbinsX();
  Double_t metBinEdges[nMetBins + 1];
  hMetShape[lessBinHistoIndex]->GetXaxis()->GetLowEdge(metBinEdges); // this doesn't include the last bin up edge (the overflow low edge)
  metBinEdges[nMetBins] = hMetShape[lessBinHistoIndex]->GetXaxis()->GetBinLowEdge(nMetBins + 1); // get low edge of overflow bin

  for (Int_t i = 0; i <= nMetBins; i++) {
    cout<<" metBinEdges["<<i<<"] = "<<metBinEdges[i]<<endl;
  }

  TCanvas *c = new TCanvas("c","");
  TLegend *leg = new TLegend(0.7,0.7,0.89,0.89); 
  c->SetLogy();
  //cout <<"Drawing histograms and ratio"<<endl;
  hMetShape[0]->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hMetShape[0]->GetXaxis()->SetTitleSize(0.04);
  hMetShape[0]->GetYaxis()->SetTitle("a.u.");
  hMetShape[0]->GetYaxis()->SetTitleSize(0.045);
  hMetShape[0]->GetYaxis()->CenterTitle();
  hMetShape[0]->Draw("HE");
  hMetShape[1]->Draw("HE SAME");
  hMetShape[2]->Draw("HE SAME");
  for (Int_t j = 0; j < nFiles; j++) {
    leg->AddEntry(hMetShape[j],Form("%s",sample[j].c_str()),"l");
  }
  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  //c->SaveAs("met_NoMu_NoEle_NuNu_comparison.root");  

  TCanvas *cZvvYieldsEstimate = new TCanvas("cZvvYieldsEstimate","");
  TLegend *legZvvYieldsEstimate = new TLegend(0.65,0.7,0.89,0.89); 
  cZvvYieldsEstimate->SetLogy();
  //cout <<"Drawing histograms and ratio"<<endl;
  hZvvYieldsEstimate[0]->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hZvvYieldsEstimate[0]->GetXaxis()->SetTitleSize(0.04);
  hZvvYieldsEstimate[0]->GetYaxis()->SetTitle("# events");
  hZvvYieldsEstimate[0]->GetYaxis()->SetTitleSize(0.045);
  hZvvYieldsEstimate[0]->GetYaxis()->CenterTitle();
  hZvvYieldsEstimate[0]->Draw("HE");
  hZvvYieldsEstimate[1]->Draw("HE SAME");
  hZvvYields->Draw("HE SAME");
  for (Int_t j = 0; j < 2; j++) {
    legZvvYieldsEstimate->AddEntry(hZvvYieldsEstimate[j],Form("%s",sample[j].c_str()),"l");
  }
  legZvvYieldsEstimate->AddEntry(hZvvYields,Form("%s",sample[2].c_str()),"l");
  gStyle->SetStatStyle(0);
  legZvvYieldsEstimate->Draw(); 
  legZvvYieldsEstimate->SetHeader("Z(#nu#nu) yields estimate"); 
  legZvvYieldsEstimate->SetMargin(0.3); 
  legZvvYieldsEstimate->SetBorderSize(0);
  //cZvvYieldsEstimate->SaveAs("met_NoMu_NoEle_NuNu_comparison.root");  

  TCanvas *cRatioMuEle = new TCanvas("cRatioMuEle","");
  TH1D *hRatioMuEle = new TH1D("hRatioMuEle","",nMetBins,metBinEdges);
  hRatioMuEle->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hRatioMuEle->GetXaxis()->SetTitleSize(0.04);
  hRatioMuEle->GetYaxis()->SetTitle("MetNoMu / MetNoEle");
  hRatioMuEle->GetYaxis()->SetTitleSize(0.045);
  // use Divide method
  hRatioMuEle->Divide(hMetShape[0],hMetShape[1]);
  hRatioMuEle->SetStats(kFALSE);
  hRatioMuEle->Draw("E");
  //cRatioMuEle->SaveAs("met_NoMu_NoEle_ratio.root");

  TCanvas *cRatioLepNu = new TCanvas("cRatioLepNu","");
  TLegend *legRatioLepNu = new TLegend(0.78,0.18,0.89,0.29);
  TH1D *hratioMuNu = new TH1D("hratioMuNu","",nMetBins,metBinEdges);
  hratioMuNu->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hratioMuNu->GetXaxis()->SetTitleSize(0.04);
  hratioMuNu->GetYaxis()->SetTitle("MetNoL(Z->ll) / MetNoMu(Z->#nu#nu)");
  hratioMuNu->GetYaxis()->SetTitleSize(0.045);
  // use Divide method
  hratioMuNu->Divide(hMetShape[0],hMetShape[2]);
  hratioMuNu->SetStats(kFALSE);
  hratioMuNu->SetLineColor(kBlue);
  hratioMuNu->GetYaxis()->SetRangeUser(0.3,1.7);
  hratioMuNu->Draw("E");
  TH1D *hratioEleNu = new TH1D("hratioEleNu","",nMetBins,metBinEdges);
  hratioEleNu->Divide(hMetShape[1],hMetShape[2]);
  hratioEleNu->SetStats(kFALSE);
  hratioEleNu->SetLineColor(kRed);
  hratioEleNu->Draw("E SAME");
  legRatioLepNu->AddEntry(hratioMuNu,"#mu / #nu","l");
  legRatioLepNu->AddEntry(hratioEleNu,"e / #nu","l");
  legRatioLepNu->Draw(); 
  legRatioLepNu->SetMargin(0.3); 
  legRatioLepNu->SetBorderSize(0);
  //cRatioLepNu->SaveAs("met_NoLep_Nu_ratio.root");

  TCanvas *cRatioZvvEst_LepNu = new TCanvas("cRatioZvvEst_LepNu","");
  TLegend *legRatioZvvEst_LepNu = new TLegend(0.78,0.18,0.89,0.29);
  TH1D *hratioZvvEst_MuNu = new TH1D("hratioZvvEst_MuNu","",nMetBins,metBinEdges);
  hratioZvvEst_MuNu->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hratioZvvEst_MuNu->GetXaxis()->SetTitleSize(0.04);
  hratioZvvEst_MuNu->GetYaxis()->SetTitle("yields ratio  Z(#nu#nu)_{CS} / Z(#nu#nu)_{MC}");
  hratioZvvEst_MuNu->GetYaxis()->SetTitleSize(0.045);
  // use Divide method
  hratioZvvEst_MuNu->Divide(hZvvYieldsEstimate[0],hZvvYields);
  hratioZvvEst_MuNu->SetStats(kFALSE);
  hratioZvvEst_MuNu->SetLineColor(kBlue);
  hratioZvvEst_MuNu->GetYaxis()->SetRangeUser(0.3,1.7);
  hratioZvvEst_MuNu->Draw("E");
  TH1D *hratioZvvEst_EleNu = new TH1D("hratioZvvEst_EleNu","",nMetBins,metBinEdges);
  hratioZvvEst_EleNu->Divide(hZvvYieldsEstimate[1],hZvvYields);
  hratioZvvEst_EleNu->SetStats(kFALSE);
  hratioZvvEst_EleNu->SetLineColor(kRed);
  hratioZvvEst_EleNu->Draw("E SAME");
  legRatioZvvEst_LepNu->AddEntry(hratioZvvEst_MuNu,"#mu / #nu","l");
  legRatioZvvEst_LepNu->AddEntry(hratioZvvEst_EleNu,"e / #nu","l");
  legRatioZvvEst_LepNu->Draw(); 
  legRatioZvvEst_LepNu->SetMargin(0.3); 
  legRatioZvvEst_LepNu->SetBorderSize(0);

}


void met_NoMu_NoEle_comparison_binDensity() 
{

  TH1::SetDefaultSumw2();   

  vector<TH1D*> harray;
  Int_t total = 0;
  vector<string> fileName;
  fileName.push_back("histZmumujetsAnaYieldsMetBinGenMu.root");
  fileName.push_back("histZeejetsAnaYieldsMetBinGenEle.root");
  fileName.push_back("histZnunujetsAnaYieldsMetBin.root");
  vector<string> sample;
  sample.push_back("Z(#mu#mu)+jets");
  sample.push_back("Z(ee)+jets");
  sample.push_back("Z(#nu#nu)+jets");
  Int_t nFiles = (Int_t)fileName.size();
  Int_t histColor[] = {kBlue,kRed,kGreen};

  // loops over all files and gets histograms
  for(Int_t i = 0; i < nFiles; i++) {
    cout<<"fileName : "<<fileName[i]<<endl;
    TKey *key;
    TFile* f = TFile::Open(fileName[i].c_str(),"READ");
    if (!f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<fileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    
    TIter next((TList *)f->GetListOfKeys());
    while (key = (TKey *)next()) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (cl->InheritsFrom("TH1")) {
	// the following line is not needed if you only want
	// to count the histograms
	TH1 *h = (TH1 *)key->ReadObj();
	//cout << "Histo found: " << h->GetName() << " - " << h->GetTitle() << endl;
	harray.push_back((TH1D*)f->Get(h->GetName()));
	total++;
      }
    }
  }
  cout << "Found " << total << " Histograms" << endl;

  TFile* f1 = TFile::Open("hist_mumu_AccEff_metBin.root","READ");
    if (!f1->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file.\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    TH1D* hmumuAcceff = ((TH1D*)f1->Get("Hacceff"));
    if (!hmumuAcceff) {
      cout << "Error: histo not read" << endl;
      exit(EXIT_FAILURE);
    }

    TFile* f2 = TFile::Open("hist_ee_AccEff_metBin.root","READ");
    if (!f2->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file.\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    TH1D* heeAcceff = ((TH1D*)f2->Get("Hacceff"));
    if (!heeAcceff) {
      cout << "Error: histo not read" << endl;
      exit(EXIT_FAILURE);
    }

  // check if 2 histogram have same number of bins (they should by definition, but better to check9
  if (harray[0]->GetNbinsX() != harray[1]->GetNbinsX()) {
    cout << " Warning: the 2 histograms have different number of bins. The lower number will be used." << endl;
  }
  // in case the number of bins is different, use the lower one (hoping at least that the first low edge is the same for both)
  Int_t lessBinHistoIndex = (harray[0]->GetNbinsX() <= harray[1]->GetNbinsX()) ? 0 : 1; 
  Int_t nMetBins = harray[lessBinHistoIndex]->GetNbinsX();
  Double_t metBinEdges[nMetBins + 1];
  harray[lessBinHistoIndex]->GetXaxis()->GetLowEdge(metBinEdges); // this doesn't include the last bin up edge (the overflow low edge)
  metBinEdges[nMetBins] = harray[lessBinHistoIndex]->GetXaxis()->GetBinLowEdge(nMetBins + 1); // get low edge of overflow bin
  for (Int_t i = 0; i <= nMetBins; i++) {
    cout<<" metBinEdges["<<i<<"] = "<<metBinEdges[i]<<endl;
  }

  //correct yields for acceptance
  if ( !harray[0]->Divide(hmumuAcceff) ) cout << "Error: in harray[0]->Divide()" <<endl;
  if ( !harray[1]->Divide(heeAcceff) ) cout << "Error: in harray[1]->Divide()" <<endl;

  for (Int_t j = 0; j < nFiles; j++) {
    harray[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    harray[j]->SetLineColor(histColor[j]);
    for (Int_t i = 1; i <= harray[j]->GetNbinsX(); i++ ) {
      harray[j]->SetBinContent(i, harray[j]->GetBinContent(i)/harray[j]->GetBinWidth(i));
      harray[j]->SetBinError(i, harray[j]->GetBinError(i)/harray[j]->GetBinWidth(i));
    }
    harray[j]->Scale(1./harray[j]->Integral(0,1 + harray[j]->GetNbinsX()));  // normalize to unity, use integral including underflow & overflow bin
    //harray[j]->SetFillColorAlpha(histColor[j], 0.80);   
  }

  
  TCanvas *c = new TCanvas("c","");
  TLegend *leg = new TLegend(0.7,0.7,0.89,0.89); 
  c->SetLogy();
  //cout <<"Drawing histograms and ratio"<<endl;
  harray[0]->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  harray[0]->GetXaxis()->SetTitleSize(0.04);
  harray[0]->GetYaxis()->SetTitle("[GeV^{-1}]");
  harray[0]->GetYaxis()->SetTitleSize(0.04);
  harray[0]->GetYaxis()->CenterTitle();
  harray[0]->Draw("HE");
  harray[1]->Draw("HE SAME");
  harray[2]->Draw("HE SAME");
  for (Int_t j = 0; j < total; j++) {
    leg->AddEntry(harray[j],Form("%s",sample[j].c_str()),"lf");
  }
  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  c->SaveAs("met_NoMu_NoEle_NuNu_comparison_binDensity.root");  

  TCanvas *cRatioMuEle = new TCanvas("cRatioMuEle","");
  TH1D *hratio = new TH1D("hratio","",nMetBins,metBinEdges);
  hratio->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hratio->GetXaxis()->SetTitleSize(0.04);
  hratio->GetYaxis()->SetTitle("MetNoMu / MetNoEle");
  hratio->GetYaxis()->SetTitleSize(0.04);
  // use Divide method
  hratio->Divide(harray[0],harray[1]);
  hratio->SetStats(kFALSE);
  hratio->Draw("E");
  cRatioMuEle->SaveAs("met_NoMu_NoEle_ratio_binDensity.root");

  TCanvas *cRatioLepNu = new TCanvas("cRatioLepNu","");
  TLegend *legRatioLepNu = new TLegend(0.7,0.7,0.89,0.89);
  TH1D *hratioMuNu = new TH1D("hratioMuNu","",nMetBins,metBinEdges);
  hratioMuNu->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hratioMuNu->GetXaxis()->SetTitleSize(0.04);
  hratioMuNu->GetYaxis()->SetTitle("MetNoL(Z->ll) / MetNoMu(Z->#nu#nu)");
  hratioMuNu->GetYaxis()->SetTitleSize(0.045);
  // use Divide method
  hratioMuNu->Divide(harray[0],harray[2]);
  hratioMuNu->SetStats(kFALSE);
  hratioMuNu->SetLineColor(kBlue);
  hratioMuNu->GetYaxis()->SetRangeUser(0.85,1.6);
  hratioMuNu->Draw("E");
  TH1D *hratioEleNu = new TH1D("hratioEleNu","",nMetBins,metBinEdges);
  hratioEleNu->Divide(harray[1],harray[2]);
  hratioEleNu->SetStats(kFALSE);
  hratioEleNu->SetLineColor(kRed);
  hratioEleNu->Draw("E SAME");
  legRatioLepNu->AddEntry(hratioMuNu,"#mu / #nu","l");
  legRatioLepNu->AddEntry(hratioEleNu,"e / #nu","l");
  legRatioLepNu->Draw(); 
  legRatioLepNu->SetMargin(0.3); 
  legRatioLepNu->SetBorderSize(0);
  cRatioLepNu->SaveAs("met_NoLep_Nu_ratio_binDensity.root");

}
