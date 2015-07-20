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

void znnEstimate_znnMC_zlljetsCS(const char* suffix = "") 
{

  TH1::SetDefaultSumw2();   

  string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/ZtoLLSamples/znunuEstimate/";
  string plotFileExtension = ".pdf";

  vector<TH1D*> hMetShape;
  
  vector<TH1D*> hAxe;                     // only for Z(mumu) and Z(ee)
  vector<TH1D*> hZvvYieldsEstimate;       // only for Z(mumu) and Z(ee)

  TH1D* hZvvYields = NULL;

  vector<string> fileName;
  fileName.push_back("zmumujetsAna.root");
  fileName.push_back("zeejetsAna.root");
  fileName.push_back("znunujetsAna.root");
  vector<string> Axe_fileName;
  Axe_fileName.push_back("zmumujets_Axe_noSkim_light.root");
  Axe_fileName.push_back("zeejets_Axe_noSkim_light.root");
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
      hZvvYieldsEstimate.push_back((TH1D*)f->Get("HzvvEstimateNoAxe"));

    }
	 
  }

  cout<<"fileName : "<<Axe_fileName[0]<<endl;

  TFile* f1 = TFile::Open(Axe_fileName[0].c_str(),"READ");
    
  if (!f1 || !f1->IsOpen()) {

    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<Axe_fileName[0]<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);

  }
    
  // use Axe computed dividing histograms, but it's the same as that computed with TEfficiency. The advantage is that i have an histogram, while TEfficiency only provides TGraphs
  hAxe.push_back((TH1D*)f1->Get("HacceffW")); 

  cout<<"fileName : "<<Axe_fileName[1]<<endl;

  TFile* f2 = TFile::Open(Axe_fileName[1].c_str(),"READ");
    
  if (!f2 || !f2->IsOpen()) {

    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<Axe_fileName[1]<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);

  }
    
  // use Axe computed dividing histograms, but it's the same as that computed with TEfficiency. The advantage is that i have an histogram, while TEfficiency only provides TGraphs
  hAxe.push_back((TH1D*)f2->Get("HacceffW")); 
	 
  cout <<"check" << endl;

  for (Int_t j = 0; j < 2; j++) {
    // correct shape dividing by Axe bin by bin
    if ( !hZvvYieldsEstimate[j]->Divide(hAxe[j]) ) cout << "Error in hZvvYieldsEstimate[" << j << "]->Divide(hAxe[" <<j << "])" <<endl;
    if ( !hMetShape[j]->Divide(hAxe[j]) ) cout << "Error in hMetShape[" << j << "]->Divide(hAxe[" <<j << "])" <<endl;
    hZvvYieldsEstimate[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hZvvYieldsEstimate[j]->SetLineColor(histColor[j]);
    hAxe[j]->SetStats(kFALSE);   // to avoid drawing statistic box
    hAxe[j]->SetLineColor(histColor[j]);

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
  c->SaveAs( (plotDirectoryPath + "metShape" + suffix + plotFileExtension).c_str()); 

  TCanvas *cZvvYieldsEstimate = new TCanvas("ZvvYieldsEstimate","");
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
  cZvvYieldsEstimate->SaveAs( (plotDirectoryPath + cZvvYieldsEstimate->GetName() + suffix + plotFileExtension).c_str()); 


  TCanvas *cRatioMuEle = new TCanvas("ratioMuEle","");
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
  cRatioMuEle->SaveAs( (plotDirectoryPath + cRatioMuEle->GetName() + suffix + plotFileExtension).c_str()); 

  TCanvas *cRatioLepNu = new TCanvas("ratioLepNu","");
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
  cRatioLepNu->SaveAs( (plotDirectoryPath + cRatioLepNu->GetName() + suffix + plotFileExtension).c_str()); 

  TCanvas *cRatioZvvEst_LepNu = new TCanvas("ratioZvvEst_LepNu","");
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
  cRatioZvvEst_LepNu->SaveAs( (plotDirectoryPath + cRatioZvvEst_LepNu->GetName() + suffix + plotFileExtension).c_str());

  TCanvas *czll_Axe = new TCanvas("zll_Axe","");
  TLegend *legzll_Axe = new TLegend(0.65,0.7,0.89,0.89); 
  hAxe[0]->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");
  hAxe[0]->GetXaxis()->SetTitleSize(0.04);
  hAxe[0]->GetYaxis()->SetTitle("A#times#epsilon");
  hAxe[0]->GetYaxis()->SetTitleSize(0.045);
  hAxe[0]->GetYaxis()->CenterTitle();
  hAxe[0]->Draw("HE");
  hAxe[1]->Draw("HE SAME");
  hZvvYields->Draw("HE SAME");
  for (Int_t j = 0; j < 2; j++) {
    legzll_Axe->AddEntry(hAxe[j],Form("%s",sample[j].c_str()),"l");
  }
  gStyle->SetStatStyle(0);
  legzll_Axe->Draw(); 
  legzll_Axe->SetHeader("A#times#epsilon for Z(ll)+jets"); 
  legzll_Axe->SetMargin(0.3); 
  legzll_Axe->SetBorderSize(0);
  //czll_Axe->SaveAs("met_NoMu_NoEle_NuNu_comparison.root");  
  czll_Axe->SaveAs( (plotDirectoryPath + czll_Axe->GetName() + suffix + plotFileExtension).c_str()); 


}
