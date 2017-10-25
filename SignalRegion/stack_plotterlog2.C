#include <fstream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "iostream"
void stack_plotterlog2()
{
  TCanvas *c = new TCanvas("c", "canvas",800,800);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  //c->SetLeftMargin(0.15);
  //c->SetLogy();
  //c->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.25,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetLogy();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  //opening the data file and adding "j1pT_10" histogram
  TFile *f_datafile_0 = new TFile("postMETdata_final.root");
  //TFile *f_datafile_1 = new TFile("postMETdata_1.root");
  TH1F *histo_j1EtaWidth_data_0 = (TH1F*)f_datafile_0->Get("j1pT_10");
  //TH1F *histo_j1EtaWidth_data_1 = (TH1F*)f_datafile_1->Get("j1pT_10");
  //histo_j1EtaWidth_data_0->Add(histo_j1EtaWidth_data_1);

  double integraldata = histo_j1EtaWidth_data_0->Integral();
  std::cout<<"integral of Data here:"<<integraldata<<std::endl;

  histo_j1EtaWidth_data_0->SetStats(0);
  histo_j1EtaWidth_data_0->SetLineWidth(2);
  histo_j1EtaWidth_data_0->SetLineColor(kWhite);
  histo_j1EtaWidth_data_0->SetTitle("");
  histo_j1EtaWidth_data_0->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_data_0->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_data_0->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_data_0->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_data_0->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_data_0->GetYaxis()->SetLabelOffset(999);
  //histo_j1EtaWidth_data_0->Draw("");

  //opening background WJets Sample file
  TFile *f_WJets_0 = new TFile("postWJets_MLM_0.root");
  TFile *f_W1Jets = new TFile("postW100to200_0.root");
  TFile *f_W2Jets = new TFile("postW200to400_0.root");
  TFile *f_W3Jets = new TFile("postW400to600_0.root");
  TFile *f_W4Jets = new TFile("postW600to800_0.root");
  TFile *f_W5Jets = new TFile("postW800to1200_0.root");
  TFile *f_W6Jets = new TFile("postW1200to2500_0.root");
  TFile *f_W7Jets = new TFile("postW2500toInf_0.root");

  TH1F *histo_j1EtaWidth_WJets_0 = (TH1F*)f_WJets_0->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W1Jets = (TH1F*)f_W1Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W2Jets = (TH1F*)f_W2Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W3Jets = (TH1F*)f_W3Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W4Jets = (TH1F*)f_W4Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W5Jets = (TH1F*)f_W5Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W6Jets = (TH1F*)f_W6Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_W7Jets = (TH1F*)f_W7Jets->Get("j1pT_10");
  
  histo_j1EtaWidth_WJets_0->SetStats(0);
  histo_j1EtaWidth_W1Jets->SetStats(0);
  histo_j1EtaWidth_W2Jets->SetStats(0);
  histo_j1EtaWidth_W3Jets->SetStats(0);
  histo_j1EtaWidth_W4Jets->SetStats(0); 
  histo_j1EtaWidth_W5Jets->SetStats(0);
  histo_j1EtaWidth_W6Jets->SetStats(0);
  histo_j1EtaWidth_W7Jets->SetStats(0);

  double rawWJets = (histo_j1EtaWidth_WJets_0->Integral())+ (histo_j1EtaWidth_W1Jets->Integral())+ (histo_j1EtaWidth_W2Jets->Integral())+ (histo_j1EtaWidth_W3Jets->Integral())+ (histo_j1EtaWidth_W4Jets->Integral())+ (histo_j1EtaWidth_W5Jets->Integral())+ (histo_j1EtaWidth_W6Jets->Integral())+ (histo_j1EtaWidth_W7Jets->Integral());
  std::cout<<"raw WJets events: "<<rawWJets<<std::endl;

   // Scaling = (1/Totalevents)*Luminosity*NNLO-cross-section
  histo_j1EtaWidth_WJets_0->Scale((1.0/57025300)*12570*50690);
  histo_j1EtaWidth_W1Jets->Scale((1.0/39617100)*12570*1345);
  histo_j1EtaWidth_W2Jets->Scale((1.0/19914100)*12570*359.7);
  histo_j1EtaWidth_W3Jets->Scale((1.0/5795950)*12570*48.91);
  histo_j1EtaWidth_W4Jets->Scale((1.0/14907300)*12570*12.05);
  histo_j1EtaWidth_W5Jets->Scale((1.0/6200210)*12570*5.501);
  histo_j1EtaWidth_W6Jets->Scale((1.0/5527630)*12570*1.329);
  histo_j1EtaWidth_W7Jets->Scale((1.0/2382300)*12570*0.03216);
  
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W1Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W2Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W3Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W4Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W5Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W6Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W7Jets);  

  double integralWJets = histo_j1EtaWidth_WJets_0->Integral();
  std::cout<<"integral of WJets bkg here:"<<integralWJets<<std::endl;

  histo_j1EtaWidth_WJets_0->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WJets_0->SetFillColor(kRed-10);

  TFile *f_Zvv_100to200 = new TFile("postZ100to200_0.root");
  TFile *f_Zvv_200to400 = new TFile("postZ200to400_0.root");
  TFile *f_Zvv_400to600 = new TFile("postZ400to600_0.root");
  TFile *f_Zvv_600to800 = new TFile("postZ600to800_0.root");
  TFile *f_Zvv_800to1200 = new TFile("postZ800to1200_0.root");
  TFile *f_Zvv_1200to2500 = new TFile("postZ1200to2500_0.root");
  TFile *f_Zvv_2500toInf = new TFile("postZ2500toInf_0.root");

  TH1F *histo_j1EtaWidth_100to200 = (TH1F*)f_Zvv_100to200->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_200to400 = (TH1F*)f_Zvv_200to400->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_400to600 = (TH1F*)f_Zvv_400to600->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_600to800 = (TH1F*)f_Zvv_600to800->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_800to1200 = (TH1F*)f_Zvv_800to1200->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_1200to2500 = (TH1F*)f_Zvv_1200to2500->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_2500toInf = (TH1F*)f_Zvv_2500toInf->Get("j1pT_10");
  histo_j1EtaWidth_100to200->SetStats(0);
  histo_j1EtaWidth_200to400->SetStats(0);
  histo_j1EtaWidth_400to600->SetStats(0);
  histo_j1EtaWidth_600to800->SetStats(0);
  histo_j1EtaWidth_800to1200->SetStats(0);
  histo_j1EtaWidth_1200to2500->SetStats(0);
  histo_j1EtaWidth_2500toInf->SetStats(0);

  double rawZvvJets = (histo_j1EtaWidth_100to200->Integral())+ (histo_j1EtaWidth_200to400->Integral())+ (histo_j1EtaWidth_400to600->Integral())+ (histo_j1EtaWidth_600to800->Integral())+ (histo_j1EtaWidth_800to1200->Integral())+ (histo_j1EtaWidth_1200to2500->Integral())+ (histo_j1EtaWidth_2500toInf->Integral()); 
  std::cout<<"raw ZvvJets bkg:"<<rawZvvJets<<std::endl;
  
  // Scaling = (1/Totalevents)*Luminosity*LO-cross-section
  histo_j1EtaWidth_100to200->Scale((1.0/19026300)*12570*280.35);
  histo_j1EtaWidth_200to400->Scale((1.0/19611900)*12570*77.67);
  histo_j1EtaWidth_400to600->Scale((1.0/8842310)*12570*10.73);
  histo_j1EtaWidth_600to800->Scale((1.0/5763180)*12570*2.559);
  histo_j1EtaWidth_800to1200->Scale((1.0/2169900)*12570*1.1796);
  histo_j1EtaWidth_1200to2500->Scale((1.0/143919)*12570*0.28833);
  histo_j1EtaWidth_2500toInf->Scale((1.0/404725)*12570*0.006945);

  //Add the ZJetsToNuNu histograms to the first one
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_200to400);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_400to600);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_600to800);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_800to1200);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_1200to2500);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_2500toInf);


  histo_j1EtaWidth_100to200->SetTitle("");
  histo_j1EtaWidth_100to200->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_100to200->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_100to200->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_100to200->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_100to200->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_100to200->GetYaxis()->SetLabelOffset(999);  
  histo_j1EtaWidth_100to200->SetFillColor(kAzure+10);


  double integralZvvJets = histo_j1EtaWidth_100to200->Integral();
  std::cout<<"integral of ZvvJets bkg here:"<<integralZvvJets<<std::endl;

  //opening background samples Gamma+jets
  TFile *f_G1Jets = new TFile("postGJets40to100.root");
  TFile *f_G2Jets = new TFile("postGJets100to200.root");
  TFile *f_G3Jets = new TFile("postGJets200to400.root");
  TFile *f_G4Jets = new TFile("postGJets400to600.root");
  TFile *f_G5Jets = new TFile("postGJets600toInf.root");
  
  TH1F *histo_j1EtaWidth_G1Jets = (TH1F*)f_G1Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_G2Jets = (TH1F*)f_G2Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_G3Jets = (TH1F*)f_G3Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_G4Jets = (TH1F*)f_G4Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_G5Jets = (TH1F*)f_G5Jets->Get("j1pT_10");

  histo_j1EtaWidth_G1Jets->SetStats(0);
  histo_j1EtaWidth_G2Jets->SetStats(0);
  histo_j1EtaWidth_G3Jets->SetStats(0);
  histo_j1EtaWidth_G4Jets->SetStats(0);
  histo_j1EtaWidth_G5Jets->SetStats(0);

  double rawGJets = (histo_j1EtaWidth_G1Jets->Integral())+ (histo_j1EtaWidth_G2Jets->Integral())+ (histo_j1EtaWidth_G3Jets->Integral())+ (histo_j1EtaWidth_G4Jets->Integral())+ (histo_j1EtaWidth_G5Jets->Integral()); 
  std::cout<<"raw GJets bkg:"<<rawGJets<<std::endl;
  //Scaling
  histo_j1EtaWidth_G1Jets->Scale((1.0/11490000)*12570*17420);
  histo_j1EtaWidth_G2Jets->Scale((1.0/14375000)*12570*5391);
  histo_j1EtaWidth_G3Jets->Scale((1.0/49455800)*12570*1168);
  histo_j1EtaWidth_G4Jets->Scale((1.0/11582100)*12570*132.5);
  histo_j1EtaWidth_G5Jets->Scale((1.0/11611900)*12570*44.05);
  
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G2Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G3Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G4Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G5Jets);

  histo_j1EtaWidth_G1Jets->SetTitle("");
  histo_j1EtaWidth_G1Jets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_G1Jets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_G1Jets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_G1Jets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_G1Jets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_G1Jets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_G1Jets->SetFillColor(kTeal-9);

  double integralGJets = histo_j1EtaWidth_G1Jets->Integral();
  std::cout<<"integral of GJets bkg here:"<<integralGJets<<std::endl;

//opening DYJetsToLL backgrounds
  TFile *f_DY1Jets = new TFile("postDY_MLM_0.root");
  TFile *f_DY2Jets = new TFile("postDY100to200.root");
  TFile *f_DY3Jets = new TFile("postDY200to400.root");
  TFile *f_DY4Jets = new TFile("postDY400to600.root");
  TFile *f_DY5Jets = new TFile("postDY600to800.root");
  TFile *f_DY6Jets = new TFile("postDY800to1200.root");
  TFile *f_DY7Jets = new TFile("postDY1200to2500.root");
  TFile *f_DY8Jets = new TFile("postDY2500toInf.root");

  TH1F *histo_j1EtaWidth_DY1Jets = (TH1F*)f_DY1Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY2Jets = (TH1F*)f_DY2Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY3Jets = (TH1F*)f_DY3Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY4Jets = (TH1F*)f_DY4Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY5Jets = (TH1F*)f_DY5Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY6Jets = (TH1F*)f_DY6Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY7Jets = (TH1F*)f_DY7Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_DY8Jets = (TH1F*)f_DY8Jets->Get("j1pT_10");

  histo_j1EtaWidth_DY1Jets->SetStats(0);
  histo_j1EtaWidth_DY2Jets->SetStats(0);
  histo_j1EtaWidth_DY3Jets->SetStats(0);
  histo_j1EtaWidth_DY4Jets->SetStats(0);
  histo_j1EtaWidth_DY5Jets->SetStats(0);
  histo_j1EtaWidth_DY6Jets->SetStats(0);
  histo_j1EtaWidth_DY7Jets->SetStats(0);
  histo_j1EtaWidth_DY8Jets->SetStats(0);

  double rawDY = (histo_j1EtaWidth_DY1Jets->Integral())+ (histo_j1EtaWidth_DY2Jets->Integral())+ (histo_j1EtaWidth_DY3Jets->Integral())+ (histo_j1EtaWidth_DY4Jets->Integral())+ (histo_j1EtaWidth_DY5Jets->Integral());  
  std::cout<<"raw DYJets bkg:"<<rawDY<<std::endl;

  histo_j1EtaWidth_DY1Jets->Scale((1.0/47644783)*12570*4895*1.23);
  histo_j1EtaWidth_DY2Jets->Scale((1.0/7855800)*12570*147.4*1.23);
  histo_j1EtaWidth_DY3Jets->Scale((1.0/8691200)*12570*40.99*1.23);
  histo_j1EtaWidth_DY4Jets->Scale((1.0/8937770)*12570*5.678*1.23);
  histo_j1EtaWidth_DY5Jets->Scale((1.0/8292150)*12570*1.367*1.23);
  histo_j1EtaWidth_DY6Jets->Scale((1.0/2668260)*12570*0.6304*1.23);
  histo_j1EtaWidth_DY7Jets->Scale((1.0/555134)*12570*0.1514*1.23);
  histo_j1EtaWidth_DY8Jets->Scale((1.0/398370)*12570*0.003565*1.23);

  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY2Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY3Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY4Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY5Jets);

  histo_j1EtaWidth_DY1Jets->SetTitle("");
  histo_j1EtaWidth_DY1Jets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_DY1Jets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_DY1Jets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_DY1Jets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_DY1Jets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_DY1Jets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_DY1Jets->SetFillColor(kGray+2);

  double integralDY = histo_j1EtaWidth_DY1Jets->Integral();
  std::cout<<"integral of DYJets bkg here:"<<integralDY<<std::endl;

  //opening background TTJets
  TFile *f_TTJets = new TFile("postTTJets_MLM.root");
  TH1F *histo_j1EtaWidth_TTJets = (TH1F*)f_TTJets->Get("j1pT_10");
  histo_j1EtaWidth_TTJets->SetStats(0);

  double rawTTJets = histo_j1EtaWidth_TTJets->Integral();
  std::cout<<"raw TTJets here:"<<rawTTJets<<std::endl;

  histo_j1EtaWidth_TTJets->Scale((1.0/10139700)*12570*502.2);
  histo_j1EtaWidth_TTJets->SetTitle("");
  histo_j1EtaWidth_TTJets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_TTJets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_TTJets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_TTJets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_TTJets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_TTJets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_TTJets->SetFillColor(kOrange-2);
 
  double integralTTJets = histo_j1EtaWidth_TTJets->Integral();
  std::cout<<"integral of integralTTJets bkg here:"<<integralTTJets<<std::endl;

 //addding some backgrounds like WW, WZ, ZZ	

  TFile *f_WW = new TFile("postWW.root");
  TH1F *histo_j1EtaWidth_WW = (TH1F*)f_WW->Get("j1pT_10");
  histo_j1EtaWidth_WW->SetStats(0);
  histo_j1EtaWidth_WW->Scale((1.0/6986980)*12570*118.7);
  /*histo_j1EtaWidth_WW->SetTitle("");
  histo_j1EtaWidth_WW->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WW->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WW->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WW->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WW->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WW->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WW->SetFillColor(kCyan+2);
*/
  double integralWW = histo_j1EtaWidth_WW->Integral();
  std::cout<<"integral of WW bkg here:"<<integralWW<<std::endl;
 
  TFile *f_WZ = new TFile("postWZ.root");
  TH1F *histo_j1EtaWidth_WZ = (TH1F*)f_WZ->Get("j1pT_10");
  histo_j1EtaWidth_WZ->SetStats(0);
  histo_j1EtaWidth_WZ->Scale((1.0/2995760)*12570*47.2);
/*  histo_j1EtaWidth_WZ->SetTitle("");
  histo_j1EtaWidth_WZ->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WZ->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WZ->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WZ->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WZ->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WZ->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WZ->SetFillColor(kAzure+10);
*/
  double integralWZ = histo_j1EtaWidth_WZ->Integral();
  std::cout<<"integral of WZ bkg here:"<<integralWZ<<std::endl;

  TFile *f_ZZ = new TFile("postZZ.root");
  TH1F *histo_j1EtaWidth_ZZ = (TH1F*)f_ZZ->Get("j1pT_10");
  histo_j1EtaWidth_ZZ->SetStats(0);
  histo_j1EtaWidth_ZZ->Scale((1.0/998014)*12570*16.6);
/*  histo_j1EtaWidth_ZZ->SetTitle("");
  histo_j1EtaWidth_ZZ->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_ZZ->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_ZZ->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_ZZ->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_ZZ->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_ZZ->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_ZZ->SetFillColor(kCyan);
*/
  double integralZZ = histo_j1EtaWidth_ZZ->Integral();
  std::cout<<"integral of ZZ bkg here:"<<integralZZ<<std::endl;
 
  //Ading WW, WZ, ZZ together
  histo_j1EtaWidth_WW->Add(histo_j1EtaWidth_WZ);
  histo_j1EtaWidth_WW->Add(histo_j1EtaWidth_ZZ);
  histo_j1EtaWidth_WW->SetTitle("");
  histo_j1EtaWidth_WW->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WW->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WW->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WW->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WW->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WW->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WW->SetFillColor(kCyan-10);

  double TotintegralDiBoson = histo_j1EtaWidth_WW->Integral();
  std::cout<<"integral of WW/WZ/ZZ bkg here:"<<TotintegralDiBoson<<std::endl;

 //opening QCD background files (METValue-binned samples)
  TFile *f_Q1Jets = new TFile("postQCD100to200_0.root");
  TFile *f_Q2Jets = new TFile("postQCD200to300_0.root");
  TFile *f_Q3Jets = new TFile("postQCD300to500_0.root");
  TFile *f_Q4Jets = new TFile("postQCD500to700_0.root");
  TFile *f_Q5Jets = new TFile("postQCD700to1000_0.root");
  TFile *f_Q6Jets = new TFile("postQCD1000to1500_0.root");
  TFile *f_Q7Jets = new TFile("postQCD1500to2000_0.root");
  TFile *f_Q8Jets = new TFile("postQCD2000toInf_0.root");

  TH1F *histo_j1EtaWidth_Q1Jets = (TH1F*)f_Q1Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q2Jets = (TH1F*)f_Q2Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q3Jets = (TH1F*)f_Q3Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q4Jets = (TH1F*)f_Q4Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q5Jets = (TH1F*)f_Q5Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q6Jets = (TH1F*)f_Q6Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q7Jets = (TH1F*)f_Q7Jets->Get("j1pT_10");
  TH1F *histo_j1EtaWidth_Q8Jets = (TH1F*)f_Q8Jets->Get("j1pT_10");

  histo_j1EtaWidth_Q1Jets->SetStats(0);
  histo_j1EtaWidth_Q2Jets->SetStats(0);
  histo_j1EtaWidth_Q3Jets->SetStats(0);
  histo_j1EtaWidth_Q4Jets->SetStats(0);
  histo_j1EtaWidth_Q5Jets->SetStats(0);
  histo_j1EtaWidth_Q6Jets->SetStats(0);
  histo_j1EtaWidth_Q7Jets->SetStats(0);
  histo_j1EtaWidth_Q8Jets->SetStats(0);

  histo_j1EtaWidth_Q1Jets->Scale((1.0/80584300)*12570*27500000);
  histo_j1EtaWidth_Q2Jets->Scale((1.0/38742400)*12570*1735000);
  histo_j1EtaWidth_Q3Jets->Scale((1.0/37501400)*12570*367000);
  histo_j1EtaWidth_Q4Jets->Scale((1.0/43210300)*12570*29370);
  histo_j1EtaWidth_Q5Jets->Scale((1.0/29695500)*12570*6524);
  histo_j1EtaWidth_Q6Jets->Scale((1.0/10358900)*12570*1064);
  histo_j1EtaWidth_Q7Jets->Scale((1.0/7853460)*12570*121.5);
  histo_j1EtaWidth_Q8Jets->Scale((1.0/4045060)*12570*25.42);  

  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q2Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q3Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q4Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q5Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q6Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q7Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q8Jets);
  
  double integral = histo_j1EtaWidth_Q1Jets->Integral();
  std::cout<<"integral of QCD bkg here:"<<integral<<std::endl;

  histo_j1EtaWidth_Q1Jets->SetTitle("");
  histo_j1EtaWidth_Q1Jets->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_Q1Jets->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_Q1Jets->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_Q1Jets->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_Q1Jets->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_Q1Jets->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_Q1Jets->SetFillColor(kGray);

  //Stack histograms using THStack
  THStack *hs_datamc = new THStack("hs_datamc","Data/MC comparison"); 
  hs_datamc->Add(histo_j1EtaWidth_Q1Jets);
  //hs_datamc->Add(histo_j1EtaWidth_WW);
  //hs_datamc->Add(histo_j1EtaWidth_100to200);
  hs_datamc->Add(histo_j1EtaWidth_DY1Jets);
  hs_datamc->Add(histo_j1EtaWidth_G1Jets);
  hs_datamc->Add(histo_j1EtaWidth_TTJets);
  //hs_datamc->Add(histo_j1EtaWidth_G1Jets);
  hs_datamc->Add(histo_j1EtaWidth_WW);
  //hs_datamc->Add(histo_j1EtaWidth_TTJets);
  hs_datamc->Add(histo_j1EtaWidth_WJets_0);
  //hs_datamc->Add(histo_j1EtaWidth_DY1Jets); 
  //hs_datamc->Add(histo_j1EtaWidth_Q1Jets);
  hs_datamc->Add(histo_j1EtaWidth_100to200);
  //hs_datamc->Add(histo_j1EtaWidth_WJets_0);
  //hs_datamc->Add(histo_j1EtaWidth_G1Jets);
  hs_datamc->SetTitle("");
  hs_datamc->Draw("HIST");
  hs_datamc->SetMinimum(0.1);
  hs_datamc->SetMaximum(pow(10,4));
  hs_datamc->Draw("HIST");
  histo_j1EtaWidth_data_0->SetLineColor(kBlack);
  histo_j1EtaWidth_data_0->SetMarkerStyle(20);
  histo_j1EtaWidth_data_0->SetMarkerSize(0.7);
  histo_j1EtaWidth_data_0->Draw("pex0same");
 
//  TFile *f_signal_1GeVfile = new TFile("postSignal_mchi1_Mzp1_MET300.root");
//  TH1F *histo_signal_1GeV = (TH1F*)f_signal_1GeVfile ->Get("j1pT_10");
//  histo_signal_1GeV->Scale((1.0/629)*12570*0.056);
//  histo_signal_1GeV->SetLineColor(kRed);
//  histo_signal_1GeV->SetLineWidth(2);
//  histo_signal_1GeV->Draw("HIST SAME");
//
//  double integralsignal_1GeV = histo_signal_1GeV->Integral(); 
//  std::cout<<"integral of signal_1GeV here:"<<integralsignal_1GeV<<std::endl;
//
  TFile *f_signal_5GeVfile = new TFile("postSignal.root");
  TH1F *histo_signal_5GeV = (TH1F*)f_signal_5GeVfile ->Get("j1pT_10");
  histo_signal_5GeV->Scale((1.0/2133)*12570*0.047);
  histo_signal_5GeV->SetLineColor(kBlue);
  histo_signal_5GeV->SetLineWidth(2);
  histo_signal_5GeV->Draw("HIST SAME");
//
  double integralsignal_5GeV = histo_signal_5GeV->Integral(); 
  std::cout<<"integral of signal_5GeV here:"<<integralsignal_5GeV<<std::endl;
//
//  TFile *f_signal_20GeVfile = new TFile("postSignal_mchi20_Mzp1_MET300.root");
//  TH1F *histo_signal_20GeV = (TH1F*)f_signal_20GeVfile ->Get("j1pT_10");
//  histo_signal_20GeV->Scale((1.0/10000)*12570*0.034);
//  histo_signal_20GeV->SetLineColor(kMagenta);
//  histo_signal_20GeV->SetLineWidth(2);
//  histo_signal_20GeV->Draw("HIST SAME");
//
//  double integralsignal_20GeV = histo_signal_20GeV->Integral(); 
//  std::cout<<"integral of signal_20GeV here:"<<integralsignal_20GeV<<std::endl;
// 
//  TFile *f_signal_50GeVfile = new TFile("postSignal_mchi50_Mzp1_MET300.root");
//  TH1F *histo_signal_50GeV = (TH1F*)f_signal_50GeVfile ->Get("j1pT_10");
//  histo_signal_50GeV->Scale((1.0/10000)*12570*0.025);
//  histo_signal_50GeV->SetLineColor(kSpring-1);
//  histo_signal_50GeV->SetLineWidth(2);
//  histo_signal_50GeV->Draw("HIST SAME");
//
//  double integralsignal_50GeV = histo_signal_50GeV->Integral(); 
//  std::cout<<"integral of signal_50GeV here:"<<integralsignal_50GeV<<std::endl;
//  
//  TFile *f_signal_10GeVfile = new TFile("postSignal_mchi10_Mzp1_MET300.root");
//  TH1F *histo_signal_10GeV = (TH1F*)f_signal_10GeVfile ->Get("j1pT_10");
//  histo_signal_10GeV->Scale((1.0/4052)*12570*0.04);
//  histo_signal_10GeV->SetLineColor(kViolet+1);
//  histo_signal_10GeV->SetLineWidth(2);
//  histo_signal_10GeV->Draw("HIST SAME");
//
//  double integralsignal_10GeV = histo_signal_10GeV->Integral(); 
//  std::cout<<"integral of signal_10GeV here:"<<integralsignal_10GeV<<std::endl;
//  
//  TFile *f_signal_100GeVfile = new TFile("postSignal_mchi100_Mzp1_MET300.root");
//  TH1F *histo_signal_100GeV = (TH1F*)f_signal_100GeVfile ->Get("j1pT_10");
//  histo_signal_100GeV->Scale((1.0/10000)*12570*0.019);
//  histo_signal_100GeV->SetLineColor(kAzure+1);
//  histo_signal_100GeV->SetLineWidth(2);
//  histo_signal_100GeV->Draw("HIST SAME");
//
//  double integralsignal_100GeV = histo_signal_100GeV->Integral(); 
//  std::cout<<"integral of signal_100GeV here:"<<integralsignal_100GeV<<std::endl;
  //TLegend *leg = new TLegend(0.181948,0.663948,0.567335,0.836868,"");
  TLegend *leg = new TLegend(0.62,0.60,0.86,0.887173,"");
  leg->AddEntry(histo_j1EtaWidth_data_0,"Data");
//  leg->AddEntry(histo_signal_1GeV, "ZprimeSignal_mchi1GeV");
  leg->AddEntry(histo_signal_5GeV, "ZprimeSignal_mchi5GeV"); 
//  leg->AddEntry(histo_signal_10GeV, "ZprimeSignal_mchi10GeV");
//  leg->AddEntry(histo_signal_20GeV, "ZprimeSignal_mchi20GeV"); 
//  leg->AddEntry(histo_signal_50GeV, "ZprimeSignal_mchi50GeV");
//  leg->AddEntry(histo_signal_100GeV, "ZprimeSignal_mchi100GeV");
  leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu","F"); 
  //leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets", "F");
  //leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD","F");
  leg->AddEntry(histo_j1EtaWidth_WJets_0,"W#rightarrowl#nu","F");
  //leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu","F"); 
  //leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD","F");
  //leg->AddEntry(histo_j1EtaWidth_TTJets, "Top Quark", "F");
  leg->AddEntry(histo_j1EtaWidth_WW,"WW/WZ/ZZ","F");
  //leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets", "F");
  leg->AddEntry(histo_j1EtaWidth_TTJets, "Top Quark", "F");
  leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets", "F");
  leg->AddEntry(histo_j1EtaWidth_DY1Jets,"DYJets#rightarrowLL","F");  
  //leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu","F");
  //leg->AddEntry(histo_j1EtaWidth_WW,"WW/WZ/ZZ","F");
  leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD","F");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.025);
  leg->Draw();  

  TLatex *texS = new TLatex(0.20,0.837173,"#sqrt{s} = 13 TeV, 12.57 fb^{-1}");
  texS->SetNDC();
  texS->SetTextFont(42);
  texS->SetTextSize(0.040);
  texS->Draw();
  TLatex *texS1 = new TLatex(0.12092,0.907173,"#bf{CMS} : #it{Preliminary}");
  texS1->SetNDC();
  texS1->SetTextFont(42);
  texS1->SetTextSize(0.040);
  texS1->Draw();

  c->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.25);
  pad2->Draw(); pad2->cd();
  pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
 
  //TStyle * style = getStyle("ZZ");
  //style->cd();
 
  int nbins = histo_j1EtaWidth_data_0->GetNbinsX();  
  TH1F* Ratio = (TH1F*)histo_j1EtaWidth_data_0->Clone("Ratio");
  TH1F *last_hist = (TH1F*)hs_datamc->GetStack()->Last();
  TH1F* last = (TH1F*)last_hist->Clone("last");
  for(int ibin=3; ibin<=nbins;ibin++) {
    double stackcontent = last->GetBinContent(ibin);
    double stackerror = last->GetBinError(ibin);
    double datacontent = histo_j1EtaWidth_data_0->GetBinContent(ibin);
    double dataerror = histo_j1EtaWidth_data_0->GetBinError(ibin);
    std::cout<<"stackcontent: "<<stackcontent<<" and data content: "<<datacontent<<std::endl;
    double ratiocontent=0;
    if(datacontent!=0){
    	ratiocontent = ( datacontent) / stackcontent ;}
    double num = 1/datacontent;
    double den = 1/stackcontent;
    double error=0;
    if(datacontent!=0){
    	 error = ratiocontent*sqrt(pow((dataerror/datacontent),2) + pow((stackerror/stackcontent),2));}
    else {error = 2.07;}
    std::cout<<"ratio content: "<<ratiocontent<<" and error: "<<error<<std::endl;
    Ratio->SetBinContent(ibin,ratiocontent);
    Ratio->SetBinError(ibin,error);
  }  
  Ratio->GetYaxis()->SetRangeUser(0.0,2.2);
  //Ratio->Draw("e1p");
  Ratio->SetStats(0);
  //Ratio->GetXaxis()->SetTitle("pfMET [GeV]");

  //Ratio->GetYaxis()->SetTitle("Data/MC");

  Ratio->GetYaxis()->CenterTitle();
  Ratio->SetMarkerStyle(20);
  Ratio->SetMarkerSize(0.7);

  TLine *line = new TLine(160, 1.,2500, 1.);
  line->SetLineStyle(8);

  Ratio->GetYaxis()->SetLabelSize(0.14);
  Ratio->GetYaxis()->SetTitleSize(0.12);
  Ratio->GetYaxis()->SetLabelFont(42);
  Ratio->GetYaxis()->SetTitleFont(42);
  Ratio->GetYaxis()->SetTitleOffset(0.25);
  Ratio->GetYaxis()->SetNdivisions(100);
  Ratio->GetYaxis()->SetTickLength(0.05);

  Ratio->GetXaxis()->SetLabelSize(0.15);
  Ratio->GetXaxis()->SetTitleSize(0.12);
  Ratio->GetXaxis()->SetLabelFont(42);
  Ratio->GetXaxis()->SetTitleFont(42);
  Ratio->GetXaxis()->SetTitleOffset(0.90);
  Ratio->GetXaxis()->SetTickLength(0.05);
  Ratio->Draw("pex0");
  line->SetLineColor(kBlack);
  line->Draw("same");

  TGaxis *xaxis = new TGaxis(160,0,2500,0,160,2500,510);
  xaxis->SetTitle("P_{T} of Leading Jet (GeV)");
  xaxis->SetLabelFont(42);
  xaxis->SetLabelSize(0.10);
  xaxis->SetTitleFont(42);
  xaxis->SetTitleSize(0.12);
  xaxis->SetTitleOffset(1.2);
  xaxis->Draw("SAME");

  TGaxis *yaxis = new TGaxis(160,0,160,2.2,0,2.2,6,"");
  yaxis->SetTitle("Data/MC");
  yaxis->SetLabelFont(42);
  yaxis->SetLabelSize(0.10);
  yaxis->SetTitleFont(42);
  yaxis->SetTitleSize(0.12);
  yaxis->SetTitleOffset(0.35);
  yaxis->Draw("SAME");

  c->Update();
  //hs_datamc->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
  hs_datamc->GetYaxis()->SetTitle("Events");
  hs_datamc->GetYaxis()->SetTitleOffset(1.0);
  //hs_datamc->GetYaxis()->SetTitleOffset(1.3);  
  
/*  double xmin = c->GetUxmin();
  double ymin = c->GetUymin();
  //double ymin = 0.1;
  double xmax = c->GetUxmax();
  double ymax = c->GetUymax();
  //double ymax = pow(10,5);
  std::cout<<"xmin: "<<xmin<<" xmax: "<<xmax<<" ymin: "<<ymin<< " ymax: "<<ymax<<std::endl;
  TGaxis *xaxis = new TGaxis(xmin,ymin,xmax,ymin,xmin,xmax,510);
  xaxis->SetTitle("pfMET [GeV]");
  xaxis->SetLabelFont(42);
  xaxis->SetLabelSize(0.030);
  xaxis->SetTitleFont(42);
  xaxis->SetTitleSize(0.035);
  xaxis->Draw("SAME");

  TGaxis *xaxis_top = new TGaxis(xmin,ymax,xmax,ymax,xmin,xmax,510,"-");
  xaxis_top->SetTitle("");
  xaxis_top->SetLabelOffset(999);
  xaxis_top->Draw("SAME");

  TGaxis *yaxis = new TGaxis(xmin,ymin,xmin,ymax,ymin,ymax,50510,"G");
  yaxis->SetTitle("Events");
  yaxis->SetLabelFont(42);
  yaxis->SetLabelSize(0.030);
  yaxis->SetTitleFont(42);
  yaxis->SetTitleSize(0.035);
  yaxis->SetTitleOffset(1.8);
  yaxis->Draw("SAME");

  TGaxis *yaxis_right = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+");
  yaxis_right->SetTitle("");
  yaxis_right->SetLabelOffset(999);
  yaxis_right->Draw("SAME");  
 
*/
  c->SaveAs("datamc_j1pT_10.pdf");
  c->SaveAs("datamc_j1pT_10.png");
}
