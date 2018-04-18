//Create: ./rootcom plotter plotter
//Usage: ./plotter variable
//Example: ./plotter h_dileptonM_8
//X axis label is contained in samplename.txt
//To add new variables add a common name (i.e. dileptonM)
//Under the common name add what you want the X axis to be labeled (i.e. Dilepton Mass (GeV))
//Cutflow plots are handled with a special statment to add an extra TLatex

#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
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
#include "TLatex.h"
#include "TGaxis.h"
#include "iostream"
#include "stdlib.h"
#include "fstream"

std::vector<float> GetTotal(std::vector<TFile*> Files)
{
  std::vector<float> Total = {};
  for (int i = 0; i < Files.size(); i++)
    {
      TH1D *cutflow = (TH1D*)Files[i]->Get("h_cutflow");
      Total.push_back(cutflow->GetBinContent(1));
    }
  return Total;
}

std::string SampleName(const char * variable)
{
  std::string line;
  std::vector<std::string> lines;
  std::string name;
  std::ifstream file ("samplenames.txt");
  if (file.is_open())
    {
      while (getline(file,line)){lines.push_back(line);}

      for (int i = 0; i < lines.size(); i++)
	{
	  if (strstr(variable,lines[i].c_str()) != NULL)
	    {
	      name=lines[i+1];
	      break;
	    }
	}
    }
  else {std::cout << "samplenames.txt could not be found" << std::endl;}
  return name;
}

void hs_save(const char * variable, std::vector<TH1F*> hs_list)
{
  const char* hs_root = "../../Plots/Histogram.root";
  const char* region = "SignalRegion";
  std::ifstream file(hs_root);
  if (file){file.close();}
  else
    {
      TFile* hs_file = TFile::Open(hs_root,"RECREATE");
      hs_file->Close();
    }
  TFile* hs_file = TFile::Open(hs_root,"UPDATE");
  if (!(hs_file->GetDirectory(region)))
    {
      hs_file->mkdir(region);
    }
  hs_file->cd(region);
  for (int i = 0; i < hs_list.size(); i++)
    {
      hs_list[i]->Write();
    }
  hs_file->Close();
}
  
std::vector<int> hs_sort(std::vector<TH1F*> hs_list)
{
  std::vector<int> hs_index;
  std::vector<Double_t> hs_order;
  for (int i = 0; i < hs_list.size(); i++)
    {
      hs_order.push_back(hs_list[i]->Integral());
    }
  std::sort(hs_order.begin(),hs_order.end());
  for (int i = 0; hs_index.size() != hs_list.size(); i++)
    {
      for (int j = 0; j < hs_list.size(); j++)
        {
          float max = std::max(hs_order[i],hs_list[j]->Integral());
          if (fabs(hs_order[i] - hs_list[j]->Integral()) <= (FLT_EPSILON*max) && hs_list[j]->Integral() != 0)
            {
              hs_index.push_back(j);
              break;
            }
        }
    }
  return hs_index;
}

void plotter(const char * variable,std::string name)
{
  double lumi_1 = 3723.664;
  double lumi_2 = 35900.;

  std::cout << name << std::endl;
  std::ifstream file("postMETdata_final.root");
  if (file)
    {
      file.close();
    }
  else
    {
      system("hadd -f postMETdata_final.root postMETdata_{0..2}.root");
    }

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
  
  //opening the data file and adding "h_dileptonM_8" histogram
  TFile *f_datafile_0 = new TFile("postMETdata_final.root");
  //TFile *f_datafile_1 = new TFile("postMETdata_1.root");
  TH1F *histo_j1EtaWidth_data_0 = (TH1F*)f_datafile_0->Get(variable);
  //TH1F *histo_j1EtaWidth_data_1 = (TH1F*)f_datafile_1->Get(variable);
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
  std::vector<const char *> WJets_FileNames = {"postWJets_MLM_0.root","postW100to200_0.root","postW200to400_0.root","postW400to600_0.root","postW600to800_0.root","postW800to1200_0.root","postW1200to2500_0.root","postW2500toInf_0.root"};
  std::vector<TFile*> WJets_Files;
  for (int i = 0; i < WJets_FileNames.size(); i++){WJets_Files.push_back(new TFile(WJets_FileNames[i]));}			    
  std::vector<float> WJets_Total = GetTotal(WJets_Files);

  TH1F *histo_j1EtaWidth_WJets_0 = (TH1F*)WJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_W1Jets = (TH1F*)WJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_W2Jets = (TH1F*)WJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_W3Jets = (TH1F*)WJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_W4Jets = (TH1F*)WJets_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_W5Jets = (TH1F*)WJets_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_W6Jets = (TH1F*)WJets_Files[6]->Get(variable);
  TH1F *histo_j1EtaWidth_W7Jets = (TH1F*)WJets_Files[7]->Get(variable);
  
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
  histo_j1EtaWidth_WJets_0->Scale((1.0/WJets_Total[0])*35900*50690);
  histo_j1EtaWidth_W1Jets->Scale((1.0/WJets_Total[1])*35900*1345);
  histo_j1EtaWidth_W2Jets->Scale((1.0/WJets_Total[2])*35900*359.7);
  histo_j1EtaWidth_W3Jets->Scale((1.0/WJets_Total[3])*35900*48.91);
  histo_j1EtaWidth_W4Jets->Scale((1.0/WJets_Total[4])*35900*12.05);
  histo_j1EtaWidth_W5Jets->Scale((1.0/WJets_Total[5])*35900*5.501);
  histo_j1EtaWidth_W6Jets->Scale((1.0/WJets_Total[6])*35900*1.329);
  histo_j1EtaWidth_W7Jets->Scale((1.0/WJets_Total[7])*35900*0.03216);
  
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W1Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W2Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W3Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W4Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W5Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W6Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W7Jets);

  double integralWJets = histo_j1EtaWidth_WJets_0->Integral();
  std::cout<<"integral of WJets bkg here:"<<integralWJets<<std::endl;

  histo_j1EtaWidth_WJets_0->SetName(((std::string(variable))+(std::string("_WJets"))).c_str());
  histo_j1EtaWidth_WJets_0->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetTickLength(0);
  histo_j1EtaWidth_WJets_0->GetXaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetTitle("");
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetTickLength(0);
  histo_j1EtaWidth_WJets_0->GetYaxis()->SetLabelOffset(999);
  histo_j1EtaWidth_WJets_0->SetFillColor(kRed-10);

  std::vector<const char *> Zvv_FileNames = {"postZ100to200_0.root","postZ200to400_0.root","postZ400to600_0.root","postZ600to800_0.root","postZ800to1200_0.root","postZ1200to2500_0.root","postZ2500toInf_0.root"};
  std::vector<TFile *> Zvv_Files;
  for (int i = 0; i < Zvv_FileNames.size(); i++) {Zvv_Files.push_back(new TFile(Zvv_FileNames[i]));}
  std::vector<float> Zvv_Total = GetTotal(Zvv_Files);

  TH1F *histo_j1EtaWidth_100to200 = (TH1F*)Zvv_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_200to400 = (TH1F*)Zvv_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_400to600 = (TH1F*)Zvv_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_600to800 = (TH1F*)Zvv_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_800to1200 = (TH1F*)Zvv_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_1200to2500 = (TH1F*)Zvv_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_2500toInf = (TH1F*)Zvv_Files[6]->Get(variable);
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
  histo_j1EtaWidth_100to200->Scale((1.0/Zvv_Total[0])*35900*280.35);
  histo_j1EtaWidth_200to400->Scale((1.0/Zvv_Total[1])*35900*77.67);
  histo_j1EtaWidth_400to600->Scale((1.0/Zvv_Total[2])*35900*10.73);
  histo_j1EtaWidth_600to800->Scale((1.0/Zvv_Total[3])*35900*2.559);
  histo_j1EtaWidth_800to1200->Scale((1.0/Zvv_Total[4])*35900*1.1796);
  histo_j1EtaWidth_1200to2500->Scale((1.0/Zvv_Total[5])*35900*0.28833);
  histo_j1EtaWidth_2500toInf->Scale((1.0/Zvv_Total[6])*35900*0.006945);

  //Add the ZJetsToNuNu histograms to the first one
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_200to400);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_400to600);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_600to800);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_800to1200);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_1200to2500);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_2500toInf);

  histo_j1EtaWidth_100to200->SetName(((std::string(variable))+(std::string("_Zvv"))).c_str());
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
  std::vector<const char *> GJets_FileNames = {"postGJets40to100.root","postGJets100to200.root","postGJets200to400.root","postGJets400to600.root","postGJets600toInf.root"};
  std::vector<TFile *> GJets_Files;
  for (int i = 0; i < GJets_FileNames.size(); i++) {GJets_Files.push_back(new TFile(GJets_FileNames[i]));}
  std::vector<float> GJets_Total = GetTotal(GJets_Files);
  
  TH1F *histo_j1EtaWidth_G1Jets = (TH1F*)GJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_G2Jets = (TH1F*)GJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_G3Jets = (TH1F*)GJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_G4Jets = (TH1F*)GJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_G5Jets = (TH1F*)GJets_Files[4]->Get(variable);

  histo_j1EtaWidth_G1Jets->SetStats(0);
  histo_j1EtaWidth_G2Jets->SetStats(0);
  histo_j1EtaWidth_G3Jets->SetStats(0);
  histo_j1EtaWidth_G4Jets->SetStats(0);
  histo_j1EtaWidth_G5Jets->SetStats(0);

  double rawGJets = (histo_j1EtaWidth_G1Jets->Integral())+ (histo_j1EtaWidth_G2Jets->Integral())+ (histo_j1EtaWidth_G3Jets->Integral())+ (histo_j1EtaWidth_G4Jets->Integral())+ (histo_j1EtaWidth_G5Jets->Integral()); 
  std::cout<<"raw GJets bkg:"<<rawGJets<<std::endl;
  //Scaling
  histo_j1EtaWidth_G1Jets->Scale((1.0/GJets_Total[0])*35900*17420);
  histo_j1EtaWidth_G2Jets->Scale((1.0/GJets_Total[1])*35900*5391);
  histo_j1EtaWidth_G3Jets->Scale((1.0/GJets_Total[2])*35900*1168);
  histo_j1EtaWidth_G4Jets->Scale((1.0/GJets_Total[3])*35900*132.5);
  histo_j1EtaWidth_G5Jets->Scale((1.0/GJets_Total[4])*35900*44.05);
  
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G2Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G3Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G4Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G5Jets);

  histo_j1EtaWidth_G1Jets->SetName(((std::string(variable))+(std::string("_GJets"))).c_str());
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
  std::vector<const char *> DYJets_FileNames = {"postDY_MLM_0.root","postDY100to200.root","postDY200to400.root","postDY400to600.root","postDY600to800.root","postDY800to1200.root","postDY1200to2500.root","postDY2500toInf.root"};
  std::vector<TFile *> DYJets_Files;
  for (int i = 0; i < DYJets_FileNames.size(); i++) {DYJets_Files.push_back(new TFile(DYJets_FileNames[i]));}
  std::vector<float> DYJets_Total = GetTotal(DYJets_Files);

  TH1F *histo_j1EtaWidth_DY1Jets = (TH1F*)DYJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_DY2Jets = (TH1F*)DYJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_DY3Jets = (TH1F*)DYJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_DY4Jets = (TH1F*)DYJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_DY5Jets = (TH1F*)DYJets_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_DY6Jets = (TH1F*)DYJets_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_DY7Jets = (TH1F*)DYJets_Files[6]->Get(variable);
  TH1F *histo_j1EtaWidth_DY8Jets = (TH1F*)DYJets_Files[7]->Get(variable);

  histo_j1EtaWidth_DY1Jets->SetStats(0);
  histo_j1EtaWidth_DY2Jets->SetStats(0);
  histo_j1EtaWidth_DY3Jets->SetStats(0);
  histo_j1EtaWidth_DY4Jets->SetStats(0);
  histo_j1EtaWidth_DY5Jets->SetStats(0);
  histo_j1EtaWidth_DY6Jets->SetStats(0);
  histo_j1EtaWidth_DY7Jets->SetStats(0);
  histo_j1EtaWidth_DY8Jets->SetStats(0);

  double rawDY = (histo_j1EtaWidth_DY1Jets->Integral()) + (histo_j1EtaWidth_DY2Jets->Integral())+ (histo_j1EtaWidth_DY3Jets->Integral())+ (histo_j1EtaWidth_DY4Jets->Integral())+ (histo_j1EtaWidth_DY5Jets->Integral()) + (histo_j1EtaWidth_DY6Jets->Integral()) + (histo_j1EtaWidth_DY7Jets->Integral()) + (histo_j1EtaWidth_DY8Jets->Integral());  
  std::cout<<"raw DYJets bkg:"<<rawDY<<std::endl;

  //histo_j1EtaWidth_DY1Jets->Scale((1.0/96657400)*35900*4895);
  histo_j1EtaWidth_DY1Jets->Scale((1.0/DYJets_Total[0])*35900*4895);
  histo_j1EtaWidth_DY2Jets->Scale((1.0/DYJets_Total[1])*35900*148);
  histo_j1EtaWidth_DY3Jets->Scale((1.0/DYJets_Total[2])*35900*40.94);
  histo_j1EtaWidth_DY4Jets->Scale((1.0/DYJets_Total[3])*35900*5.497);
  histo_j1EtaWidth_DY5Jets->Scale((1.0/DYJets_Total[4])*35900*1.354);
  histo_j1EtaWidth_DY6Jets->Scale((1.0/DYJets_Total[5])*35900*0.6250);
  histo_j1EtaWidth_DY7Jets->Scale((1.0/DYJets_Total[6])*35900*0.1511);
  histo_j1EtaWidth_DY8Jets->Scale((1.0/DYJets_Total[7])*35900*0.003647);

  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY2Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY3Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY4Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY5Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY6Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY7Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY8Jets);
  
  histo_j1EtaWidth_DY1Jets->SetName(((std::string(variable))+(std::string("_DYJets"))).c_str());
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
  std::vector<const char *> TTJets_FileNames = {"postTTJets_MLM.root"};
  std::vector<TFile *> TTJets_Files;
  for(int i = 0; i < TTJets_FileNames.size(); i++) {TTJets_Files.push_back(new TFile(TTJets_FileNames[i]));}
  std::vector<float> TTJets_Total = GetTotal(TTJets_Files);

  TH1F *histo_j1EtaWidth_TTJets = (TH1F*)TTJets_Files[0]->Get(variable);
  histo_j1EtaWidth_TTJets->SetStats(0);

  double rawTTJets = histo_j1EtaWidth_TTJets->Integral();
  std::cout<<"raw TTJets here:"<<rawTTJets<<std::endl;

  histo_j1EtaWidth_TTJets->Scale((1.0/TTJets_Total[0])*35900*502.2);

  histo_j1EtaWidth_TTJets->SetName(((std::string(variable))+(std::string("_TTJets"))).c_str());
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

  std::vector<const char *> WW_FileNames = {"postWW.root"};
  std::vector<TFile *> WW_Files;
  for (int i = 0; i < WW_FileNames.size(); i++) {WW_Files.push_back(new TFile(WW_FileNames[i]));}
  std::vector<float> WW_Total = GetTotal(WW_Files);

  TH1F *histo_j1EtaWidth_WW = (TH1F*)WW_Files[0]->Get(variable);
  histo_j1EtaWidth_WW->SetStats(0);
  histo_j1EtaWidth_WW->Scale((1.0/WW_Total[0])*35900*118.7);
  histo_j1EtaWidth_WW->SetName(((std::string(variable))+(std::string("_WW"))).c_str());
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
 
  std::vector<const char *> WZ_FileNames = {"postWZ.root"};
  std::vector<TFile *> WZ_Files;
  for (int i = 0; i < WZ_FileNames.size(); i++) {WZ_Files.push_back(new TFile(WZ_FileNames[i]));}
  std::vector<float> WZ_Total = GetTotal(WW_Files);
  
  TH1F *histo_j1EtaWidth_WZ = (TH1F*)WZ_Files[0]->Get(variable);
  histo_j1EtaWidth_WZ->SetStats(0);
  histo_j1EtaWidth_WZ->Scale((1.0/WZ_Total[0])*35900*47.2);

  histo_j1EtaWidth_WZ->SetName(((std::string(variable))+(std::string("_WZ"))).c_str());
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

  std::vector<const char *> ZZ_FileNames = {"postZZ.root"};
  std::vector<TFile *> ZZ_Files;
  for (int i = 0; i < ZZ_FileNames.size(); i++) {ZZ_Files.push_back(new TFile(ZZ_FileNames[i]));}
  std::vector<float> ZZ_Total = GetTotal(ZZ_Files);

  TH1F *histo_j1EtaWidth_ZZ = (TH1F*)ZZ_Files[0]->Get(variable);
  histo_j1EtaWidth_ZZ->SetStats(0);
  histo_j1EtaWidth_ZZ->Scale((1.0/ZZ_Total[0])*35900*16.6);

  histo_j1EtaWidth_ZZ->SetName(((std::string(variable))+(std::string("_ZZ"))).c_str());
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
  
  histo_j1EtaWidth_WW->SetName(((std::string(variable))+(std::string("_WW_WZ_ZZ"))).c_str());
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

 //opening QCD background files (HT-binned samples)
  std::vector<const char *> QJets_FileNames = {"postQCD100to200_0.root","postQCD200to300_0.root","postQCD300to500_0.root","postQCD500to700_0.root","postQCD700to1000_0.root","postQCD1000to1500_0.root","postQCD1500to2000_0.root","postQCD2000toInf_0.root"};
  std::vector<TFile *> QJets_Files;
  for (int i = 0; i < QJets_FileNames.size(); i++) {QJets_Files.push_back(new TFile(QJets_FileNames[i]));}
  std::vector<float> QJets_Total = GetTotal(QJets_Files);

  TH1F *histo_j1EtaWidth_Q1Jets = (TH1F*)QJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_Q2Jets = (TH1F*)QJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_Q3Jets = (TH1F*)QJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_Q4Jets = (TH1F*)QJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_Q5Jets = (TH1F*)QJets_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_Q6Jets = (TH1F*)QJets_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_Q7Jets = (TH1F*)QJets_Files[6]->Get(variable);
  TH1F *histo_j1EtaWidth_Q8Jets = (TH1F*)QJets_Files[7]->Get(variable);

  histo_j1EtaWidth_Q1Jets->SetStats(0);
  histo_j1EtaWidth_Q2Jets->SetStats(0);
  histo_j1EtaWidth_Q3Jets->SetStats(0);
  histo_j1EtaWidth_Q4Jets->SetStats(0);
  histo_j1EtaWidth_Q5Jets->SetStats(0);
  histo_j1EtaWidth_Q6Jets->SetStats(0);
  histo_j1EtaWidth_Q7Jets->SetStats(0);
  histo_j1EtaWidth_Q8Jets->SetStats(0);

  histo_j1EtaWidth_Q1Jets->Scale((1.0/QJets_Total[0])*35900*27500000);
  histo_j1EtaWidth_Q2Jets->Scale((1.0/QJets_Total[1])*35900*1735000);
  histo_j1EtaWidth_Q3Jets->Scale((1.0/QJets_Total[2])*35900*367000);
  histo_j1EtaWidth_Q4Jets->Scale((1.0/QJets_Total[3])*35900*29370);
  histo_j1EtaWidth_Q5Jets->Scale((1.0/QJets_Total[4])*35900*6524);
  histo_j1EtaWidth_Q6Jets->Scale((1.0/QJets_Total[5])*35900*1064);
  histo_j1EtaWidth_Q7Jets->Scale((1.0/QJets_Total[6])*35900*121.5);
  histo_j1EtaWidth_Q8Jets->Scale((1.0/QJets_Total[7])*35900*25.42);

  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q2Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q3Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q4Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q5Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q6Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q7Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q8Jets);
  
  double integral = histo_j1EtaWidth_Q1Jets->Integral();
  std::cout<<"integral of QCD bkg here:"<<integral<<std::endl;

  histo_j1EtaWidth_Q1Jets->SetName(((std::string(variable))+(std::string("_QJets"))).c_str());
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
  std::vector<TH1F*> hs_list = {histo_j1EtaWidth_100to200,histo_j1EtaWidth_G1Jets,histo_j1EtaWidth_TTJets,histo_j1EtaWidth_Q1Jets,histo_j1EtaWidth_WW,histo_j1EtaWidth_DY1Jets,histo_j1EtaWidth_WJets_0};
  hs_save(variable,hs_list);
  std::vector<int> hs_index = hs_sort(hs_list);
  //hs_datamc->Add(histo_j1EtaWidth_WW);
  //hs_datamc->Add(histo_j1EtaWidth_100to200);
  //hs_datamc->Add(histo_j1EtaWidth_DY1Jets);
  //hs_datamc->Add(histo_j1EtaWidth_G1Jets);
  //hs_datamc->Add(histo_j1EtaWidth_TTJets);
  //hs_datamc->Add(histo_j1EtaWidth_WJets_0);
  //hs_datamc->Add(histo_j1EtaWidth_DY1Jets); 
  //hs_datamc->Add(histo_j1EtaWidth_Q1Jets);
  for (int i = 0; i < hs_list.size(); i++)
    {
      hs_datamc->Add(hs_list[hs_index[i]]);
    }
  //hs_datamc->Add(histo_j1EtaWidth_WJets_0);
  //hs_datamc->Add(histo_j1EtaWidth_G1Jets);
  hs_datamc->SetTitle("");
  hs_datamc->Draw("HIST");
  //hs_datamc->SetMinimum(0);
  //hs_datamc->SetMaximum(500);
  hs_datamc->SetMinimum(0.1);
  hs_datamc->SetMaximum(hs_datamc->GetMaximum()*pow(10,2.5));
  hs_datamc->Draw("HIST");
  histo_j1EtaWidth_data_0->SetLineColor(kBlack);
  histo_j1EtaWidth_data_0->SetMarkerStyle(20);
  histo_j1EtaWidth_data_0->SetMarkerSize(0.7);
  histo_j1EtaWidth_data_0->Draw("pex0same");
 
  //TFile *f_signal_1GeVfile = new TFile("postSignal_mchi1GeV.root");
  //TH1F *histo_signal_1GeV = (TH1F*)f_signal_1GeVfile ->Get(variable);
  //histo_signal_1GeV->Scale((1.0/629)*35900*0.056);
  //histo_signal_1GeV->SetLineColor(kRed);
  //histo_signal_1GeV->SetLineWidth(2);
  //histo_signal_1GeV->Draw("HIST SAME");

  //double integralsignal_1GeV = histo_signal_1GeV->Integral(); 
  //std::cout<<"integral of signal_1GeV here:"<<integralsignal_1GeV<<std::endl;

  TFile *f_signal_5GeVfile = new TFile("postSignal.root");
  TH1F *histo_signal_5GeV = (TH1F*)f_signal_5GeVfile ->Get(variable);
  histo_signal_5GeV->Scale((1.0/2133)*35900*0.047);
  histo_signal_5GeV->SetLineColor(kBlue);
  histo_signal_5GeV->SetLineWidth(2);
  histo_signal_5GeV->Draw("HIST SAME");

  double integralsignal_5GeV = histo_signal_5GeV->Integral(); 
  std::cout<<"integral of signal_5GeV here:"<<integralsignal_5GeV<<std::endl;
  std::vector<TH1F*> signal = {histo_signal_5GeV};

  //TFile *f_signal_10GeVfile = new TFile("postSignal_mchi10GeV.root");
  //TH1F *histo_signal_10GeV = (TH1F*)f_signal_10GeVfile ->Get(variable);
  //histo_signal_10GeV->Scale((1.0/4052)*35900*0.04);
  //histo_signal_10GeV->SetLineColor(kViolet+1);
  //histo_signal_10GeV->SetLineWidth(2);
  //histo_signal_10GeV->Draw("HIST SAME");

  //double integralsignal_10GeV = histo_signal_10GeV->Integral(); 
  //std::cout<<"integral of signal_10GeV here:"<<integralsignal_10GeV<<std::endl;

  //TFile *f_signal_20GeVfile = new TFile("postSignal_mchi20GeV.root");
  //TH1F *histo_signal_20GeV = (TH1F*)f_signal_20GeVfile ->Get(variable);
  //histo_signal_20GeV->Scale((1.0/9998)*35900*0.034);
  //histo_signal_20GeV->SetLineColor(kMagenta);
  //histo_signal_20GeV->SetLineWidth(2);
  //histo_signal_20GeV->Draw("HIST SAME");

  //double integralsignal_20GeV = histo_signal_20GeV->Integral(); 
  //std::cout<<"integral of signal_20GeV here:"<<integralsignal_20GeV<<std::endl;
 
  //TFile *f_signal_50GeVfile = new TFile("postSignal_mchi50GeV.root");
  //TH1F *histo_signal_50GeV = (TH1F*)f_signal_50GeVfile ->Get(variable);
  //histo_signal_50GeV->Scale((1.0/9999)*35900*0.025);
  //histo_signal_50GeV->SetLineColor(kSpring-1);
  //histo_signal_50GeV->SetLineWidth(2);
  //histo_signal_50GeV->Draw("HIST SAME");

  //double integralsignal_50GeV = histo_signal_50GeV->Integral(); 
  //std::cout<<"integral of signal_50GeV here:"<<integralsignal_50GeV<<std::endl;
////  
////  
  //TFile *f_signal_100GeVfile = new TFile("postSignal_mchi100GeV.root");
  //TH1F *histo_signal_100GeV = (TH1F*)f_signal_100GeVfile ->Get(variable);
  //histo_signal_100GeV->Scale((1.0/9994)*35900*0.019);
  //histo_signal_100GeV->SetLineColor(kAzure+1);
  //histo_signal_100GeV->SetLineWidth(2);
  //histo_signal_100GeV->Draw("HIST SAME");

  //double integralsignal_100GeV = histo_signal_100GeV->Integral(); 
  //std::cout<<"integral of signal_100GeV here:"<<integralsignal_100GeV<<std::endl;

  //TLegend *leg = new TLegend(0.181948,0.663948,0.567335,0.836868,"");
  TLegend *leg = new TLegend(0.62,0.60,0.86,0.887173,"");
  //leg->AddEntry(histo_j1EtaWidth_data_0,"Data");
  //leg->AddEntry(histo_signal_1GeV, "ZprimeSignal_mchi1GeV");
  leg->AddEntry(histo_signal_5GeV, "ZprimeSignal_mchi5GeV"); 
  //leg->AddEntry(histo_signal_10GeV, "ZprimeSignal_mchi10GeV");
  //leg->AddEntry(histo_signal_20GeV, "ZprimeSignal_mchi20GeV"); 
  //leg->AddEntry(histo_signal_50GeV, "ZprimeSignal_mchi50GeV");
  //leg->AddEntry(histo_signal_100GeV, "ZprimeSignal_mchi100GeV");  
  leg->AddEntry(histo_j1EtaWidth_WJets_0,"W#rightarrowl#nu","f");
  leg->AddEntry(histo_j1EtaWidth_DY1Jets,"Z#rightarrow ll","F"); 
  //leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu","F"); 
  //leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets", "F");
  //leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD","F");
  //leg->addentry(histo_j1etawidth_wjets_0,"w#rightarrowl#nu","f");
  //leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu","F"); 
  //leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD","F");
  //leg->AddEntry(histo_j1EtaWidth_TTJets, "Top Quark", "F");
  leg->AddEntry(histo_j1EtaWidth_WW,"WW/WZ/ZZ","F");
  leg->AddEntry(histo_j1EtaWidth_Q1Jets, "QCD","F");
  //leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets", "F");
  leg->AddEntry(histo_j1EtaWidth_TTJets, "Top Quark", "F"); 
  leg->AddEntry(histo_j1EtaWidth_G1Jets,"#gamma+jets", "F");
  //leg->AddEntry(histo_j1EtaWidth_DY1Jets,"DYJets#rightarrowLL","F");  
  leg->AddEntry(histo_j1EtaWidth_100to200,"Z#rightarrow#nu#nu","F");
  //leg->AddEntry(histo_j1EtaWidth_WW,"WW/WZ/ZZ","F");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.025);
  leg->Draw();  

  TLatex *texS = new TLatex(0.20,0.837173,"#sqrt{s} = 13 TeV, 35.9 fb^{-1}");
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
  for(int ibin=0; ibin<=nbins;ibin++) {
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

  TLine *line = new TLine(hs_datamc->GetXaxis()->GetXmin(), 1.,hs_datamc->GetXaxis()->GetXmax(), 1.);
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
  Ratio->GetXaxis()->SetTitleOffset(0.9);
  Ratio->GetXaxis()->SetTickLength(0.05);
  Ratio->Draw("pex0");
  line->SetLineColor(kBlack);
  line->Draw("same");

  double xmin = hs_datamc->GetXaxis()->GetXmin();
  double xmax = hs_datamc->GetXaxis()->GetXmax();
  double xwmin = xmin;
  double xwmax = xmax;

  TGaxis *xaxis = new TGaxis(xmin,0,xmax,0,xwmin,xwmax,510);
  xaxis->SetTitle(name.c_str());
  xaxis->SetLabelFont(42);
  xaxis->SetLabelSize(0.10);
  xaxis->SetTitleFont(42);
  xaxis->SetTitleSize(0.12);
  xaxis->SetTitleOffset(1.2);
  xaxis->Draw("SAME");

  if (name == "Cutflow")
    {
      xaxis->SetTitle("");
      for (int i = 1; i <= nbins; i++)
	{
	  TLatex *label = new TLatex(i-0.5,-0.3,hs_datamc->GetXaxis()->GetBinLabel(i));
	  label->SetTextSize(0.065);
	  label->SetTextAngle(-30.);
	  label->Draw("SAME");
	}
    }
      

  TGaxis *yaxis = new TGaxis(xmin,0,xmin,2.2,0,2.2,6,"");
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
  hs_datamc->GetYaxis()->SetTitleOffset(1.25);
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
  c->SaveAs((std::string("../../Plots/SignalRegionPlots_EWK/datamc_")+std::string(variable)+std::string(".pdf")).c_str());
  c->SaveAs((std::string("../../Plots/SignalRegionPlots_EWK/datamc_")+std::string(variable)+std::string(".png")).c_str());
}

int main(int argc, const char *argv[])
{
  std::vector<const char*> variable;
  for (int i = 1; i < argc; i++)
    {
      variable.push_back(argv[i]);
    }
  for (int i = 0; i < variable.size(); i++)
    {
      std::string name = SampleName(variable[i]);
      plotter(variable[i],name);
    } 
  return 0;
}
