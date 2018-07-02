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

TH1F* GetHistogram(std::vector<const char*> Sample_FileNames,std::vector<double> Sample_Xsec,const char* variable,std::string SampleName,double lumi)
{
  std::vector<TFile*> Sample_Files;
  for(int i = 0; i < Sample_FileNames.size(); i++){Sample_Files.push_back(new TFile(Sample_FileNames[i]));}
  std::vector<float> Sample_Total = GetTotal(Sample_Files);
  std::vector<TH1F*> Sample_Histo;
  double rawEvents = 0;
  for(int i = 0; i < Sample_FileNames.size(); i++)
    {
      Sample_Histo.push_back((TH1F*)Sample_Files[i]->Get(variable));
      Sample_Histo[i]->SetStats(0);
      rawEvents+=Sample_Histo[i]->Integral();
      //Scaling = (1/Totalevents)*Luminosity*NNLO-cross-section
      if (Sample_Xsec[i] > 0){Sample_Histo[i]->Scale((1.0/Sample_Total[i])*lumi*Sample_Xsec[i]);}
    }
  for(int i = 1; i < Sample_FileNames.size(); i++){Sample_Histo[0]->Add(Sample_Histo[i]);}
  Sample_Histo[0]->SetName(((std::string(variable)+std::string("_")+SampleName)).c_str());

  return Sample_Histo[0];
}

THStack* StackHistogram(std::vector<TH1F*> hs_list)
{
  THStack* hs_datamc = new THStack("hs_datamc","Data/MC comparison");
  std::vector<int> hs_index = hs_sort(hs_list);
  for (int i = 1; i < hs_list.size(); i ++){hs_datamc->Add(hs_list[hs_index[i]]);}
  hs_datamc->SetTitle("");
  hs_datamc->Draw("HIST");
  hs_datamc->SetMinimum(0.1);
  hs_datamc->SetMaximum(hs_datamc->GetMaximum()*pow(10,2.5));
  hs_datamc->Draw("HIST");
    
  return hs_datamc;
}

TH1F* GetSignal(int mchi,const char* variable, double lumi)
{
  if (mchi > 0)
    {
      if (mchi == 1)
	{
	  TFile *f_signal_1GeVfile = new TFile("postSignal_mchi1GeV.root");
	  TH1F *histo_signal_1GeV = (TH1F*)f_signal_1GeVfile ->Get(variable);
	  std::vector<TFile *> SignalFile = {f_signal_1GeVfile};
	  std::vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_1GeV->Scale((1.0/Signal_Totals[0])*lumi*0.056);
	  histo_signal_1GeV->SetName("ZprimeSignal_mchi1GeV");
	  histo_signal_1GeV->SetLineColor(kRed);
	  histo_signal_1GeV->SetLineWidth(2);
	  
	  return histo_signal_1GeV;
	}
      else if (mchi == 5)
	{
	  TFile *f_signal_5GeVfile = new TFile("postSignal.root");
	  TH1F *histo_signal_5GeV = (TH1F*)f_signal_5GeVfile ->Get(variable);
	  std::vector<TFile *> SignalFile = {f_signal_5GeVfile};
	  std::vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_5GeV->Scale((1.0/Signal_Totals[0])*lumi*0.047);
	  histo_signal_5GeV->SetName("ZprimeSignal_mchi5GeV");
	  histo_signal_5GeV->SetLineColor(kBlue);
	  histo_signal_5GeV->SetLineWidth(2);
	  
	  return histo_signal_5GeV;
	}
      else if (mchi == 10)
	{
	  TFile *f_signal_10GeVfile = new TFile("postSignal_mchi10GeV.root");
	  TH1F *histo_signal_10GeV = (TH1F*)f_signal_10GeVfile ->Get(variable);
	  std::vector<TFile *> SignalFile = {f_signal_10GeVfile};
	  std::vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_10GeV->Scale((1.0/Signal_Totals[0])*lumi*0.04);
	  histo_signal_10GeV->SetName("ZprimeSignal_mchi10GeV");
	  histo_signal_10GeV->SetLineColor(kViolet+1);
	  histo_signal_10GeV->SetLineWidth(2);
	  
	  return histo_signal_10GeV;
	}
      else if (mchi == 20)
	{
	  TFile *f_signal_20GeVfile = new TFile("postSignal_mchi20GeV.root");
	  TH1F *histo_signal_20GeV = (TH1F*)f_signal_20GeVfile ->Get(variable);
	  std::vector<TFile *> SignalFile = {f_signal_20GeVfile};
	  std::vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_20GeV->Scale((1.0/Signal_Totals[0])*lumi*0.034);
	  histo_signal_20GeV->SetName("ZprimeSignal_mchi20GeV");
	  histo_signal_20GeV->SetLineColor(kMagenta);
	  histo_signal_20GeV->SetLineWidth(2);
	  
	  return histo_signal_20GeV;
	}
      else if (mchi == 50)
	{
	  TFile *f_signal_50GeVfile = new TFile("postSignal_mchi50GeV.root");
	  TH1F *histo_signal_50GeV = (TH1F*)f_signal_50GeVfile ->Get(variable);
	  std::vector<TFile *> SignalFile = {f_signal_50GeVfile};
	  std::vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_50GeV->Scale((1.0/Signal_Totals[0])*lumi*0.025);
	  histo_signal_50GeV->SetName("ZprimeSignal_mchi50GeV");
	  histo_signal_50GeV->SetLineColor(kSpring-1);
	  histo_signal_50GeV->SetLineWidth(2);
	  
	  return histo_signal_50GeV;  
	}
      else if (mchi == 100)
	{
	  TFile *f_signal_100GeVfile = new TFile("postSignal_mchi100GeV.root");
	  TH1F *histo_signal_100GeV = (TH1F*)f_signal_100GeVfile ->Get(variable);
	  std::vector<TFile *> SignalFile = {f_signal_100GeVfile};
	  std::vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_100GeV->Scale((1.0/Signal_Totals[0])*lumi*0.019);
	  histo_signal_100GeV->SetName("ZprimeSignal_mchi100GeV");
	  histo_signal_100GeV->SetLineColor(kAzure+1);
	  histo_signal_100GeV->SetLineWidth(2);
	  
	  return histo_signal_100GeV;
	}
    }
  else
    {
      TH1F* noSignal = new TH1F();
      noSignal->SetName("NoSignal");
      return noSignal;
    }
}

void DrawRatio(TH1F* histo_Data, THStack* hs_datamc,TCanvas* c, std::string name)
{
  int nbins = histo_Data->GetNbinsX();  
  TH1F* Ratio = (TH1F*)histo_Data->Clone("Ratio");
  TH1F *last_hist = (TH1F*)hs_datamc->GetStack()->Last();
  TH1F* last = (TH1F*)last_hist->Clone("last");
  for(int ibin=0; ibin<=nbins;ibin++) {
    double stackcontent = last->GetBinContent(ibin);
    double stackerror = last->GetBinError(ibin);
    double datacontent = histo_Data->GetBinContent(ibin);
    double dataerror = histo_Data->GetBinError(ibin);
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
  Ratio->SetStats(0);
  
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
  hs_datamc->GetYaxis()->SetTitle("Events");
  hs_datamc->GetYaxis()->SetTitleOffset(1.25);
}
void plotter(const char * variable,std::string name,int mchi)
{
  double lumi_1 = 3723.664;
  double lumi_2 = 35900.;
  double lumi_3 = 1885.122;

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
  
  std::vector<const char*> Data_FileNames = {"postMETdata_final.root"};
  std::vector<double> Data_Xsec = {-1};
  TH1F* histo_Data = GetHistogram(Data_FileNames,Data_Xsec,variable,"data",lumi_3);
  histo_Data->SetLineWidth(2);
  histo_Data->SetLineColor(kWhite);
  histo_Data->SetTitle("");
  histo_Data->GetXaxis()->SetTitle("");
  histo_Data->GetXaxis()->SetTickLength(0);
  histo_Data->GetXaxis()->SetLabelOffset(999);
  histo_Data->GetYaxis()->SetTitle("");
  histo_Data->GetYaxis()->SetTickLength(0);
  histo_Data->GetYaxis()->SetLabelOffset(999);

  //opening background WJets Sample file
  std::vector<const char *> WJets_FileNames = {"postWJets_MLM_0.root","postW100to200_0.root","postW200to400_0.root","postW400to600_0.root","postW600to800_0.root","postW800to1200_0.root","postW1200to2500_0.root","postW2500toInf_0.root"};
  std::vector<double> WJets_Xsec = {50690,1345,359.7,48.91,12.05,5.501,1.329,0.03216};
  TH1F* histo_WJets = GetHistogram(WJets_FileNames,WJets_Xsec,variable,"WJets",lumi_3);

  histo_WJets->SetTitle("");
  histo_WJets->GetXaxis()->SetTitle("");
  histo_WJets->GetXaxis()->SetTickLength(0);
  histo_WJets->GetXaxis()->SetLabelOffset(999);
  histo_WJets->GetYaxis()->SetTitle("");
  histo_WJets->GetYaxis()->SetTickLength(0);
  histo_WJets->GetYaxis()->SetLabelOffset(999);
  histo_WJets->SetFillColor(kRed-10);

  std::cout<<"integral of WJets bkg here:"<<histo_WJets->Integral()<<std::endl;

  std::vector<const char *> ZJets_FileNames = {"postZ100to200_0.root","postZ200to400_0.root","postZ400to600_0.root","postZ600to800_0.root","postZ800to1200_0.root","postZ1200to2500_0.root","postZ2500toInf_0.root"};
  std::vector<double> ZJets_Xsec = {280.35,77.67,10.73,2.559,1.1796,0.28833,0.006945};
  TH1F* histo_ZJets = GetHistogram(ZJets_FileNames,ZJets_Xsec,variable,"ZJets",lumi_3);

  histo_ZJets->SetTitle("");
  histo_ZJets->GetXaxis()->SetTitle("");
  histo_ZJets->GetXaxis()->SetTickLength(0);
  histo_ZJets->GetXaxis()->SetLabelOffset(999);
  histo_ZJets->GetYaxis()->SetTitle("");
  histo_ZJets->GetYaxis()->SetTickLength(0);
  histo_ZJets->GetYaxis()->SetLabelOffset(999);  
  histo_ZJets->SetFillColor(kAzure+10);

  std::cout<<"integral of ZvvJets bkg here:"<<histo_ZJets->Integral()<<std::endl;

  //opening background samples Gamma+jets
  std::vector<const char *> GJets_FileNames = {"postGJets40to100.root","postGJets100to200.root","postGJets200to400.root","postGJets400to600.root","postGJets600toInf.root"};
  std::vector<double> GJets_Xsec = {17420,5391,1168,132.5,44.05};
  TH1F* histo_GJets = GetHistogram(GJets_FileNames,GJets_Xsec,variable,"GJets",lumi_3);

  histo_GJets->SetTitle("");
  histo_GJets->GetXaxis()->SetTitle("");
  histo_GJets->GetXaxis()->SetTickLength(0);
  histo_GJets->GetXaxis()->SetLabelOffset(999);
  histo_GJets->GetYaxis()->SetTitle("");
  histo_GJets->GetYaxis()->SetTickLength(0);
  histo_GJets->GetYaxis()->SetLabelOffset(999);
  histo_GJets->SetFillColor(kTeal-9);

  std::cout<<"integral of GJets bkg here:"<<histo_GJets->Integral()<<std::endl;

//opening DYJetsToLL backgrounds
  std::vector<const char *> DYJets_FileNames = {"postDY_MLM_0.root","postDY100to200.root","postDY200to400.root","postDY400to600.root","postDY600to800.root","postDY800to1200.root","postDY1200to2500.root","postDY2500toInf.root"};
  std::vector<double> DYJets_Xsec = {4895,148,40.94,5.497,1.354,0.6250,0.1511,0.003647};
  TH1F* histo_DYJets = GetHistogram(DYJets_FileNames,DYJets_Xsec,variable,"DYJets",lumi_3);

  histo_DYJets->SetTitle("");
  histo_DYJets->GetXaxis()->SetTitle("");
  histo_DYJets->GetXaxis()->SetTickLength(0);
  histo_DYJets->GetXaxis()->SetLabelOffset(999);
  histo_DYJets->GetYaxis()->SetTitle("");
  histo_DYJets->GetYaxis()->SetTickLength(0);
  histo_DYJets->GetYaxis()->SetLabelOffset(999);
  histo_DYJets->SetFillColor(kGray+2);

  std::cout<<"integral of DYJets bkg here:"<<histo_DYJets->Integral()<<std::endl;

  //opening background TTJets
  std::vector<const char *> TTJets_FileNames = {"postTTJets_MLM.root"};
  std::vector<double> TTJets_Xsec = {502.2};
  TH1F* histo_TTJets = GetHistogram(TTJets_FileNames,TTJets_Xsec,variable,"TTJets",lumi_3);

  histo_TTJets->SetTitle("");
  histo_TTJets->GetXaxis()->SetTitle("");
  histo_TTJets->GetXaxis()->SetTickLength(0);
  histo_TTJets->GetXaxis()->SetLabelOffset(999);
  histo_TTJets->GetYaxis()->SetTitle("");
  histo_TTJets->GetYaxis()->SetTickLength(0);
  histo_TTJets->GetYaxis()->SetLabelOffset(999);
  histo_TTJets->SetFillColor(kOrange-2);
 
  double integralTTJets = histo_TTJets->Integral();
  std::cout<<"integral of integralTTJets bkg here:"<<integralTTJets<<std::endl;

 //addding some backgrounds like WW, WZ, ZZ	

  std::vector<const char *> DiBoson_FileNames = {"postWW.root","postWZ.root","postZZ.root"};
  std::vector<double> DiBoson_Xsec = {118.7,47.2,16.6};
  TH1F* histo_DiBoson = GetHistogram(DiBoson_FileNames,DiBoson_Xsec,variable,"Diboson",lumi_3);

  histo_DiBoson->SetTitle("");
  histo_DiBoson->GetXaxis()->SetTitle("");
  histo_DiBoson->GetXaxis()->SetTickLength(0);
  histo_DiBoson->GetXaxis()->SetLabelOffset(999);
  histo_DiBoson->GetYaxis()->SetTitle("");
  histo_DiBoson->GetYaxis()->SetTickLength(0);
  histo_DiBoson->GetYaxis()->SetLabelOffset(999);
  histo_DiBoson->SetFillColor(kCyan-10);

  std::cout<<"integral of DiBoson bkg here:"<<histo_DiBoson->Integral()<<std::endl;

 //opening QCD background files (HT-binned samples)
  std::vector<const char *> QCD_FileNames = {"postQCD100to200_0.root","postQCD200to300_0.root","postQCD300to500_0.root","postQCD500to700_0.root","postQCD700to1000_0.root","postQCD1000to1500_0.root","postQCD1500to2000_0.root","postQCD2000toInf_0.root"};
  std::vector<double> QCD_Xsec = {27500000,1735000,367000,29370,6524,1064,121.5,25.42};
  std::vector<TFile *> QJets_Files;
  TH1F* histo_QCD = GetHistogram(QCD_FileNames,QCD_Xsec,variable,"QCD",lumi_3);

  histo_QCD->SetTitle("");
  histo_QCD->GetXaxis()->SetTitle("");
  histo_QCD->GetXaxis()->SetTickLength(0);
  histo_QCD->GetXaxis()->SetLabelOffset(999);
  histo_QCD->GetYaxis()->SetTitle("");
  histo_QCD->GetYaxis()->SetTickLength(0);
  histo_QCD->GetYaxis()->SetLabelOffset(999);
  histo_QCD->SetFillColor(kGray);

  //Stack histograms using THStack
  std::vector<TH1F*> hs_list = {histo_ZJets,histo_GJets,histo_TTJets,histo_QCD,histo_DiBoson,histo_DYJets,histo_WJets};
  hs_save(variable,hs_list);
  std::vector<int> hs_index = hs_sort(hs_list);
  THStack* hs_datamc = StackHistogram(hs_list);
  histo_Data->SetLineColor(kBlack);
  histo_Data->SetMarkerStyle(20);
  histo_Data->SetMarkerSize(0.7);
  histo_Data->Draw("pex0same");

  TH1F* histo_Signal = GetSignal(mchi,variable,lumi_3);
  if (mchi > 0){histo_Signal->Draw("HIST");}
 
  TLegend *leg = new TLegend(0.62,0.60,0.86,0.887173,"");
  if (mchi > 0){leg->AddEntry(histo_Signal, histo_Signal->GetName());};   
  leg->AddEntry(histo_WJets,"W#rightarrowl#nu","f");
  leg->AddEntry(histo_DYJets,"Z#rightarrow ll","F"); 
  leg->AddEntry(histo_DiBoson,"WW/WZ/ZZ","F");
  leg->AddEntry(histo_QCD, "QCD","F");
  leg->AddEntry(histo_TTJets, "Top Quark", "F"); 
  leg->AddEntry(histo_GJets,"#gamma+jets", "F"); 
  leg->AddEntry(histo_ZJets,"Z#rightarrow#nu#nu","F");
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
 
  DrawRatio(histo_Data,hs_datamc,c,name);

  c->SaveAs((std::string("../../Plots/SignalRegionPlots_EWK/datamc_")+std::string(variable)+std::string(".pdf")).c_str());
  c->SaveAs((std::string("../../Plots/SignalRegionPlots_EWK/datamc_")+std::string(variable)+std::string(".png")).c_str());
}

int main(int argc, const char *argv[])
{
  int mchi = atof(argv[1]);
  std::vector<const char*> variable;
  for (int i = 2; i < argc; i++)
    {
      variable.push_back(argv[i]);
    }
  for (int i = 0; i < variable.size(); i++)
    {
      std::string name = SampleName(variable[i]);
      plotter(variable[i],name,mchi);
    } 
  return 0;
}
