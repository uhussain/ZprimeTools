//Create: ./rootcom plotter plotter
//Usage: ./plotter variable
//Example: ./plotter h_dileptonM_8
//X axis label is contained in samplename.txt
//To add new variables add a common name (i.e. dileptonM)
//Under the common name add what you want the X axis to be labeled (i.e. Dilepton Mass (GeV))
//Cutflow plots are handled with a special statment to add an extra TLatex
#define PlotTool_cxx
#include "PlotTool.h"
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
#include "TColor.h"
#include "iostream"
#include "stdlib.h"
#include "fstream"
#include "string"
using namespace std;

int main(int argc, const char* argv[])
{
  if (string(argv[1]).compare("-h") == 0)
    {
      cout<<"usage: ./PlotTool option variable\n";
      cout<<"|Options|\n";
      cout<<"\tnone : Plots variables specified in command line argument (i.e. ./PlotTool h_dileptonM_8 h_cutflow)\n\t\tWhen plotting in Signal Region specify the mass of the Z' sample you would like to use as the first command line argument\n\t\t(i.e. ./PlotTool 5GeV h_cutflow pfMET_8)\n";
      cout<<"\n\t-i : Prints the events in all the categories (NoCat, Cat1, Cat2, and Cat3)\n\t\tUse secondary option 0 to print only Data and Monte Carlo events\n\t\tUse secondary option 1 to print all Z' mass events\n\t\tUse secondary option -1 to print all Data, Monte Carlo, and all Z' mass events\n\t\tUsage: ./PlotTool -i 1\n";
      cout<<"\n\t-s : Saves plots specified in command line argument to a root file for use in shape analysis\n\t\tSpecify the Z' mass and the plots to save to root file (i.e. ./PlotTool 1GeV pfMET_15 pfMET_31 pfMET_44) example given for NoCat\n\t\tUse the secondary option pdf to save the PDF uncertainty plots to all the created root files from the previous option (i.e. ./PlotTool -s pdf)\n"; 
    }
  else
    {
      PlotTool t(argc,argv);
      t.Options(argc,argv);
    }
  return 0;
}

void PlotTool::Options(int argc,const char* argv[])
{
  double lumi_1 = 3723.664;
  double lumi_2 = 35900.;
  double lumi_3 = 1885.122;

  lumi=lumi_2;
  
  cout<<"Running in "<<Region<<endl;
  if (Region == "SignalRegion"){lumi=lumi_3;}
  regionFile = haddRegion();
  vector<const char*> realargv;
  if (string(argv[1]).compare("-i") == 0)
    {
      cout<<"Integrals"<<endl;
      integralOption(atof(argv[2]));
    }
  else if (string(argv[1]).compare("-s") == 0)
    {
      if (Region != "SignalRegion")
	{
	  cout<<"Please only run ./PlotTool -s in SignalRegion"<<endl;
	}
      else
	{
	  cout<<"Saving Uncertainty Plots"<<endl;
	  for (int j = 2; j < argc; j++){ realargv.push_back(argv[j]);}
	  saveplotOption(realargv.size(),realargv);
	}
    }
  else if (string(argv[1]).compare("-th2f") == 0)
    {
      if (Region != "SignalRegion")
	{
	  cout<<"Please only run ./PlotTool -th2f in PFUncertainty"<<endl;
	}
      cout<<"Making TH2F Uncertianty Plots"<<endl;
      for (int j = 1; j < argc; j++){ realargv.push_back(argv[j]);}
      plotTH2FOption(realargv.size(),realargv);
    }
  else
    {
      cout<<"Plotting at "<<lumi<<" pb^{-1}"<<endl;
      for (int j = 1; j < argc; j++){realargv.push_back(argv[j]);}
      plotterOption(realargv.size(),realargv);
    }
}
  
void PlotTool::plotTH2FOption(int argc,vector<const char*> argv)
{
  string mchi = "-1";
  vector<const char*> variable;
  for (int i = 1; i < argc; i++){variable.push_back(argv[i]);}
  for (int i = 0; i < variable.size(); i++)
    {
      string name = SampleName(variable[i]);
      plotTH2F(variable[i],name,mchi);
    }
}

void PlotTool::saveplotOption(int argc,vector<const char*> argv)
{
  string mchi = string(argv[1]);
  vector<const char*> variable;
  if (mchi != "pdf")
    {
      for (int i = 2; i < argc; i++)
	{
	  variable.push_back(argv[i]);
	}
      for (int i = 0; i < variable.size(); i++)
	{
	  string name = SampleName(variable[i]);
	  vector<string> label = GetName(variable[i]);
	  string varname = label[0];
	  string cat = label[1];
	  saveplot(variable[i],name,varname,cat,mchi);
	} 
    }
  else if (mchi == "pdf")
    {
      cout<<"Saving PDF"<<endl;
      savepdf();
    }
}

void PlotTool::integralOption(int print)
{
  vector<const char*> variable;
  if (regionFile == "postMETdata_final"){variable = {"pfMET_15","pfMET_10","pfMET_12","pfMET_14"};}
  else{variable = {"nJets_8"};}//,"nJets_10","nJets_12","nJets_14"};}
  vector<string> type = {"NoCat","Cat1","Cat2","Cat3"};
  for (int i = 0; i < variable.size(); i++)
    {
      string name = type[i];
      integral(variable[i],name,"-1",print);
    } 
}

void PlotTool::plotterOption(int argc, vector<const char*> argv)
{
  string mchi;
  int argIter;
  if (Region == "SignalRegion")
    {
      mchi = string(argv[0]);
      argIter = 1;
    }
  else
    {
      mchi = "0";
      argIter = 0;
    }
  vector<const char*> variable;
  for (argIter; argIter < argc; argIter++){variable.push_back(argv[argIter]);}
  for (int i = 0; i < variable.size(); i++)
    {
      string name = SampleName(variable[i]);
      plotter(variable[i],name,mchi);
    } 
}

void PlotTool::haddAll(vector<const char*> Data_FileNames)
{
  vector<vector<const char*>> MC_FileNames{Data_FileNames,WJets_FileNames,ZJets_FileNames,GJets_FileNames,DYJets_FileNames,DiBoson_FileNames,TTJets_FileNames,QCD_FileNames};
  for (int i = 0; i < MC_FileNames.size(); i++)
    {
      for (int j = 0; j< MC_FileNames[i].size(); j++)
	{
	  ifstream curFile(string(MC_FileNames[i][j])+".root");
	  if (curFile.is_open())
	    {
	      curFile.close();
	    }
	  else
	    {
	      int nfile;
	      system(("ls -f .output/"+string(MC_FileNames[i][j])+"_* | wc -l > .filelist/"+string(MC_FileNames[i][j])).c_str());
	      ifstream mcFile(".filelist/"+string(MC_FileNames[i][j]));
	      if (mcFile.is_open())
		{
		  string line;
		  while(getline(mcFile,line)) nfile = stoi(line);
		  system(("hadd -f "+string(MC_FileNames[i][j])+".root .output/"+string(MC_FileNames[i][j])+"_{1.."+to_string(nfile)+"}.root").c_str());
		}
	      else cout<<"Unable to find"+string(MC_FileNames[i][j])<<endl;
	    }
	}
    }
}

const char* PlotTool::haddRegion()
{
  const char* regionFile;
  if (Region == "SignalRegion")
    {
      haddAll(SignalData_FileNames);
      ifstream file("postMETdata_final.root");
      if (file)
	{
	  file.close();
	}
      else
	{
	  system("hadd -f postMETdata_final.root postMETdata_{0..2}.root");
	}
      regionFile="postMETdata_final";
    }
  else if (Region == "SingleEleCR")
    {
      haddAll(SingleEleData_FileNames);
      ifstream file("postSingleEle_final.root");
      if (file)
	{
	  file.close();
	}
      else
	{
	  system("hadd -f postSingleEle_final.root postSingleEle_{0..19}.root");
	}
      regionFile="postSingleEle_final";
    }
  else if (Region == "SingleMuCR")
    {
      haddAll(SingleMuData_FileNames);
      ifstream file("postSingleMu_final.root");
      if (file)
	{
	  file.close();
	}
      else
	{
	  system("hadd -f postSingleMu_final.root postSingleMu_{0..19}.root");
	}
      regionFile="postSingleMu_final";
    }
  else if (Region == "DoubleEleCR")
    {
      haddAll(DoubleEleData_FileNames);
      ifstream file("postDoubleEle_final.root");
      if (file)
	{
	  file.close();
	} 
      else
	{
	  system("hadd -f postDoubleEle_final.root postDoubleEle_{0..19}.root");
	}
      regionFile="postDoubleEle_final";
    }
  else if (Region == "DoubleMuCR")
    {
      haddAll(DoubleMuData_FileNames);
      ifstream file("postDoubleMu_final.root");
      if (file)
	{
	  file.close();
	}
      else
	{
	  system("hadd -f postDoubleMu_final.root postDoubleMu_{0..19}.root");
	}
      regionFile="postDoubleMu_final";
    }
  
  return regionFile;
}

vector<float> PlotTool::GetTotal(vector<TFile*> Files)
{
  vector<float> Total = {};
  for (int i = 0; i < Files.size(); i++)
    {
      TH1D *cutflow = (TH1D*)Files[i]->Get("h_cutflow");
      Total.push_back(cutflow->GetBinContent(1));
    }
  return Total;
}

string PlotTool::SampleName(const char * variable)
{
  string line;
  vector<string> lines;
  string name;
  ifstream file ("../samplenames.txt");
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
  else {cout << "../samplenames.txt could not be found" << endl;}
  return name;
}

string PlotTool::GetCategory(int hs_num)
{
  string cat;
  if (hs_num == 15) cat = "_NoCat";
  else if (hs_num == 10) cat = "_Cat1";
  else if (hs_num == 12) cat = "_Cat2";
  else if (hs_num == 14) cat = "_Cat3";

  return cat;
}

vector<string> PlotTool::GetName(const char * variable)
{
  string name = string(variable);
  name.erase(name.begin(),name.begin()+6);
  int n = stoi(name);
  string cat;
  if (8<=n && n<=15)
    {
      cat = GetCategory(n);
      name="";
    }
  else if (24<=n && n<=31)
    {
      cat = GetCategory(n-16);
      name="_jesUp";
    }
  else if (37<=n && n<=44)
    {
      cat = GetCategory(n-29);
      name="_jesDown";
    }
  vector<string> label = {name,cat};
  return label;
}

int PlotTool::hs_save(string mchi,string cat, const char * variable, vector<TH1F*> histo)
{
  cout<<"Category"<<cat<<endl;
  const char* hs_root = (string("../../Systematics")+mchi+cat+string(".root")).c_str();
  cout<<"Saving to "<<hs_root<<endl;
  TFile* hs_file = TFile::Open(hs_root,"UPDATE");
  for (int i = 0; i < histo.size(); i++)
    {
      if (gDirectory->GetListOfKeys()->Contains(histo[i]->GetName()))
	{
	  gDirectory->Delete((string(histo[i]->GetName())+";1").c_str());
	}
      histo[i]->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
      histo[i]->GetYaxis()->SetTitle("Events");
      histo[i]->SetTitle("");
      histo[i]->Write();
    }
  cout<<"Saving to "<<hs_file->GetName()<<endl;
  hs_file->Close();
}
  
vector<int> PlotTool::hs_sort(vector<TH1F*> hs_list)
{
  vector<int> hs_index;
  vector<Double_t> hs_order;
  for (int i = 0; i < hs_list.size(); i++)
    {
      hs_order.push_back(hs_list[i]->Integral());
    }
  sort(hs_order.begin(),hs_order.end());
  for (int i = 0; hs_index.size() != hs_list.size(); i++)
    {
      for (int j = 0; j < hs_list.size(); j++)
        {
          float maxval = max(hs_order[i],hs_list[j]->Integral());
          if (fabs(hs_order[i] - hs_list[j]->Integral()) <= (FLT_EPSILON*maxval) && hs_list[j]->Integral() != 0)
            {
              hs_index.push_back(j);
              break;
            }
        }
    }
  return hs_index;
}

TH1F* PlotTool::GetHistogram(vector<const char*> Sample_FileNames,vector<double> Sample_Xsec,const char* variable,string SampleName)
{
  vector<TFile*> Sample_Files;
  for(int i = 0; i < Sample_FileNames.size(); i++){Sample_Files.push_back(new TFile((string(Sample_FileNames[i])+string(".root")).c_str()));}
  vector<float> Sample_Total = GetTotal(Sample_Files);
  vector<TH1F*> Sample_Histo;
  double rawEvents = 0;
  for(int i = 0; i < Sample_FileNames.size(); i++)
    {
      Sample_Histo.push_back((TH1F*)Sample_Files[i]->Get(variable));
      Sample_Histo[i]->SetStats(0);
      rawEvents+=Sample_Histo[i]->Integral();
      //Scaling = (1/Totalevents)*Luminosity*NNLO-cross-section
      if (Sample_Xsec[i] > 0)
	{
	  double scaleEvents = (Sample_Histo[i]->Integral()/Sample_Total[i])*lumi*Sample_Xsec[i];
	  //cout<<lumi<<"\t"<<Sample_Files[i]->GetName()<<" scaled events: "<<scaleEvents<<endl;
	  Sample_Histo[i]->Scale((1.0/Sample_Total[i])*lumi*Sample_Xsec[i]);
	}
    }
  for(int i = 1; i < Sample_FileNames.size(); i++){Sample_Histo[0]->Add(Sample_Histo[i]);}
  Sample_Histo[0]->SetName(((string(variable)+string("_")+SampleName)).c_str());

  return Sample_Histo[0];
}

TH2F* PlotTool::GetTH2F(vector<const char*> Sample_FileNames,vector<double> Sample_Xsec,const char* variable,string SampleName)
{
  vector<TFile*> Sample_Files;
  for(int i = 0; i < Sample_FileNames.size(); i++){Sample_Files.push_back(new TFile((string(Sample_FileNames[i])+string(".root")).c_str()));}
  vector<float> Sample_Total = GetTotal(Sample_Files);
  vector<TH2F*> Sample_Histo;
  double rawEvents = 0;
  for(int i = 0; i < Sample_FileNames.size(); i++)
    {
      Sample_Histo.push_back((TH2F*)Sample_Files[i]->Get(variable));
      Sample_Histo[i]->SetStats(0);
      rawEvents+=Sample_Histo[i]->Integral();
      //Scaling = (1/Totalevents)*Luminosity*NNLO-cross-section
      if (Sample_Xsec[i] > 0)
	{
	  double scaleEvents = (Sample_Histo[i]->Integral()/Sample_Total[i])*lumi*Sample_Xsec[i];
	  //cout<<lumi<<"\t"<<Sample_Files[i]->GetName()<<" scaled events: "<<scaleEvents<<endl;
	  Sample_Histo[i]->Scale((1.0/Sample_Total[i])*lumi*Sample_Xsec[i]);
	}
    }
  for(int i = 1; i < Sample_FileNames.size(); i++){Sample_Histo[0]->Add(Sample_Histo[i]);}
  Sample_Histo[0]->SetName(((string(variable)+string("_")+SampleName)).c_str());

  return Sample_Histo[0];
}

THStack* PlotTool::StackHistogram(vector<TH1F*> hs_list,string name)
{
  THStack* hs_datamc = new THStack("hs_datamc","Data/MC comparison");
  vector<int> hs_index = hs_sort(hs_list);
  for (int i = 0; i < hs_list.size(); i ++){hs_datamc->Add(hs_list[hs_index[i]]);}
  hs_datamc->SetTitle("");
  double min,max;
  if (name != "Z Mass (GeV)")
    {
      min=0.1;
      max=pow(10,2.5);
    }
  else
    {
      min = 0;
      max = 1.3;
    }
  hs_datamc->SetMinimum(min);
  hs_datamc->SetMaximum(hs_datamc->GetMaximum()*max);
    
  return hs_datamc;
}

vector<TH1F*> PlotTool::GetSignal(string mchi,const char* variable)
{
  vector<TH1F*> histo_Signal;
  if (mchi != "0")
    {
      if (mchi == "-1" || mchi == "1GeV")
	{
	  TFile *f_signal_1GeVfile = new TFile("postSignal_mchi1GeV.root");
	  TH1F *histo_signal_1GeV = (TH1F*)f_signal_1GeVfile ->Get(variable);
	  vector<TFile *> SignalFile = {f_signal_1GeVfile};
	  vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_1GeV->Scale((1.0/Signal_Totals[0])*lumi*0.056);
	  histo_signal_1GeV->SetName("ZprimeSignal_mchi1GeV");
	  histo_signal_1GeV->SetLineColor(kRed);
	  histo_signal_1GeV->SetLineWidth(2);
	  
	  histo_Signal.push_back(histo_signal_1GeV);
	}
       if (mchi == "-1" || mchi == "5GeV")
	{
	  TFile *f_signal_5GeVfile = new TFile("postSignal.root");
	  TH1F *histo_signal_5GeV = (TH1F*)f_signal_5GeVfile ->Get(variable);
	  vector<TFile *> SignalFile = {f_signal_5GeVfile};
	  vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_5GeV->Scale((1.0/Signal_Totals[0])*lumi*0.047);
	  histo_signal_5GeV->SetName("ZprimeSignal_mchi5GeV");
	  histo_signal_5GeV->SetLineColor(kBlue);
	  histo_signal_5GeV->SetLineWidth(2);
	  
	  histo_Signal.push_back(histo_signal_5GeV);
	}
       if (mchi == "-1" || mchi == "10GeV")
	{
	  TFile *f_signal_10GeVfile = new TFile("postSignal_mchi10GeV.root");
	  TH1F *histo_signal_10GeV = (TH1F*)f_signal_10GeVfile ->Get(variable);
	  vector<TFile *> SignalFile = {f_signal_10GeVfile};
	  vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_10GeV->Scale((1.0/Signal_Totals[0])*lumi*0.04);
	  histo_signal_10GeV->SetName("ZprimeSignal_mchi10GeV");
	  histo_signal_10GeV->SetLineColor(kViolet+1);
	  histo_signal_10GeV->SetLineWidth(2);
	  
	  histo_Signal.push_back(histo_signal_10GeV);
	}
       if (mchi == "-1" || mchi == "20GeV")
	{
	  TFile *f_signal_20GeVfile = new TFile("postSignal_mchi20GeV.root");
	  TH1F *histo_signal_20GeV = (TH1F*)f_signal_20GeVfile ->Get(variable);
	  vector<TFile *> SignalFile = {f_signal_20GeVfile};
	  vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_20GeV->Scale((1.0/Signal_Totals[0])*lumi*0.034);
	  histo_signal_20GeV->SetName("ZprimeSignal_mchi20GeV");
	  histo_signal_20GeV->SetLineColor(kMagenta);
	  histo_signal_20GeV->SetLineWidth(2);
	  
	  histo_Signal.push_back(histo_signal_20GeV);
	}
       if (mchi == "-1" || mchi == "50GeV")
	{
	  TFile *f_signal_50GeVfile = new TFile("postSignal_mchi50GeV.root");
	  TH1F *histo_signal_50GeV = (TH1F*)f_signal_50GeVfile ->Get(variable);
	  vector<TFile *> SignalFile = {f_signal_50GeVfile};
	  vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_50GeV->Scale((1.0/Signal_Totals[0])*lumi*0.025);
	  histo_signal_50GeV->SetName("ZprimeSignal_mchi50GeV");
	  histo_signal_50GeV->SetLineColor(kSpring-1);
	  histo_signal_50GeV->SetLineWidth(2);
	  
	  histo_Signal.push_back(histo_signal_50GeV);  
	}
       if (mchi == "-1" || mchi == "100GeV")
	{
	  TFile *f_signal_100GeVfile = new TFile("postSignal_mchi100GeV.root");
	  TH1F *histo_signal_100GeV = (TH1F*)f_signal_100GeVfile ->Get(variable);
	  vector<TFile *> SignalFile = {f_signal_100GeVfile};
	  vector<float> Signal_Totals = GetTotal(SignalFile);
	  histo_signal_100GeV->Scale((1.0/Signal_Totals[0])*lumi*0.019);
	  histo_signal_100GeV->SetName("ZprimeSignal_mchi100GeV");
	  histo_signal_100GeV->SetLineColor(kAzure+1);
	  histo_signal_100GeV->SetLineWidth(2);
	  
	  histo_Signal.push_back(histo_signal_100GeV);
	}
    }
  return histo_Signal;
}

void PlotTool::DrawRatio(TH1F* histo_Data, THStack* hs_datamc,TCanvas* c)
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
    //cout<<"stackcontent: "<<stackcontent<<" and data content: "<<datacontent<<endl;
    double ratiocontent=0;
    if(datacontent!=0){
      ratiocontent = ( datacontent) / stackcontent ;}
    double num = 1/datacontent;
    double den = 1/stackcontent;
    double error=0;
    if(datacontent!=0){
      error = ratiocontent*sqrt(pow((dataerror/datacontent),2) + pow((stackerror/stackcontent),2));}
    else {error = 2.07;}
    //cout<<"ratio content: "<<ratiocontent<<" and error: "<<error<<endl;
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
 
  c->Update();
  hs_datamc->GetYaxis()->SetTitle("Events");
  hs_datamc->GetYaxis()->SetTitleOffset(1.5);
}

void PlotTool::DrawAxis(TH1F* histo_Data,THStack* hs_datamc, TCanvas* c,string name)
{
  int nbins = histo_Data->GetNbinsX();
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
}

vector<TH1F*> PlotTool::GetPDF(vector<string> filenames)
{
  TFile* main = TFile::Open(filenames[0].c_str());
  string name = string(filenames[0].c_str());
  name.erase(name.begin()+1,name.begin()+10);
  name.erase(name.end()-5,name.end());
  TH1F* pdfUp = (TH1F*)main->Get((name+"_pdfUp").c_str());
  TH1F* pdfDo = (TH1F*)main->Get((name+"_pdfDown").c_str());

  for (int i = 1; i < filenames.size(); i++)
    {
      TFile* file = TFile::Open(filenames[i].c_str());
      name = string(filenames[i].c_str());
      name.erase(name.begin()+1,name.begin()+10);
      name.erase(name.end()-5,name.end());
      TH1F* UpTemp = (TH1F*)file->Get((name+"_pdfUp").c_str());
      TH1F* DoTemp = (TH1F*)file->Get((name+"_pdfDown").c_str());
      pdfUp->Add(UpTemp);
      pdfDo->Add(DoTemp);
    }

  vector<TH1F*> pdfSys = {pdfUp,pdfDo};
  return pdfSys;
}

void PlotTool::pdf_save(vector< vector<TH1F*> > hs_list)
{
  vector<const char*> mchi = {"1GeV","5GeV","10GeV","20GeV","50GeV","100GeV"};
  vector<const char*> cat = {"_NoCat","_Cat1","_Cat2","_Cat3"};
  for (int k = 0; k < mchi.size(); k++)
    {
      for (int j = 0; j < cat.size(); j++)
	{
	  const char* hs_root = (string("../../Systematics")+string(mchi[k])+string(cat[j])+".root").c_str();
	  TFile* hs_file = TFile::Open(hs_root,"UPDATE");
	  for (int i = 0; i < hs_list.size(); i++)
	    {
	      hs_list[i][0]->SetStats(1);
	      hs_list[i][0]->Write();
	    }
	  for (int i = 0; i < hs_list.size(); i++)
	    {
	      hs_list[i][1]->SetStats(1);
	      hs_list[i][1]->Write();
	    }
	  hs_file->Close();
	}
    }
}
void PlotTool::savepdf()
{
  vector<string> DYFiles = {"histos_MET_postDY100to200.root","histos_MET_postDY1200to2500.root","histos_MET_postDY200to400.root","histos_MET_postDY2500toInf.root","histos_MET_postDY400to600.root","histos_MET_postDY600to800.root","histos_MET_postDY800to1200.root","histos_MET_postDY_MLM_0.root"};
  
  vector<TH1F*> DYpdfSys = GetPDF(DYFiles);

  DYpdfSys[0]->SetName("DYJets_pdfUp");
  DYpdfSys[1]->SetName("DYJets_pdfDown");

  vector<string> GJetFiles = {"histos_MET_postGJets100to200.root","histos_MET_postGJets200to400.root","histos_MET_postGJets400to600.root","histos_MET_postGJets400to600.root","histos_MET_postGJets600toInf.root"};

  vector<TH1F*> GJetpdfSys = GetPDF(GJetFiles);
  
  GJetpdfSys[0]->SetName("GJets_pdfUp");
  GJetpdfSys[1]->SetName("GJets_pdfDown");

  vector<string> QCDFiles = {"histos_MET_postQCD1000to1500_0.root","histos_MET_postQCD100to200_0.root","histos_MET_postQCD1500to2000_0.root","histos_MET_postQCD2000toInf_0.root","histos_MET_postQCD200to300_0.root","histos_MET_postQCD300to500_0.root","histos_MET_postQCD500to700_0.root","histos_MET_postQCD700to1000_0.root"};

  vector<TH1F*> QCDpdfSys = GetPDF(QCDFiles);

  QCDpdfSys[0]->SetName("QCD_pdfUp");
  QCDpdfSys[1]->SetName("QCD_pdfDown");

  vector<string> TTJetFiles = {"histos_MET_postTTJets_MLM.root"};

  vector<TH1F*> TTJetpdfSys = GetPDF(TTJetFiles);

  TTJetpdfSys[0]->SetName("TTJets_pdfUp");
  TTJetpdfSys[1]->SetName("TTJets_pdfDown");
  
  vector<string> WJetFiles = {"histos_MET_postW100to200_0.root","histos_MET_postW1200to2500_0.root","histos_MET_postW200to400_0.root","histos_MET_postW2500toInf_0.root","histos_MET_postW400to600_0.root","histos_MET_postW600to800_0.root","histos_MET_postW800to1200_0.root","histos_MET_postWJets_MLM_0.root"};

  vector<TH1F*> WJetpdfSys = GetPDF(WJetFiles);

  WJetpdfSys[0]->SetName("WJets_pdfUp");
  WJetpdfSys[1]->SetName("WJets_pdfDown");

  vector<string> ZJetFiles = {"histos_MET_postZ100to200_0.root","histos_MET_postZ1200to2500_0.root","histos_MET_postZ200to400_0.root","histos_MET_postZ2500toInf_0.root","histos_MET_postZ400to600_0.root","histos_MET_postZ600to800_0.root","histos_MET_postZ800to1200_0.root"};

  vector<TH1F*> ZJetpdfSys = GetPDF(ZJetFiles);

  ZJetpdfSys[0]->SetName("ZJets_pdfUp");
  ZJetpdfSys[1]->SetName("ZJets_pdfDown");

  vector< vector<TH1F*> > PdfSys = {DYpdfSys,GJetpdfSys,QCDpdfSys,TTJetpdfSys,WJetpdfSys,ZJetpdfSys};

  pdf_save(PdfSys);
}

void PlotTool::saveplot(const char* variable,string name,string varname,string cat,string mchi)
{
  vector<TH1F*> hs_list;
  
  vector<const char*> Data_FileNames = {"postMETdata_final.root"};
  vector<double> Data_Xsec = {-1};
  TH1F* histo_Data = GetHistogram(Data_FileNames,Data_Xsec,variable,"data");
  histo_Data->SetName(((string("data_obs")+varname)).c_str());
  
  cout<<"integral of Data here:"<<histo_Data->Integral()<<endl;
  
  if (varname.find("jes") == string::npos){hs_list.push_back(histo_Data);}

  for(int i = 0; i < MC_FileNames.size(); i++)
    {
      hs_list.push_back(GetHistogram(MC_FileNames[i],MC_Xsec[i],variable,MC_Label[i]));
      hs_list[i]->SetName((MC_Label[i]+varname).c_str());
      cout<<"integral of "<<MC_Label[i]<<" here:"<<hs_list[i]->Integral()<<endl;
    }
  
  vector<TH1F*> histo_Signal = GetSignal(mchi,variable);
  histo_Signal[0]->SetName((string("DM")+varname).c_str());
  if (varname.find("jes") == string::npos){hs_list.push_back(histo_Signal[0]);}

  hs_save(mchi,cat,variable,hs_list);
}

void PlotTool::integral(const char* variable,string name,string mchi,int print)
{
  cout<<name<<endl;
  if (print == 0 || print == -1)
    {
      vector<const char*> Data_FileNames = {regionFile};
      vector<double> Data_Xsec = {-1};

      TH1F* histo_Data = GetHistogram(Data_FileNames,Data_Xsec,variable,"data");

      cout<<"integral of Data here:"<<histo_Data->Integral()<<endl;

      vector<TH1F*> hs_list;
      for(int i = 0; i < MC_FileNames.size(); i++)
	{
	  hs_list.push_back(GetHistogram(MC_FileNames[i],MC_Xsec[i],variable,MC_Label[i]));
	  cout<<"integral of "<<MC_Label[i]<<" here:"<<hs_list[i]->Integral()<<endl;
	}
    }
  if (print == 1 || print == -1)
    {
      vector<TH1F*> histo_Signal = GetSignal(mchi,variable);
      for (int i = 0; i < histo_Signal.size(); i++){cout<<"integral of "<<histo_Signal[i]->GetName()<<" here:"<<histo_Signal[i]->Integral()<<endl;}
    }
}

void PlotTool::plotter(const char * variable,string name,string mchi)
{
  cout << name <<endl;

  TCanvas *c = new TCanvas("c", "canvas",800,800);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  //c->SetLeftMargin(0.15);
  //c->SetLogy();
  //c->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.25,0.99,0.99);
  pad1->Draw(); pad1->cd();
  if (name != "Z Mass (GeV)") pad1->SetLogy();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  //opening the data file and adding "h_dileptonM_8" histogram
  
  vector<const char*> Data_FileNames = {regionFile};
  vector<double> Data_Xsec = {-1};
  TH1F* histo_Data = GetHistogram(Data_FileNames,Data_Xsec,variable,"data");
  histo_Data->SetLineWidth(2);
  histo_Data->SetLineColor(kWhite);
  histo_Data->SetTitle("");
  histo_Data->GetXaxis()->SetTitle("");
  histo_Data->GetXaxis()->SetTickLength(0);
  histo_Data->GetXaxis()->SetLabelOffset(999);
  histo_Data->GetYaxis()->SetTitle("");
  histo_Data->GetYaxis()->SetTickLength(0);
  histo_Data->GetYaxis()->SetLabelOffset(999);

  cout<<"integral of Data here:"<<histo_Data->Integral()<<endl;

  vector<TH1F*> hs_list;
  for(int i = 0; i < MC_FileNames.size(); i++)
    {
      hs_list.push_back(GetHistogram(MC_FileNames[i],MC_Xsec[i],variable,MC_Label[i]));
      hs_list[i]->SetTitle("");
      hs_list[i]->GetXaxis()->SetTitle("");
      hs_list[i]->GetXaxis()->SetTickLength(0);
      hs_list[i]->GetXaxis()->SetLabelOffset(999);
      hs_list[i]->GetYaxis()->SetTitle("");
      hs_list[i]->GetYaxis()->SetTickLength(0);
      hs_list[i]->GetYaxis()->SetLabelOffset(999);
      hs_list[i]->SetFillColor(MC_Color[i]);
      cout<<"integral of "<<MC_Label[i]<<" here:"<<hs_list[i]->Integral()<<endl;
    }

  //Stack histograms using THStack
  THStack* hs_datamc = StackHistogram(hs_list,name);
  hs_datamc->Draw("HIST");

  histo_Data->SetLineColor(kBlack);
  histo_Data->SetMarkerStyle(20);
  histo_Data->SetMarkerSize(0.7);
  histo_Data->Draw("pex0same");
  vector<TH1F*> histo_Signal = GetSignal(mchi,variable);
  if (mchi != "0"){histo_Signal[0]->Draw("HIST same");}
  TLegend *leg = new TLegend(0.62,0.60,0.86,0.887173,"");
  if (mchi != "0"){leg->AddEntry(histo_Signal[0], histo_Signal[0]->GetName());};   
  leg->AddEntry(hs_list[0],"W#rightarrowl#nu","f");
  leg->AddEntry(hs_list[3],"Z#rightarrow ll","F"); 
  leg->AddEntry(hs_list[5],"WW/WZ/ZZ","F");
  leg->AddEntry(hs_list[6], "QCD","F");
  leg->AddEntry(hs_list[4], "Top Quark", "F"); 
  leg->AddEntry(hs_list[2],"#gamma+jets", "F"); 
  leg->AddEntry(hs_list[1],"Z#rightarrow#nu#nu","F");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.025);
  leg->Draw();
  TLatex *texS = new TLatex(0.20,0.837173,"#sqrt{s} = 13 TeV, 1.89 fb^{-1}");
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
 
  DrawRatio(histo_Data,hs_datamc,c);
  
  DrawAxis(histo_Data,hs_datamc,c,name);

  system("echo '${PWD##*/}' > .filelist/dir.txt");
  ifstream dirfile(".filelist/dir.txt");
  string dir;
  if (dirfile.is_open())
    {
      string line;
      while(getline(dirfile,line)) dir=line;
    }
  //c->SaveAs((string(variable)+string(".png")).c_str());
  c->SaveAs((string(variable)+string(".pdf")).c_str());
  c->SaveAs((string(variable)+string(".png")).c_str());
  system((string("mv ")+string(variable)+string(".pdf ")+string("/afs/hep.wisc.edu/home/ekoenig4/public_html/Plots/")+dir+string("Plots_EWK/datamc_")+string(variable)+string(".pdf")).c_str());
  system((string("mv ")+string(variable)+string(".png ")+string("/afs/hep.wisc.edu/home/ekoenig4/public_html/Plots/")+dir+string("Plots_EWK/datamc_")+string(variable)+string(".png")).c_str());
}

void PlotTool::plotTH2F(const char * variable,string name,string mchi)
{
  cout << name <<endl;

  TCanvas *c = new TCanvas("c", "canvas",800,800);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  //c->SetLeftMargin(0.15);
  //c->SetLogy();
  //c->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.25,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  //opening the data file and adding "h_dileptonM_8" histogram

  vector<TH2F*> hs_list;

  TH2F* mainBkg;
  TH2F* allBkg;
  for(int i = 0; i < MC_FileNames.size(); i++)
    {
      hs_list.push_back(GetTH2F(MC_FileNames[i],MC_Xsec[i],variable,MC_Label[i]));
      if (i == 0)
	{
	  mainBkg=hs_list[0];
	  allBkg=hs_list[0];
	}
      if (i == 1)mainBkg->Add(hs_list[i]);
      allBkg->Add(hs_list[i]);
      cout<<"integral of "<<MC_Label[i]<<" here:"<<hs_list[i]->Integral()<<endl;
    }

  mainBkg->GetYaxis()->SetTitle("Uncertainty");
  allBkg->GetYaxis()->SetTitle("Uncertainty");
  
  system("echo '${PWD##*/}' > .filelist/dir.txt");
  ifstream dirfile(".filelist/dir.txt");
  string dir;
  if (dirfile.is_open())
    {
      string line;
      while(getline(dirfile,line)) dir=line;
    }

  mainBkg->Draw();
  //c->SaveAs((string(variable)+string(".png")).c_str());
  c->SaveAs((string(variable)+string("MainBkg.pdf")).c_str());
  c->SaveAs((string(variable)+string("MainBkg.png")).c_str());
  system((string("mv ")+string(variable)+string("MainBkg.pdf ")+string("/afs/hep.wisc.edu/home/ekoenig4/public_html/Plots/")+dir+string("Plots_EWK/datamc_")+string(variable)+string("MainBkg.pdf")).c_str());
  system((string("mv ")+string(variable)+string("MainBkg.png ")+string("/afs/hep.wisc.edu/home/ekoenig4/public_html/Plots/")+dir+string("Plots_EWK/datamc_")+string(variable)+string("MainBkg.png")).c_str());

  allBkg->Draw();
  c->SaveAs((string(variable)+string("AllBkg.pdf")).c_str());
  c->SaveAs((string(variable)+string("AllBkg.png")).c_str());
  system((string("mv ")+string(variable)+string("AllBkg.pdf ")+string("/afs/hep.wisc.edu/home/ekoenig4/public_html/Plots/")+dir+string("Plots_EWK/datamc_")+string(variable)+string("AllBkg.pdf")).c_str());
  system((string("mv ")+string(variable)+string("AllBkg.png ")+string("/afs/hep.wisc.edu/home/ekoenig4/public_html/Plots/")+dir+string("Plots_EWK/datamc_")+string(variable)+string("AllBkg.png")).c_str());    
}
