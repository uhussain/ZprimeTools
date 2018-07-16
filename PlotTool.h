#ifndef PlotTool_h
#define PlotTool_h

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

class PlotTool{
 public:

  double lumi;
  string Region;
  const char* regionFile;

  //List of Signal Files and Xsec
  vector<const char*> Mx10_Mv100 = {"postMx10_Mv100"};vector<const char*> Mx50_Mv100 = {"postMx50_Mv100"};vector<const char*> Mx100_Mv100 = {"postMx100_Mv100"};
  vector<double> Mx10_Mv100_Xsec = {0.006065};        vector<double> Mx50_Mv100_Xsec = {0.002986};        vector<double> Mx100_Mv100_Xsec = {0.001705};
  
  vector<const char*> Mx10_Mv200 = {"postMx10_Mv200"};vector<const char*> Mx50_Mv200 = {"postMx50_Mv200"};vector<const char*> Mx100_Mv200 = {"postMx100_Mv200"};
  vector<double> Mx10_Mv200_Xsec = {0.006522};        vector<double> Mx50_Mv200_Xsec = {0.003184};        vector<double> Mx100_Mv200_Xsec = {0.001786};
  
  vector<const char*> Mx10_Mv500 = {"postMx10_Mv500"};vector<const char*> Mx50_Mv500 = {"postMx50_Mv500"};vector<const char*> Mx100_Mv500 = {"postMx100_Mv500"};
  vector<double> Mx10_Mv500_Xsec = {0.01165};         vector<double> Mx50_Mv500_Xsec = {0.00515};         vector<double> Mx100_Mv500_Xsec = {0.002607};
  
  vector<const char*> Mx10_Mv1000 = {"postMx10_Mv1000"};vector<const char*> Mx50_Mv1000 = {"postMx50_Mv1000"};vector<const char*> Mx100_Mv1000 = {"postMx100_Mv1000"};
  vector<double> Mx10_Mv1000_Xsec = {0.3};              vector<double> Mx50_Mv1000_Xsec = {0.1473};           vector<double> Mx100_Mv1000_Xsec = {0.08242};
  
  vector<const char*> Mx10_Mv1500 = {"postMx10_Mv1500"};vector<const char*> Mx50_Mv1500 = {"postMx50_Mv1500"};vector<const char*> Mx100_Mv1500 = {"postMx100_Mv1500"};
  vector<double> Mx10_Mv1500_Xsec = {0.1559};           vector<double> Mx50_Mv1500_Xsec = {0.08896};           vector<double> Mx100_Mv1500_Xsec = {0.06155};
  
  vector<const char*> Mx10_Mv1800 = {"postMx10_Mv1800"};vector<const char*> Mx50_Mv1800 = {"postMx50_Mv1800"};vector<const char*> Mx100_Mv1800 = {"postMx100_Mv1800"};
  vector<double> Mx10_Mv1800_Xsec = {0.09052};          vector<double> Mx50_Mv1800_Xsec = {0.05555};           vector<double> Mx100_Mv1800_Xsec = {0.04012};
  
  vector<const char*> Mx10_Mv2000 = {"postMx10_Mv2000"};vector<const char*> Mx50_Mv2000 = {"postMx50_Mv2000"};vector<const char*> Mx100_Mv2000 = {"postMx100_Mv2000"};
  vector<double> Mx10_Mv2000_Xsec = {0.06387};          vector<double> Mx50_Mv2000_Xsec = {0.03968};          vector<double> Mx100_Mv2000_Xsec = {0.02939};
  
  vector<const char*> Mx10_Mv2500 = {"postMx10_Mv2500"};vector<const char*> Mx50_Mv2500 = {"postMx50_Mv2500"};vector<const char*> Mx100_Mv2500 = {"postMx100_Mv2500"};
  vector<double> Mx10_Mv2500_Xsec = {0.02553};          vector<double> Mx50_Mv2500_Xsec = {0.01682};          vector<double> Mx100_Mv2500_Xsec = {0.01288};
  
  vector<const char*> Mx10_Mv3500 = {"postMx10_Mv3500"};vector<const char*> Mx50_Mv3500 = {"postMx50_Mv3500"};vector<const char*> Mx100_Mv3500 = {"postMx100_Mv3500"};
  vector<double> Mx10_Mv3500_Xsec = {0.00454};          vector<double> Mx50_Mv3500_Xsec = {0.003093};          vector<double> Mx100_Mv3500_Xsec = {0.00243};

  vector<vector<const char*>> Mx10_Mv = {Mx10_Mv100,Mx10_Mv200,Mx10_Mv1000,Mx10_Mv1500,Mx10_Mv1800,Mx10_Mv2000,Mx10_Mv2500,Mx10_Mv3500};
  vector<vector<double>> Mx10_Mv_Xsec = {Mx10_Mv100_Xsec,Mx10_Mv200_Xsec,Mx10_Mv1000_Xsec,Mx10_Mv1500_Xsec,Mx10_Mv1800_Xsec,Mx10_Mv2000_Xsec,Mx10_Mv2500_Xsec,Mx10_Mv3500_Xsec};

  vector<vector<const char*>> Mx50_Mv = {Mx50_Mv100,Mx50_Mv200,Mx50_Mv1000,Mx50_Mv1500,Mx50_Mv1800,Mx50_Mv2000,Mx50_Mv2500,Mx50_Mv3500};
  vector<vector<double>> Mx50_Mv_Xsec = {Mx50_Mv100_Xsec,Mx50_Mv200_Xsec,Mx50_Mv1000_Xsec,Mx50_Mv1500_Xsec,Mx50_Mv1800_Xsec,Mx50_Mv2000_Xsec,Mx50_Mv2500_Xsec,Mx50_Mv3500_Xsec};

  vector<vector<const char*>> Mx100_Mv = {Mx100_Mv100,Mx100_Mv200,Mx100_Mv1000,Mx100_Mv1500,Mx100_Mv1800,Mx100_Mv2000,Mx100_Mv2500,Mx100_Mv3500};
  vector<vector<double>> Mx100_Mv_Xsec = {Mx100_Mv100_Xsec,Mx100_Mv200_Xsec,Mx100_Mv1000_Xsec,Mx100_Mv1500_Xsec,Mx100_Mv1800_Xsec,Mx100_Mv2000_Xsec,Mx100_Mv2500_Xsec,Mx100_Mv3500_Xsec};

  vector<vector<vector<const char*>>> Mx_Mv = {Mx10_Mv,Mx50_Mv,Mx100_Mv};
  vector<vector<vector<double>>> Mx_Mv_Xsec = {Mx10_Mv_Xsec,Mx50_Mv_Xsec,Mx100_Mv_Xsec};

  vector<string> Mx_Label = {"10","50","100"};
  vector<string> Mv_Label = {"100","200","1000","1500","1800","2000","2500","3500"};
  
  //List of Region Data Files
  vector<const char*> SignalData_FileNames = {"postMETdata_0","postMETdata_1","postMETdata_2"};
  vector<const char*> SingleEleData_FileNames = {"postSingleEle_0","postSingleEle_1","postSingleEle_2","postSingleEle_3","postSingleEle_4","postSingleEle_5","postSingleEle_6","postSingleEle_7","postSingleEle_8","postSingleEle_9","postSingleEle_10","postSingleEle_11","postSingleEle_12","postSingleEle_13","postSingleEle_14","postSingleEle_15","postSingleEle_16","postSingleEle_17","postSingleEle_18","postSingleEle_19"};
  vector<const char*> SingleMuData_FileNames = {"postSingleMu_0","postSingleMu_1","postSingleMu_2","postSingleMu_3","postSingleMu_4","postSingleMu_5","postSingleMu_6","postSingleMu_7","postSingleMu_8","postSingleMu_9","postSingleMu_10","postSingleMu_11","postSingleMu_12","postSingleMu_13","postSingleMu_14","postSingleMu_15","postSingleMu_16","postSingleMu_17","postSingleMu_18","postSingleMu_19"};
  vector<const char*> DoubleEleData_FileNames = {"postDoubleEle_0","postDoubleEle_1","postDoubleEle_2","postDoubleEle_3","postDoubleEle_4","postDoubleEle_5","postDoubleEle_6","postDoubleEle_7","postDoubleEle_8","postDoubleEle_9","postDoubleEle_10","postDoubleEle_11","postDoubleEle_12","postDoubleEle_13","postDoubleEle_14","postDoubleEle_15","postDoubleEle_16","postDoubleEle_17","postDoubleEle_18","postDoubleEle_19"};
  vector<const char*> DoubleMuData_FileNames = {"postDoubleMu_0","postDoubleMu_1","postDoubleMu_2","postDoubleMu_3","postDoubleMu_4","postDoubleMu_5","postDoubleMu_6","postDoubleMu_7","postDoubleMu_8","postDoubleMu_9","postDoubleMu_10","postDoubleMu_11","postDoubleMu_12","postDoubleMu_13","postDoubleMu_14","postDoubleMu_15","postDoubleMu_16","postDoubleMu_17","postDoubleMu_18","postDoubleMu_19"};

  //List of Sample Files and Xsec
  vector<const char *> WJets_FileNames = {"postWJets_MLM_0","postW100to200_0","postW200to400_0","postW400to600_0","postW600to800_0","postW800to1200_0","postW1200to2500_0","postW2500toInf_0"};
  vector<double> WJets_Xsec =            {50690            ,1345             ,359.7            ,48.91            ,12.05            ,5.501             ,1.329              ,0.03216};
  
  vector<const char *> ZJets_FileNames = {"postZ100to200_0","postZ200to400_0","postZ400to600_0","postZ600to800_0","postZ800to1200_0","postZ1200to2500_0","postZ2500toInf_0"};
  vector<double> ZJets_Xsec =            {280.35           ,77.67            ,10.73            ,2.559            ,1.1796            ,0.28833            ,0.006945};
  
  vector<const char *> GJets_FileNames = {"postGJets40to100","postGJets100to200","postGJets200to400","postGJets400to600","postGJets600toInf"};
  vector<double> GJets_Xsec =            {17420             ,5391               ,1168               ,132.5              ,44.05};
  
  vector<const char *> DYJets_FileNames = {"postDY_MLM_0","postDY100to200","postDY200to400","postDY400to600","postDY600to800","postDY800to1200","postDY1200to2500","postDY2500toInf"};
  vector<double> DYJets_Xsec =            {4895          ,148             ,40.94           ,5.497           ,1.354           ,0.6250           ,0.1511            ,0.003647};
  
  vector<const char *> TTJets_FileNames = {"postTTJets_MLM"};
  vector<double> TTJets_Xsec =            {502.2};
  
  vector<const char *> DiBoson_FileNames = {"postWW","postWZ","postZZ"};
  vector<double> DiBoson_Xsec =            {118.7   ,47.2    ,16.6};
  
  vector<const char *> QCD_FileNames = {"postQCD100to200_0","postQCD200to300_0","postQCD300to500_0","postQCD500to700_0","postQCD700to1000_0","postQCD1000to1500_0","postQCD1500to2000_0","postQCD2000toInf_0"};
  vector<double> QCD_Xsec =            {27500000           ,1735000            ,367000             ,29370              ,6524                ,1064                 ,121.5                ,25.42};

  vector<vector<const char*>> MC_FileNames = {WJets_FileNames,ZJets_FileNames,GJets_FileNames,DYJets_FileNames,TTJets_FileNames,DiBoson_FileNames,QCD_FileNames};
  vector<vector<double>> MC_Xsec =           {WJets_Xsec     ,ZJets_Xsec     ,GJets_Xsec     ,DYJets_Xsec     ,TTJets_Xsec     ,DiBoson_Xsec     ,QCD_Xsec};
  vector<Color_t> MC_Color =                 {kRed-10        ,kAzure+10      ,kTeal-9        ,kGray+2         ,kOrange-2       ,kCyan-10         ,kGray};
  vector<string> MC_Label =                  {"WJets"        ,"ZJets"        ,"GJets"        ,"DYJets"        ,"TTJets"        ,"DiBoson"        ,"QCD"};
  
  PlotTool(int argc,const char* argv[]);
  virtual void Options(int argc,const char* argv[]);
  virtual void haddAll(vector<const char*> Data_FileNames);
  virtual const char* haddRegion();
  virtual vector<float> GetTotal(vector<TFile*> Files);
  virtual string SampleName(const char * variable);
  virtual string GetCategory(int hs_num);
  virtual vector<string> GetName(const char * variable);
  virtual int hs_save(string mchi,string cat, const char * variable, vector<TH1F*> histo);
  virtual vector<int> hs_sort(vector<TH1F*> hs_list);
  virtual TH1F* GetHistogram(vector<const char*> Sample_FileNames,vector<double> Sample_Xsec,const char* variable,string SampleName);
  virtual TH2F* GetTH2F(vector<const char*> Sample_FileNames,vector<double> Sample_Xsec,const char* variable,string SampleName);
  virtual THStack* StackHistogram(vector<TH1F*> hs_list,string name);
  virtual vector<TH1F*> GetMx_Mv(string mchi,const char* variable);
  virtual vector<TH1F*> GetSignal(string mchi,const char* variable);
  virtual void DrawRatio(TH1F* histo_Data, THStack* hs_datamc,TCanvas* c);
  virtual void DrawAxis(TH1F* histo_Data,THStack* hs_datamc, TCanvas* c,string name);
  virtual vector<TH1F*> GetPDF(vector<string> filenames);
  virtual void pdf_save(vector< vector<TH1F*> > hs_list);
  virtual void savepdf();
  virtual void saveplot(const char* variable,string name,string varname,string cat,string mchi);
  virtual void integral(const char* variable,string name,string mchi,int print);
  virtual void plotter(const char * variable,string name,string mchi);
  virtual void plotTH2F(const char * variable,string name,string mchi);
  virtual void saveplotOption(int argc,vector<const char*> argv);
  virtual void integralOption(int print);
  virtual void plotterOption(int argc, vector<const char*> argv);
  virtual void plotTH2FOption(int argc,vector<const char*> argv);
  
};

#endif

#ifdef PlotTool_cxx

PlotTool::PlotTool(int argc,const char* argv[])
{
  vector<const char*> preRegionData = {".output/postMETdata_0_1.root",".output/postSingleEle_0_1.root",".output/postSingleMu_0_1.root",".output/postDoubleEle_0_1.root",".output/postDoubleMu_0_1.root"};
  vector<const char*> postRegionData ={"postMETdata_0.root","postSingleEle_0.root","postSingleMu_0.root","postDoubleEle_0.root","postDoubleMu_0.root"}; 
  vector<string> RegionName = {"SignalRegion","SingleEleCR","SingleMuCR","DoubleEleCR","DoubleMuCR"};
  for (int i = 0; i < preRegionData.size(); i++)
    {
      ifstream prefile(preRegionData[i]);
      ifstream postfile(postRegionData[i]);
      if (prefile || postfile)
	{
	  Region=RegionName[i];
	  break;
	}
    }
  if (Region.size() == 0)
    {
      cout<<"Error: No Region Files Detected"<<endl;
    }
}

#endif
