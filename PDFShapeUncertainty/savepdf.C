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
#include "TKey.h"
#include "iostream"
#include "stdlib.h"

std::vector<TH1F*> GetHisto(std::vector<std::string> filenames)
{
  TFile* main = TFile::Open(filenames[0].c_str());
  std::string name = std::string(filenames[0].c_str());
  name.erase(name.begin()+1,name.begin()+10);
  name.erase(name.end()-5,name.end());
  TH1F* pdfUp = (TH1F*)main->Get((name+"_pdfUp").c_str());
  TH1F* pdfDo = (TH1F*)main->Get((name+"_pdfDown").c_str());

  for (int i = 1; i < filenames.size(); i++)
    {
      TFile* file = TFile::Open(filenames[i].c_str());
      name = std::string(filenames[i].c_str());
      name.erase(name.begin()+1,name.begin()+10);
      name.erase(name.end()-5,name.end());
      TH1F* UpTemp = (TH1F*)file->Get((name+"_pdfUp").c_str());
      TH1F* DoTemp = (TH1F*)file->Get((name+"_pdfDown").c_str());
      pdfUp->Add(UpTemp);
      pdfDo->Add(DoTemp);
    }

  std::vector<TH1F*> pdfSys = {pdfUp,pdfDo};
  return pdfSys;
}

void SaveHisto(std::vector< std::vector<TH1F*> > hs_list)
{
  const char* hs_root = "../../LimitSyst/Systematics_etaWidth.root";
  const char* region = "SignalRegion";
  std::ifstream file(hs_root);
  if (file){file.close();}
  else
    {
      TFile* hs_file = TFile::Open(hs_root,"RECREATE");
      hs_file->Close();
    }
  TFile* hs_file = TFile::Open(hs_root,"UPDATE");
  /*if (!(hs_file->GetDirectory(region)))
    {
      hs_file->mkdir(region);
    }
    hs_file->cd(region);*/
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

void savepdf()
{
  std::vector<std::string> DYFiles = {"histos_MET_postDY100to200.root","histos_MET_postDY1200to2500.root","histos_MET_postDY200to400.root","histos_MET_postDY2500toInf.root","histos_MET_postDY400to600.root","histos_MET_postDY600to800.root","histos_MET_postDY800to1200.root","histos_MET_postDY_MLM_0.root"};
  
  std::vector<TH1F*> DYpdfSys = GetHisto(DYFiles);

  DYpdfSys[0]->SetName("DYJets_PDFUp");
  DYpdfSys[1]->SetName("DYJets_PDFDown");

  std::vector<std::string> GJetFiles = {"histos_MET_postGJets100to200.root","histos_MET_postGJets200to400.root","histos_MET_postGJets400to600.root","histos_MET_postGJets400to600.root","histos_MET_postGJets600toInf.root"};

  std::vector<TH1F*> GJetpdfSys = GetHisto(GJetFiles);
  
  GJetpdfSys[0]->SetName("GJets_PDFUp");
  GJetpdfSys[1]->SetName("GJets_PDFDown");

  std::vector<std::string> QCDFiles = {"histos_MET_postQCD1000to1500_0.root","histos_MET_postQCD100to200_0.root","histos_MET_postQCD1500to2000_0.root","histos_MET_postQCD2000toInf_0.root","histos_MET_postQCD200to300_0.root","histos_MET_postQCD300to500_0.root","histos_MET_postQCD500to700_0.root","histos_MET_postQCD700to1000_0.root"};

  std::vector<TH1F*> QCDpdfSys = GetHisto(QCDFiles);

  QCDpdfSys[0]->SetName("QCD_PDFUp");
  QCDpdfSys[1]->SetName("QCD_PDFDown");

  std::vector<std::string> TTJetFiles = {"histos_MET_postTTJets_MLM.root"};

  std::vector<TH1F*> TTJetpdfSys = GetHisto(TTJetFiles);

  TTJetpdfSys[0]->SetName("TTJets_PDFUp");
  TTJetpdfSys[1]->SetName("TTJets_PDFDown");
  
  std::vector<std::string> WJetFiles = {"histos_MET_postW100to200_0.root","histos_MET_postW1200to2500_0.root","histos_MET_postW200to400_0.root","histos_MET_postW2500toInf_0.root","histos_MET_postW400to600_0.root","histos_MET_postW600to800_0.root","histos_MET_postW800to1200_0.root","histos_MET_postWJets_MLM_0.root"};

  std::vector<TH1F*> WJetpdfSys = GetHisto(WJetFiles);

  WJetpdfSys[0]->SetName("WJets_PDFUp");
  WJetpdfSys[1]->SetName("WJets_PDFDown");

  std::vector<std::string> ZJetFiles = {"histos_MET_postZ100to200_0.root","histos_MET_postZ1200to2500_0.root","histos_MET_postZ200to400_0.root","histos_MET_postZ2500toInf_0.root","histos_MET_postZ400to600_0.root","histos_MET_postZ600to800_0.root","histos_MET_postZ800to1200_0.root"};

  std::vector<TH1F*> ZJetpdfSys = GetHisto(ZJetFiles);

  ZJetpdfSys[0]->SetName("ZJets_PDFUp");
  ZJetpdfSys[1]->SetName("ZJets_PDFDown");

  std::vector< std::vector<TH1F*> > PdfSys = {DYpdfSys,GJetpdfSys,QCDpdfSys,TTJetpdfSys,WJetpdfSys,ZJetpdfSys};

  SaveHisto(PdfSys);
}
