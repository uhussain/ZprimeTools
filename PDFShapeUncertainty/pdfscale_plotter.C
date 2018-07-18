#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include <bitset>
#include <climits>
#include <string>
#include <cstring>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
using namespace std;

// xaxis_variable must be one of "Pt", "MET", or "Mt"
// int_lumi measured in pb-1
void plot(string xaxis_variable, string filename, Float_t int_lumi, Float_t nevents, Float_t xsec)
{
  // Set things according to the variable specified in xaxis_variable
  TString xaxis_title = TString("TITLE NOT SET");
  int nBins_temp = -1;
  string histname_input = "HISTNAME_INPUT NOT SET";
  if(xaxis_variable == "j1Pt"){
    xaxis_title = TString("p_{T} of Leading Jet (GeV)");
    nBins_temp = 48;
    histname_input = "j1pT";
  }
  else if(xaxis_variable == "MET"){
    //string filename_prefix = filename.substr(0,5);
    //if(filename_prefix == "postZ" || filename_prefix == "postW"){
      xaxis_title = TString("E_{T}^{miss} (GeV)");
      nBins_temp = 48;
      histname_input = "pfMET";
    //}
    //else{
    //  std::cout<<"ERROR: Invalid prefix for filename "<<filename<<std::endl;
    //}
  }
  else if(xaxis_variable == "Mt"){
    string filename_prefix = filename.substr(0,4);
    if(filename_prefix == "ZnnG" || filename_prefix == "GJet"){
      xaxis_title = TString("Photon-MET #it{M}_{T} [GeV]");
      nBins_temp = 9;
      histname_input = "Mt";
    }
    else if(filename_prefix == "WenG" || filename_prefix == "WmnG" || filename_prefix == "ZeeG" || filename_prefix == "ZmmG"){
      xaxis_title = TString("Photon-MET #it{M}_{T} [GeV]");
      nBins_temp = 9;
      histname_input = "h_phoRecoilMt";
    }
  }
  else{
    std::cout<<"ERROR: Invalid value for xaxis_variable"<<std::endl;
  }
  const int nBins = nBins_temp;
  
  // Other useful constants
  TString yaxis_title = TString("Events / bin");
  //Float_t scale_factor = 0.984*1.002;
  Float_t scale_factor = 1.0;
  const int nPDFreplicas = 100;
  const int nPDFStart = 16; //Thats where my PDFScale histos start : hist_16 -> hist_116

  // The names of the final output histograms as they should appear in the data card
  TString outputName_default_unscaled = TString("h_"+filename+"_default_unscaled");
  TString outputName_facUp_unscaled = TString("h_"+filename+"_facUp_unscaled");
  TString outputName_facDown_unscaled = TString("h_"+filename+"_facDown_unscaled");
  TString outputName_renUp_unscaled = TString("h_"+filename+"_renUp_unscaled");
  TString outputName_renDown_unscaled = TString("h_"+filename+"_renDown_unscaled");
  TString outputName_bothUp_unscaled = TString("h_"+filename+"_bothUp_unscaled");
  TString outputName_bothDown_unscaled = TString("h_"+filename+"_bothDown_unscaled");
  TString outputName_pdfUp_unscaled = TString("h_"+filename+"_pdfUp_unscaled");
  TString outputName_pdfDown_unscaled = TString("h_"+filename+"_pdfDown_unscaled");
  TString outputName_uncorrected_unscaled = TString("h_"+filename+"_uncorrected_unscaled");
  TString outputName_onlyEWK_unscaled = TString("h_"+filename+"_onlyEWK_unscaled");
  TString outputName_onlyNNLO_unscaled = TString("h_"+filename+"_onlyNNLO_unscaled");
  
  TString outputName_default = TString("h_"+filename+"_default");
  TString outputName_facUp = TString("h_"+filename+"_facUp");
  TString outputName_facDown = TString("h_"+filename+"_facDown");
  TString outputName_renUp = TString("h_"+filename+"_renUp");
  TString outputName_renDown = TString("h_"+filename+"_renDown");
  TString outputName_bothUp = TString("h_"+filename+"_bothUp");
  TString outputName_bothDown = TString("h_"+filename+"_bothDown");
  TString outputName_pdfUp = TString("h_"+filename+"_pdfUp");
  TString outputName_pdfDown = TString("h_"+filename+"_pdfDown");
  TString outputName_uncorrected = TString("h_"+filename+"_uncorrected");
  TString outputName_onlyEWK = TString("h_"+filename+"_onlyEWK");
  TString outputName_onlyNNLO = TString("h_"+filename+"_onlyNNLO");

  // The names of the histograms produced by the postAnalyzers, to be read by this script
  TString inputName_default = TString(histname_input+"_8");
  TString inputName_facUp = TString(histname_input+"_118");
  TString inputName_facDown = TString(histname_input+"_119");
  TString inputName_renUp = TString(histname_input+"_120");
  TString inputName_renDown = TString(histname_input+"_123");
  TString inputName_bothUp = TString(histname_input+"_121");
  TString inputName_bothDown = TString(histname_input+"_125");
  TString inputName_uncorrected = TString(histname_input+"_111");
  TString inputName_onlyEWK = TString(histname_input+"_112");
  TString inputName_onlyNNLO = TString(histname_input+"_113");
  std::vector<TString> histnames_pdf;
  histnames_pdf.clear();
  //my pdf histos are from hist_16 to hist_116
  //histnames_pdf starts from 0 though :)
  for(int i = nPDFStart; i < nPDFStart+nPDFreplicas+1; i++){
    char histindex_char[100];
    sprintf(histindex_char, "_%d", i);
    std::string histindex(histindex_char);
    TString histname = TString(histname_input+histindex);
    histnames_pdf.push_back(histname);
  }

  TFile* inputfile = new TFile(TString(filename+".root"));

  // Default selection
  TH1F *histo_default_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_default))->Clone(outputName_default_unscaled);
  histo_default_unscaled->SetBinContent(nBins, histo_default_unscaled->GetBinContent(nBins)+histo_default_unscaled->GetBinContent(nBins+1));
  histo_default_unscaled->ClearUnderflowAndOverflow();
  TH1F *histo_default = (TH1F*)histo_default_unscaled->Clone(outputName_default);
  histo_default->SetStats(0);
  histo_default->Scale(xsec*int_lumi*scale_factor/nevents);
  Float_t int_histo_default = histo_default->Integral();
  cout<<filename<<": "<<int_histo_default<<endl;
  histo_default->SetTitle("");
  histo_default->GetXaxis()->SetTitle(xaxis_title);
  histo_default->GetYaxis()->SetTitle(yaxis_title);

  // Apply only one, or none, of the NLO EWK + NNLO QCD corrections
  TH1F *histo_uncorrected_unscaled;
  TH1F *histo_uncorrected;
  TH1F *histo_onlyEWK_unscaled;
  TH1F *histo_onlyEWK;
  TH1F *histo_onlyNNLO_unscaled;
  TH1F *histo_onlyNNLO;
  if(filename == "ZnnG_pdfscale_ZNuNuGJets" || filename == "ZnnG_pdfscale_ZNuNuGJets_ext" || filename == "ZnnG_pdfscale_WGJets" || filename == "ZnnG_pdfscale_WGJets_ext" || filename == "ZnnG_pdfscale_WGJets_mitExt" || filename == "WenG_pdfscale_WGJets" || filename == "WenG_pdfscale_WGJets_ext" || filename == "WenG_pdfscale_WGJets_mitExt" || filename == "WmnG_pdfscale_WGJets" || filename == "WmnG_pdfscale_WGJets_ext" || filename == "WmnG_pdfscale_WGJets_mitExt"){
    histo_uncorrected_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_uncorrected))->Clone(outputName_uncorrected_unscaled);
    histo_uncorrected_unscaled->SetBinContent(nBins, histo_uncorrected_unscaled->GetBinContent(nBins)+histo_uncorrected_unscaled->GetBinContent(nBins+1));
    histo_uncorrected_unscaled->ClearUnderflowAndOverflow();
    histo_uncorrected = (TH1F*)histo_uncorrected_unscaled->Clone(outputName_uncorrected);
    histo_uncorrected->SetStats(0);
    histo_uncorrected->Scale(xsec*int_lumi*scale_factor/nevents);
    histo_uncorrected->SetTitle("");
    histo_uncorrected->GetXaxis()->SetTitle(xaxis_title);
    histo_uncorrected->GetYaxis()->SetTitle(yaxis_title);
    histo_onlyEWK_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_onlyEWK))->Clone(outputName_onlyEWK_unscaled);
    histo_onlyEWK_unscaled->SetBinContent(nBins, histo_onlyEWK_unscaled->GetBinContent(nBins)+histo_onlyEWK_unscaled->GetBinContent(nBins+1));
    histo_onlyEWK_unscaled->ClearUnderflowAndOverflow();
    histo_onlyEWK = (TH1F*)histo_onlyEWK_unscaled->Clone(outputName_onlyEWK);
    histo_onlyEWK->SetStats(0);
    histo_onlyEWK->Scale(xsec*int_lumi*scale_factor/nevents);
    histo_onlyEWK->SetTitle("");
    histo_onlyEWK->GetXaxis()->SetTitle(xaxis_title);
    histo_onlyEWK->GetYaxis()->SetTitle(yaxis_title);
    histo_onlyNNLO_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_onlyNNLO))->Clone(outputName_onlyNNLO_unscaled);
    histo_onlyNNLO_unscaled->SetBinContent(nBins, histo_onlyNNLO_unscaled->GetBinContent(nBins)+histo_onlyNNLO_unscaled->GetBinContent(nBins+1));
    histo_onlyNNLO_unscaled->ClearUnderflowAndOverflow();
    histo_onlyNNLO = (TH1F*)histo_onlyNNLO_unscaled->Clone(outputName_onlyNNLO);
    histo_onlyNNLO->SetStats(0);
    histo_onlyNNLO->Scale(xsec*int_lumi*scale_factor/nevents);
    histo_onlyNNLO->SetTitle("");
    histo_onlyNNLO->GetXaxis()->SetTitle(xaxis_title);
    histo_onlyNNLO->GetYaxis()->SetTitle(yaxis_title);
  }

  // Scale variations
  TH1F* histo_facUp_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_facUp))->Clone(outputName_facUp_unscaled);
  histo_facUp_unscaled->SetBinContent(nBins, histo_facUp_unscaled->GetBinContent(nBins)+histo_facUp_unscaled->GetBinContent(nBins+1));
  histo_facUp_unscaled->ClearUnderflowAndOverflow();
  TH1F* histo_facUp = (TH1F*)histo_facUp_unscaled->Clone(outputName_facUp);
  histo_facUp->SetStats(0);
  histo_facUp->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_facUp->SetTitle("");
  histo_facUp->GetXaxis()->SetTitle(xaxis_title);
  histo_facUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_facDown_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_facDown))->Clone(outputName_facDown_unscaled);
  histo_facDown_unscaled->SetBinContent(nBins, histo_facDown_unscaled->GetBinContent(nBins)+histo_facDown_unscaled->GetBinContent(nBins+1));
  histo_facDown_unscaled->ClearUnderflowAndOverflow();
  TH1F* histo_facDown = (TH1F*)histo_facDown_unscaled->Clone(outputName_facDown);
  histo_facDown->SetStats(0);
  histo_facDown->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_facDown->SetTitle("");
  histo_facDown->GetXaxis()->SetTitle(xaxis_title);
  histo_facDown->GetYaxis()->SetTitle(yaxis_title);
  
  TH1F* histo_renUp_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_renUp))->Clone(outputName_renUp_unscaled);
  histo_renUp_unscaled->SetBinContent(nBins, histo_renUp_unscaled->GetBinContent(nBins)+histo_renUp_unscaled->GetBinContent(nBins+1));
  histo_renUp_unscaled->ClearUnderflowAndOverflow();
  TH1F* histo_renUp = (TH1F*)histo_renUp_unscaled->Clone(outputName_renUp);
  histo_renUp->SetStats(0);
  histo_renUp->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_renUp->SetTitle("");
  histo_renUp->GetXaxis()->SetTitle(xaxis_title);
  histo_renUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_renDown_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_renDown))->Clone(outputName_renDown_unscaled);
  histo_renDown_unscaled->SetBinContent(nBins, histo_renDown_unscaled->GetBinContent(nBins)+histo_renDown_unscaled->GetBinContent(nBins+1));
  histo_renDown_unscaled->ClearUnderflowAndOverflow();
  TH1F* histo_renDown = (TH1F*)histo_renDown_unscaled->Clone(outputName_renDown);
  histo_renDown->SetStats(0);
  histo_renDown->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_renDown->SetTitle("");
  histo_renDown->GetXaxis()->SetTitle(xaxis_title);
  histo_renDown->GetYaxis()->SetTitle(yaxis_title);

  TH1F* histo_bothUp_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_bothUp))->Clone(outputName_bothUp_unscaled);
  histo_bothUp_unscaled->SetBinContent(nBins, histo_bothUp_unscaled->GetBinContent(nBins)+histo_bothUp_unscaled->GetBinContent(nBins+1));
  histo_bothUp_unscaled->ClearUnderflowAndOverflow();
  TH1F* histo_bothUp = (TH1F*)histo_bothUp_unscaled->Clone(outputName_bothUp);
  histo_bothUp->SetStats(0);
  histo_bothUp->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_bothUp->SetTitle("");
  histo_bothUp->GetXaxis()->SetTitle(xaxis_title);
  histo_bothUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_bothDown_unscaled = (TH1F*)((TH1F*)inputfile->Get(inputName_bothDown))->Clone(outputName_bothDown_unscaled);
  histo_bothDown_unscaled->SetBinContent(nBins, histo_bothDown_unscaled->GetBinContent(nBins)+histo_bothDown_unscaled->GetBinContent(nBins+1));
  histo_bothDown_unscaled->ClearUnderflowAndOverflow();
  TH1F* histo_bothDown = (TH1F*)histo_bothDown_unscaled->Clone(outputName_bothDown);
  histo_bothDown->SetStats(0);
  histo_bothDown->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_bothDown->SetTitle("");
  histo_bothDown->GetXaxis()->SetTitle(xaxis_title);
  histo_bothDown->GetYaxis()->SetTitle(yaxis_title);
  
  // PDF variations
  std::vector<std::vector<Float_t>> pdf_sums;
  pdf_sums.clear();
  for(int j = 1; j <= nBins; j++){ // j = bin number
    std::vector<Float_t> pdf_sum_bin;
    pdf_sum_bin.clear();
    pdf_sums.push_back(pdf_sum_bin);
  }
  //histnames_pdf is a vector defined in L111 containing the names of the hists from the postAnalyzer
  //100 hists in it
  for(int i = 0; i <= nPDFreplicas; i++){ // i = pdf replica number
    TH1F* hist_pdf = (TH1F*)inputfile->Get(histnames_pdf[i]);
    for(int j = 1; j <= nBins; j++){ // j = bin number
      pdf_sums[j-1].push_back(hist_pdf->GetBinContent(j));
    }
  }
  
  std::vector<double> sum_of_sums_eachBin;
  sum_of_sums_eachBin.clear();
  for(int i = 1; i <= nBins; i++){
    sum_of_sums_eachBin.push_back(0.0);
  }
  for(int i = 1; i <= nPDFreplicas; i++){ // i = pdf replica number
    for(int j = 1; j <= nBins; j++){ // j = bin number
      sum_of_sums_eachBin[j-1] += pdf_sums[j-1][i];
    }
  }
  
  std::vector<double> mean_Npassing_eachBin;
  mean_Npassing_eachBin.clear();
  for(int j = 1; j <= nBins; j++){ // j = bin number
    mean_Npassing_eachBin.push_back(sum_of_sums_eachBin[j-1]/nPDFreplicas);
  }
  
  std::vector<double> sum_of_squared_residuals_eachBin;
  sum_of_squared_residuals_eachBin.clear();
  for(int j = 1; j <= nBins; j++){
    sum_of_squared_residuals_eachBin.push_back(0.0);
  }
  for(int i = 1; i <= nPDFreplicas; i++){ // i = pdf replica number
    for(int j = 1; j <= nBins; j++){ // j = bin number
      sum_of_squared_residuals_eachBin[j-1] += pow((pdf_sums[j-1][i] - mean_Npassing_eachBin[j-1]), 2.0);
    }
  }
  
  std::vector<Float_t> rms_error_eachBin;
  rms_error_eachBin.clear();
  for(int j = 1; j <= nBins; j++){ // j = bin number
    rms_error_eachBin.push_back(sqrt(sum_of_squared_residuals_eachBin[j-1]/(nPDFreplicas-1))); // RMS error in event weight sum for this bin
  }
  
  TH1F* histo_pdfUp_unscaled = (TH1F*)((TH1F*)inputfile->Get(histnames_pdf[0]))->Clone(outputName_pdfUp_unscaled);
  histo_pdfUp_unscaled->SetBinContent(nBins, histo_pdfUp_unscaled->GetBinContent(nBins)+histo_pdfUp_unscaled->GetBinContent(nBins+1));
  histo_pdfUp_unscaled->ClearUnderflowAndOverflow();
  for(int j = 1; j <= nBins; j++){
    histo_pdfUp_unscaled->SetBinContent(j,histo_pdfUp_unscaled->GetBinContent(j) + rms_error_eachBin[j-1]);
  }
  TH1F* histo_pdfUp = (TH1F*)histo_pdfUp_unscaled->Clone(outputName_pdfUp);
  histo_pdfUp->SetStats(0);
  histo_pdfUp->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_pdfUp->SetTitle("");
  histo_pdfUp->GetXaxis()->SetTitle(xaxis_title);
  histo_pdfUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_pdfDown_unscaled = (TH1F*)((TH1F*)inputfile->Get(histnames_pdf[0]))->Clone(outputName_pdfDown_unscaled);
  histo_pdfDown_unscaled->SetBinContent(nBins, histo_pdfDown_unscaled->GetBinContent(nBins)+histo_pdfDown_unscaled->GetBinContent(nBins+1));
  histo_pdfDown_unscaled->ClearUnderflowAndOverflow();
  for(int j = 1; j <= nBins; j++){
    histo_pdfDown_unscaled->SetBinContent(j,TMath::Max(histo_pdfDown_unscaled->GetBinContent(j) - rms_error_eachBin[j-1], 0.0));
  }
  TH1F* histo_pdfDown = (TH1F*)histo_pdfDown_unscaled->Clone(outputName_pdfDown);
  histo_pdfDown->SetStats(0);
  histo_pdfDown->Scale(xsec*int_lumi*scale_factor/nevents);
  histo_pdfDown->SetTitle("");
  histo_pdfDown->GetXaxis()->SetTitle(xaxis_title);
  histo_pdfDown->GetYaxis()->SetTitle(yaxis_title);
  
  Float_t sum_scale_difference = 0.0;
  Float_t sum_pdf_difference = 0.0;
  for(int j = 1; j <= nBins; j++){ // j = bin number
    Float_t int_bin_nominal = histo_default->GetBinContent(j);
    
    Float_t int_bin_facUp = histo_facUp->GetBinContent(j);
    Float_t int_bin_facDown = histo_facDown->GetBinContent(j);
    Float_t int_bin_renUp = histo_renUp->GetBinContent(j);
    Float_t int_bin_renDown = histo_renDown->GetBinContent(j);
    Float_t int_bin_bothUp = histo_bothUp->GetBinContent(j);
    Float_t int_bin_bothDown = histo_bothDown->GetBinContent(j);
    Float_t int_bin_variations[] = {int_bin_facUp, int_bin_facDown, int_bin_renUp, int_bin_renDown, int_bin_bothUp, int_bin_bothDown};
    Float_t max_int_bin = *max_element(int_bin_variations, int_bin_variations+6);
    Float_t min_int_bin = *min_element(int_bin_variations, int_bin_variations+6);
    sum_scale_difference += (fabs(max_int_bin - int_bin_nominal) + fabs(int_bin_nominal - min_int_bin))/2.0;
    
    Float_t int_bin_pdfUp = histo_pdfUp->GetBinContent(j);
    Float_t int_bin_pdfDown = histo_pdfDown->GetBinContent(j);
    Float_t avg_pdf_difference = (fabs(int_bin_pdfUp - int_bin_nominal) + fabs(int_bin_pdfDown - int_bin_nominal))/2.0;
    // cout<<"Fractional PDF uncertainty bin "<<j<<": "<<(avg_pdf_difference/int_bin_nominal)<<endl;
    sum_pdf_difference += avg_pdf_difference;
  }
  cout<<"Scale uncertainty: "<<sum_scale_difference<<" net, "<<(sum_scale_difference/int_histo_default)<<" fractional"<<endl;
  cout<<"Overall PDF uncertainty: "<<sum_pdf_difference<<" net, "<<(sum_pdf_difference/int_histo_default)<<" fractional"<<endl;

  // Set zero bins to a small value to avoid potential fitting issues
  for(int i = 1; i <= nBins; i++){
    if(histo_default->GetBinContent(i) == 0.0)
      histo_default->SetBinContent(i,1e-6);
    if(histo_facUp->GetBinContent(i) == 0.0)
      histo_facUp->SetBinContent(i,0.97e-6);
    if(histo_facDown->GetBinContent(i) == 0.0)
      histo_facDown->SetBinContent(i,1.03e-6);
    if(histo_renUp->GetBinContent(i) == 0.0)
      histo_renUp->SetBinContent(i,0.97e-6);
    if(histo_renDown->GetBinContent(i) == 0.0)
      histo_renDown->SetBinContent(i,1.03e-6);
    if(histo_pdfUp->GetBinContent(i) == 0.0)
      histo_pdfUp->SetBinContent(i,1.03e-6);
    if(histo_pdfDown->GetBinContent(i) == 0.0)
      histo_pdfDown->SetBinContent(i,0.97e-6);
  }
  
  // Write output histograms
  TFile *outputFile = new TFile(TString("histos_"+xaxis_variable+"_"+filename+".root"),"RECREATE");
  outputFile->cd();
  histo_default_unscaled->Write();
  histo_facUp_unscaled->Write();
  histo_facDown_unscaled->Write();
  histo_renUp_unscaled->Write();
  histo_renDown_unscaled->Write();
  histo_bothUp_unscaled->Write();
  histo_bothDown_unscaled->Write();
  histo_pdfUp_unscaled->Write();
  histo_pdfDown_unscaled->Write();
  histo_default->Write();
  histo_facUp->Write();
  histo_facDown->Write();
  histo_renUp->Write();
  histo_renDown->Write();
  histo_bothUp->Write();
  histo_bothDown->Write();
  histo_pdfUp->Write();
  histo_pdfDown->Write();
  if(filename == "ZnnG_pdfscale_ZNuNuGJets" || filename == "ZnnG_pdfscale_ZNuNuGJets_ext" || filename == "ZnnG_pdfscale_WGJets" || filename == "ZnnG_pdfscale_WGJets_ext" || filename == "ZnnG_pdfscale_WGJets_mitExt" || filename == "WenG_pdfscale_WGJets" || filename == "WenG_pdfscale_WGJets_ext" || filename == "WenG_pdfscale_WGJets_mitExt" || filename == "WmnG_pdfscale_WGJets" || filename == "WmnG_pdfscale_WGJets_ext" || filename == "WmnG_pdfscale_WGJets_mitExt"){
    histo_uncorrected_unscaled->Write();
    histo_uncorrected->Write();
    histo_onlyEWK_unscaled->Write();
    histo_onlyEWK->Write();
    histo_onlyNNLO_unscaled->Write();
    histo_onlyNNLO->Write();
  }
  
  cout<<"File "<<("histos_"+xaxis_variable+"_"+filename+".root")<<" created"<<endl;

  inputfile->Close();
  outputFile->Close();
}
//GetTotal number of Events
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
void pdfscale_plotter()
{
  std::vector<string> filenames;
  filenames.clear();
  std::vector<Float_t> int_lumis;
  int_lumis.clear();
  std::vector<Float_t> nevents;
  nevents.clear();
  std::vector<Float_t> xsecs;
  xsecs.clear();

  //WJets Samples
  std::vector<string> WJets_FileNames = {"postWJets_MLM_0","postW100to200_0","postW200to400_0","postW400to600_0","postW600to800_0","postW800to1200_0","postW1200to2500_0","postW2500toInf_0"};
  std::vector<TFile*> WJets_Files;
  std::vector<Float_t> WJets_xsecs={50690,1345,359.7,48.91,12.05,5.501,1.329,0.03216};
  for (int i = 0; i < WJets_FileNames.size(); i++){WJets_Files.push_back(new TFile(TString(WJets_FileNames[i]+".root")));}			    
  std::vector<float> WJets_Total = GetTotal(WJets_Files);

  //Loop over the WJets Files
  for(int i=0; i < WJets_FileNames.size();i++){
    filenames.push_back(WJets_FileNames[i]);
    int_lumis.push_back(35900);
    nevents.push_back(WJets_Total[i]);
    xsecs.push_back(WJets_xsecs[i]);
  }
  
  //ZvvJets Samples
  std::vector<string> Zvv_FileNames = {"postZ100to200_0","postZ200to400_0","postZ400to600_0","postZ600to800_0","postZ800to1200_0","postZ1200to2500_0","postZ2500toInf_0"};
  std::vector<TFile *> Zvv_Files;
  std::vector<Float_t> Zvv_xsecs={280.35,77.67,10.73,2.559,1.1796,0.28833,0.006945};
  for (int i = 0; i < Zvv_FileNames.size(); i++) {Zvv_Files.push_back(new TFile(TString(Zvv_FileNames[i]+".root")));}
  std::vector<float> Zvv_Total = GetTotal(Zvv_Files);

  //Loop over the Zvv Files
  for(int i=0; i < Zvv_FileNames.size();i++){
    filenames.push_back(Zvv_FileNames[i]);
    int_lumis.push_back(35900);
    nevents.push_back(Zvv_Total[i]);
    xsecs.push_back(Zvv_xsecs[i]);
  }

  //DiBoson Samples

//  std::vector<string> DiBoson_FileNames = {"postWW","postWZ","postZZ"};
//  std::vector<TFile *> DiBoson_Files;
//  std::vector<Float_t> DiBoson_xsecs={118.7,47.2,16.6};
//  for (int i = 0; i < DiBoson_FileNames.size(); i++) {DiBoson_Files.push_back(new TFile(TString(DiBoson_FileNames[i]+".root")));}
//  std::vector<float> DiBoson_Total = GetTotal(DiBoson_Files);
//
//  //Loop over the DiBoson Files
//  for(int i=0; i < DiBoson_FileNames.size();i++){
//    std::cout<<"DiBosonSample: "<<DiBoson_FileNames[i]<<std::endl;
//    filenames.push_back(DiBoson_FileNames[i]);
//    int_lumis.push_back(35900);
//    nevents.push_back(DiBoson_Total[i]);
//    xsecs.push_back(DiBoson_xsecs[i]);
//  }

  //Top Quark Samples

  std::vector<string> TopQuark_FileNames = {"postTTJets_MLM"};
  std::vector<TFile *> TopQuark_Files;
  std::vector<Float_t> TopQuark_xsecs={502.2};
  for (int i = 0; i < TopQuark_FileNames.size(); i++) {TopQuark_Files.push_back(new TFile(TString(TopQuark_FileNames[i]+".root")));}
  std::vector<float> TopQuark_Total = GetTotal(TopQuark_Files);

  //Loop over the TopQuark Files
  for(int i=0; i < TopQuark_FileNames.size();i++){
    std::cout<<"TopQuarkSample: "<<TopQuark_FileNames[i]<<std::endl;
    filenames.push_back(TopQuark_FileNames[i]);
    int_lumis.push_back(35900);
    nevents.push_back(TopQuark_Total[i]);
    xsecs.push_back(TopQuark_xsecs[i]);
  }

  //G+Jets Samples

  std::vector<string> GJets_FileNames = {"postGJets40to100","postGJets100to200","postGJets200to400","postGJets400to600","postGJets600toInf"};
  std::vector<TFile *> GJets_Files;
  std::vector<Float_t> GJets_xsecs={17420,5391,1168,132.5,44.05};
  for (int i = 0; i < GJets_FileNames.size(); i++) {GJets_Files.push_back(new TFile(TString(GJets_FileNames[i]+".root")));}
  std::vector<float> GJets_Total = GetTotal(GJets_Files);

  //Loop over the GJets Files
  for(int i=0; i < GJets_FileNames.size();i++){
    std::cout<<"GJetsSample: "<<GJets_FileNames[i]<<std::endl;
    filenames.push_back(GJets_FileNames[i]);
    int_lumis.push_back(35900);
    nevents.push_back(GJets_Total[i]);
    xsecs.push_back(GJets_xsecs[i]);
  }

  //DY Samples

  std::vector<string> DYJets_FileNames = {"postDY_MLM_0","postDY100to200","postDY200to400","postDY400to600","postDY600to800","postDY800to1200","postDY1200to2500","postDY2500toInf"};
  std::vector<TFile *> DYJets_Files;
  std::vector<Float_t> DYJets_xsecs={4895,148,40.94,5.497,1.354,0.6250,0.1511,0.003647};
  for (int i = 0; i < DYJets_FileNames.size(); i++) {DYJets_Files.push_back(new TFile(TString(DYJets_FileNames[i]+".root")));}
  std::vector<float> DYJets_Total = GetTotal(DYJets_Files);

  //Loop over the DYJets Files
  for(int i=0; i < DYJets_FileNames.size();i++){
    std::cout<<"DYJetsSample: "<<DYJets_FileNames[i]<<std::endl;
    filenames.push_back(DYJets_FileNames[i]);
    int_lumis.push_back(35900);
    nevents.push_back(DYJets_Total[i]);
    xsecs.push_back(DYJets_xsecs[i]);
  }

  //QCD Samples

  std::vector<string> QCD_FileNames = {"postQCD100to200_0","postQCD200to300_0","postQCD300to500_0","postQCD500to700_0","postQCD700to1000_0","postQCD1000to1500_0","postQCD1500to2000_0","postQCD2000toInf_0"};
  std::vector<TFile *> QCD_Files;
  std::vector<Float_t> QCD_xsecs={27500000,1735000,367000,29370,6524,1064,121.5,25.42};
  for (int i = 0; i < QCD_FileNames.size(); i++) {QCD_Files.push_back(new TFile(TString(QCD_FileNames[i]+".root")));}
  std::vector<float> QCD_Total = GetTotal(QCD_Files);

  //Loop over the QCD Files
  for(int i=0; i < QCD_FileNames.size();i++){
    std::cout<<"QCDSample: "<<QCD_FileNames[i]<<std::endl;
    filenames.push_back(QCD_FileNames[i]);
    int_lumis.push_back(35900);
    nevents.push_back(QCD_Total[i]);
    xsecs.push_back(QCD_xsecs[i]);
  }

//  std::cout<<"-----------------------------"<<std::endl;
//  std::cout<<"j1Pt:"<<std::endl;
//  for(int i = 0; i < filenames.size(); i++){
//    plot("j1Pt",filenames[i],int_lumis[i],nevents[i],xsecs[i]);
//    // DEBUG
//    cout<<filenames[i]<<" finished"<<endl;
//  }
  std::cout<<"-----------------------------"<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<"-----------------------------"<<std::endl;
  std::cout<<"MET:"<<std::endl;
  for(int i = 0; i < filenames.size(); i++){
    plot("MET",filenames[i],int_lumis[i],nevents[i],xsecs[i]);
  }
  std::cout<<"-----------------------------"<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<"-----------------------------"<<std::endl;
}
