//For use with Ntuples made from JetAnalyzer
////Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
////
////To compile using rootcom to an executable named 'analyze':
////$ ./rootcom ZprimeJetsClass_MC_WJets analyze
////
////To run, assuming this is compiled to an executable named 'analyze':
////$ ./analyze /hdfs/store/user/uhussain/Zprime_Ntuples/ /cms/uhussain/MonoZprimeJet/CMSSW_8_0_8/src/LightZPrimeAnalysis/JetAnalyzer/test/output.root -1 10000
////Runs over every event in the folder Zprime_Ntuples, reporting progress every 10000 events
////and storing the resulting histograms in the file output.root.
////
//
#define ZprimeJetsClass_MC_WJets_cxx
#include "ZprimeJetsClass_MC_WJets.h"
#include <TH2.h>
#include<TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
#include <stdio.h> 
using namespace std;
using std::vector;

int main(int argc, const char* argv[])
{ 
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
      return 1;
    }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
      return 1;
    }
  ZprimeJetsClass_MC_WJets t(argv[1],argv[2],atoi(argv[6]),atoi(argv[7]));
  
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void ZprimeJetsClass_MC_WJets::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;   
  // Book histograms for recording properties of leading jet that passes dR and MET cut
  // TFile* histFile = new TFile(file2, "RECREATE");
  // h_deltar = new TH1F("j1deltaR","j1deltaR; #DeltaR of Leading Jet",50,0,0.51);h_deltar->Sumw2();
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Coming in: "<<std::endl;
  std::cout<<"nentries:"<<nentries<<std::endl;
  Long64_t nentriesToCheck = nentries;   

  int initialIndex = -1;//Should end up being 9.
  bool initialIndexNotSet = true;
  int nMCreplicas = 101;
  string initialID = "111"; //111 for basically everything: NNPDF30_lo_as_0130_nf_4 (LHAID 263400)
  std::vector<int> vecIndices;//The vector indices of all the MC replicas. Should just increase sequentially after the initial value.
  vecIndices.clear();
  std::vector<float> sum_passing;//Sum of weights of passing events for each MC replica.
  sum_passing.clear();
  std::vector<int> jetveto;

  double nTotalEvents,nFilters, nHLT, nMET200, nMETcut,nLeptonIDs,nbtagVeto, nDphiJetMET,nJetSelection,Norm,JESUp,JESDo;
  nTotalEvents = nFilters = nHLT = nMET200 = nMETcut = nLeptonIDs = nDphiJetMET = nbtagVeto = nJetSelection = Norm = JESUp = JESDo = 0;

  //getPFCandidates
  std::vector<int>PFCandidates;
 
  //This is the PU histogram obtained from Nick's recipe
  TFile *weights = TFile::Open("PU_Central.root");
  TH1D* PU = (TH1D*)weights->Get("pileup");

  //This is the root file with EWK Corrections
  TFile *file = new TFile("kfactors.root");
  TH1D *ewkCorrection = (TH1D*)file->Get("EWKcorr/W");
  TH1D *NNLOCorrection = (TH1D*)file->Get("WJets_LO/inv_pt");
  float dphimin=-99;
  //Event is rejected if it contains a HighPtMuon faking MET

  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;
  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++) {
    
    //Clear from previous event
    jetveto.clear();  
    jetCand.clear();
    j1PFConsPt.clear();
    j1PFConsEta.clear();
    j1PFConsPhi.clear();
    j1PFConsPID.clear();
    j1PFConsPtUnc.clear();
    EcalCand.clear();
    TrackerCand.clear();
    HcalCand.clear();
    PtFraction_to_use.clear();
    Pt123Fraction_to_use.clear();

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double event_weight=1.0;
    double EWK_corrected_weight=1.0;
    double NNLO_weight = 1.0;
    double kfactor = 1.0;
    
    //For each event we find the bin in the PU histogram that corresponds to puTrue->at(0) and store
    //binContent as event_weight
    int bin = PU->GetXaxis()->FindBin(puTrue->at(0));
    event_weight = PU->GetBinContent(bin);
    //std::cout<<"event_weight: "<<event_weight<<std::endl;

    //PDF/QCD Scale Uncertainity Setup

    if(initialIndexNotSet)
    {
      //Match the initial vector index to the specified MC replica ID.
      for(int i = 0; i < lheWeightIDs->size(); i++)
      {
        if(lheWeightIDs->at(i) == initialID)
        {
          initialIndex = i;
          break;
        }
      }
      //DEBUG
      cout<<"initialIndex = "<<initialIndex<<endl;
      //Set up all the vector indices, and sum_passing while we're at it.
      for(int i = 0; i < nMCreplicas; i++)
      {
        vecIndices.push_back(initialIndex + i);
        sum_passing.push_back(0.0);
      }

      cout<<"lheWeightIDsSize: "<<lheWeightIDs->size()<<endl;
      //Now vecIndices.size() and sum_passing.size() should both == nMCreplicas == 101
      //DEBUG
      for(int i = 0; i < vecIndices.size(); i++)
        cout<<"vector index: "<<vecIndices[i]<<", ID: "<<lheWeightIDs->at(vecIndices[i])<<endl;
      cout<<endl;

      cout<<"vecIndicesSize: "<<vecIndices.size()<<endl;
      //Don't initialize more than once.
      initialIndexNotSet = false;
    }
    int bosonPID;
    double bosonPt;
    bool Wfound = false;
    //check which mc particle is W boson
    for(int i=0; i<nMC;i++){
      if((*mcPID)[i] == 24){
        Wfound=true;
        bosonPID = (*mcPID)[i];
        bosonPt = (*mcPt)[i];
      }
    }
    //if(Wfound){
    ////std::cout<<"bosonPID: "<<bosonPID<<std::endl;
    ////std::cout<<"bosonPt: "<<bosonPt<<std::endl;
    //} 
    jetveto = JetVetoDecision(0);
    jetCand = getJetCand(jetveto,200,2.4,0.8,0.1,0);
    AllPFCand(jetCand,PFCandidates);
    Pt123Fraction_to_use=getPt123Frac(0);
    PtFraction_to_use=getPtFrac(0);
    fillHistos(jetCand,0,event_weight);
    float metcut= 0.0;
    metcut = (fabs(pfMET-caloMET))/pfMET;
    MET_to_use = pfMET;
    METPhi_to_use = pfMETPhi;
    //std::cout<<"|caloMET-pfMET|/pfMET: "<<metcut<<std::endl;
    nTotalEvents+=event_weight;
    if (metFilters==0)
      {
        EWK_corrected_weight = 1.0*(ewkCorrection->GetBinContent(ewkCorrection->GetXaxis()->FindBin(bosonPt)));
        NNLO_weight = 1.0*(NNLOCorrection->GetBinContent(NNLOCorrection->GetXaxis()->FindBin(bosonPt)));
        if(EWK_corrected_weight!=0 && NNLO_weight!=0){
          kfactor = (EWK_corrected_weight/NNLO_weight);}
        else{kfactor=1.21;}
        event_weight*=kfactor;
        nFilters+=event_weight;
	fillHistos(jetCand,1,event_weight); 
	if (true) 
	  {
	    nHLT+=event_weight;
	    fillHistos(jetCand,2,event_weight);
	    if (jetCand.size()>0)
	      {
		nJetSelection+=event_weight;
		fillHistos(jetCand,3,event_weight);
		Norm+=event_weight;
		if (pfMET>250)
		  {
		    nMET200+=event_weight;
		    fillHistos(jetCand,4,event_weight);
		    h_metcut->Fill(metcut,event_weight);
		    if(metcut<0.5)
		      {
			nMETcut+=event_weight;
			fillHistos(jetCand,5,event_weight);
			if(electron_veto_looseID(jetCand[0].first,10) &&  muon_veto_looseID(jetCand[0].first,10))
			  {
			    nLeptonIDs+=event_weight;
			    fillHistos(jetCand,6,event_weight);
			    if(btagVeto(0))
			      {
				nbtagVeto+=event_weight;
				fillHistos(jetCand,7,event_weight);
				double minDPhiJetMET = TMath::Pi();
				double minDPhiJetMET_first4 = TMath::Pi();
				for(int j = 0; j < jetveto.size(); j++)
				  {
				    if(DeltaPhi(jetPhi->at(jetveto[j]),pfMETPhi) < minDPhiJetMET)
				      {
					minDPhiJetMET = DeltaPhi(jetPhi->at(jetveto[j]),pfMETPhi);
					if(j < 4){
					  minDPhiJetMET_first4 = DeltaPhi(jetPhi->at(jetveto[j]),pfMETPhi);}
				      } 
				  }
				h_dphimin->Fill(minDPhiJetMET_first4,event_weight);
				if(dPhiJetMETcut(jetveto,METPhi_to_use))
				  {
				    nDphiJetMET+=event_weight;
				    fillHistos(jetCand,8,event_weight);    
				    if (Pt123Fraction_to_use[0]>0.75)
				      {
					fillHistos(jetCand,9,event_weight);
	        //These histos for pdf uncert calculations is for Pt123Fraction 
          //Start with FillHistos(10)
          //nMCreplicas=101
          for(unsigned int i = 0; i < nMCreplicas; i++)
          {
            fillHistos(jetCand,i+10,event_weight*lheNormalizedWeights->at(vecIndices[i]));
          }
          //Scale variations
          fillHistos(jetCand,nMCreplicas+10,event_weight*genWeight_QCDscale_muR1_muF1); //Default
          fillHistos(jetCand,nMCreplicas+11,event_weight*genWeight_QCDscale_muR1_muF2); //Fac. up
          fillHistos(jetCand,nMCreplicas+12,event_weight*genWeight_QCDscale_muR1_muF0p5); //Fac. down
          fillHistos(jetCand,nMCreplicas+13,event_weight*genWeight_QCDscale_muR2_muF1); //Ren. up
          fillHistos(jetCand,nMCreplicas+14,event_weight*genWeight_QCDscale_muR2_muF2); //Both up
          fillHistos(jetCand,nMCreplicas+15,event_weight*genWeight_QCDscale_muR2_muF0p5);
          fillHistos(jetCand,nMCreplicas+16,event_weight*genWeight_QCDscale_muR0p5_muF1); //Ren. down
          fillHistos(jetCand,nMCreplicas+17,event_weight*genWeight_QCDscale_muR0p5_muF2);
          fillHistos(jetCand,nMCreplicas+18,event_weight*genWeight_QCDscale_muR0p5_muF0p5); //Both down
          //fillHistos(jetCand,101+18=119)
				      }
				  }
			      }   
			  }	
		      }
		  }
	      }
	  }
      }

    tree->Fill();

    if (jentry%reportEvery == 0)
      {
	std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
      }
  
  }
   
  h_cutflow->SetBinContent(1,nTotalEvents); 
  h_cutflow->SetBinContent(2,nFilters);
  h_cutflow->SetBinContent(3,nHLT);
  h_cutflow->SetBinContent(4,nJetSelection);
  h_cutflow->SetBinContent(5,nMET200);
  h_cutflow->SetBinContent(6,nMETcut);
  h_cutflow->SetBinContent(7,nLeptonIDs);
  h_cutflow->SetBinContent(8,nbtagVeto);
  h_cutflow->SetBinContent(9,nDphiJetMET);
   
  //save the histograms
  //   histFile->Write();
  //   histFile->Close();
   
}//Closing the Loop function

void ZprimeJetsClass_MC_WJets::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ZprimeJet","ZprimeJet");
  fileName->cd();

  h_cutflow = new TH1D("h_cutflow","h_cutflow",9,0,9);h_cutflow->Sumw2();
  h_cutflow->GetXaxis()->SetBinLabel(1,"Total Events");
  h_cutflow->GetXaxis()->SetBinLabel(2,"metFilters");
  h_cutflow->GetXaxis()->SetBinLabel(3,"Trigger");
  h_cutflow->GetXaxis()->SetBinLabel(4,"GoodJet");
  h_cutflow->GetXaxis()->SetBinLabel(5,"MetCut");
  h_cutflow->GetXaxis()->SetBinLabel(6,"caloMET cut");
  h_cutflow->GetXaxis()->SetBinLabel(7,"LeptonIDs");
  h_cutflow->GetXaxis()->SetBinLabel(8,"B-JetVeto");
  h_cutflow->GetXaxis()->SetBinLabel(9,"DeltaPhiCut");

  float MtBins[51]={180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		    780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.,3000.};
  
  float MetBins[49]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};

  float PtBins[49]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		    780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};

  double UncBins[59]={0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		    780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};
  h_dphimin = new TH1F("h_dphimin","h_dphimin; Minimum dPhiJetMET",50,0,3.2);h_dphimin->Sumw2();
  h_metcut  = new TH1F("h_metcut","h_metcut; |pfMET-caloMET|/pfMET", 50,0,1.2);h_metcut->Sumw2();
  for(int i=0; i<125; i++){

    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    std::string histname(ptbins);
    h_nJets[i]   = new TH1F(("nJets"+histname).c_str(), "nJets;Number of Jets", 50, 0, 100);h_nJets[i]->Sumw2(); 
    h_pfMETall[i] =  new TH1F(("pfMETall"+histname).c_str(), "pfMET",50,0,2000);h_pfMETall[i] ->Sumw2(); 
    h_pfMET200[i] = new TH1F(("pfMET200"+histname).c_str(), "pfMET",50,170,1500);h_pfMET200[i] ->Sumw2(); 
    h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "E_{T}^{miss} (GeV)",48,MetBins);h_pfMET[i] ->Sumw2();
    h_pfMETPhi[i] = new TH1F(("pfMETPhi"+histname).c_str(), "pfMETPhi",50,-4,4);h_pfMETPhi[i]->Sumw2();
    h_j1Pt[i]  = new TH1F(("j1pT"+histname).c_str(), "j1pT;p_{T} of Leading Jet (GeV)", 48,PtBins);h_j1Pt[i]->Sumw2();
    h_j1Eta[i] = new TH1F(("j1Eta"+histname).c_str(), "j1Eta; #eta of Leading Jet", 50, -2.5, 2.5);h_j1Eta[i]->Sumw2();
    h_j1Phi[i] = new TH1F(("j1Phi"+histname).c_str(), "j1Phi; #phi of Leading Jet", 50, -3.0, 3.0);h_j1Phi[i]->Sumw2();     
    h_j1etaWidth[i] = new TH1F(("j1etaWidth"+histname).c_str(),"j1etaWidh; #eta width of Leading Jet", 50,0,0.25);h_j1etaWidth[i] ->Sumw2();
    h_j1phiWidth[i] = new TH1F(("j1phiWidth"+histname).c_str(),"j1phiWidth; #phi width of Leading Jet", 50, 0,0.5);h_j1phiWidth[i]->Sumw2();
    h_j1nCons[i] = new TH1F (("j1nCons"+histname).c_str(),"j1NConstituents; Number of Constituents of Leading Jet",25, 0, 50);h_j1nCons[i]->Sumw2();
    h_j1nCategory1[i] = new TH1F(("j1nCategory1"+histname).c_str(),"j1nCategory1: Number of events with exactly two charged Hadrons",2,-0.5,1.5);h_j1nCategory1[i]->Sumw2();
    h_j1nCategory2[i] = new TH1F(("j1nCategory2"+histname).c_str(),"j1nCategory2: Number of events with two charged Hadrons and one Photon",2,-0.5,1.5);h_j1nCategory2[i]->Sumw2(); 
    h_PF123PtFraction[i]= new TH1F(("PF123PtFraction"+histname).c_str(), "PF123PtFraction;P_{T} fraction carried by 3 leading daughters of the Pencil Jet" ,50,0,1.0);h_PF123PtFraction[i]->Sumw2();   
    h_j1PF12PtFrac_ID_1[i]= new TH1F(("j1PF12PtFrac_ID_1"+histname).c_str(), "j1PF12PtFrac_ID_1;P_{T} fraction carried by charged hadrons of the leading Jet: Category 1" ,50,0,1.0);h_j1PF12PtFrac_ID_1[i]->Sumw2();   
    h_j1dRPF12_ID_1[i] = new TH1F(("j1dRPF12_ID_1"+histname).c_str(),"j1dRPF12_ID_1; deltaR between oppositely charged hadrons of the leading Jet: Category 1",50,0,0.15);h_j1dRPF12_ID_1[i]->Sumw2();
    h_j1PF12PtFrac_ID_2[i]= new TH1F(("j1PF12PtFrac_ID_2"+histname).c_str(), "j1PF12PtFrac_ID_2;P_{T} fraction carried by charged hadrons of the leading Jet: Category 2" ,50,0,1.0);h_j1PF12PtFrac_ID_2[i]->Sumw2();
    h_j1dRPF12_ID_2[i] = new TH1F(("j1dRPF12_ID_2"+histname).c_str(),"j1dRPF12_ID_2; deltaR between oppositely charged hadrons of the leading Jet: Category 2",50,0,0.15);h_j1dRPF12_ID_2[i]->Sumw2();
    h_j1PFPtFrac_ID_2[i] = new TH1F(("j1PFPtFrac_ID_2"+histname).c_str(),"j1PFPtFrac_ID_2;P_{T} fraction carried by charged hadrons+Photon of the leading Jet: Category 2" ,50,0,1.0);h_j1PFPtFrac_ID_2[i]->Sumw2();  
    h_j1TotPFCands[i] = new TH1F(("j1TotPFCands"+histname).c_str(),"j1TotPFCands;# of all PF candidates in Leading Jet",25,0,50);h_j1TotPFCands[i]->Sumw2();
    h_j1ChPFCands[i] = new TH1F(("j1ChPFCands"+histname).c_str(),"j1ChPFCands;# of PF charged hadrons in Leading Jet",25,0,50);h_j1ChPFCands[i]->Sumw2();
    h_j1NeutPFCands[i] = new TH1F(("j1NeutPFCands"+histname).c_str(),"j1NeutPFCands;# of PF neutral hadrons in Leading Jet",15,0,15);h_j1NeutPFCands[i]->Sumw2();
    h_j1GammaPFCands[i] = new TH1F(("j1GammaPFCands"+histname).c_str(),"j1GammaPFCands;# of PF gammas in Leading Jet",20,0,40);h_j1GammaPFCands[i]->Sumw2();
    h_j1CHF[i] = new TH1F(("j1CHF"+histname).c_str(),"j1CHF;Charged Hadron Energy Fraction in Leading Jet",50,0,1.0);h_j1CHF[i]->Sumw2(); 
    h_j1NHF[i] = new TH1F(("j1NHF"+histname).c_str(),"j1NHF;Neutral Hadron Energy Fraction in Leading Jet",50,0,1.0);h_j1NHF[i]->Sumw2(); 
    h_j1ChMultiplicity[i] = new TH1F(("j1ChMultiplicity"+histname).c_str(),"j1ChMultiplicity;Charged Multiplicity of Leading Jet",25,0,50);h_j1ChMultiplicity[i]->Sumw2();
    h_j1NeutMultiplicity[i] = new TH1F(("j1NeutMultiplicity"+histname).c_str(),"j1NeutMultiplicity;Neutral Multiplicity of Leading Jet",25,0,50);h_j1NeutMultiplicity[i]->Sumw2(); 
    h_j1Mt[i]  = new TH1F(("j1Mt"+histname).c_str(), "j1Mt;M_{T} of Leading Jet (GeV)", 50,MtBins);h_j1Mt[i]->Sumw2(); 
    h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(),"nVtx;nVtx",70,0,70);h_nVtx[i]->Sumw2();
    //Category3 variables
    h_ChPionPt[i]=new TH1F(("ChPionPt"+histname).c_str(),"ChPionPt;p_{T} of Charged Pion in 3rd Signal Category",50,0,2000);h_ChPionPt[i]->Sumw2();
    h_PhotonPt[i]=new TH1F(("PhotonPt"+histname).c_str(),"PhotonPt;p_{T} of Photon in 3rd Signal Category",50,0,2000);h_PhotonPt[i]->Sumw2();
    h_dRPionPhoton[i]=new TH1F(("dRPionPhoton"+histname).c_str(),"dRPionPhoton;deltaR between ChPion and Photon 3rd Signal Category",50,0,0.5);h_dRPionPhoton[i]->Sumw2(); 
    h_EcalPtUnc[i]=new TH2F(("EcalPtUnc"+histname).c_str(),"ECAL P_{T} Uncertainty;Photon P_{T} (GeV);Uncertainty",58,UncBins,50,0.,1.);
    h_TrackerPtUnc[i]=new TH2F(("TrackerPtUnc"+histname).c_str(),"Tracker P_{T} Uncertainty;Charged Hadrons P_{T} (GeV);Uncertainty",58,UncBins,50,0.,1.);
    h_HcalPtUnc[i]=new TH2F(("HcalPtUnc"+histname).c_str(),"HCAL P_{T} Uncertainty;Neutral Hadron P_{T} (GeV);Uncertainty",58,UncBins,50,0.,1.);
    h_TrackerPtFrac[i]=new TH1F(("TrackerPtFraction"+histname).c_str(), "TrackerPtFraction;P_{T}^{123} fraction with TrackerUnc" ,50,0,1.0);h_TrackerPtFrac[i]->Sumw2();
    h_EcalPtFrac[i]=new TH1F(("EcalPtFraction"+histname).c_str(), "EcalPtFraction;P_{T}^{123} fraction with EcalUnc" ,50,0,1.0);h_EcalPtFrac[i]->Sumw2();
    h_HcalPtFrac[i]=new TH1F(("HcalPtFraction"+histname).c_str(), "HcalPtFraction;P_{T}^{123} fraction with HcalUnc" ,50,0,1.0);h_HcalPtFrac[i]->Sumw2();
    //PFCons Fractions
    h_j1PFCons1PtFraction[i] = new TH1F(("h_j1PFCons1PtFraction"+histname).c_str(),"h_j1PFCons1PtFraction; p_{T} Fraction of PFConsNo.1 of Pencil Jet",50,0,1.0);h_j1PFCons1PtFraction[i]->Sumw2();
    h_j1PFCons2PtFraction[i] = new TH1F(("h_j1PFCons2PtFraction"+histname).c_str(),"h_j1PFCons2PtFraction; p_{T} Fraction of PFConsNo.2 of Pencil Jet",50,0,1.0);h_j1PFCons2PtFraction[i]->Sumw2();
    h_j1PFCons3PtFraction[i] = new TH1F(("h_j1PFCons3PtFraction"+histname).c_str(),"h_j1PFCons3PtFraction; p_{T} Fraction of PFConsNo.3 of Pencil Jet",50,0,1.0);h_j1PFCons3PtFraction[i]->Sumw2();
  }
}

//double ZprimeJetsClass_MC_WJets::dR(double jetetaWidth, double jetphiWidth)
//{
//  double deltar = sqrt(jetetaWidth*jetetaWidth + jetphiWidth*jetphiWidth);
//  return deltar;
//}


void ZprimeJetsClass_MC_WJets::fillHistos(std::vector<std::pair<int,double> > jetCand_to_use,int histoNumber,double event_weight)
{
  h_nVtx[histoNumber]->Fill(nVtx,event_weight);
  h_nJets[histoNumber]->Fill(nJet,event_weight);
  h_pfMETall[histoNumber]->Fill(pfMET,event_weight);
  h_pfMET200[histoNumber]->Fill(pfMET,event_weight);
  h_pfMET[histoNumber]->Fill(pfMET,event_weight);
  h_pfMETPhi[histoNumber]->Fill(pfMETPhi,event_weight);
  if(jetCand_to_use.size()>0){
    h_j1Pt[histoNumber]->Fill(jetCand_to_use[0].second,event_weight);
    h_j1Eta[histoNumber]->Fill(jetEta->at(jetCand_to_use[0].first),event_weight);
    h_j1Phi[histoNumber]->Fill(jetPhi->at(jetCand_to_use[0].first),event_weight); 
    h_j1nCategory1[histoNumber]->Fill(TwoChPFCons,event_weight); 
    h_j1nCategory2[histoNumber]->Fill(TwoChPFConsPlusPho,event_weight); 
    h_PF123PtFraction[histoNumber]->Fill(Pt123Fraction_to_use[0],event_weight);
    h_j1PF12PtFrac_ID_1[histoNumber]->Fill(PF12PtFrac_ID_1,event_weight);
    h_j1dRPF12_ID_1[histoNumber]->Fill(dR_PF12_ID_1,event_weight);
    h_j1PF12PtFrac_ID_2[histoNumber]->Fill(PF12PtFrac_ID_2,event_weight);
    h_j1dRPF12_ID_2[histoNumber]->Fill(dR_PF12_ID_2,event_weight);
    h_j1PFPtFrac_ID_2[histoNumber]->Fill(PF123PtFrac_ID_2,event_weight);
    h_j1TotPFCands[histoNumber]->Fill(TotalPFCandidates,event_weight);
    h_j1ChPFCands[histoNumber]->Fill(ChargedPFCandidates,event_weight);
    h_j1NeutPFCands[histoNumber]->Fill(NeutralPFCandidates,event_weight);
    h_j1GammaPFCands[histoNumber]->Fill(GammaPFCandidates,event_weight); 
    h_j1CHF[histoNumber]->Fill(jetCHF->at(jetCand_to_use[0].first),event_weight);
    h_j1NHF[histoNumber]->Fill(jetNHF->at(jetCand_to_use[0].first),event_weight);
    h_j1ChMultiplicity[histoNumber]->Fill(jetNCH->at(jetCand_to_use[0].first),event_weight);
    h_j1NeutMultiplicity[histoNumber]->Fill(jetNNP->at(jetCand_to_use[0].first),event_weight);
    h_j1Mt[histoNumber]->Fill(jetMt->at(jetCand_to_use[0].first),event_weight);
    h_j1etaWidth[histoNumber]->Fill(jetetaWidth->at(jetCand_to_use[0].first),event_weight);
    h_j1phiWidth[histoNumber]->Fill(jetphiWidth->at(jetCand_to_use[0].first),event_weight);
    h_j1nCons[histoNumber]->Fill((jetnPhotons->at(jetCand_to_use[0].first)+jetnCHPions->at(jetCand_to_use[0].first)+jetnMisc->at(jetCand_to_use[0].first)),event_weight);
    //These 3 histograms store Pt123Fraction but with each Unc applied separately
    h_TrackerPtFrac[histoNumber]->Fill(Pt123Fraction_to_use[1],event_weight);
    h_EcalPtFrac[histoNumber]->Fill(Pt123Fraction_to_use[2],event_weight);
    h_HcalPtFrac[histoNumber]->Fill(Pt123Fraction_to_use[3],event_weight);
    //And for these 2D plots we only fill with one of the 3 leading Constituents based on which category it falls in
    for(int i=0;i<TrackerCand.size();i++)
      {
        if (j1PFConsPt.at(TrackerCand[i]) > 1. && TrackerCand[i] < 3)
            {
              h_TrackerPtUnc[histoNumber]->Fill(j1PFConsPt.at(TrackerCand[i]),j1PFConsPtUnc.at(TrackerCand[i]),event_weight);
            }
      }
    for(int i=0;i<EcalCand.size();i++)
      {
        if (j1PFConsPt.at(EcalCand[i]) > 1. && EcalCand[i] < 3)
            {
              h_EcalPtUnc[histoNumber]->Fill(j1PFConsPt.at(EcalCand[i]),j1PFConsPtUnc.at(EcalCand[i]),event_weight);
            }
      }
    for(int i=0;i<HcalCand.size();i++)
      {
        if (j1PFConsPt.at(HcalCand[i]) > 1. && HcalCand[i] < 3)
            {
              h_HcalPtUnc[histoNumber]->Fill(j1PFConsPt.at(HcalCand[i]),j1PFConsPtUnc.at(HcalCand[i]),event_weight);
            }
      }
    //PtFraction of the first 3 Const.
    h_j1PFCons1PtFraction[histoNumber]->Fill(PtFraction_to_use[1],event_weight);
    h_j1PFCons2PtFraction[histoNumber]->Fill(PtFraction_to_use[2],event_weight);
    h_j1PFCons3PtFraction[histoNumber]->Fill(PtFraction_to_use[3],event_weight);
  }
}

std::vector<double>ZprimeJetsClass_MC_WJets::getPt123Frac(int UncType)
{
    std::vector<double> Pt123FractionToUse;
    Pt123FractionToUse.clear();
    vector<double> Pt123 = {0.0,0.0,0.0,0.0};
    vector<double> jetPtAll = {0.0,0.0,0.0,0.0};
    vector<int> index(j1PFConsPID.size());
    iota(begin(index),end(index),0);
    vector<vector<int>> ConsIndex = {index,TrackerCand,EcalCand,HcalCand};
    //vector<string> ConsName = {"All","Tracker","Ecal","Hcal"};
    for (int i = 0; i < ConsIndex.size(); i++)
        {
        for (int j = 0; j < j1PFConsPID.size(); j++)
          {
          if (find(ConsIndex[i].begin(),ConsIndex[i].end(),j) != ConsIndex[i].end())
              {
                jetPtAll[i]+=j1PFConsPt.at(j)+UncType*j1PFConsPtUnc.at(j);
                if (j < 3) Pt123[i]+=j1PFConsPt.at(j)+UncType*j1PFConsPtUnc.at(j);
              }
          else
            {
              jetPtAll[i]+=j1PFConsPt.at(j);
              if (j < 3) Pt123[i]+=j1PFConsPt.at(j);
            }
          }
        Pt123FractionToUse.push_back((Pt123[i]/jetPtAll[i]));
      }
    return Pt123FractionToUse;
}
std::vector<double>ZprimeJetsClass_MC_WJets::getPtFrac(int UncType)
{

  std::vector<double> PtFractionToUse;
  PtFractionToUse.clear();
  double Pt1Cons,Pt2Cons,Pt3Cons,Pt123Cons;
  Pt1Cons=Pt2Cons=Pt3Cons=Pt123Cons=0;
  double jetPtAll=0.0;
  for (int i = 0; i < j1PFConsPID.size(); i++)
    {
      jetPtAll+=j1PFConsPt.at(i)+UncType*j1PFConsPtUnc.at(i);
      if (i < 3) Pt123Cons+=j1PFConsPt.at(i)+UncType*j1PFConsPtUnc.at(i);
      if (i==0) Pt1Cons=j1PFConsPt.at(0)+UncType*j1PFConsPtUnc.at(0);
      if (i==1) Pt2Cons+=j1PFConsPt.at(1)+UncType*j1PFConsPtUnc.at(1);
      if (i==2) Pt3Cons+=j1PFConsPt.at(2)+UncType*j1PFConsPtUnc.at(2);
    }
  PtFractionToUse.push_back((Pt123Cons/jetPtAll));//Index 0 of the PtFractionToUse Vector
  PtFractionToUse.push_back((Pt1Cons/jetPtAll));//Index 1 of the PtFractionToUse Vector
  PtFractionToUse.push_back((Pt2Cons/jetPtAll));//Index 2 of the PtFractionToUse Vector
  PtFractionToUse.push_back((Pt3Cons/jetPtAll));//Index 3 of the PtFractionToUse Vector
  
  return PtFractionToUse;
}

void ZprimeJetsClass_MC_WJets::AllPFCand(std::vector<std::pair<int,double>> jetCand,std::vector<int>PFCandidates)
{
  //getPFCandidatesMethod for the Pencil Jet -> jetCand[0].first
    TotalPFCandidates=ChargedPFCandidates=NeutralPFCandidates=GammaPFCandidates=0;
    PFCandidates = getPFCandidates();
    //std::cout<<"Vector of Pairs should have size 4: "<<PFCandidates.size()<<std::endl;
    if(PFCandidates.size()>0){
      TotalPFCandidates=PFCandidates.at(0);}
    //std::cout<<"TotalPFCandidates: "<<TotalPFCandidates<<std::endl;}

    if(PFCandidates.size()>1){
      ChargedPFCandidates=PFCandidates.at(1);}
    //std::cout<<"TotalChargedPFCandidates: "<<ChargedPFCandidates<<std::endl;}
    
    if(PFCandidates.size()>2){
      GammaPFCandidates=PFCandidates.at(2);}
    //std::cout<<"TotalGammaPFCandidates: "<<GammaPFCandidates<<std::endl;}

    if(PFCandidates.size()>3){
      NeutralPFCandidates=PFCandidates.at(3);}
    //std::cout<<"TotalNeutralPFCandidates: "<<NeutralPFCandidates<<std::endl;}
    //We are using these conditions so we only calculate the following quantities for the signal we are interested in
    //This will also make it faster to process the events
    if(pfMET>250 && jetCand.size()>0){
      j1PFConsPt=JetsPFConsPt->at(jetCand[0].first);
      j1PFConsEta=JetsPFConsEta->at(jetCand[0].first);
      j1PFConsPhi=JetsPFConsPhi->at(jetCand[0].first);
      j1PFConsPID=JetsPFConsPID->at(jetCand[0].first);
      for(int i=0;i<j1PFConsPID.size();i++)
	{
	  if (abs(j1PFConsPID.at(i)) == 211 || abs(j1PFConsPID.at(i)) == 13)
	    {
	      //Tracker Uncertainty
	      //deltaPt=(1/100)*sqrt((0.015*Pt)^2+(0.5)^2)
	      j1PFConsPtUnc.push_back((1/100.)*sqrt(pow(0.015*j1PFConsPt.at(i),2)+pow(0.5,2)));
	      TrackerCand.push_back(i);
	    }
	  else if (abs(j1PFConsPID.at(i)) == 22 || abs(j1PFConsPID.at(i)) == 11)
	    {
	      //ECAL Uncertainty
	      //deltaPt=(1/100)*sqrt((2.8)^2/Pt+(12.8/Pt)^2+(0.3)^2)
	      j1PFConsPtUnc.push_back((1/100.)*sqrt((pow(2.8,2)/j1PFConsPt.at(i))+(pow(12.8,2)/pow(j1PFConsPt.at(i),2))+pow(0.3,2)));
	      EcalCand.push_back(i);
	    }
	  else if (abs(j1PFConsPID.at(i)) == 130)
	    {
	      //HCAL Uncertainty
	      //deltaPt=(1/100)*sqrt((115)^2/Pt+(5.5)^2)
	      j1PFConsPtUnc.push_back((1/100.)*sqrt(pow(115,2)/j1PFConsPt.at(i)+pow(5.5,2)));
	      HcalCand.push_back(i);
	    }
	  else
	    {
	      j1PFConsPtUnc.push_back(0);
	    }
	}
    }
}

//Function to calculate regular deltaR separate from jet width variable 'dR'
double ZprimeJetsClass_MC_WJets::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//Gives the (minimum) separation in phi between the specified phi values
////Must return a positive value
float ZprimeJetsClass_MC_WJets::DeltaPhi(float phi1, float phi2)
{
  float pi = TMath::Pi();
  float dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}

float ZprimeJetsClass_MC_WJets::dPhiJetMETmin(std::vector<int> jets)
{
  float dPhimin=TMath::Pi();
  int njetsMax = jets.size();
  if(njetsMax > 4)
    njetsMax = 4; 
  for(int j=0;j< njetsMax; j++)
    {
      float dPhi = DeltaPhi((*jetPhi)[j],METPhi_to_use);
      //std::cout<<"DeltaPhi: "<<dPhi<<std::endl;
      if(dPhi < dPhimin){
        dPhimin = dPhi;
      }
    }
  return dPhimin;
}

std::vector<std::pair<int,double>> ZprimeJetsClass_MC_WJets::getJetCand(std::vector<int> jets,double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut,int UncType){

  //save the Pt of the jetCand as well, whether normal, shiftedUp or shiftedDown 
  std::vector<std::pair<int,double>> tmpCand;
  tmpCand.clear();
  //So only check if leading candidate will pass these cuts?
  int njets = jets.size();
  //for(int p=0;p<njets;p++)
    //{
      //UncType = +1(up), or -1(down) or 0(normal)
  if(njets>0){
    int p=jets[0];
    Float_t jetPt_to_use;
    if (UncType == 0){jetPt_to_use = (*jetPt)[p];}
    else if (UncType == 1){jetPt_to_use = (*jetPt)[p]*(1.+(*jetJECUnc)[p]);}
    else if (UncType == -1){jetPt_to_use = (*jetPt)[p]*(1.-(*jetJECUnc)[p]);}
      
    bool kinematic = jetPt_to_use > jetPtCut && (*jetNHF)[p] < jetNHFCut && (*jetCHF)[p] > jetCHFCut && fabs((*jetEta)[p])<jetEtaCut;

    if((*jetPFLooseId)[p]==1 && kinematic){
      tmpCand.push_back(std::make_pair(p,jetPt_to_use));
      }
    }

  return tmpCand;

}
std::vector<int> ZprimeJetsClass_MC_WJets::JetVetoDecision(int UncType) {

  bool jetVeto=true;
  std::vector<int> jetindex;

  for(int i = 0; i < nJet; i++)
    {
      Float_t jetPt_to_use;
      if (UncType == 0){jetPt_to_use = (*jetPt)[i];}
      else if (UncType == 1){jetPt_to_use = (*jetPt)[i]*(1+(*jetJECUnc)[i]);}
      else if (UncType == -1){jetPt_to_use =(*jetPt)[i]*(1-(*jetJECUnc)[i]);}
      if(jetPt_to_use >30.0 && jetPFLooseId->at(i)==1)
        {
          jetindex.push_back(i);
        }


    }

  return jetindex;

}

//Return a vector of pairs. "0" = #pfCands, "1"=#chargedPFCands , "3"=#neutralPFCands,"2"=#photonPFCands
std::vector<int>ZprimeJetsClass_MC_WJets::getPFCandidates(){
  std::vector<int>PFCands;
  for(int i=0;i<nJet;i++)
    {
      int TotPFCands;
      if(i==0){
	TotPFCands = j1PFConsPID.size();
	//std::cout<<TotPFCands<<std::endl;
	PFCands.push_back(TotPFCands);
	int ChPFCands,NeuPFCands,GammaPFCands;
	ChPFCands=NeuPFCands=GammaPFCands=0;
	for(int j=0;j<TotPFCands;j++){
	  if(abs(j1PFConsPID.at(j))==211){
	    ChPFCands++;
	  }
	  if(j1PFConsPID.at(j)==22){
	    GammaPFCands++;
	  }
	  if(j1PFConsPID.at(j)==130){
	    NeuPFCands++;
	  }
	}
	PFCands.push_back(ChPFCands);
	PFCands.push_back(GammaPFCands);
	PFCands.push_back(NeuPFCands);
      }
    }
  return PFCands;
}
bool ZprimeJetsClass_MC_WJets::btagVeto(int UncType) {

  bool btagVeto = true;
  for(int i = 0; i < nJet; i++)
    {
      Float_t jetPt_to_use;
      if (UncType == 0){jetPt_to_use = (*jetPt)[i];}
      else if (UncType == 1){jetPt_to_use = (*jetPt)[i]*(1+(*jetJECUnc)[i]);}
      else if (UncType == -1){jetPt_to_use =(*jetPt)[i]*(1-(*jetJECUnc)[i]);}

      if(jetPt_to_use >20.0 && jetEta->at(i) < 2.4 && jetCSV2BJetTags->at(i) > 0.8)
        btagVeto = false;
    }
  return btagVeto;
}

bool ZprimeJetsClass_MC_WJets::dPhiJetMETcut(std::vector<int> jets, float METPhi)
{
  //reject jet if it is found within DeltaPhi(jet,MET) < 0.5                                                                                              \
  
  bool passes = false;

  int njetsMax = jets.size();
  //std::cout<<"njets: "<<njetsMax<<std::endl;
  //Only look at first four jets (because that's what monojet analysis do)
  if(njetsMax > 4)
    njetsMax = 4;
  int j=0;
  for(;j< njetsMax; j++){
    //std::cout<<"DeltaPhi b/w Jet and MET"<<std::endl;
    //std::cout<<"jet "<<j<<":"<<DeltaPhi((*jetPhi)[j],pfMETPhi)<<std::endl;
    if(DeltaPhi((*jetPhi)[j],METPhi) < 0.5)
      break;
  }

  if(j==njetsMax)
    passes = true;

  return passes;

}
bool ZprimeJetsClass_MC_WJets::electron_veto_looseID(int jet_index, float elePtCut)
{
  bool veto_passed = true; //pass veto if no good electron found
  
  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false; 
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue                                                                                               
  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
    {
      pass_SigmaIEtaIEtaFull5x5 = false;
      pass_dEtaIn = false;
      pass_dPhiIn = false;
      pass_HoverE = false;
      pass_iso = false;
      pass_ooEmooP = false;
      pass_d0 = false;
      pass_dz = false;
      pass_missingHits = false;
      pass_convVeto = false;
      //Find EA for corrected relative iso.     
      if(abs(eleSCEta->at(i)) <= 1.0)
        EA = 0.1752;
      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
        EA = 0.1862;
      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
        EA = 0.1411;
      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
        EA = 0.1534;
      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
        EA = 0.1903;
      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
        EA = 0.2243;
      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
        EA = 0.2687;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));
      if(abs(eleSCEta->at(i)) <= 1.479)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
          pass_HoverE = eleHoverE->at(i) < 0.104;
          pass_iso = EAcorrIso < 0.0893;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
          pass_d0 = abs(eleD0->at(i)) < 0.0261;
          pass_dz = abs(eleDz->at(i)) < 0.41;
          pass_missingHits = eleMissHits->at(i) <= 2;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
          pass_HoverE = eleHoverE->at(i) < 0.0897;
          pass_iso = EAcorrIso < 0.121;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
          pass_d0 = abs(eleD0->at(i)) < 0.118;
          pass_dz = abs(eleDz->at(i)) < 0.822;
          pass_missingHits = eleMissHits->at(i) <= 1;
          pass_convVeto = eleConvVeto->at(i) == 1;
        } 
      //Electron passes Loose Electron ID cuts  
      if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
        { 
          //Electron passes pt cut                                                                                                                                
          if(elePt->at(i) > elePtCut)
            {
              //Electron does not overlap with jet
              if(deltaR(eleSCEta->at(i),eleSCPhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.5)
                {
                  veto_passed = false;
                  break;
                }
            }   
        }         
    }             
  return veto_passed;
}                  
bool ZprimeJetsClass_MC_WJets::muon_veto_looseID(int jet_index, float muPtCut)
{
  bool veto_passed = true; //pass veto if no good muon found 
  bool pass_iso = false;  
                                                                                                                                                        
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
    
  for(int i = 0; i < nMu; i++)
    {
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;

      if(muPt->at(i) > muPtCut)
        {
          if(pass_iso)
            {
	      //muon does not overlap with jet
              if(deltaR(muEta->at(i),muPhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.5)
                {
                  veto_passed = false;
                  break;
                }
            }
        }
    }      
  return veto_passed;
}                
  
