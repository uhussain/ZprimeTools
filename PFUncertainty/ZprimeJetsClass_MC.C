//For use with Ntuples made from JetAnalyzer
////Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
////
////To compile using rootcom to an executable named 'analyze':
////$ ./rootcom ZprimeJetsClass_MC analyze
////
////To run, assuming this is compiled to an executable named 'analyze':
////$ ./analyze /hdfs/store/user/uhussain/Zprime_Ntuples/ /cms/uhussain/MonoZprimeJet/CMSSW_8_0_8/src/LightZPrimeAnalysis/JetAnalyzer/test/output.root -1 10000
////Runs over every event in the folder Zprime_Ntuples, reporting progress every 10000 events
////and storing the resulting histograms in the file output.root.
////
//
#define ZprimeJetsClass_MC_cxx
#include "ZprimeJetsClass_MC.h"
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
  ZprimeJetsClass_MC t(argv[1],argv[2],atoi(argv[6]),atoi(argv[7]));
  
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void ZprimeJetsClass_MC::Loop(Long64_t maxEvents, int reportEvery)
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

  std::vector<int> jetveto;

  double nTotalEvents,nFilters, nHLT, nMET200, nMETcut,nLeptonIDs,nbtagVeto, nDphiJetMET,nJetSelection;
  nTotalEvents = nFilters = nHLT = nMET200 = nMETcut = nLeptonIDs = nDphiJetMET = nbtagVeto = nJetSelection = 0;

  //getPFCandidates
  std::vector<int>PFCandidates;
   
  //This is the PU histogram obtained from Nick's recipe
  TFile *weights = TFile::Open("PU_Central.root");
  TH1D* PU = (TH1D*)weights->Get("pileup");


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

    Pt123Fraction_to_use = {-1,-1,-1,-1};
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //double kfactor = 1.0;

    double event_weight=1.0;
    //For each event we find the bin in the PU histogram that corresponds to puTrue->at(0) and store
    //binContent as event_weight
    int bin = PU->GetXaxis()->FindBin(puTrue->at(0));
    event_weight = PU->GetBinContent(bin);
    //std::cout<<"event_weight: "<<event_weight<<std::endl;
    jetveto = JetVetoDecision(0);
    jetCand = getJetCand(jetveto,200,2.4,0.8,0.1,0);
    AllPFCand(jetCand,PFCandidates);
    getPt123Frac(0);
    fillHistos(jetCand,0,event_weight);
    float metcut= 0.0;
    metcut = (fabs(pfMET-caloMET))/pfMET;
    MET_to_use = pfMET;
    METPhi_to_use = pfMETPhi;
    //std::cout<<"|caloMET-pfMET|/pfMET: "<<metcut<<std::endl;
    nTotalEvents+=event_weight;
    if (metFilters==0)
      { 
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
				    
				    //Category 1: Exactly Two Charged Hadrons
				    if(TwoChPFCons==1)
				      {
					fillHistos(jetCand,9,event_weight);
					//Effectiveness of this cut in Category 1 Events
					if(PF12PtFrac_ID_1>0.7)
					  {
					    fillHistos(jetCand,10,event_weight);}
				      } 
				    //Category 2: Two charged Hadrons + Photon
				    if(TwoChPFConsPlusPho==1)
				      {
					fillHistos(jetCand,11,event_weight);
					//Effectiveness of this cut in Category 2 Events
					if(PF123PtFrac_ID_2>0.7)
					  {
					    fillHistos(jetCand,12,event_weight);}
				      }
				    //Category of events with < 2 charged Hadrons
				    if(TwoChPFCons==0 && TwoChPFConsPlusPho==0)
				      {
					fillHistos(jetCand,13,event_weight);
					//Calculating the effectiveness of this cut in only events with < 2 oppositely charged Hadrons
					if(jetetaWidth->at(jetCand[0].first)<0.04)
					  {
					    fillHistos(jetCand,14,event_weight);
					  }}
				    //This is for comparison with previous results (for all events)
				    if (jetetaWidth->at(jetCand[0].first)<0.04)
				      {
					fillHistos(jetCand,15,event_weight);
				      }
				    if (Pt123Fraction_to_use[0]>0.8)
				      {
					fillHistos(jetCand,16,event_weight);
				      }
				    if (Pt123Fraction_to_use[0]>0.85)
				      {
					fillHistos(jetCand,17,event_weight);
				      }
				    if (Pt123Fraction_to_use[0]>0.9)
				      {
					fillHistos(jetCand,18,event_weight);
				      }
				    

				    getPt123Frac(1);
				    fillHistos(jetCand,19,event_weight);
				    //Category 1: Exactly Two Charged Hadrons
				    if(TwoChPFCons==1)
				      {
					fillHistos(jetCand,20,event_weight);
					//Effectiveness of this cut in Category 1 Events
					if(PF12PtFrac_ID_1>0.7)
					  {
					    fillHistos(jetCand,21,event_weight);}
				      } 
				    //Category 2: Two charged Hadrons + Photon
				    if(TwoChPFConsPlusPho==1)
				      {
					fillHistos(jetCand,22,event_weight);
					//Effectiveness of this cut in Category 2 Events
					if(PF123PtFrac_ID_2>0.7)
					  {
					    fillHistos(jetCand,23,event_weight);}
				      }
				    //Category of events with < 2 charged Hadrons
				    if(TwoChPFCons==0 && TwoChPFConsPlusPho==0)
				      {
					fillHistos(jetCand,24,event_weight);
					//Calculating the effectiveness of this cut in only events with < 2 oppositely charged Hadrons
					if(jetetaWidth->at(jetCand[0].first)<0.04)
					  {
					    fillHistos(jetCand,25,event_weight);
					  }}
				    //This is for comparison with previous results (for all events)
				    if (jetetaWidth->at(jetCand[0].first)<0.04)
				      {
					fillHistos(jetCand,26,event_weight);
				      }

				    if (Pt123Fraction_to_use[0] > 0.8)
				      {//Correlated Uncertainty
					fillHistos(jetCand,27,event_weight);
				      }
				    if (Pt123Fraction_to_use[1] > 0.8)
				      {//Tracker Uncertainty
					fillHistos(jetCand,28,event_weight);
				      }
				    if (Pt123Fraction_to_use[2] > 0.8)
				      {//Ecal Uncertainty
					fillHistos(jetCand,29,event_weight);
				      }
				    if (Pt123Fraction_to_use[3] > 0.8)
				      {//Hcal Uncertainty
					fillHistos(jetCand,30,event_weight);
				      }

				    if (Pt123Fraction_to_use[0] > 0.85)
				      {//Correlated Uncertainty
					fillHistos(jetCand,31,event_weight);
				      }
				    if (Pt123Fraction_to_use[1] > 0.85)
				      {//Tracker Uncertainty
					fillHistos(jetCand,32,event_weight);
				      }
				    if (Pt123Fraction_to_use[2] > 0.85)
				      {//Ecal Uncertainty
					fillHistos(jetCand,33,event_weight);
				      }
				    if (Pt123Fraction_to_use[3] > 0.85)
				      {//Hcal Uncertainty
					fillHistos(jetCand,34,event_weight);
				      }
				    if (Pt123Fraction_to_use[0] > 0.9)
				      {//Correlated Uncertainty
					fillHistos(jetCand,35,event_weight);
				      }
				    if (Pt123Fraction_to_use[1] > 0.9)
				      {//Tracker Uncertainty
					fillHistos(jetCand,36,event_weight);
				      }
				    if (Pt123Fraction_to_use[2] > 0.9)
				      {//Ecal Uncertainty
					fillHistos(jetCand,37,event_weight);
				      }
				    if (Pt123Fraction_to_use[3] > 0.9)
				      {//Hcal Uncertainty
					fillHistos(jetCand,38,event_weight);
				      }

				    getPt123Frac(-1);
				    fillHistos(jetCand,39,event_weight);
				    //Category 1: Exactly Two Charged Hadrons
				    if(TwoChPFCons==1)
				      {
					fillHistos(jetCand,40,event_weight);
					//Effectiveness of this cut in Category 1 Events
					if(PF12PtFrac_ID_1>0.7)
					  {
					    fillHistos(jetCand,41,event_weight);}
				      } 
				    //Category 2: Two charged Hadrons + Photon
				    if(TwoChPFConsPlusPho==1)
				      {
					fillHistos(jetCand,42,event_weight);
					//Effectiveness of this cut in Category 2 Events
					if(PF123PtFrac_ID_2>0.7)
					  {
					    fillHistos(jetCand,43,event_weight);}
				      }
				    //Category of events with < 2 charged Hadrons
				    if(TwoChPFCons==0 && TwoChPFConsPlusPho==0)
				      {
					fillHistos(jetCand,44,event_weight);
					//Calculating the effectiveness of this cut in only events with < 2 oppositely charged Hadrons
					if(jetetaWidth->at(jetCand[0].first)<0.04)
					  {
					    fillHistos(jetCand,46,event_weight);
					  }}
				    //This is for comparison with previous results (for all events)
				    if (jetetaWidth->at(jetCand[0].first)<0.04)
				      {
					fillHistos(jetCand,47,event_weight);
				      }

				    if (Pt123Fraction_to_use[0] > 0.8)
				      {//Correlated Uncertainty
					fillHistos(jetCand,48,event_weight);
				      }
				    if (Pt123Fraction_to_use[1] > 0.8)
				      {//Tracker Uncertainty
					fillHistos(jetCand,49,event_weight);
				      }
				    if (Pt123Fraction_to_use[2] > 0.8)
				      {//Ecal Uncertainty
					fillHistos(jetCand,50,event_weight);
				      }
				    if (Pt123Fraction_to_use[3] > 0.8)
				      {//Hcal Uncertainty
					fillHistos(jetCand,51,event_weight);
				      }

				    if (Pt123Fraction_to_use[0] > 0.85)
				      {//Correlated Uncertainty
					fillHistos(jetCand,52,event_weight);
				      }
				    if (Pt123Fraction_to_use[1] > 0.85)
				      {//Tracker Uncertainty
					fillHistos(jetCand,53,event_weight);
				      }
				    if (Pt123Fraction_to_use[2] > 0.85)
				      {//Ecal Uncertainty
					fillHistos(jetCand,54,event_weight);
				      }
				    if (Pt123Fraction_to_use[3] > 0.85)
				      {//Hcal Uncertainty
					fillHistos(jetCand,55,event_weight);
				      }
				    if (Pt123Fraction_to_use[0] > 0.9)
				      {//Correlated Uncertainty
					fillHistos(jetCand,56,event_weight);
				      }
				    if (Pt123Fraction_to_use[1] > 0.9)
				      {//Tracker Uncertainty
					fillHistos(jetCand,57,event_weight);
				      }
				    if (Pt123Fraction_to_use[2] > 0.9)
				      {//Ecal Uncertainty
					fillHistos(jetCand,58,event_weight);
				      }
				    if (Pt123Fraction_to_use[3] > 0.9)
				      {//Hcal Uncertainty
					fillHistos(jetCand,59,event_weight);
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

void ZprimeJetsClass_MC::BookHistos(const char* file2)
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

  double PtBins[49]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};

  h_dphimin = new TH1F("h_dphimin","h_dphimin; Minimum dPhiJetMET",50,0,3.2);h_dphimin->Sumw2();
  h_metcut  = new TH1F("h_metcut","h_metcut; |pfMET-caloMET|/pfMET", 50,0,1.2);h_metcut->Sumw2();
  for(int i=0; i<60; i++){

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
    h_PF123PtFraction[i]= new TH1F(("PF123PtFraction"+histname).c_str(), "PF123PtFraction;P_{T} fraction carried by 3 leading daughters of the Pencil Jet" ,50,0,1.1);h_PF123PtFraction[i]->Sumw2();   
    h_j1PF12PtFrac_ID_1[i]= new TH1F(("j1PF12PtFrac_ID_1"+histname).c_str(), "j1PF12PtFrac_ID_1;P_{T} fraction carried by charged hadrons of the leading Jet: Category 1" ,50,0,1.1);h_j1PF12PtFrac_ID_1[i]->Sumw2();   
    h_j1dRPF12_ID_1[i] = new TH1F(("j1dRPF12_ID_1"+histname).c_str(),"j1dRPF12_ID_1; deltaR between oppositely charged hadrons of the leading Jet: Category 1",50,0,0.15);h_j1dRPF12_ID_1[i]->Sumw2();
    h_j1PF12PtFrac_ID_2[i]= new TH1F(("j1PF12PtFrac_ID_2"+histname).c_str(), "j1PF12PtFrac_ID_2;P_{T} fraction carried by charged hadrons of the leading Jet: Category 2" ,50,0,1.1);h_j1PF12PtFrac_ID_2[i]->Sumw2();
    h_j1dRPF12_ID_2[i] = new TH1F(("j1dRPF12_ID_2"+histname).c_str(),"j1dRPF12_ID_2; deltaR between oppositely charged hadrons of the leading Jet: Category 2",50,0,0.15);h_j1dRPF12_ID_2[i]->Sumw2();
    h_j1PFPtFrac_ID_2[i] = new TH1F(("j1PFPtFrac_ID_2"+histname).c_str(),"j1PFPtFrac_ID_2;P_{T} fraction carried by charged hadrons+Photon of the leading Jet: Category 2" ,50,0,1.1);h_j1PFPtFrac_ID_2[i]->Sumw2();  
    h_j1TotPFCands[i] = new TH1F(("j1TotPFCands"+histname).c_str(),"j1TotPFCands;# of all PF candidates in Leading Jet",25,0,50);h_j1TotPFCands[i]->Sumw2();
    h_j1ChPFCands[i] = new TH1F(("j1ChPFCands"+histname).c_str(),"j1ChPFCands;# of PF charged hadrons in Leading Jet",25,0,50);h_j1ChPFCands[i]->Sumw2();
    h_j1NeutPFCands[i] = new TH1F(("j1NeutPFCands"+histname).c_str(),"j1NeutPFCands;# of PF neutral hadrons in Leading Jet",15,0,15);h_j1NeutPFCands[i]->Sumw2();
    h_j1GammaPFCands[i] = new TH1F(("j1GammaPFCands"+histname).c_str(),"j1GammaPFCands;# of PF gammas in Leading Jet",20,0,40);h_j1GammaPFCands[i]->Sumw2();
    h_j1CHF[i] = new TH1F(("j1CHF"+histname).c_str(),"j1CHF;Charged Hadron Energy Fraction in Leading Jet",50,0,1.1);h_j1CHF[i]->Sumw2(); 
    h_j1NHF[i] = new TH1F(("j1NHF"+histname).c_str(),"j1NHF;Neutral Hadron Energy Fraction in Leading Jet",50,0,1.1);h_j1NHF[i]->Sumw2(); 
    h_j1ChMultiplicity[i] = new TH1F(("j1ChMultiplicity"+histname).c_str(),"j1ChMultiplicity;Charged Multiplicity of Leading Jet",25,0,50);h_j1ChMultiplicity[i]->Sumw2();
    h_j1NeutMultiplicity[i] = new TH1F(("j1NeutMultiplicity"+histname).c_str(),"j1NeutMultiplicity;Neutral Multiplicity of Leading Jet",25,0,50);h_j1NeutMultiplicity[i]->Sumw2(); 
    h_j1Mt[i]  = new TH1F(("j1Mt"+histname).c_str(), "j1Mt;M_{T} of Leading Jet (GeV)", 50,MtBins);h_j1Mt[i]->Sumw2(); 
    h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(),"nVtx;nVtx",70,0,70);h_nVtx[i]->Sumw2();
    //Category3 variables
    h_ChPionPt[i]=new TH1F(("ChPionPt"+histname).c_str(),"ChPionPt;p_{T} of Charged Pion in 3rd Signal Category",50,0,2000);h_ChPionPt[i]->Sumw2();
    h_PhotonPt[i]=new TH1F(("PhotonPt"+histname).c_str(),"PhotonPt;p_{T} of Photon in 3rd Signal Category",50,0,2000);h_PhotonPt[i]->Sumw2();
    h_dRPionPhoton[i]=new TH1F(("dRPionPhoton"+histname).c_str(),"dRPionPhoton;deltaR between ChPion and Photon 3rd Signal Category",50,0,0.5);h_dRPionPhoton[i]->Sumw2();
    h_EcalPtUnc[i]=new TH2F(("EcalPtUnc"+histname).c_str(),"ECAL P_{T} Uncertainty;Photon P_{T} (GeV);Uncertainty",50,0.,2500.,50,0.,1.);
    h_TrackerPtUnc[i]=new TH2F(("TrackerPtUnc"+histname).c_str(),"Tracker P_{T} Uncertainty;Charged Hadrons P_{T} (GeV);Uncertainty",50,0.,2500.,50,0.,1.);
    h_HcalPtUnc[i]=new TH2F(("HcalPtUnc"+histname).c_str(),"HCAL P_{T} Uncertainty;Neutral Hadron P_{T} (GeV);Uncertainty",50,0.,2500.,50,0.,1.);
    h_TrackerPtFrac[i]=new TH1F(("TrackerPtFraction"+histname).c_str(), "TrackerPtFraction;P_{T} fraction carried by Charged Hadrons of the Pencil Jet" ,50,0,1.1);h_TrackerPtFrac[i]->Sumw2();
    h_EcalPtFrac[i]=new TH1F(("EcalPtFraction"+histname).c_str(), "EcalPtFraction;P_{T} fraction carried by Photons of the Pencil Jet" ,50,0,1.1);h_EcalPtFrac[i]->Sumw2();
    h_HcalPtFrac[i]=new TH1F(("HcalPtFraction"+histname).c_str(), "HcalPtFraction;P_{T} fraction carried by Neutral Hadrons of the Pencil Jet" ,50,0,1.1);h_HcalPtFrac[i]->Sumw2();
  }
}

//double ZprimeJetsClass_MC::dR(double jetetaWidth, double jetphiWidth)
//{
//  double deltar = sqrt(jetetaWidth*jetetaWidth + jetphiWidth*jetphiWidth);
//  return deltar;
//}


void ZprimeJetsClass_MC::fillHistos(std::vector<std::pair<int,double>> jetCand_to_use,int histoNumber,double event_weight)
{
  h_nVtx[histoNumber]->Fill(nVtx,event_weight);
  h_nJets[histoNumber]->Fill(nJet,event_weight);
  h_pfMETall[histoNumber]->Fill(MET_to_use,event_weight);
  h_pfMET200[histoNumber]->Fill(MET_to_use,event_weight);
  h_pfMET[histoNumber]->Fill(MET_to_use,event_weight);
  h_pfMETPhi[histoNumber]->Fill(METPhi_to_use,event_weight);
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
    h_TrackerPtFrac[histoNumber]->Fill(Pt123Fraction_to_use[1],event_weight);
    h_EcalPtFrac[histoNumber]->Fill(Pt123Fraction_to_use[2],event_weight);
    h_HcalPtFrac[histoNumber]->Fill(Pt123Fraction_to_use[3],event_weight);
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
  }
}

void ZprimeJetsClass_MC::getPt123Frac(int UncType)
{
  vector<double> Pt123 = {0.0,0.0,0.0,0.0};
  vector<double> jetPtAll = {0.0,0.0,0.0,0.0};
  vector<int> index(j1PFConsPID.size());
  iota(begin(index),end(index),0);
  vector<vector<int>> ConsIndex = {index,TrackerCand,EcalCand,HcalCand};
  vector<string> ConsName = {"All","Tracker","Ecal","Hcal"};
  for (int i = 0; i < ConsIndex.size(); i++)
    {
      //cout<<"Running on "<<ConsName[i]<<": ";
      for (int j = 0; j < j1PFConsPID.size(); j++)
	{
	  if (find(ConsIndex[i].begin(),ConsIndex[i].end(),j) != ConsIndex[i].end())
	    {
	      //cout<<j<<":"<<j1PFConsPtUnc.at(j)<<" ";
	      jetPtAll[i]+=j1PFConsPt.at(j)*(1+UncType*j1PFConsPtUnc.at(j));
	      if (j < 3) Pt123[i]+=j1PFConsPt.at(j)*(1+UncType*j1PFConsPtUnc.at(j));
	    }
	  else
	    {
	      //cout<<j<<":"<<0<<" ";
	      jetPtAll[i]+=j1PFConsPt.at(j);
	      if (j < 3) Pt123[i]+=j1PFConsPt.at(j);
	    }
	}
      //cout<<endl;
      Pt123Fraction_to_use[i]=(Pt123[i]/jetPtAll[i]);
    }
}

void ZprimeJetsClass_MC::AllPFCand(std::vector<std::pair<int,double>> jetCand, std::vector<int> PFCandidates)
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
    
  TwoChPFCons=TwoChPFConsPlusPho=0; 
  PF12PtFrac_ID_1=PF12PtFrac_ID_2=dR_PF12_ID_1=dR_PF12_ID_2=PF123PtFrac_ID_2=0.0;
  NoPosPFCons=NoNegPFCons=NoPhoPFCons=0;
  j1PFPosConsPt= j1PFPosConsEta=j1PFPosConsPhi=j1PFNegConsPt=j1PFNegConsEta=j1PFNegConsPhi=j1PFPhoConsPt=j1PFPhoConsEta=j1PFPhoConsPhi=0.0;
  //Category 3 variables
  dR_PionPhoton_3=Cat3_ChPionPt=Cat3_PhotonPt=Cat3_ChPionEta=Cat3_PhotonEta=Cat3_ChPionPhi=Cat3_PhotonPhi=0.0;
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
	    //cout<<"T"<<TrackerCand[TrackerCand.size()-1]<<":"<<j1PFConsPtUnc[TrackerCand[TrackerCand.size()-1]]<<" ";
	  }
	else if (abs(j1PFConsPID.at(i)) == 22 || abs(j1PFConsPID.at(i)) == 11)
	  {
	    //ECAL Uncertainty
	    //deltaPt=(1/100)*sqrt((2.8)^2/Pt+(12.8/Pt)^2+(0.3)^2)
	    j1PFConsPtUnc.push_back((1/100.)*sqrt(pow(2.8,2)/j1PFConsPt.at(i)+pow(12.8/j1PFConsPt.at(i),2)+pow(0.3,2)));
	    EcalCand.push_back(i);
	    //cout<<"E"<<EcalCand[EcalCand.size()-1]<<":"<<j1PFConsPtUnc[EcalCand[EcalCand.size()-1]]<<" ";
	  }
	else if (abs(j1PFConsPID.at(i)) == 130)
	  {
	    //HCAL Uncertainty
	    //deltaPt=(1/100)*sqrt((115)^2/Pt+(5.5)^2)
	    j1PFConsPtUnc.push_back((1/100.)*sqrt(pow(115,2)/j1PFConsPt.at(i)+pow(5.5,2)));
	    HcalCand.push_back(i);
	    //cout<<"H"<<HcalCand[HcalCand.size()-1]<<":"<<j1PFConsPtUnc[HcalCand[HcalCand.size()-1]]<<" ";
	  }
	else
	  {
	    j1PFConsPtUnc.push_back(0);
	    //cout<<"N"<<i<<":"<<j1PFConsPtUnc[i]<<" ";
	  }
      }
    //cout<<endl;
    //Positively charged hadron Cons of the Pencil Jet
    if(j1PFConsPID.size()>0 && j1PFConsPID.at(0)==+211)
      {
	j1PFPosConsPt = j1PFConsPt.at(0);
	j1PFPosConsEta = j1PFConsEta.at(0);
	j1PFPosConsPhi = j1PFConsPhi.at(0);    
      }
    else if(j1PFConsPID.size()>1 && j1PFConsPID.at(1)==+211)
      {
	j1PFPosConsPt = j1PFConsPt.at(1);
	j1PFPosConsEta = j1PFConsEta.at(1);
	j1PFPosConsPhi = j1PFConsPhi.at(1);    
      }
    else if(j1PFConsPID.size()>2 && j1PFConsPID.at(2)==+211)
      {
	j1PFPosConsPt = j1PFConsPt.at(2);
	j1PFPosConsEta = j1PFConsEta.at(2);
	j1PFPosConsPhi = j1PFConsPhi.at(2);    
      }
    else{NoPosPFCons=1;}
    //Negatively charged hadron Cons of the Pencil Jet
    if(j1PFConsPID.size()>0 && j1PFConsPID.at(0)==-211)
      {
	j1PFNegConsPt = j1PFConsPt.at(0);
	j1PFNegConsEta = j1PFConsEta.at(0);
	j1PFNegConsPhi = j1PFConsPhi.at(0);    
      }
    else if(j1PFConsPID.size()>1 && j1PFConsPID.at(1)==-211)
      {
	j1PFNegConsPt = j1PFConsPt.at(1);
	j1PFNegConsEta = j1PFConsEta.at(1);
	j1PFNegConsPhi = j1PFConsPhi.at(1);    
      }
    else if(j1PFConsPID.size()>2 && j1PFConsPID.at(2)==-211)
      {
	j1PFNegConsPt = j1PFConsPt.at(2);
	j1PFNegConsEta = j1PFConsEta.at(2);
	j1PFNegConsPhi = j1PFConsPhi.at(2);    
      }
    else{
      //std::cout<<"Where is the error:"<<std::endl;
      NoNegPFCons=1;}
    //Photon PFCons of the Pencil Jet
    if(j1PFConsPID.size()>0 && j1PFConsPID.at(0)==22)
      {
	j1PFPhoConsPt = j1PFConsPt.at(0);
	j1PFPhoConsEta = j1PFConsEta.at(0);
	j1PFPhoConsPhi = j1PFConsPhi.at(0);    
      }
    else if(j1PFConsPID.size()>1 && j1PFConsPID.at(1)==22)
      {
	j1PFPhoConsPt = j1PFConsPt.at(1);
	j1PFPhoConsEta = j1PFConsEta.at(1);
	j1PFPhoConsPhi = j1PFConsPhi.at(1);    
      }
    else if(j1PFConsPID.size()>2 && j1PFConsPID.at(2)==22)
      {
	j1PFPhoConsPt = j1PFConsPt.at(2);
	j1PFPhoConsEta = j1PFConsEta.at(2);
	j1PFPhoConsPhi = j1PFConsPhi.at(2);    
      }
    else{NoPhoPFCons=1;}
      
    //Category I: Exactly Two Charged Hadrons/Tracks
    if(NoPosPFCons==0 && NoNegPFCons==0 && NoPhoPFCons==1){
      TwoChPFCons=1;
      PF12PtFrac_ID_1 =(j1PFPosConsPt+j1PFNegConsPt)/(jetCand[0].second);
      dR_PF12_ID_1 = deltaR(j1PFPosConsEta,j1PFPosConsPhi,j1PFNegConsEta,j1PFNegConsPhi);
    }
    //Category II: Exactly Two Charged Hadrons/Tracks + One Photon
    if(NoPosPFCons==0 && NoNegPFCons==0 && NoPhoPFCons==0){
      TwoChPFConsPlusPho=1;
      PF12PtFrac_ID_2 =(j1PFPosConsPt+j1PFNegConsPt)/(jetCand[0].second);
      dR_PF12_ID_2 = deltaR(j1PFPosConsEta,j1PFPosConsPhi,j1PFNegConsEta,j1PFNegConsPhi);
      PF123PtFrac_ID_2 = (j1PFPosConsPt+j1PFNegConsPt+j1PFPhoConsPt)/(jetCand[0].second);
    }
    //Category3
    if(TwoChPFCons==0 && TwoChPFConsPlusPho==0){
      if(j1PFConsPID.size()>0){
	if(abs(j1PFConsPID.at(0))==211){
	  Cat3_ChPionPt=j1PFConsPt.at(0); 
	  Cat3_ChPionEta=j1PFConsEta.at(0);
	  Cat3_ChPionPhi=j1PFConsPhi.at(0);}
	else if(abs(j1PFConsPID.at(0))==22){
	  Cat3_PhotonPt=j1PFConsPt.at(0); 
	  Cat3_PhotonEta=j1PFConsEta.at(0);
	  Cat3_PhotonPhi=j1PFConsPhi.at(0);}
      }
      if(j1PFConsPID.size()>1){
	if(abs(j1PFConsPID.at(1))==211){
	  //Confirm that it does not get overwritten with smaller value
	  if(j1PFConsPt.at(1)>Cat3_ChPionPt){
	    Cat3_ChPionPt=j1PFConsPt.at(1); 
	    Cat3_ChPionEta=j1PFConsEta.at(1);
	    Cat3_ChPionPhi=j1PFConsPhi.at(1);}}
	else if(abs(j1PFConsPID.at(1))==22){
	  if(j1PFConsPt.at(1)>Cat3_PhotonPt){
	    Cat3_PhotonPt=j1PFConsPt.at(1); 
	    Cat3_PhotonEta=j1PFConsEta.at(1);
	    Cat3_PhotonPhi=j1PFConsPhi.at(1);}}
      }
      if(j1PFConsPID.size()>2){
	if(abs(j1PFConsPID.at(2))==211){
	  //Confirm that it does not get overwritten with smaller value
	  if(j1PFConsPt.at(2)>Cat3_ChPionPt){
	    Cat3_ChPionPt=j1PFConsPt.at(2); 
	    Cat3_ChPionEta=j1PFConsEta.at(2);
	    Cat3_ChPionPhi=j1PFConsPhi.at(2);}}
	else if(abs(j1PFConsPID.at(2))==22){
	  if(j1PFConsPt.at(2)>Cat3_PhotonPt){
	    Cat3_PhotonPt=j1PFConsPt.at(2); 
	    Cat3_PhotonEta=j1PFConsEta.at(2);
	    Cat3_PhotonPhi=j1PFConsPhi.at(2);}}
      }
      if(Cat3_ChPionPt>0 && Cat3_PhotonPt>0){
	dR_PionPhoton_3 = deltaR(Cat3_ChPionEta,Cat3_ChPionPhi,Cat3_PhotonEta,Cat3_PhotonPhi);
      }
    }
  }
}

//Function to calculate regular deltaR separate from jet width variable 'dR'
double ZprimeJetsClass_MC::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//Gives the (minimum) separation in phi between the specified phi values
////Must return a positive value
float ZprimeJetsClass_MC::DeltaPhi(float phi1, float phi2)
{
  float pi = TMath::Pi();
  float dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}

float ZprimeJetsClass_MC::dPhiJetMETmin(std::vector<int> jets)
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
std::vector<std::pair<int,double>> ZprimeJetsClass_MC::getJetCand(std::vector<int> jets,double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut,int UncType){

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

std::vector<int> ZprimeJetsClass_MC::JetVetoDecision(int UncType) {

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
std::vector<int>ZprimeJetsClass_MC::getPFCandidates(){
  std::vector<int>PFCands;
  for(int i=0;i<nJet;i++)
    {
      int TotPFCands;
      if(i==0){
	TotPFCands = j1PFConsPID.size();
	//std::cout<<"Total PFCands: "<<TotPFCands<<std::endl;
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
bool ZprimeJetsClass_MC::btagVeto(int UncType) {

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

bool ZprimeJetsClass_MC::dPhiJetMETcut(std::vector<int> jets, float METPhi)
{
  //reject jet if it is found within DeltaPhi(jet,MET) < 0.5                                                                                              \
  
  bool passes = false;

  int njetsMax = jets.size();
  //std::cout<<"METPhi_to_use(JetMETCut): "<<METPhi<<std::endl;
  //Only look at first four jets (because that's what monojet analysis do)
  if(njetsMax > 4)
    njetsMax = 4;
  int j=0;
  for(;j< njetsMax; j++){
    //std::cout<<"DeltaPhi b/w Jet and MET"<<std::endl;
    //std::cout<<"jet "<<j<<":"<<DeltaPhi((*jetPhi)[j],METPhi_to_use)<<std::endl;
    if(DeltaPhi((*jetPhi)[j],METPhi) < 0.5)
      break;
  }

  if(j==njetsMax)
    passes = true;

  return passes;

}
bool ZprimeJetsClass_MC::electron_veto_looseID(int jet_index, float elePtCut)
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
bool ZprimeJetsClass_MC::muon_veto_looseID(int jet_index, float muPtCut)
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
