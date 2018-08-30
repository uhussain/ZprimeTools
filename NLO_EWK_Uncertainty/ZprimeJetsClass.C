//For use with Ntuples made from JetAnalyzer
////Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
////
////To compile using rootcom to an executable named 'analyze':
////$ ./rootcom ZprimeJetsClass analyze
////
////To run, assuming this is compiled to an executable named 'analyze':
////$ ./analyze /hdfs/store/user/uhussain/Zprime_Ntuples/ /cms/uhussain/MonoZprimeJet/CMSSW_8_0_8/src/LightZPrimeAnalysis/JetAnalyzer/test/output.root -1 10000
////Runs over every event in the folder Zprime_Ntuples, reporting progress every 10000 events
////and storing the resulting histograms in the file output.root.
////
//
#define ZprimeJetsClass_cxx
#include "ZprimeJetsClass.h"
#include <TH2.h>
#include<TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TH2D.h"
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
  //const char* file2 = argv[2];
  ZprimeJetsClass t(argv[1],argv[2],atoi(argv[6]),atoi(argv[7]));
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void ZprimeJetsClass::Loop(Long64_t maxEvents, int reportEvery)
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

  //double ntest, nbtagVeto;
  //ntest=nbtagVeto=0;
  std::vector<int> jetveto;

  double nTotalEvents,nFilters, nHLT, nMET200, nMETcut,nLeptonIDs,nbtagVeto, nDphiJetMET,nJetSelection;
  nTotalEvents = nFilters = nHLT = nMET200 = nMETcut = nLeptonIDs = nDphiJetMET = nbtagVeto = nJetSelection = 0;

  //getPFCandidates
  std::vector<int>PFCandidates;

  float dphimin=-99;
  //Event is rejected if it contains a HighPtMuon faking MET

  
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;
  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry+=4) {
    
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
    
    fillHistos(jetCand,0);
    float metcut= 0.0;
    MET_to_use = pfMET;
    METPhi_to_use = pfMETPhi;
    jetveto = JetVetoDecision(0);
    jetCand = getJetCand(jetveto,200,2.4,0.8,0.1,0);
    AllPFCand(jetCand,PFCandidates);
    nTotalEvents++;
    if (metFilters==1536)
      {    
        nFilters++;
        fillHistos(jetCand,1);
        if ((HLTJet>>4&1 == 1) || (HLTJet>>5&1 == 1) || (HLTJet>>6&1 == 1) || (HLTJet>>8&1 == 1) )//Mono-jet triggers
	  {
	    nHLT++;
	    fillHistos(jetCand,2);
	    if(jetCand.size()>0)
	      {
		nJetSelection++;
		fillHistos(jetCand,3);
		if (pfMET>250)
		  {
		    nMET200++;
		    fillHistos(jetCand,4);
		    metcut = (fabs(MET_to_use-caloMET))/MET_to_use;
		    h_metcut->Fill(metcut);
		    if(metcut<0.5)
		      {
			nMETcut++;
			fillHistos(jetCand,5);
			if(electron_veto_looseID(jetCand[0].first,10) &&  muon_veto_looseID(jetCand[0].first,10))
			  {
			    nLeptonIDs++;
			    fillHistos(jetCand,6);
			    if(btagVeto(0))
			      {
				nbtagVeto++;
				fillHistos(jetCand,7);
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
				h_dphimin->Fill(minDPhiJetMET_first4);
				if(dPhiJetMETcut(jetveto,METPhi_to_use))
				  {
				    nDphiJetMET++;
				    for (int i = 0; i < 7; i++) fillHistos(jetCand,8+2*i);
				    for (int i = 0; i < 7; i++) if (Pt123Fraction_to_use[0]>0.5) fillHistos(jetCand,9+2*i);
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
  
    // std::cout<<"ntest: "<<ntest<<std::endl;
    // std::cout<<"nbtagVeto: "<<nbtagVeto<<std::endl;
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

void ZprimeJetsClass::BookHistos(const char* file2)
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
  
  float MetBins[45]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1400.,1800.,2000.,2500.};

  double PtBins[49]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};
  
  h_metcut  = new TH1F("h_metcut","h_metcut; |pfMET-caloMET|/pfMET", 50,0,1.2);h_metcut->Sumw2();
  h_dphimin = new TH1F("h_dphimin","h_dphimin; Minimum dPhiJetMET",50,0,3.2);h_dphimin->Sumw2();
  for(int i=0; i<22; i++){

    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    std::string histname(ptbins);
    h_metFilters[i] = new TH1F(("metFilters"+histname).c_str(),"metFilters",50,0,3000); h_metFilters[i]->Sumw2();
    h_nJets[i]   = new TH1F(("nJets"+histname).c_str(), "nJets;Number of Jets", 50, 0, 100);h_nJets[i]->Sumw2();
    h_pfMETall[i] =  new TH1F(("pfMETall"+histname).c_str(), "pfMET",50,0,2000);h_pfMETall[i] ->Sumw2(); 
    h_pfMET200[i] = new TH1F(("pfMET200"+histname).c_str(), "pfMET",50,170,1500);h_pfMET200[i] ->Sumw2(); 
    h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "E_{T}^{miss} (GeV)",44,MetBins);h_pfMET[i] ->Sumw2();
    h_pfMETPhi[i] = new TH1F(("pfMETPhi"+histname).c_str(), "pfMETPhi",50,-4,4);h_pfMETPhi[i]->Sumw2();
    h_j1Pt[i]  = new TH1F(("j1pT"+histname).c_str(), "j1pT;p_{T} of Leading Jet (GeV)", 48,PtBins);h_j1Pt[i]->Sumw2();
    h_j1Eta[i] = new TH1F(("j1Eta"+histname).c_str(), "j1Eta; #eta of Leading Jet", 50, -2.5, 2.5);h_j1Eta[i]->Sumw2();
    h_j1Phi[i] = new TH1F(("j1Phi"+histname).c_str(), "j1Phi; #phi of Leading Jet", 50, -3.0, 3.0);h_j1Phi[i]->Sumw2();
    h_j1etaWidth[i] = new TH1F(("j1etaWidth"+histname).c_str(),"j1etaWidh; #eta width of Leading Jet", 50,0,0.25);h_j1etaWidth[i] ->Sumw2();
    h_j1phiWidth[i] = new TH1F(("j1phiWidth"+histname).c_str(),"j1phiWidth; #phi width of Leading Jet", 50, 0,0.5);h_j1phiWidth[i]->Sumw2();
    h_j1nCons[i] = new TH1F (("j1nCons"+histname).c_str(),"j1nCons; Number of Constituents of Leading Jet",25, 0, 50);h_j1nCons[i]->Sumw2();
    h_j1nCategory1[i] = new TH1F(("j1nCategory1"+histname).c_str(),"j1nCategory1: Number of events with exactly two charged Hadrons",2,-0.5,1.5);h_j1nCategory1[i]->Sumw2();
    h_j1nCategory2[i] = new TH1F(("j1nCategory2"+histname).c_str(),"j1nCategory2: Number of events with two charged Hadrons and one Photon",2,-0.5,1.5);h_j1nCategory2[i]->Sumw2(); 
    h_PF123PtFraction[i]= new TH1F(("PF123PtFraction"+histname).c_str(), "PF123PtFraction;P_{T} fraction carried by 3 leading daughters of the Pencil Jet" ,50,0,1);h_PF123PtFraction[i]->Sumw2();   
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
    h_dnKUnc[i]=new TH2F(("NLO_EWK"+histname).c_str(),"NLO EWK Uncertainty;Pt123Fraction;Uncertainty",50,0,1.,50,-1.,1.);
  }
}

//double ZprimeJetsClass::dR(double jetetaWidth, double jetphiWidth)
//{
//  double deltar = sqrt(jetetaWidth*jetetaWidth + jetphiWidth*jetphiWidth);
//  return deltar;
//}


void ZprimeJetsClass::fillHistos(std::vector<std::pair<int,double>> jetCand_to_use,int histoNumber)
{
  h_nVtx[histoNumber]->Fill(nVtx);
  h_nJets[histoNumber]->Fill(nJet);
  h_pfMETall[histoNumber]->Fill(MET_to_use);
  h_pfMET200[histoNumber]->Fill(MET_to_use);
  h_pfMET[histoNumber]->Fill(MET_to_use);
  h_pfMETPhi[histoNumber]->Fill(METPhi_to_use);
  if(jetCand_to_use.size()>0){
    h_j1Pt[histoNumber]->Fill(jetCand_to_use[0].second);
    h_j1Eta[histoNumber]->Fill(jetEta->at(jetCand_to_use[0].first));
    h_j1Phi[histoNumber]->Fill(jetPhi->at(jetCand_to_use[0].first));
    h_j1nCategory1[histoNumber]->Fill(TwoChPFCons); 
    h_j1nCategory2[histoNumber]->Fill(TwoChPFConsPlusPho);
    h_PF123PtFraction[histoNumber]->Fill(Pt123Fraction_to_use[0]);
    h_j1PF12PtFrac_ID_1[histoNumber]->Fill(PF12PtFrac_ID_1);
    h_j1dRPF12_ID_1[histoNumber]->Fill(dR_PF12_ID_1);
    h_j1PF12PtFrac_ID_2[histoNumber]->Fill(PF12PtFrac_ID_2);
    h_j1dRPF12_ID_2[histoNumber]->Fill(dR_PF12_ID_2);
    h_j1PFPtFrac_ID_2[histoNumber]->Fill(PF123PtFrac_ID_2);
    h_j1TotPFCands[histoNumber]->Fill(TotalPFCandidates);
    h_j1ChPFCands[histoNumber]->Fill(ChargedPFCandidates);
    h_j1NeutPFCands[histoNumber]->Fill(NeutralPFCandidates);
    h_j1GammaPFCands[histoNumber]->Fill(GammaPFCandidates);
    h_j1CHF[histoNumber]->Fill(jetCHF->at(jetCand_to_use[0].first));
    h_j1NHF[histoNumber]->Fill(jetNHF->at(jetCand_to_use[0].first));
    h_j1ChMultiplicity[histoNumber]->Fill(jetNCH->at(jetCand_to_use[0].first));
    h_j1NeutMultiplicity[histoNumber]->Fill(jetNNP->at(jetCand_to_use[0].first));
    h_j1Mt[histoNumber]->Fill(jetMt->at(jetCand_to_use[0].first));
    //Category3
    h_ChPionPt[histoNumber]->Fill(Cat3_ChPionPt);
    h_PhotonPt[histoNumber]->Fill(Cat3_PhotonPt);
    if(Cat3_ChPionPt>0 && Cat3_PhotonPt>0){
      h_dRPionPhoton[histoNumber]->Fill(dR_PionPhoton_3);}
    h_j1etaWidth[histoNumber]->Fill(jetetaWidth->at(jetCand_to_use[0].first));
    h_j1phiWidth[histoNumber]->Fill(jetphiWidth->at(jetCand_to_use[0].first));
    h_j1nCons[histoNumber]->Fill(jetnPhotons->at(jetCand_to_use[0].first)+jetnCHPions->at(jetCand_to_use[0].first)+jetnMisc->at(jetCand_to_use[0].first));
    h_dnKUnc[histoNumber]->Fill(Pt123Fraction_to_use[0],0);
  }
}

void ZprimeJetsClass::getPt123Frac(std::vector<std::pair<int,double>> jetCand,int UncType)
{
  vector<double> Pt123 = {0.0,0.0,0.0,0.0};
  vector<double> jetPtAll = {0.0,0.0,0.0,0.0};
  vector<int> index(j1PFConsPID.size());
  iota(begin(index),end(index),0);
  vector<vector<int>> ConsIndex = {index,TrackerCand,EcalCand,HcalCand};
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
      jetPtAll[0]=jetPt->at(jetCand[0].first);
      Pt123Fraction_to_use[i]=(Pt123[i]/jetPtAll[i]);
    }
}

void ZprimeJetsClass::AllPFCand(std::vector<std::pair<int,double>> jetCand, std::vector<int> PFCandidates)
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
	if (abs(j1PFConsPID.at(i)) == 221 || abs(j1PFConsPID.at(i)) == 13)
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
	    j1PFConsPtUnc.push_back((1/100.)*sqrt(pow(2.8,2)/j1PFConsPt.at(i)+pow(12.8/j1PFConsPt.at(i),2)+pow(0.3,2)));
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
    getPt123Frac(jetCand,0);
  }
}
    
//Function to calculate regular deltaR separate from jet width variable 'dR'
double ZprimeJetsClass::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
float ZprimeJetsClass::DeltaPhi(float phi1, float phi2)
{
  float pi = TMath::Pi();
  float dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}
float ZprimeJetsClass::dPhiJetMETmin(std::vector<int> jets)
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
//getJetCand function should be acting on the output of JetVetoDecision function
//This function should check if leading Jet passes the cuts for the leading Jet
std::vector<std::pair<int,double>> ZprimeJetsClass::getJetCand(std::vector<int> jets,double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut,int UncType){

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
std::vector<int> ZprimeJetsClass::JetVetoDecision(int UncType) {

  std::vector<int> jetindex;
  //let's only look at the jets selected from our getJetCand cut
  //int njets = jets.size();
  for(int i = 0; i < nJet; i++)
    {
      Float_t jetPt_to_use; // = jets[i].second;
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
//get PF Candidates of the selected Jet ->jetCand[0].first
std::vector<int>ZprimeJetsClass::getPFCandidates(){
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
bool ZprimeJetsClass::btagVeto(int UncType) {

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
bool ZprimeJetsClass::dPhiJetMETcut(std::vector<int> jets, float METPhi)
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

bool ZprimeJetsClass::electron_veto_looseID(int jet_index, float elePtCut)
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
      //std::cout<<"electronveto goes into loop"<<std::endl;
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
      //std::cout<<"electronveto works until here?"<<std::endl;      
      //std::cout<<"electronPt: "<<elePt->at(i)<<std::endl;
      //std::cout<<"elePFPhoIso: "<<elePFPhoIso->at(i)<<std::endl;
      //std::cout<<"elePFNeuIso: "<<elePFNeuIso->at(i)<<std::endl;
      //std::cout<<"elePFChIso: "<<elePFChIso->at(i)<<std::endl;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));
      //std::cout<<"EAcorrIso: "<<EAcorrIso<<std::endl;
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
      //Electron passes LoosonvVeto
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

bool ZprimeJetsClass::muon_veto_looseID(int jet_index, float muPtCut)
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
	      //Muon does not overlap jet  
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
