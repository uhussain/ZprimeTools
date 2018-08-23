//For use with Ntuples made from JetAnalyzer
////Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all, events), 4 is reportEvery
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
  //jetCandidate that passes the basic pt,eta, NHF, CHF cuts
  jetCand.clear();

  std::vector<int> jetveto;
  jetveto.clear();

  double nTotalEvents,nFilters, nHLT, nCRSelection, nMET200, pfMET50, nNoMuons, nMETcut,nbtagVeto, nDphiJetMET,nJetSelection;
  nTotalEvents = nFilters = nHLT = nCRSelection = nMET200 = pfMET50 = nNoMuons = nMETcut = nDphiJetMET = nbtagVeto = nJetSelection = 0;
  
  //getPFCandidates
  std::vector<int>PFCandidates;
  
  //jetCandidate that passes dPhiJetMET cut out of the above jetCand
  std::vector<int> jetCand1;
  jetCand1.clear();
   
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
    jetCand = getJetCand(200,2.4,0.8,0.1);

    float metcut= 0.0;
    //metcut = (fabs(pfMET-caloMET))/pfMET;
    AllPFCand(jetCand,PFCandidates);
    //closing the pfMET>300 and goodJets condition
    
    //CR Variables
    lepindex = -1;
    nTotalEvents+=event_weight;
    if (metFilters==0)
      {   
        nFilters+=event_weight;
        fillHistos(0,event_weight);
	//if ((HLTEleMuX>>4&1 == 1) || (HLTEleMuX>>38&1 == 1))//"HLT_Ele27_WPTight_Gsf_v or HLT_Ele115_CaloIdVT_GsfTrkIdT_v"
        if (true)
	  {
	    nHLT+=event_weight;
	    fillHistos(1,event_weight);
	    if(jetCand.size()>0)
	      {
		nJetSelection+=event_weight;
		//CR code
		//At least one of the two electrons passes the tight selection
		std::vector<int> elelist = electron_veto_tightID(jetCand[0],30.0);
		std::vector<int> looseEles = electron_veto_looseID(jetCand[0],0,10.0);
		std::vector<int> mulist;
		mulist.clear();  
    
		if(elelist.size() == 1 && looseEles.size() == 1)
                  { 
		    lepindex = elelist[0];
                    mulist = muon_veto_looseID(jetCand[0],elelist[0],10.0);
                    jetveto = JetVetoDecision(jetCand[0],elelist[0]);

                    TLorentzVector lep_4vec;
                    lep_4vec.SetPtEtaPhiE(elePt->at(elelist[0]),eleEta->at(elelist[0]),elePhi->at(elelist[0]),eleEn->at(elelist[0]));

                    lepton_pt = lep_4vec.Pt();         
		    TLorentzVector met_4vec;
		    met_4vec.SetPtEtaPhiE(pfMET,0.,pfMETPhi,pfMET);
		    TLorentzVector leptoMET_4vec = lep_4vec+met_4vec;
		    Double_t leptoMET = fabs(leptoMET_4vec.Pt());
		    Double_t leptoMET_phi = leptoMET_4vec.Phi();
		    nCRSelection+=event_weight;
		    Recoil = leptoMET;
		    metcut = (fabs(pfMET-caloMET))/Recoil;
		    fillHistos(2,event_weight);
		    if (leptoMET>250)
		      {
			//leptoMET_phi_to_use = leptoMET_phi;
			nMET200+=event_weight;
			fillHistos(3,event_weight);
			if(mulist.size() == 0)
			  {
			    nNoMuons+=event_weight;
			    fillHistos(4,event_weight);
			    Float_t dPhi_lepMET = DeltaPhi(elePhi->at(lepindex),pfMETPhi);
			    Float_t lepMET_MT = sqrt(2*elePt->at(lepindex)*pfMET*(1-TMath::Cos(dPhi_lepMET)));
			    h_lepMET_MT->Fill(lepMET_MT,event_weight);
			    if(pfMET > 50)
			      {
				pfMET50+=event_weight;
				fillHistos(5,event_weight);
				h_metcut->Fill(metcut,event_weight);
				if(metcut<0.5)
				  {
				    nMETcut+=event_weight;
				    fillHistos(6,event_weight);
				    if(btagVeto())
				      {
					nbtagVeto+=event_weight;
					fillHistos(7,event_weight);
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
					if(dPhiJetMETcut(jetveto))
					  {
					    nDphiJetMET+=event_weight;
					    fillHistos(8,event_weight);
					    //Category 1: Exactly Two Charged Hadrons
					    if(TwoChPFCons==1)
					      {
						fillHistos(9,event_weight);
						//Effectiveness of this cut in Category 1 Events
						if(PF12PtFrac_ID_1>0.7)
						  {
						    fillHistos(10,event_weight);}
					      } 
					    //Category 2: Two charged Hadrons + Photon
					    if(TwoChPFConsPlusPho==1)
					      {
						fillHistos(11,event_weight);
						//Effectiveness of this cut in Category 2 Events
						if(PF123PtFrac_ID_2>0.7)
						  {
						    fillHistos(12,event_weight);}
					      }
					    //Category of events with < 2 charged Hadrons
					    if(TwoChPFCons==0 && TwoChPFConsPlusPho==0)
					      {
						fillHistos(13,event_weight);
						//Calculating the effectiveness of this cut in only events with < 2 oppositely charged Hadrons
						if(jetetaWidth->at(jetCand[0])<0.04)
						  {
						    fillHistos(14,event_weight);
						  }}
					    //This is for comparison with previous results (for all events)
					    if (jetetaWidth->at(jetCand[0])<0.04)
					      {
						fillHistos(15,event_weight);
					      }
					  }
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
  h_cutflow->SetBinContent(5,nCRSelection);
  h_cutflow->SetBinContent(6,nMET200);
  h_cutflow->SetBinContent(7,nNoMuons);
  h_cutflow->SetBinContent(8,pfMET50);
  h_cutflow->SetBinContent(9,nMETcut);
  h_cutflow->SetBinContent(10,nbtagVeto);
  h_cutflow->SetBinContent(11,nDphiJetMET);
  //save the histograms
  //   histFile->Write();
  //   histFile->Close();
  
}//Closing the Loop function

void ZprimeJetsClass_MC::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ZprimeJet","ZprimeJet");
  fileName->cd();

  h_cutflow = new TH1D("h_cutflow","h_cutflow",11,0,11);h_cutflow->Sumw2();
  h_cutflow->GetXaxis()->SetBinLabel(1,"Total Events");
  h_cutflow->GetXaxis()->SetBinLabel(2,"metFilters");
  h_cutflow->GetXaxis()->SetBinLabel(3,"Trigger");
  h_cutflow->GetXaxis()->SetBinLabel(4,"GoodJet");
  h_cutflow->GetXaxis()->SetBinLabel(5,"CRSelection"); 
  h_cutflow->GetXaxis()->SetBinLabel(6,"leptoMetCut");
  h_cutflow->GetXaxis()->SetBinLabel(7,"NoMuons");
  h_cutflow->GetXaxis()->SetBinLabel(8,"pfMET50");
  h_cutflow->GetXaxis()->SetBinLabel(9,"caloMET cut");
  h_cutflow->GetXaxis()->SetBinLabel(10,"B-JetVeto");
  h_cutflow->GetXaxis()->SetBinLabel(11,"DeltaPhiCut");

  float MtBins[51]={180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		    780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.,3000.};
  
  float MetBins[45]={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		     780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1400.,1800.,2000.,2500.};

  float PtBins[51]={160,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,
		    780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.,1000.,1050.,1100.,1200.,1300.,1400.,1500.,2000.,2500.};

  h_dphimin = new TH1F("h_dphimin","h_dphimin; Minimum dPhiJetMET",50,0,3.2);h_dphimin->Sumw2();
  h_metcut  = new TH1F("h_metcut","h_metcut; |pfMET-caloMET|/pfMET", 50,0,1.2);h_metcut->Sumw2();
  h_lepMET_MT = new TH1F("h_lepMET_MT","h_lepMET_MT; transverse mass of the lepton-Emiss system",40,0,400);h_lepMET_MT->Sumw2();
  for(int i=0; i<16; i++){

    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    std::string histname(ptbins);
    h_nJets[i]   = new TH1F(("nJets"+histname).c_str(), "nJets;Number of Jets", 10, 0, 10);h_nJets[i]->Sumw2();
    h_pfMETall[i] =  new TH1F(("pfMETall"+histname).c_str(), "pfMET",50,0,2000);h_pfMETall[i] ->Sumw2(); 
    h_pfMET200[i] = new TH1F(("pfMET200"+histname).c_str(), "pfMET",50,170,1500);h_pfMET200[i] ->Sumw2(); 
    h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "E_{T}^{miss} (GeV)",44,MetBins);h_pfMET[i] ->Sumw2();
    h_pfMETPhi[i] = new TH1F(("pfMETPhi"+histname).c_str(), "pfMETPhi",50,-4,4);h_pfMETPhi[i]->Sumw2();
    h_j1Pt[i]  = new TH1F(("j1pT"+histname).c_str(), "j1pT;p_{T} of Leading Jet (GeV)", 50,PtBins);h_j1Pt[i]->Sumw2();
    h_j1Eta[i] = new TH1F(("j1Eta"+histname).c_str(), "j1Eta; #eta of Leading Jet", 50, -2.5, 2.5);h_j1Eta[i]->Sumw2();
    h_j1Phi[i] = new TH1F(("j1Phi"+histname).c_str(), "j1Phi; #phi of Leading Jet", 50, -3.0, 3.0);h_j1Phi[i]->Sumw2();     
    h_j1etaWidth[i] = new TH1F(("j1etaWidth"+histname).c_str(),"j1etaWidh; #eta width of Leading Jet", 50,0,0.25);h_j1etaWidth[i] ->Sumw2();
    h_j1phiWidth[i] = new TH1F(("j1phiWidth"+histname).c_str(),"j1phiWidth; #phi width of Leading Jet", 50, 0,0.5);h_j1phiWidth[i]->Sumw2();
    h_j1nCons[i] = new TH1F (("j1nCons"+histname).c_str(),"j1NConstituents; Number of Constituents of Leading Jet",25, 0, 50);h_j1nCons[i]->Sumw2();
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
    //CR Histograms
    h_LeptonPt[i] = new TH1F(("h_LeptonPt"+histname).c_str(),"h_LeptonPt",30,0.,1500.);h_LeptonPt[i]->Sumw2();
    h_LeptonEta[i] = new TH1F(("h_LeptonEta"+histname).c_str(),"h_LeptonEta",30,-3.,3.);h_LeptonEta[i]->Sumw2();
    h_LeptonPhi[i] = new TH1F(("h_LeptonPhi"+histname).c_str(),"h_LeptonPhi",30,0.,3.1416);h_LeptonPhi[i]->Sumw2();
    h_recoil[i] = new TH1F(("h_recoil"+histname).c_str(), "Recoil (GeV)",44,MetBins);h_recoil[i] ->Sumw2();
  }
}

//double ZprimeJetsClass_MC::dR(double jetetaWidth, double jetphiWidth)
//{
//  double deltar = sqrt(jetetaWidth*jetetaWidth + jetphiWidth*jetphiWidth);
//  return deltar;
//}


void ZprimeJetsClass_MC::fillHistos(int histoNumber,double event_weight)
{
  h_nVtx[histoNumber]->Fill(nVtx,event_weight);
  h_nJets[histoNumber]->Fill(nJet,event_weight);
  h_pfMETall[histoNumber]->Fill(pfMET,event_weight);
  h_pfMET200[histoNumber]->Fill(pfMET,event_weight);
  h_pfMET[histoNumber]->Fill(pfMET,event_weight);
  h_pfMETPhi[histoNumber]->Fill(pfMETPhi,event_weight);
  if(jetCand.size()>0){
    h_j1Pt[histoNumber]->Fill(jetPt->at(jetCand[0]),event_weight);
    h_j1Eta[histoNumber]->Fill(jetEta->at(jetCand[0]),event_weight);
    h_j1Phi[histoNumber]->Fill(jetPhi->at(jetCand[0]),event_weight); 
    h_j1nCategory1[histoNumber]->Fill(TwoChPFCons,event_weight); 
    h_j1nCategory2[histoNumber]->Fill(TwoChPFConsPlusPho,event_weight);
    h_PF123PtFraction[histoNumber]->Fill(Pt123Fraction,event_weight);
    h_j1PF12PtFrac_ID_1[histoNumber]->Fill(PF12PtFrac_ID_1,event_weight);
    h_j1dRPF12_ID_1[histoNumber]->Fill(dR_PF12_ID_1,event_weight);
    h_j1PF12PtFrac_ID_2[histoNumber]->Fill(PF12PtFrac_ID_2,event_weight);
    h_j1dRPF12_ID_2[histoNumber]->Fill(dR_PF12_ID_2,event_weight);
    h_j1PFPtFrac_ID_2[histoNumber]->Fill(PF123PtFrac_ID_2,event_weight);
    h_j1TotPFCands[histoNumber]->Fill(TotalPFCandidates,event_weight);
    h_j1ChPFCands[histoNumber]->Fill(ChargedPFCandidates,event_weight);
    h_j1NeutPFCands[histoNumber]->Fill(NeutralPFCandidates,event_weight);
    h_j1GammaPFCands[histoNumber]->Fill(GammaPFCandidates,event_weight); 
    h_j1CHF[histoNumber]->Fill(jetCHF->at(jetCand[0]),event_weight);
    h_j1NHF[histoNumber]->Fill(jetNHF->at(jetCand[0]),event_weight);
    h_j1ChMultiplicity[histoNumber]->Fill(jetNCH->at(jetCand[0]),event_weight);
    h_j1NeutMultiplicity[histoNumber]->Fill(jetNNP->at(jetCand[0]),event_weight);
    h_j1Mt[histoNumber]->Fill(jetMt->at(jetCand[0]),event_weight);
    h_j1etaWidth[histoNumber]->Fill(jetetaWidth->at(jetCand[0]),event_weight);
    h_j1phiWidth[histoNumber]->Fill(jetphiWidth->at(jetCand[0]),event_weight);
    h_j1nCons[histoNumber]->Fill((jetnPhotons->at(jetCand[0])+jetnCHPions->at(jetCand[0])+jetnMisc->at(jetCand[0])),event_weight);
  }
  //CR Histograms
  if(lepindex >= 0){ 
    h_LeptonPt[histoNumber]->Fill(elePt->at(lepindex),event_weight);
    h_LeptonEta[histoNumber]->Fill(eleEta->at(lepindex),event_weight);
    h_LeptonPhi[histoNumber]->Fill(elePhi->at(lepindex),event_weight);
  }    
  if(lepton_pt > 0){  
    h_recoil[histoNumber]->Fill(Recoil,event_weight);}
}

void ZprimeJetsClass_MC::getPt123Frac()
{
  double Pt123=0.;
  double jetPtAll=0.;
  for (int j = 0; j < j1PFConsPID.size(); j++)
    {
      jetPtAll+=j1PFConsPt.at(j);
      if (j < 3) Pt123+=j1PFConsPt.at(j);
    }
  Pt123Fraction=(Pt123/jetPtAll);
}

void ZprimeJetsClass_MC::AllPFCand(std::vector<int> jetCand, std::vector<int> PFCandidates)
{
  //getPFCandidatesMethod
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
    
  TwoChPFCons=TwoChPFConsPlusPho=0; 
  PF12PtFrac_ID_1=PF12PtFrac_ID_2=dR_PF12_ID_1=dR_PF12_ID_2=PF123PtFrac_ID_2=0.0;
  NoPosPFCons=NoNegPFCons=NoPhoPFCons=0;
  j1PFPosConsPt= j1PFPosConsEta=j1PFPosConsPhi=j1PFNegConsPt=j1PFNegConsEta=j1PFNegConsPhi=j1PFPhoConsPt=j1PFPhoConsEta=j1PFPhoConsPhi=0.0;
    
  Pt123=Pt123Fraction=0.0;
  //Category 3 variables
  dR_PionPhoton_3=Cat3_ChPionPt=Cat3_PhotonPt=Cat3_ChPionEta=Cat3_PhotonEta=Cat3_ChPionPhi=Cat3_PhotonPhi=0.0;
  //We are using these conditions so we only calculate the following quantities for the signal we are interested in
  //This will also make it faster to process the events
  if(jetCand.size()>0){
    j1PFConsPt=JetsPFConsPt->at(jetCand[0]);
    j1PFConsEta=JetsPFConsEta->at(jetCand[0]);
    j1PFConsPhi=JetsPFConsPhi->at(jetCand[0]);
    j1PFConsPID=JetsPFConsPID->at(jetCand[0]);

    getPt123Frac();
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
      PF12PtFrac_ID_1 =(j1PFPosConsPt+j1PFNegConsPt)/(jetPt->at(jetCand[0]));
      dR_PF12_ID_1 = deltaR(j1PFPosConsEta,j1PFPosConsPhi,j1PFNegConsEta,j1PFNegConsPhi);
    }
    //Category II: Exactly Two Charged Hadrons/Tracks + One Photon
    if(NoPosPFCons==0 && NoNegPFCons==0 && NoPhoPFCons==0){
      TwoChPFConsPlusPho=1;
      PF12PtFrac_ID_2 =(j1PFPosConsPt+j1PFNegConsPt)/(jetPt->at(jetCand[0]));
      dR_PF12_ID_2 = deltaR(j1PFPosConsEta,j1PFPosConsPhi,j1PFNegConsEta,j1PFNegConsPhi);
      PF123PtFrac_ID_2 = (j1PFPosConsPt+j1PFNegConsPt+j1PFPhoConsPt)/(jetPt->at(jetCand[0]));
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
      float dPhi = DeltaPhi((*jetPhi)[j],pfMETPhi);
      //std::cout<<"DeltaPhi: "<<dPhi<<std::endl;
      if(dPhi < dPhimin){
        dPhimin = dPhi;
      }
    }
  return dPhimin;
}
std::vector<int> ZprimeJetsClass_MC::getJetCand(double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut){

  std::vector<int> tmpCand;
  tmpCand.clear();

  for(int p=0;p<nJet;p++)
    {

      bool kinematic = (*jetPt)[p] > jetPtCut && (*jetNHF)[p] < jetNHFCut && (*jetCHF)[p] > jetCHFCut && fabs((*jetEta)[p])<jetEtaCut;

      if((*jetPFLooseId)[p]==1 && kinematic){
        tmpCand.push_back(p);
      }
    }

  return tmpCand;

}

std::vector<int> ZprimeJetsClass_MC::JetVetoDecision(int jet_index, int ele_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;

  for(int i = 0; i < nJet; i++)
    {
      double deltar_ele = 0.0;
      double deltar_jet = 0.0;
      deltar_ele = deltaR(jetEta->at(i),jetPhi->at(i),eleEta->at(ele_index),elePhi->at(ele_index));
      deltar_jet = deltaR(jetEta->at(i),jetPhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index));
      if(deltar_ele>0.4 && deltar_jet>0.4 && jetPt->at(i) >30.0 && jetPFLooseId->at(i)==1)
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
bool ZprimeJetsClass_MC::btagVeto() {

  bool btagVeto = true;
  for(int i = 0; i < nJet; i++)
    {
      if(jetPt->at(i) >20.0 && jetEta->at(i) < 2.4 && jetCSV2BJetTags->at(i) > 0.8)
        btagVeto = false;
    }
  return btagVeto;
}

bool ZprimeJetsClass_MC::dPhiJetMETcut(std::vector<int> jets)
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
    if(DeltaPhi((*jetPhi)[j],pfMETPhi) < 0.5)
      break;
  }

  if(j==njetsMax)
    passes = true;

  return passes;

}

std::vector<int> ZprimeJetsClass_MC::electron_veto_tightID(int jet_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();

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
      //Make sure these get reset for every electron
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
	  pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0101;
	  pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00926;
	  pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.0336;
	  pass_HoverE = eleHoverE->at(i) < 0.0597;
	  pass_iso = EAcorrIso < 0.0354;
	  pass_ooEmooP = eleEoverPInv->at(i) < 0.012;
	  pass_d0 = abs(eleD0->at(i)) < 0.0111;
	  pass_dz = abs(eleDz->at(i)) < 0.0466;
	  pass_missingHits = eleMissHits->at(i) <= 2;
	  pass_convVeto = eleConvVeto->at(i) == 1;
	}
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
	{
	  pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0279;
	  pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00724;
	  pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.0918;
	  pass_HoverE = eleHoverE->at(i) < 0.0615;
	  pass_iso = EAcorrIso < 0.0646;
	  pass_ooEmooP = eleEoverPInv->at(i) < 0.00999;
	  pass_d0 = abs(eleD0->at(i)) < 0.0351;
	  pass_dz = abs(eleDz->at(i)) < 0.417;
	  pass_missingHits = eleMissHits->at(i) <= 1;
	  pass_convVeto = eleConvVeto->at(i) == 1;
	}
      //Electron passes Loose Electron ID cuts
      if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
	{
	  //Electron passes pt cut
	  if(elePt->at(i) > elePtCut)
	    {
	      //Electron does not overlap photon
	      if(deltaR(eleEta->at(i),elePhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.5)
		{
		  ele_cands.push_back(i);
		}
	    }
	}
    }
  return ele_cands;
}

std::vector<int> ZprimeJetsClass_MC::muon_veto_tightID(int jet_index, float muPtCut)
{
  // bool veto_passed = true; //pass veto if no good muon found
  std::vector<int> mu_cands;
  mu_cands.clear();

  bool pass_PFMuon = true;
  bool pass_globalMuon = true;
  // bool pass_trackerMuon = true;
  bool pass_chi2ndf = false;
  bool pass_chamberHit = false;
  bool pass_matchedStations = false;
  bool pass_dxy = false;
  bool pass_dz = false;
  bool pass_pixelHits = false;
  bool pass_trackLayers = false;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      // pass_globalMuon = muIsGlobalMuon->at(i);
      // pass_PFMuon = muIsPFMuon->at(i);
      // pass_trackerMuon = muIsTrackerMuon->at(i);
      pass_chi2ndf = muChi2NDF->at(i) < 10.0;
      pass_chamberHit = muMuonHits->at(i) > 0;
      pass_matchedStations = muStations->at(i) > 1;
      pass_dxy = fabs(muInnerD0->at(i)) < 0.2;
      pass_dz = fabs(muInnerDz->at(i)) < 0.5;
      pass_pixelHits = muPixelHits->at(i) > 0;
      pass_trackLayers = muTrkLayers->at(i) > 5;

      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.15;
      //Muon passes Tight Muon ID
      // if(pass_iso && pass_globalMuon && pass_PFMuon && pass_chi2ndf && pass_chamberHit && pass_matchedStations && pass_dxy && pass_dz && pass_pixelHits && pass_trackLayers)
      if(pass_iso && muIDbit->at(i)>>2&1==1)
	{
	  //Muon passes pt cut
	  if(muPt->at(i) > muPtCut)
	    {
	      //Muon does not overlap photon
	      if(deltaR(muEta->at(i),muPhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.5)
		{
		  mu_cands.push_back(i);
		}
	    }
	}
    }
  return mu_cands;
}

std::vector<int> ZprimeJetsClass_MC::electron_veto_looseID(int jet_index,int mu_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();

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
      //Make sure these get reset for every electron
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
	      //Electron does not overlap photon
	      if(deltaR(eleEta->at(i),elePhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.5)
		{
		  ele_cands.push_back(i);
		}
	    }
	}
    }
  return ele_cands;
}



//Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5
std::vector<int> ZprimeJetsClass_MC::muon_veto_looseID(int jet_index, int ele_index, float muPtCut)
{
  std::vector<int> mu_cands;
  mu_cands.clear();

  bool veto_passed = true; //pass veto if no good muon found
  bool pass_PFMuon = true;
  bool pass_globalMuon = true;
  bool pass_trackerMuon = true;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      // pass_PFMuon = muIsPFMuon->at(i);
      // pass_globalMuon = muIsGlobalMuon->at(i);
      // pass_trackerMuon = muIsTrackerMuon->at(i);
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      //Muon passes Loose Muon ID and PF-based combined relative, dBeta-corrected Loose Muon Isolation cuts
      //      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
      // if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon))
      if(muIDbit->at(i)>>0&1==1)
	{
	  //Muon passes pt cut
	  if(muPt->at(i) > muPtCut)
	    {
	      //Muon does not overlap photon
	      if(deltaR(muEta->at(i),muPhi->at(i),jetEta->at(jet_index),jetPhi->at(jet_index)) > 0.5 && deltaR(muEta->at(i),muPhi->at(i),eleEta->at(ele_index),elePhi->at(ele_index)) > 0.5)
		{
		  mu_cands.push_back(i);
		}
	    }
	}
    }
  return mu_cands;
}

