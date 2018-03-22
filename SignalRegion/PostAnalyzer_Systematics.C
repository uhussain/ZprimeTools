
////postAnalyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom postAnalyzer analyze
//
//To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jjbuch/LatestNtuples/ /afs/hep.wisc.edu/user/jjbuchanan/private/CMSSW_7_4_9/src/output.root -1 10000
//Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
//and storing the resulting histograms in the file output.root.
//
//To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
//root[0] TFile *f = new TFile("output.root");
//root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
//root[2] efficiency->Draw("AP")
//root[3] efficiency->SetTitle("Single photon trigger efficiency")
//root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
//root[5] efficiency->Draw("AP")
//

#define postAnalyzer_cxx
#include "postAnalyzer_mc_systematics.h"
#include <TH2.h>
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
  postAnalyzer t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void postAnalyzer::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  int nInspected;
  nInspected = 0;
  double nInspected_genWeighted;
  nInspected_genWeighted = 0.0;
  int nPhoCand, nMET170, nDphiPhoMET, nLeptonVeto, nDphiJetMET;
  nPhoCand = nMET170 = nDphiPhoMET = nLeptonVeto = nDphiJetMET = 0;
  double nPhoCand_weighted, nMET170_weighted, nDphiPhoMET_weighted, nLeptonVeto_weighted, nDphiJetMET_weighted;
  nPhoCand_weighted = nMET170_weighted = nDphiPhoMET_weighted = nLeptonVeto_weighted = nDphiJetMET_weighted = 0.0;

  std::vector<int> phoCand;
  phoCand.clear();
  std::vector<int> phoCandUp;
  phoCandUp.clear();
  std::vector<int> phoCandDown;
  phoCandDown.clear();

  std::vector<int> qcdden;
  qcdden.clear();
   
  std::vector<int> jetveto;
  jetveto.clear();
   
  TFile *file = new TFile("ewk_corr.root");
  TH1D *ewkCorrection = (TH1D*)file->Get("wnlg-130-o_p");
  TH1D *ewkCorrection_straightUp = (TH1D*)file->Get("wnlg-130-o_p_straightUp");
  TH1D *ewkCorrection_straightDown = (TH1D*)file->Get("wnlg-130-o_p_straightDown");
  TH1D *ewkCorrection_twistedUp = (TH1D*)file->Get("wnlg-130-o_p_twistedUp");
  TH1D *ewkCorrection_twistedDown = (TH1D*)file->Get("wnlg-130-o_p_twistedDown");
  TH1D *ewkCorrection_gammaUp = (TH1D*)file->Get("wnlg-130-o_p_gammaUp");
  TH1D *ewkCorrection_gammaDown = (TH1D*)file->Get("wnlg-130-o_p_gammaDown");
  TH1D *ewkCorrection_m = (TH1D*)file->Get("wnlg-130-o_m");
  TH1D *ewkCorrection_straightUp_m = (TH1D*)file->Get("wnlg-130-o_m_straightUp");
  TH1D *ewkCorrection_straightDown_m = (TH1D*)file->Get("wnlg-130-o_m_straightDown");
  TH1D *ewkCorrection_twistedUp_m = (TH1D*)file->Get("wnlg-130-o_m_twistedUp");
  TH1D *ewkCorrection_twistedDown_m = (TH1D*)file->Get("wnlg-130-o_m_twistedDown");
  TH1D *ewkCorrection_gammaUp_m = (TH1D*)file->Get("wnlg-130-o_m_gammaUp");
  TH1D *ewkCorrection_gammaDown_m = (TH1D*)file->Get("wnlg-130-o_m_gammaDown");
  ewkCorrection->Add(ewkCorrection_m);
  ewkCorrection_straightUp->Add(ewkCorrection_straightUp_m);
  ewkCorrection_straightDown->Add(ewkCorrection_straightDown_m);
  ewkCorrection_twistedUp->Add(ewkCorrection_twistedUp_m);
  ewkCorrection_twistedDown->Add(ewkCorrection_twistedDown_m);
  ewkCorrection_gammaUp->Add(ewkCorrection_gammaUp_m);
  ewkCorrection_gammaDown->Add(ewkCorrection_gammaDown_m);
  ewkCorrection->Scale(0.5);
  ewkCorrection_straightUp->Scale(0.5);
  ewkCorrection_straightDown->Scale(0.5);
  ewkCorrection_twistedUp->Scale(0.5);
  ewkCorrection_twistedDown->Scale(0.5);
  ewkCorrection_gammaUp->Scale(0.5);
  ewkCorrection_gammaDown->Scale(0.5);

  bool debug=true;
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Coming in: "<<std::endl;
  std::cout<<"nentries:"<<nentries<<std::endl;
  //Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;

  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
    {
    
      event_.clear();
      event_info.clear();
    
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      double inspected_event_weight = 1.0;
      fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0; //Generator may have given event negative weight
      nInspected_genWeighted += inspected_event_weight;
      nInspected += 1;
      //=1.0 for real data
      int ewkCorrectionBin = min(max(ewkCorrection->GetXaxis()->FindBin(LHephoET), 1), ewkCorrection->GetXaxis()->GetNbins());
      double event_weight=1.0;
      double event_weight_straightUp=1.0;
      double event_weight_straightDown=1.0;
      double event_weight_twistedUp=1.0;
      double event_weight_twistedDown=1.0;
      double event_weight_gammaUp=1.0;
      double event_weight_gammaDown=1.0;
      double uncorrected_weight=1.0;
    
      phoCand   = getPhoCand(175,1.4442,1);
      phoCandUp   = getPhoCandUp(175,1.4442,1);
      phoCandDown   = getPhoCandDown(175,1.4442,1);

      if(true || metFilters==1536)
	{
	  //  nMETFiltersPassed++;
	  //  if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1) || (HLTPho>>22&1 == 1))
	  //  {
	  //    nHLTPassed++;
	  if(phoCand.size() >0)
	    {
	      // if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCand[0]] - rho*EAchargedworst((*phoSCEta)[phoCand[0]]) ), 0.0) < 1.37 )
	      // {
	      Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCand[0]]/TMath::CosH((*phoSCEta)[phoCand[0]]));
	      event_weight = 1.0;
	      Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_straightUp = ewkCorrection_straightUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_straightDown = ewkCorrection_straightDown->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_twistedUp = ewkCorrection_twistedUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_twistedDown = ewkCorrection_twistedDown->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_gammaUp = ewkCorrection_gammaUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_gammaDown = ewkCorrection_gammaDown->GetBinContent(ewkCorrectionBin);
	      event_weight *= EWK_percent_adjustment;
	      event_weight *= NNLOCorrection(LHephoET);
	      event_weight*=(1.002 - 0.00004395*phoEt->at(phoCand[0])); //Trigger inefficiency correction
	      fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0; //Generator may have given event negative weight
	      event_weight_straightUp = event_weight*EWK_percent_adjustment_straightUp/EWK_percent_adjustment;
	      event_weight_straightDown = event_weight*EWK_percent_adjustment_straightDown/EWK_percent_adjustment;
	      event_weight_twistedUp = event_weight*EWK_percent_adjustment_twistedUp/EWK_percent_adjustment;
	      event_weight_twistedDown = event_weight*EWK_percent_adjustment_twistedDown/EWK_percent_adjustment;
	      event_weight_gammaUp = event_weight*EWK_percent_adjustment_gammaUp/EWK_percent_adjustment;
	      event_weight_gammaDown = event_weight*EWK_percent_adjustment_gammaDown/EWK_percent_adjustment;
	      uncorrected_weight *= (1.002 - 0.00004395*phoEt->at(phoCand[0]));
	      fabs(genWeight) > 0.0 ? uncorrected_weight *= genWeight/fabs(genWeight) : uncorrected_weight = 0;
	      nPhoCand++;
	      nPhoCand_weighted += event_weight;
        
	      jetveto = JetVetoDecision(phoCand[0]);
                
	      if(pfMET>170)
		{
		  nMET170++;
		  nMET170_weighted += event_weight;
          
		  Float_t MET_to_use = pfMET;
		  Float_t METPhi_to_use = pfMETPhi;
          
		  fillHistos(0,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  if(DeltaPhi(phoPhi->at(phoCand[0]),METPhi_to_use)>0.5)
		    {
		      nDphiPhoMET++;
		      nDphiPhoMET_weighted += event_weight;
		      fillHistos(1,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		      if(electron_veto_looseID(phoCand[0],10) && muon_veto_looseID(phoCand[0],10))
			{
			  nLeptonVeto++;
			  nLeptonVeto_weighted += event_weight;
			  fillHistos(2,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			  if(dPhiJetMET_veto(jetveto,METPhi_to_use))
			    {
			      nDphiJetMET++;
			      nDphiJetMET_weighted += event_weight;
			      fillHistos(3,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      fillHistos(30,event_weight*NNLOCorrection_err(LHephoET),phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      if(uncorrectedPhoEt/MET_to_use < 1.4){
				fillHistos(24,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(21,event_weight_straightUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(22,event_weight_straightDown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(31,event_weight_twistedUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(32,event_weight_twistedDown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(33,event_weight_gammaUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(34,event_weight_gammaDown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(20,uncorrected_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
				fillHistos(23,event_weight*NNLOCorrection_err(LHephoET),phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      }
                
			      //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
			    }
			}
		    }
		}
	      if(pfMET_T1JESUp > 170)
		{

		  Float_t MET_to_use = pfMET_T1JESUp;
		  Float_t METPhi_to_use = pfMETPhi_T1JESUp;
          
		  fillHistos(4,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  if(DeltaPhi(phoPhi->at(phoCand[0]),METPhi_to_use)>0.5)
		    {
		      fillHistos(5,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		      if(electron_veto_looseID(phoCand[0],10) && muon_veto_looseID(phoCand[0],10))
			{
			  fillHistos(6,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			  if(dPhiJetMET_veto(jetveto,METPhi_to_use))
			    {
			      fillHistos(7,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      if(uncorrectedPhoEt/MET_to_use < 1.4){
				fillHistos(25,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      }
			      //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
			    }
			}
		    }
		}
	      if(pfMET_T1JESDo > 170)
		{
		  Float_t MET_to_use = pfMET_T1JESDo;
		  Float_t METPhi_to_use = pfMETPhi_T1JESDo;
          
		  fillHistos(8,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  if(DeltaPhi(phoPhi->at(phoCand[0]),METPhi_to_use)>0.5)
		    {
		      fillHistos(9,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		      if(electron_veto_looseID(phoCand[0],10) && muon_veto_looseID(phoCand[0],10))
			{
			  fillHistos(10,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			  if(dPhiJetMET_veto(jetveto,METPhi_to_use))
			    {
			      fillHistos(11,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      if(uncorrectedPhoEt/MET_to_use < 1.4){
				fillHistos(26,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      }
			      //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
			    }
			}
		    }
		}
	      // }
	    }
	  if(phoCandUp.size() > 0)
	    {
	      // if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCandUp[0]] - rho*EAchargedworst((*phoSCEta)[phoCandUp[0]]) ), 0.0) < 1.37 )
	      // {
	      event_weight=1.0;
	      Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCandUp[0]]/TMath::CosH((*phoSCEta)[phoCandUp[0]]));
	      uncorrectedPhoEt += 0.006*uncorrectedPhoEt;
	      Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_straightUp = ewkCorrection_straightUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_straightDown = ewkCorrection_straightDown->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_twistedUp = ewkCorrection_twistedUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_twistedDown = ewkCorrection_twistedDown->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_gammaUp = ewkCorrection_gammaUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_gammaDown = ewkCorrection_gammaDown->GetBinContent(ewkCorrectionBin);
	      event_weight *= EWK_percent_adjustment;
	      event_weight *= NNLOCorrection(LHephoET);
	      event_weight *= (1.002 - 0.00004395*phoEt->at(phoCandUp[0])); //Trigger inefficiency correction
	      fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0; //Generator may have given event negative weight
	      event_weight_straightUp = event_weight*EWK_percent_adjustment_straightUp/EWK_percent_adjustment;
	      event_weight_straightDown = event_weight*EWK_percent_adjustment_straightDown/EWK_percent_adjustment;
	      event_weight_twistedUp = event_weight*EWK_percent_adjustment_twistedUp/EWK_percent_adjustment;
	      event_weight_twistedDown = event_weight*EWK_percent_adjustment_twistedDown/EWK_percent_adjustment;
	      event_weight_gammaUp = event_weight*EWK_percent_adjustment_gammaUp/EWK_percent_adjustment;
	      event_weight_gammaDown = event_weight*EWK_percent_adjustment_gammaDown/EWK_percent_adjustment;
        
	      jetveto = JetVetoDecision(phoCandUp[0]);
        
                
	      if(pfMET_T1PESUp>170)
		{
		  Float_t MET_to_use = pfMET_T1PESUp;
		  Float_t METPhi_to_use = pfMETPhi_T1PESUp;
          
		  fillHistos(12,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  if(DeltaPhi(phoPhi->at(phoCandUp[0]),METPhi_to_use)>0.5)
		    {
		      fillHistos(13,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		      if(electron_veto_looseID(phoCandUp[0],10) && muon_veto_looseID(phoCandUp[0],10))
			{
			  fillHistos(14,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			  if(dPhiJetMET_veto(jetveto,METPhi_to_use))
			    {
			      fillHistos(15,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      if(uncorrectedPhoEt/MET_to_use < 1.4){
				fillHistos(27,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      }
			      //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
			    }
			}
		    }
		}
	      // }
	    }
	  if(phoCandDown.size() > 0)
	    {
	      // if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCandDown[0]] - rho*EAchargedworst((*phoSCEta)[phoCandDown[0]]) ), 0.0) < 1.37 )
	      // {
	      event_weight=1.0;
	      Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCandDown[0]]/TMath::CosH((*phoSCEta)[phoCandDown[0]]));
	      uncorrectedPhoEt -= 0.006*uncorrectedPhoEt;
	      Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_straightUp = ewkCorrection_straightUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_straightDown = ewkCorrection_straightDown->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_twistedUp = ewkCorrection_twistedUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_twistedDown = ewkCorrection_twistedDown->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_gammaUp = ewkCorrection_gammaUp->GetBinContent(ewkCorrectionBin);
	      Double_t EWK_percent_adjustment_gammaDown = ewkCorrection_gammaDown->GetBinContent(ewkCorrectionBin);
	      event_weight *= EWK_percent_adjustment;
	      event_weight *= NNLOCorrection(LHephoET);
	      event_weight *= (1.002 - 0.00004395*phoEt->at(phoCandDown[0])); //Trigger inefficiency correction
	      fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0; //Generator may have given event negative weight
	      event_weight_straightUp = event_weight*EWK_percent_adjustment_straightUp/EWK_percent_adjustment;
	      event_weight_straightDown = event_weight*EWK_percent_adjustment_straightDown/EWK_percent_adjustment;
	      event_weight_twistedUp = event_weight*EWK_percent_adjustment_twistedUp/EWK_percent_adjustment;
	      event_weight_twistedDown = event_weight*EWK_percent_adjustment_twistedDown/EWK_percent_adjustment;
	      event_weight_gammaUp = event_weight*EWK_percent_adjustment_gammaUp/EWK_percent_adjustment;
	      event_weight_gammaDown = event_weight*EWK_percent_adjustment_gammaDown/EWK_percent_adjustment;
        
	      jetveto = JetVetoDecision(phoCandDown[0]);
        
                
	      if(pfMET_T1PESDo>170)
		{
		  Float_t MET_to_use = pfMET_T1PESDo;
		  Float_t METPhi_to_use = pfMETPhi_T1PESDo;
          
		  fillHistos(16,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  if(DeltaPhi(phoPhi->at(phoCandDown[0]),METPhi_to_use)>0.5)
		    {
		      fillHistos(17,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		      if(electron_veto_looseID(phoCandDown[0],10) && muon_veto_looseID(phoCandDown[0],10))
			{
			  fillHistos(18,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			  if(dPhiJetMET_veto(jetveto,METPhi_to_use))
			    {
			      fillHistos(19,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      if(uncorrectedPhoEt/MET_to_use < 1.4){
				fillHistos(28,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
			      }
			      //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
			    }
			}
		    }
		}
	      // }
	    }
	}
    
      tree->Fill();
    
      if (jentry%reportEvery == 0)
	{
	  std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
	}
    }

  if((nentriesToCheck-1)%reportEvery != 0)
    std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop();
  std::cout<<"All events checked."<<std::endl;
  //Report
  std::cout << "RealTime : " << sw.RealTime() / 60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << sw.CpuTime()  / 60.0 << " minutes" << std::endl;
  std::cout << std::endl;
  std::cout << "Number of events inspected: " << nInspected << std::endl;
  std::cout << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl;
  std::cout<<std::endl;
  cout<<"Unweighted"<<endl;
  cout<<"----------"<<endl;
  cout<<"nPhoCand: "<<nPhoCand<<", "<<nPhoCand/nInspected<<" of previous"<<endl;
  cout<<"nMET170: "<<nMET170<<", "<<nMET170/nPhoCand<<" of previous"<<endl;
  cout<<"nDphiPhoMET: "<<nDphiPhoMET<<", "<<nDphiPhoMET/nMET170<<" of previous"<<endl;
  cout<<"nLeptonVeto: "<<nLeptonVeto<<", "<<nLeptonVeto/nDphiPhoMET<<" of previous"<<endl;
  cout<<"nDphiJetMET: "<<nDphiJetMET<<", "<<nDphiJetMET/nLeptonVeto<<" of previous"<<endl;
  cout<<endl;
  cout<<"Weighted"<<endl;
  cout<<"--------"<<endl;
  cout<<"nPhoCand_weighted: "<<nPhoCand_weighted<<", "<<nPhoCand_weighted/nInspected<<" of previous"<<endl;
  cout<<"nMET170_weighted: "<<nMET170_weighted<<", "<<nMET170_weighted/nPhoCand_weighted<<" of previous"<<endl;
  cout<<"nDphiPhoMET_weighted: "<<nDphiPhoMET_weighted<<", "<<nDphiPhoMET_weighted/nMET170_weighted<<" of previous"<<endl;
  cout<<"nLeptonVeto_weighted: "<<nLeptonVeto_weighted<<", "<<nLeptonVeto_weighted/nDphiPhoMET_weighted<<" of previous"<<endl;
  cout<<"nDphiJetMET_weighted: "<<nDphiJetMET_weighted<<", "<<nDphiJetMET_weighted/nLeptonVeto_weighted<<" of previous"<<endl;
}

void postAnalyzer::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ADD","ADD");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);
  fileName->cd();
  
  Float_t PtBins[7]={175.,200.,250., 300., 400., 600., 1000.0};
  Float_t MetBins[7]={170.,200.,250., 300., 400., 600., 1000.0};
  Float_t MTBins[10]={0.,200.,300.,400.,500.,600.,700.,800.,1000.,1200.};
  Float_t dPhiJetMETBins[14]={0.0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.1416};
  //h_phoIEtaIPhi = new TH2F("h_phoIEtaIPhi","Photon p_{T} > 175 GeV, E^{miss}_{T} > 140 GeV",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi->Sumw2();

  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<40; i++)
    {
      char ptbins[100];
      sprintf(ptbins, "_%d", i);
      std::string histname(ptbins);
      h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
      h_photon_Et[i] = new TH1F(("Photon_Et"+histname).c_str(), "Photon_Et",6,PtBins);h_photon_Et[i]->Sumw2();
      h_photon_Et_range[i] = new TH1F(("Photon_Et_range"+histname).c_str(), "Photon_Et",6,PtBins);h_photon_Et_range[i]->Sumw2();
      h_photon_eta[i] = new TH1F(("Photon_eta"+histname).c_str(), "Photon_eta",20,-1.4442,1.4442);h_photon_eta[i]->Sumw2();
      h_photon_SCEta[i] = new TH1F(("Photon_SCeta"+histname).c_str(), "Photon_SCeta",20,-1.4442,1.4442);h_photon_SCEta[i]->Sumw2();
      h_photon_phi[i] = new TH1F(("Photon_phi"+histname).c_str(), "Photon_phi", 64,0,3.2);h_photon_phi[i]->Sumw2();
      h_photon_SCPhi[i] = new TH1F(("Photon_SCphi"+histname).c_str(), "Photon_SCphi", 20,0,3.1416);h_photon_SCPhi[i]->Sumw2();
      h_photon_IDbit[i] = new TH1F(("Photon_ID_bit"+histname).c_str(), "Photon_ID_bit",8,0,8);h_photon_IDbit[i]->Sumw2();
      h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",6,MetBins);h_pfMET[i]->Sumw2();
      h_pfMET_300[i] = new TH1F(("h_pfMET_300"+histname).c_str(), "pfMET",25,0,300);h_pfMET_300[i]->Sumw2();
      h_dPhi[i] = new TH1F(("h_dPhi"+histname).c_str(),"h_dPhi",40,0,3.2);h_dPhi[i]->Sumw2();
      h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",20,0,20);h_nJet[i]->Sumw2();
      h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",30,20,1000);h_leadingJetPt[i]->Sumw2();
      h_leadingJetPt_300[i] = new TH1F(("leadingJetPt_300"+histname).c_str(),"leadingJetPt_300",25,0,300);h_leadingJetPt_300[i]->Sumw2();
      h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",40,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();
      // h_phoIEtaIPhi[i] = new TH2F(("h_phoIEtaIPhi"+histname).c_str(),"Photon p_{T} > 175 GeV, E^{miss}_{T} > 140 GeV",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi[i]->Sumw2();
      h_PTMET[i] = new TH1F(("PTMET"+histname).c_str(),"P_{T}/Missing E_{T}",50,0,3);h_PTMET[i]->Sumw2();
      h_Mt[i]= new TH1F(("Mt"+histname).c_str(),"MT",9,MTBins);h_Mt[i]->Sumw2();
      h_min_dphijetmet[i] = new TH1F(("h_min_dphijetmet"+histname).c_str(),"h_min_dphijetmet",13,dPhiJetMETBins);h_min_dphijetmet[i]->Sumw2();
      h_pfMETsumEt[i] = new TH1F(("pfMETsumEt"+histname).c_str(),"pfMETsumEt",6,MetBins);
      h_METoverSqrtSumEt_extended[i] = new TH1F(("METoverSqrtSumEt_extended"+histname).c_str(),"METoverSqrtSumEt",30,0,30);
      h_METoverSqrtSumEt[i] = new TH1F(("METoverSqrtSumEt"+histname).c_str(),"METoverSqrtSumEt",30,0,10);
    }
}

//Fill the sequential histos at a particular spot in the sequence
void postAnalyzer::fillHistos(int histoNumber, double event_weight,int index,std::vector<int> jets,Float_t MET_to_use,Float_t uncorrectedPhoEt,Float_t METPhi_to_use)
{

  h_photon_Et[histoNumber]->Fill(uncorrectedPhoEt,event_weight);
  h_photon_Et_range[histoNumber]->Fill(uncorrectedPhoEt,event_weight);
  h_photon_eta[histoNumber]->Fill(phoEta->at(index),event_weight);
  h_photon_SCEta[histoNumber]->Fill(phoSCEta->at(index),event_weight);
  h_photon_phi[histoNumber]->Fill(phoPhi->at(index),event_weight);
  h_photon_SCPhi[histoNumber]->Fill(phoSCPhi->at(index),event_weight);
  h_pfMET[histoNumber]->Fill(MET_to_use,event_weight);
  h_pfMET_300[histoNumber]->Fill(MET_to_use,event_weight);
  double dPhi = DeltaPhi(phoPhi->at(index),METPhi_to_use);
  h_dPhi[histoNumber]->Fill(dPhi,event_weight);
  h_nJet[histoNumber]->Fill(jets.size(),event_weight);
  // h_phoIEtaIPhi[histoNumber]->Fill(phoIPhi->at(index),phoIEta->at(index),event_weight);
  //h_PTMET[histoNumber]->Fill(phoEt->at(index)/pfMET,event_weight);
  h_PTMET[histoNumber]->Fill((phoSCRawE->at(index)/TMath::CosH(phoSCEta->at(index)))/MET_to_use,event_weight);
  h_Mt[histoNumber]->Fill(sqrt(2*uncorrectedPhoEt*MET_to_use*(1-TMath::Cos(dPhi))),event_weight);
  h_nVtx[histoNumber]->Fill(nVtx);
  if(jets.size()>0){
    h_leadingJetPt[histoNumber]->Fill(jetPt->at(jets[0]),event_weight);
    h_leadingJetEta[histoNumber]->Fill(jetEta->at(jets[0]),event_weight);
    int max_njets = jets.size();
    if(jets.size() > 4)
      max_njets = 4;
    double min_dphijetmet = TMath::Pi();
    for(int i = 0; i < max_njets; i++)
      {
        double dphijetmet = DeltaPhi(jetPhi->at(jets[i]),METPhi_to_use);
        if(dphijetmet < min_dphijetmet)
          min_dphijetmet = dphijetmet;
      }
    h_min_dphijetmet[histoNumber]->Fill(min_dphijetmet,event_weight);
  }
  h_pfMETsumEt[histoNumber]->Fill(pfMETsumEt,event_weight);
  h_METoverSqrtSumEt_extended[histoNumber]->Fill(MET_to_use/sqrt(pfMETsumEt),event_weight);
  h_METoverSqrtSumEt[histoNumber]->Fill(MET_to_use/sqrt(pfMETsumEt),event_weight);
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
double postAnalyzer::DeltaPhi(double phi1, double phi2)
{
  double pi = TMath::Pi();
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}


//---------------------------------------------------                                                                                                                                
// get a photon candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> postAnalyzer::getPhoCand(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {
      double  uncorrectedPhoEt = (phoSCRawE->at(p)/TMath::CosH(phoSCEta->at(p)));
      bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
		       ((*phoHoverE)[p]                <  0.0260   ) &&
		       ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.01040 ) &&
		       ((*phohasPixelSeed)[p]              ==  0      ) &&
		       ( TMath::Max( ( (*phoYuPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 1.146 )  &&
		       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (2.792 + (0.0112 * uncorrectedPhoEt) + (0.000028 * pow(uncorrectedPhoEt, 2.0))) )  &&
		       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (2.176 + (0.0043 * uncorrectedPhoEt)) ) 
		       );
      
      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001 && (*phoR9)[p] < 1;
      //      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9;

      if(photonId && kinematic && noncoll){
	tmpCand.push_back(p);
      }                                                                                                                                                              
    }                                                                                                                                                                

  return tmpCand;

}

std::vector<int> postAnalyzer::getPhoCandUp(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {
      double  uncorrectedPhoEt = (phoSCRawE->at(p)/TMath::CosH(phoSCEta->at(p)));
      uncorrectedPhoEt += 0.006*uncorrectedPhoEt;
      bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
		       ((*phoHoverE)[p]                <  0.0260   ) &&
		       ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.01040 ) &&
		       ((*phohasPixelSeed)[p]              ==  0      ) &&
		       ( TMath::Max( ( (*phoYuPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 1.146 )  &&
		       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (2.792 + (0.0112 * uncorrectedPhoEt) + (0.000028 * pow(uncorrectedPhoEt, 2.0))) )  &&
		       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (2.176 + (0.0043 * uncorrectedPhoEt)) ) 
		       );
      
      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001 && (*phoR9)[p] < 1;
      //      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9;

      if(photonId && kinematic && noncoll){
	tmpCand.push_back(p);
      }                                                                                                                                                              
    }                                                                                                                                                                

  return tmpCand;

}

std::vector<int> postAnalyzer::getPhoCandDown(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {
      double  uncorrectedPhoEt = (phoSCRawE->at(p)/TMath::CosH(phoSCEta->at(p)));
      uncorrectedPhoEt -= 0.006*uncorrectedPhoEt;
      bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
		       ((*phoHoverE)[p]                <  0.0260   ) &&
		       ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.01040 ) &&
		       ((*phohasPixelSeed)[p]              ==  0      ) &&
		       ( TMath::Max( ( (*phoYuPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 1.146 )  &&
		       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (2.792 + (0.0112 * uncorrectedPhoEt) + (0.000028 * pow(uncorrectedPhoEt, 2.0))) )  &&
		       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (2.176 + (0.0043 * uncorrectedPhoEt)) ) 
		       );
      
      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001 && (*phoR9)[p] < 1;
      //      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9;

      if(photonId && kinematic && noncoll){
	tmpCand.push_back(p);
      }                                                                                                                                                              
    }                                                                                                                                                                

  return tmpCand;

}






std::vector<int> postAnalyzer::getPhoCand1(double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {

      bool kinematic = fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
		       ((*phoHoverE)[p]                <  0.05   ) &&
		       ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0101 ) &&
		       ((*phohasPixelSeed)[p]              ==  0      ) &&
		       ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.21 )  &&
		       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (0.65 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) )  &&
		       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (0.18 + (0.0053 * (*phoEt)[p])) ) );
      
      if(photonId && kinematic){
	tmpCand.push_back(p);
      }                                                                                                                                                              
    }                                                                                                                                                               

  return tmpCand;

}


std::vector<int> postAnalyzer::getQcdden(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {
      bool upperBound=false;
      bool lowerBound =false;

      double  maxPFCharged= TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]);
      double  maxPFPhoton = TMath::Min(5.0*(0.81 + (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]);
      double  maxPFNeutral= TMath::Min(5.0*(1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]);
      
      bool kinematic = (*phoEt)[p] > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      upperBound = (
		    ((*phoHoverE)[p]                <  0.05   ) &&
		    ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0106 ) &&
		    ((*phohasPixelSeed)[p]              ==  0      ) &&
		    ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < maxPFCharged )  &&
		    ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < maxPFNeutral )  &&
		    ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < maxPFPhoton ));

      lowerBound = (
		    ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) > 3.32 )  ||
		    ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) > (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  ||
		    ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) > (0.81 + (0.0053 * (*phoEt)[p])) ));

      
      if(upperBound && lowerBound && kinematic){
	tmpCand.push_back(p);
      }                                                                                                                                                              
    }                                                                                                                                                                

  return tmpCand;

}

// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5
Double_t postAnalyzer::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0;

  return EffectiveArea;
}

Double_t postAnalyzer::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1064;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1026;
  return EffectiveArea;
}

Double_t postAnalyzer::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0597;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0807;
  return EffectiveArea;
}

Double_t postAnalyzer::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1210;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1107;
  return EffectiveArea;
}

//Returns true if veto passed                                                                                                                                                                                                                                                                                                
//Veto failed if an electron is found that passes Loose Electron ID and elePtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                                 
//Always true if no electrons with |SC eta| < 2.5, since ID always fails for |SC eta| > 2.

bool postAnalyzer::electron_veto_looseID(int pho_index, float elePtCut)
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
      // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
      if(eleIDbit->at(i)>>1&1==1)
	{
	  //Electron passes pt cut                                                                                                                                                                                                                                                                             
	  if(elePt->at(i) > elePtCut)
	    {
	      //Electron does not overlap photon                                                                                                                                                                                                                                                           
	      if(dR(eleSCEta->at(i),eleSCPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
		{
		  veto_passed = false;
		  break;
		}
	    }
	}
    }
  return veto_passed;
}



//Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                    
bool postAnalyzer::muon_veto_looseID(int pho_index, float muPtCut)
{
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
      // if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon))
      if(muIDbit->at(i)>>0&1==1)
	{
	  //Muon passes pt cut                                                                                                                                                                                                                                                                                 
	  if(muPt->at(i) > muPtCut)
	    {
	      //Muon does not overlap photon                                                                                                                                                                                                                                                               
	      if(dR(muEta->at(i),muPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
		{
		  veto_passed = false;
		  break;
		}
	    }
	}
    }
  return veto_passed;
}


std::vector<int> postAnalyzer::JetVetoDecision(int pho_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;
  float value =0.0;

  for(int i = 0; i < nJet; i++)
    {

      if(0.0 < abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.5)
        value =-0.8;
      else if(2.5 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.75)
        value =-0.95;
      else if(2.75 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <3.0)
        value =-0.97;
      else if(3.00 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <5.0)
        value =-0.99;
      else
        continue;



      //      std::cout<<"Jet size: "<<nJet<<std::endl;
      double deltar = 0.0 ;
      //      std::cout<<"Jet no:"<<i<<"coming here pujetid: "<<pfJet_pt[i]<<std::endl;
      //if(OverlapWithMuon(jetEta->at(i),jetPhi->at(i)))     continue;
      //      std::cout<<"Jet no:"<<i<<"coming here OverlapWithMuon: "<<pfJet_pt[i]<<std::endl;
      //      if(OverlapWithElectron(jetEta->at(i),jetPhi->at(i)))   continue;
      if(pho_index>=0){
        deltar= dR(jetEta->at(i),jetPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index));
        //std::cout<<"Delta R between photon and jet="<<dR<< "jet pt"<<pfJet_pt[i]<<std::endl;                                                                                       
      }
      if(deltar>0.4 && jetPt->at(i) >30.0 && jetPFLooseId->at(i)==1 && jetPUID->at(i)>value)
        {
          jetindex.push_back(i);
        }

      
    }


  //  std::cout<<"Jet size: "<< jetindex.size()<<std::endl;
  //if(jetindex.size()>1)jetVeto = false;                                                                                                                                            
  return jetindex;

}



bool postAnalyzer::OverlapWithMuon(double eta, double phi){
  
  bool overlap = false;
  //  std::cout<<"No of muon:"<<Muon_n<<std::endl;
  bool pass_PFMuon = false;
  bool pass_globalMuon = false;
  bool pass_trackerMuon = false;
  bool pass_iso = false;

  Float_t zero1 = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      pass_PFMuon = muIsPFMuon->at(i);
      pass_globalMuon = muIsGlobalMuon->at(i);
      pass_trackerMuon = muIsTrackerMuon->at(i);
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero1,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
        {
          if(muPt->at(i) > 10.)
            {
              if(dR(muEta->at(i),muPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }

  return overlap;


}



bool postAnalyzer::OverlapWithElectron(double eta, double phi){
  bool overlap = false;

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
      // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
      if(eleIDbit->at(i)>>1&1==1)
        {
          //Electron passes pt cut                                                                                                                                                                                                                                                                                          
          if(elePt->at(i) > 10.)
            {
              //Electron does not overlap photon                                                                                                                                                                                                                                                                            
              if(dR(eleSCEta->at(i),eleSCPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }





}


double postAnalyzer::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}





bool postAnalyzer::dPhiJetMET_veto(std::vector<int> jets, Float_t METPhi_to_use)
{
  //pass veto if no jet is found within DeltaPhi(jet,MET) < 0.5                                                                                         
  bool passes = false;
  int njetsMax = jets.size();
  //Only look at first four jets                                                                                                                        
  if(njetsMax > 4)
    njetsMax = 4;
  int j = 0;
  for(; j < njetsMax; j++)
    {
      //fail veto if a jet is found with DeltaPhi(jet,MET) < 0.5                                                                                    
      if(DeltaPhi(jetPhi->at(jets[j]),METPhi_to_use) < 0.5)
	break;
    }
  if(j==njetsMax)
    passes = true;

  return passes;
}

//WG only
Double_t postAnalyzer::NNLOCorrection(Double_t pt){
  Float_t corr = 1.40131;
  if(pt >= 175.0   && pt <190.0   ) corr = 1.40131;
  if(pt >= 190.0   && pt <250.0   ) corr = 1.37097;
  if(pt >= 250.0   && pt <400.0   ) corr = 1.3102;
  if(pt >= 400.0   && pt <700.0   ) corr = 1.261;
  if(pt >= 700.0  ) corr = 1.14949;
  return corr;
}

//WG only
Double_t postAnalyzer::NNLOCorrection_err(Double_t pt){
  Float_t corr = 1.0662;
  if(pt >= 175.0   && pt <190.0   ) corr = 1.0662 ; 
  if(pt >= 190.0   && pt <250.0   ) corr = 1.0665 ;
  if(pt >= 250.0   && pt <400.0   ) corr = 1.0701 ;
  if(pt >= 400.0   && pt <700.0   ) corr = 1.0762  ;
  if(pt >= 700.0  ) corr = 1.0864 ;
  return corr;
}

double postAnalyzer::phoSF(Float_t phoET, int sysbin)
{
  // int kinematic_bin = -1;
  // if (phoET > 175.0 && phoET <= 200.0) kinematic_bin = 0;
  // else if (phoET > 200.0 && phoET <= 250.0) kinematic_bin = 1;
  // else if (phoET > 250.0 && phoET <= 300.0) kinematic_bin = 2;
  // else if (phoET > 300.0 && phoET <= 350.0) kinematic_bin = 3;
  // else if (phoET > 350.0 && phoET <= 400.0) kinematic_bin = 4;
  // else if (phoET > 400.0) kinematic_bin = 5;
  
  // double central[6] = {1.00865, 0.999165, 1.01629, 0.997199, 0.986974, 0.998906};
  // double errUp[6] = {1.02474, 1.01356, 1.03537, 1.01968, 1.00945, 1.01444};
  // double errDown[6] = {0.992552, 0.984771, 0.997207, 0.974716, 0.964501, 0.983372};
  
  // double return_value = -999999.9;
  // if (kinematic_bin >= 0 && kinematic_bin <= 5){
  //   if (sysbin == 0) return_value = central[kinematic_bin];
  //   else if (sysbin == 1) return_value = errUp[kinematic_bin];
  //   else if (sysbin == 2) return_value = errDown[kinematic_bin];
  // }
  
  return 1.0;
}
