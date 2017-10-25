# Turn Zprime ntuples into Plots

```
This repository contains pacakges to analyze ntuples for the signal region and the control region for the Mono-Z' Jet + MET analysis.
Jobs are submitted to condor through the MakeCondorFiles*.csh scripts.
1) MakeCondorFiles_data.csh for data
2) MakeCondorFiles.csh for MC
3) MakeCondorFiles_WZ.csh for Wlv+Jets and Zvv+Jets MC as EWK+NNLO corrections need to be applied using the Monojet recipe via kfactors.root
Note: In the Control Region directory, MakeCondorFiles_WZ.csh is also used for the DYJetsToLL samples as we apply the same EWK+NNLO corrections
to the DYSamples for Control Region Plots.

The submit*.sh scripts can be changed as necessary and are just submitting jobs en masse.

Instructions:

```bash
cmsrel CMSSW_8_0_26_patch1
cd $CMSSW_BASE/src
cmsenv
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
voms-proxy-init --voms=cms --valid=192:00
./submit.sh
```

Stitching DY inclusive sample and WJets inclusive sample withe the HT-binned samples. Currently done by hand.

Edit the corresponding ZprimeJetsClass_MC*.C files to add the genHT < 100 cut in the beginning and then process the inclusive samples as follows:

Signal Region:

```bash
./rootcom ZprimeJetsClass_MC analyze1
./MakeCondorFiles.csh analyze1 /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL/170810_152817/0000/ postDY_MLM_0.root -1 10000 DYMLM_0 PU_Central.root
./rootcom ZprimeJetsClass_MC_WJets analyzeWJets1
./MakeCondorFiles_WZ.csh analyzeWJets1 /hdfs/store/user/uhussain/Zprime_Ntuples_Aug10/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJets/170810_152529/0000/ postWJets_MLM_0.root -1 10000 W_0 kfactors.root PU_Central.root
```
Control Region:

```bash
./rootcom ZprimeJetsClass_MC_ZJets analyze1
./MakeCondorFiles_WZ.csh analyze1 /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL/170810_152817/0000/ postDY_MLM_0.root -1 10000 DYMLM_0 kfactors.root PU_Central.root
./rootcom ZprimeJetsClass_MC_WJets analyzeWJets1
./MakeCondorFiles_WZ.csh analyzeWJets1 /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJets/170810_152529/0000/ postWJets_MLM_0.root -1 10000 W_0 kfactors.root PU_Central.root
```
After your condor jobs are finished and you have post*.root files, make plots using stackplotter*.C scripts and edit/automate the script to your needs

```bash
root -l stackplotter.C
```
