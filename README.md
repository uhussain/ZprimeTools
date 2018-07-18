# Mono-Z' Jet + MET analysis

This repository contains pacakges to analyze ntuples for the signal region and the control region 
for the Mono-Z' Jet + MET analysis.
Jobs are submitted to condor through the SubmitCondor*.sh scripts written by E.Koenig.
1) SubmitCondor.sh for data
2) SubmitCondor_signal.sh for MC
For Wlv+Jets and Zvv+Jets MC as EWK+NNLO corrections need to be applied 
using the Monojet recipe via kfactors.root
4) kfactors.root can be found in /nfs_scratch/uhussain/MonoZprimeJet_postanalyzer_jobsubmission/ 

The submit*.sh scripts can be changed as necessary and are just submitting jobs en masse.

Update instructions regarding Control Regions. Currently up to date in Evan's branch.

Instructions:

```bash
cmsrel CMSSW_8_0_26_patch1
cd $CMSSW_BASE/src
cmsenv
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
voms-proxy-init --voms=cms --valid=192:00
./submit.sh
./submitSignal.sh
```


After your condor jobs are finished and you have post*.root files, you can merge them using merge.sh script or plot directly.
      
Then you can save histos using saveplot.C script or use the plotter.C and edit/automate the script to your needs.
Add more variables that you want to plot in samplenames.txt file.

```bash
./rootcom saveplot save
./save pfMET_8 pfMET_9 pfMET_10
./rootcom plotter plot
./plot h_cutflowow
```

PileupCorrection Recipe adopted from Nick Smith (U.Wisconsin)
https://gitlab.cern.ch/ncsmith/PileupWeights
Get your "processedLumis.json" from crab --report
```bash
./getDataDistribution.sh
python calculatePileupCorrections.py
```
You will need to edit the calculatePileupCorrections.py file in case you are using a different mc profile than what is included
