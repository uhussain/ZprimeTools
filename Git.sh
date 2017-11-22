#!/bin/bash
gitdir=~/WORKINGAREA/ZprimeTools

cp changes.txt ${gitdir}/changes.txt
cp README.md ${gitdir}/README.md

cd ControlRegion/
cp allsubmit.sh ${gitdir}/ControlRegion/allsubmit.sh
cp kfactors.root ${gitdir}/ControlRegion/kfactors.root
cp MakeCondorFiles.csh ${gitdir}/ControlRegion/MakeCondorFiles.csh
cp MakeCondorFiles_data.csh ${gitdir}/ControlRegion/MakeCondorFiles_data.csh
cp MakeCondorFiles_WZ.csh ${gitdir}/ControlRegion/MakeCondorFiles_WZ.csh
cp PU_Central.root ${gitdir}/ControlRegion/PU_Central.root
cp PU_minBiasDOWN.root ${gitdir}/ControlRegion/PU_minBiasDOWN.root
cp PU_minBiasUP.root ${gitdir}/ControlRegion/PU_minBiasUP.root
cp reset.sh ${gitdir}/ControlRegion/reset.sh
cp rootcom ${gitdir}/ControlRegion/rootcom
cp stack_plotterlog.C ${gitdir}/ControlRegion/stack_plotterlog.C
cp stack_plotterlogRebin.C ${gitdir}/ControlRegion/stack_plotterlogRebin.C
cp submitDY.sh ${gitdir}/ControlRegion/submitDY.sh
cp submitGH.sh ${gitdir}/ControlRegion/submitGH.sh
cp submit.sh ${gitdir}/ControlRegion/submit.sh
cp ZprimeJetsClass.C ${gitdir}/ControlRegion/ZprimeJetsClass.C
cp ZprimeJetsClass.h ${gitdir}/ControlRegion/ZprimeJetsClass.h
cp ZprimeJetsClass_MC.C ${gitdir}/ControlRegion/ZprimeJetsClass_MC.C
cp ZprimeJetsClass_MC.h ${gitdir}/ControlRegion/ZprimeJetsClass_MC.h
cp ZprimeJetsClass_MC_WJets.C ${gitdir}/ControlRegion/ZprimeJetsClass_MC_WJets.C
cp ZprimeJetsClass_MC_WJets.h ${gitdir}/ControlRegion/ZprimeJetsClass_MC_WJets.h
cp ZprimeJetsClass_MC_ZJets.C ${gitdir}/ControlRegion/ZprimeJetsClass_MC_ZJets.C
cp ZprimeJetsClass_MC_ZJets.h ${gitdir}/ControlRegion/ZprimeJetsClass_MC_ZJets.h

cd ../SignalRegion/
cp allsubmit.sh ${gitdir}/SignalRegion/allsubmit.sh
cp kfactors.root ${gitdir}/SignalRegion/kfactors.root
cp MakeCondorFiles.csh ${gitdir}/SignalRegion/MakeCondorFiles.csh
cp MakeCondorFiles_data.csh ${gitdir}/SignalRegion/MakeCondorFiles_data.csh
cp MakeCondorFiles_WZ.csh ${gitdir}/SignalRegion/MakeCondorFiles_WZ.csh
cp PU_Central.root ${gitdir}/SignalRegion/PU_Central.root
cp PU_minBiasDOWN.root ${gitdir}/SignalRegion/PU_minBiasDOWN.root
cp PU_minBiasUP.root ${gitdir}/SignalRegion/PU_minBiasUP.root
cp reset.sh ${gitdir}/SignalRegion/reset.sh
cp rootcom ${gitdir}/SignalRegion/rootcom
cp stack_plotterlog2.C ${gitdir}/SignalRegion/stack_plotterlog2.C
cp stack_plotterlog.C ${gitdir}/SignalRegion/stack_plotterlog.C
cp stack_plotterlogRebin.C ${gitdir}/SignalRegion/stack_plotterlogRebin.C
cp submitDY.sh ${gitdir}/SignalRegion/submitDY.sh
cp submitGH.sh ${gitdir}/SignalRegion/submitGH.sh
cp submit.sh ${gitdir}/SignalRegion/submit.sh
cp ZprimeJetsClass.C ${gitdir}/SignalRegion/ZprimeJetsClass.C
cp ZprimeJetsClass.h ${gitdir}/SignalRegion/ZprimeJetsClass.h
cp ZprimeJetsClass_MC.C ${gitdir}/SignalRegion/ZprimeJetsClass_MC.C
cp ZprimeJetsClass_MC.h ${gitdir}/SignalRegion/ZprimeJetsClass_MC.h
cp ZprimeJetsClass_MC_WJets.C ${gitdir}/SignalRegion/ZprimeJetsClass_MC_WJets.C
cp ZprimeJetsClass_MC_WJets.h ${gitdir}/SignalRegion/ZprimeJetsClass_MC_WJets.h
cp ZprimeJetsClass_MC_ZJets.C ${gitdir}/SignalRegion/ZprimeJetsClass_MC_ZJets.C
cp ZprimeJetsClass_MC_ZJets.h ${gitdir}/SignalRegion/ZprimeJetsClass_MC_ZJets.h