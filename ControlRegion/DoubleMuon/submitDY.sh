DATE=Oct24

echo "Do the DY MC samples"

./rootcom ZprimeJetsClass_MC_ZJets analyzeDY

./MakeCondorFiles.csh analyzeDY /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT100to200/171024_*/0000/ postDY100to200.root -1 10000 DY100 kfactors.root 

