DATE=Mar15

echo "Do the MC samples"

./rootcom ZprimeJetsClass_MC analyze


./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT100to200/170309_*/0000/ postDY100to200.root -1 10000 DY100

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT200to400/170309_*/0000/ postDY200to400.root -1 10000 DY200

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT400to600/170309_*/0000/ postDY400to600.root -1 10000 DY400

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT600to800/170309_*/0000/ postDY600to800.root -1 10000 DY600

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT800to1200/170309_*/0000/ postDY800to1200.root -1 10000 DY800

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT1200to2500/170309_*/0000/ postDY1200to2500.root -1 10000 DY1200

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT2500toInf/170309_*/0000/ postDY2500toInf.root -1 10000 DY2500

