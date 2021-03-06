DATE=May2018
DATE_2=29May

echo "Do the Data samples"

./rootcom ZprimeJetsClass analyzedata

./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/gomber/MonoZprime_Ntuples_${DATE_2}_SR/MET/crab_dataset6/180529_*/0000/ postMETdata_0.root -1 10000 MET_0

./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/gomber/MonoZprime_Ntuples_${DATE_2}_SR/MET/crab_dataset6/180529_*/0001/ postMETdata_1.root -1 10000 MET_1

./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/gomber/MonoZprime_Ntuples_${DATE_2}_SR/MET/crab_dataset6/180529_*/0002/ postMETdata_2.root -1 10000 MET_2


echo "Do the MC samples"

./rootcom ZprimeJetsClass_MC_inclusive analyze1

./MakeCondorFiles.csh analyze1 /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL/180528_*/0000/ postDY_MLM_0.root -1 10000 DYMLM_0 PU_Central.root

./rootcom ZprimeJetsClass_MC analyze

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT100to200/180528_*/0000/ postDY100to200.root -1 10000 DY100 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT200to400/180528_*/0000/ postDY200to400.root -1 10000 DY200 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT400to600/180528_*/0000/ postDY400to600.root -1 10000 DY400 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT600to800/180528_*/0000/ postDY600to800.root -1 10000 DY600 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT800to1200/180528_*/0000/ postDY800to1200.root -1 10000 DY800 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT1200to2500/180528_*/0000/ postDY1200to2500.root -1 10000 DY1200 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJets_HT2500toInf/180528_*/0000/ postDY2500toInf.root -1 10000 DY2500 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT100to200/180528_*/0000/ postQCD100to200_0.root -1 10000 QCD100_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT200to300/180528_*/0000/ postQCD200to300_0.root -1 10000 QCD200_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT300to500/180528_*/0000/ postQCD300to500_0.root -1 10000 QCD300_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT500to700/180528_*/0000/ postQCD500to700_0.root -1 10000 QCD500_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT700to1000/180528_*/0000/ postQCD700to1000_0.root -1 10000 QCD700_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500/180528_*/0000/ postQCD1000to1500_0.root -1 10000 QCD1000_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1500to2000/180528_*/0000/ postQCD1500to2000_0.root -1 10000 QCD1500_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT2000toInf/180528_*/0000/ postQCD2000toInf_0.root -1 10000 QCD2000_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets/180528_*/0000/ postTTJets_MLM.root -1 10000 TTJets_MLM PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/180528_*/0000/ postGJets40to100.root -1 10000 GJets40_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/180528_*/0000/ postGJets100to200.root -1 10000 GJets100_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/180528_*/0000/ postGJets200to400.root -1 10000 GJets200_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/180528_*/0000/ postGJets400to600.root -1 10000 GJets400_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/180528_*/0000/ postGJets600toInf.root -1 10000 GJets600_0 PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/180528_*/0000/ postWW.root -1 10000 WW PU_Central.root 

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/180528_*/0000/ postWZ.root -1 10000 WZ PU_Central.root

./MakeCondorFiles.csh analyze /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/180528_*/0000/ postZZ.root -1 10000 ZZ PU_Central.root


#WJets
./rootcom ZprimeJetsClass_MC_WJets_inclusive analyzeWJets1

./MakeCondorFiles_WZ.csh analyzeWJets1 /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJets/180528_*/0000/ postWJets_MLM_0.root -1 10000 W_0 kfactors.root PU_Central.root

./rootcom ZprimeJetsClass_MC_WJets analyzeWJets

./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W1Jets/180528_*/0000/ postW100to200_0.root -1 10000 W100_0 kfactors.root PU_Central.root 
                                                                                                                                                                                                                  
./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W2Jets/180528_*/0000/ postW200to400_0.root -1 10000 W200_0 kfactors.root PU_Central.root

./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W3Jets/180528_*/0000/ postW400to600_0.root -1 10000 W400_0 kfactors.root PU_Central.root
                                                                                                                                                                                                                   
./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W4Jets/180528_*/0000/ postW600to800_0.root -1 10000 W600_0 kfactors.root PU_Central.root

./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W5Jets/180528_*/0000/ postW800to1200_0.root -1 10000 W800_0 kfactors.root PU_Central.root
                                                                                                                                                                                                                      
./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W6Jets/180528_*/0000/ postW1200to2500_0.root -1 10000 W1200_0 kfactors.root PU_Central.root
                                                                                                                                                                                                                      
./MakeCondorFiles_WZ.csh analyzeWJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_W7Jets/180528_*/0000/ postW2500toInf_0.root -1 10000 W2500_0 kfactors.root PU_Central.root

#ZJets

./rootcom ZprimeJetsClass_MC_ZJets analyzeZJets

./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/crab_Zvv2500toInf/180528_*/0000/ postZ2500toInf_0.root -1 10000 Z2500_0 kfactors.root PU_Central.root  
                                                                                                                                                                                                      
./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/crab_Zvv1200to2500/180528_*/0000/ postZ1200to2500_0.root -1 10000 Z1200_0 kfactors.root PU_Central.root
                                                                                                                                                                                                      
./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/crab_Zvv800to1200/180528_*/0000/ postZ800to1200_0.root -1 10000 Z800_0 kfactors.root PU_Central.root

./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-600To800_13TeV-madgraph/crab_Zvv600to800/180528_*/0000/ postZ600to800_0.root -1 10000 Z600_0 kfactors.root PU_Central.root
                                                                                                                                                                                                
./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-400To600_13TeV-madgraph/crab_Zvv400to600/180528_*/0000/ postZ400to600_0.root -1 10000 Z400_0 kfactors.root PU_Central.root
                                                                                                                                                                                                
./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-200To400_13TeV-madgraph/crab_Zvv200to400/180528_*/0000/ postZ200to400_0.root -1 10000 Z200_0 kfactors.root PU_Central.root
                                                                                                                                                                                                
./MakeCondorFiles_WZ.csh analyzeZJets /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/ZJetsToNuNu_HT-100To200_13TeV-madgraph/crab_Zvv100to200/180528_*/0000/ postZ100to200_0.root -1 10000 Z100_0 kfactors.root PU_Central.root
