DATE=Aug10

echo "Do the GH Data samples"

./rootcom ZprimeJetsClass analyzedata

./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset4/170908_*/0000/ postMETdata_8.root -1 10000 MET_8

./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset4/170908_*/0001/ postMETdata_9.root -1 10000 MET_9

./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset4/170908_*/0002/ postMETdata_10.root -1 10000 MET_10

#./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset5/170908_*/0000/ postMETdata_11.root -1 10000 MET_11
#
#./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset5/170908_*/0001/ postMETdata_12.root -1 10000 MET_12
#
#./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset5/170908_*/0002/ postMETdata_13.root -1 10000 MET_13
#
#./MakeCondorFiles_data.csh analyzedata /hdfs/store/user/uhussain/Zprime_Ntuples_${DATE}/MET/crab_dataset5/170908_*/0003/ postMETdata_14.root -1 10000 MET_14

