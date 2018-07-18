DATE=July10

echo "Do the Signal samples"

./rootcom ZprimeJetsClass_MC_signal analyze

Mx="10 50 100"
Mv="100 200 500 1000 1500 1800 2000 2500 3500"

for dm in $Mx;do
  for Med in $Mv;do
    dir="Mx${dm}_Mv${Med}_MINIAOD"
    echo " "
    echo "----------------------"   ${dir}  "-------------------------"
    ./SubmitCondor_signal analyze /hdfs/store/user/uhussain/${dir}-ZprimeSignalJobs_${DATE}/ postSignal_Mx${dm}_Mv${Med}.root -1 1000 signal_Mx${dm}_Mv${Med} split_1
  done
done


