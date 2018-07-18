#!/bin/sh
#Author: E.Koenig
#Create Directories to put condor files in
if [ ! -d .output ]; then
    #Where executable and output files go
    mkdir .output
fi
if [ ! -d .filelist ]; then
    #Where PlotTool.C puts number of each type of sample file
    mkdir .filelist
fi
if [ ! -d .status ]; then
    #Where all condor output, log, and error files go
    mkdir .status
fi
fileNum=$( ls -f ${2}*.root | wc -l) #Total number of root files in directory
dirsplit=0
#Set step size for splitting up directories
for arg in "$@"; do
    prefix="split_"
    exec="analyze"
    outfile="post"
    if echo $arg | grep -q $exec; then
	#Assure executable file is in .output/
	if [ -f $arg ]; then
	    mv $arg .output/$arg
	fi
    elif echo $arg | grep -q $outfile; then
	#Remove hadd file so that PlotTool.C does not get confused
	if [ -f $arg ]; then
	    rm $arg
	    if echo $arg | grep -q "postMETdata_"; then
		if [ -f "postMETdata_final.root" ]; then
		    rm postMETdata_final.root
		fi
	    fi
	fi
    elif echo $arg | grep -q $prefix; then
	dirsplit=${arg#$prefix}
	echo $displit
	break
    fi
done
#If split_-1 is used program will set custom split for each directory so that there are nfile of files in each batch
if [ $dirsplit -eq -1 ]; then
    nfile=60
    dirsplit=$(echo $(($fileNum/$nfile)))
    #Dealing with some edge cases
    if [ $dirsplit -eq 0 ]; then
	dirsplit=1
    elif [ $dirsplit -eq 1 ]; then
	dirsplit=2
    fi
fi

#Writing Job file
export pwd=/cms/uhussain/CMSSW_8_0_18_patch1/src
cat>.output/Job_${6}.sh<<EOF
#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/uhussain/CMSSW_8_0_26_patch1/src
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
./${1} \${1} \${2} \${3} \${4} \${5} \${6} \${7}
EOF

chmod 775 .output/Job_${6}.sh

#Beginning to write condor_submit file
cat>.output/condor_${6}<<EOF
x509userproxy = /tmp/x509up_u4318
universe = vanilla
Executable = Job_${6}.sh
Notification         = never
WhenToTransferOutput = On_Exit
ShouldTransferFiles  = yes
Requirements = (TARGET.UidDomain == "hep.wisc.edu" && TARGET.HAS_CMS_HDFS) 
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv = true
request_memory       = 1992
request_disk         = 2048000
Transfer_Input_Files = ${1},../kfactors.root,../PU_Central.root
output               = ../.status/\$(Process)_${6}.out
error                = ../.status/\$(Process)_${6}.err
Log                  = ../.status/\$(Process)_${6}.log
EOF

#Get how many files are in each batch
binsize=$(($fileNum/$dirsplit))

#If directory ends in */0001/ or higher, files in directory won't start at 1.root
#shift is used to deal with that
shift=0
if ! echo ${2} | grep -q "*/0000/"; then
    chld=${2#*/000}
    chld=${chld%/}
    for ((i=0;$(($i < $chld));i++)); do
	parntdir=${2%*000${chld}/}000${i}/
	shift=$(($shift+$( ls -f $parntdir*.root | wc -l)))
    done
fi

#Determine range of files that will be used in each batch run
start=$shift
for ((i=1;$(($i <= $dirsplit));i++)); do
    end=$(($start+$binsize))
    if [ $i -eq $dirsplit ]; then
	end=$(($fileNum+$shift))
    fi
    argv=""
    j=0
    for arg in "$@"; do
	if echo $arg | grep -q "split_"; then
	    break #split_ option not necessary for actual condor job
	elif [[ $j -eq 2 ]]; then
	    output=${arg%.root} #${2} command line argument is always output file name 
	    output="${output}_$i.root"
	    argv="$argv $output"
	elif ! echo $arg | grep -q "analyze";then
	    argv="$argv $arg" #omit executable argument since it is already in the Job file
	fi
	j=$(($j+1))
    done
    
    if echo $1 | grep -q "analyze" ; then
	echo Running $i $start $end ${6}
	#./debug.sh ${1} $argv $start $end

	#Append argument lines to condor_submit file adding the start and end file numbers for this batch
	cat>>.output/condor_${6}<<EOF
Arguments = $argv $start $end
Queue
EOF
    fi

    start=$end
done

#Move into .output/ and run newly made condor_submit file
cd .output
condor_submit condor_${6}
			   
