#!/usr/bin/python

from sys import argv
from os import path, system, mkdir, listdir, rename, remove, chdir

#Create Directories to put condor files in
#Where executable and output files go
if not path.isdir(".output/"): mkdir(".output/")
#Where PlotTool.C puts number of each type of sample file
if not path.isdir(".filelist/"): mkdir(".filelist/")
#Where all condor output, log, and error files go
if not path.isdir(".status/"): mkdir(".status/")

rootFiles = [fn.replace(".root","") for fn in listdir(argv[2]) if fn.endswith(".root")];
#Getting file number instead of the entire filename
dataset="null"
for i in range(len(rootFiles[0]) - 1, -1,-1): #Loop through filename backwards
    if rootFiles[0][i] == "_" or rootFiles[0][i] == "-":
        dataset = rootFiles[0][:i+1]
        break
for i in range(len(rootFiles)):rootFiles[i] = rootFiles[i].replace(dataset,"")
rootFiles.sort(key=int)

#Setting Command Line arguments
#Assure executable file is in .output/
executable = argv[1]
if path.isfile(executable): rename(executable,".output/"+executable)

directory = argv[2]

#Remove hadd file so that PlotTool.C does not get confused
output = argv[3]
maxEvents = argv[4]
reportEvery = argv[5]
label = argv[6]
#Remove any old condor files
# system("rm .status/*"+label+"* >/dev/null 2>&1" )
# system("rm .output/*"+label+"* >/dev/null 2>&1" )
# system("rm .output/"+output.replace(".root","")+"* >/dev/null 2>&1")

if (len(argv) == 8):nBatches=int(argv[7].replace("split_",""))
else:nBatches = 1

#If split_-1 is used program will set custom split for each directory so that there are nfile of files in each batch
if nBatches == -1:
    nfile_per_batch = 60
    nBatches = len(rootFiles)/nfile_per_batch
    #Dealing with some edge cases
    if nBatches == 0: nBatches = 1
    elif nBatches == 1: nBatches = 2

#Writing Job File
with open(".output/Job_"+label+".sh","w") as jobfile:
    jobfile.write("#!/bin/sh\n"
                + "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
                + "cd /cms/uhussain/CMSSW_8_0_26_patch1/src\n"
                + "cmsenv\n"
                + "cd ${_CONDOR_SCRATCH_DIR}\n"
                + "./"+argv[1]+" ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8}\n")

system("chmod 775 .output/Job_"+label+".sh")

files_to_transfer=argv[1]+",../kfactors.root,../PU_Central.root"

#If NLO EWK files in directory, transfer them
if path.isfile("WJets_NLO_EWK.root"): files_to_transfer=argv[1]+",../kfactors.root,../PU_Central.root,../WJets_NLO_EWK.root,../ZJets_NLO_EWK.root"
if path.isfile("egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"): files_to_transfer=argv[1]+",../kfactors.root,../PU_Central.root,../egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root,../egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root"

#Beginning to write condor_submit file
with open(".output/condor_"+label,"w") as condor:
    condor.write("x509userproxy = /tmp/x509up_u23216\n"
                + "universe = vanilla\n"
                + " Executable = Job_"+label+".sh\n"
                + " Notification         = never\n"
                + " WhenToTransferOutput = On_Exit\n"
                + " ShouldTransferFiles  = yes\n"
                + " Requirements = (TARGET.UidDomain == \"hep.wisc.edu\" && TARGET.HAS_CMS_HDFS)\n"
                + " on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))\n"
                + " +IsFastQueueJob      = True\n"
                + " getenv = true\n"
                + " request_memory       = 1992\n"
                + " request_disk         = 2048000\n"
                + " Transfer_Input_Files = "+files_to_transfer+"\n"
                + " output               = ../.status/\\$(Process)_"+label+".out\n"
                + " error                = ../.status/\\$(Process)_"+label+".err\n"
                + " Log                  = ../.status/\\$(Process)_"+label+".log\n")

    #Get how many files are in each batch
    binsize = len(rootFiles)/nBatches
    print binsize
    for i in range(nBatches):
        if nBatches == 1: fileRange = "-1"
        else:
            #0-100/200-250/300-300
            fileRange = rootFiles[0]
            if len(rootFiles)/binsize == 1:
                binsize = len(rootFiles)
            for j in range(1,binsize):
                if (int(rootFiles[j]) - int(rootFiles[j-1]) != 1):
                    fileRange += "-"+rootFiles[j-1]+"/"+rootFiles[j]

            fileRange += "-"+rootFiles[binsize - 1]
            for j in range(binsize - 1, -1, -1):rootFiles.pop(j) #Remove files already accounted for

        print "Running",i+1,fileRange,label

        #Append argument lines to condor_submit file adding the file range for this batch
        condor.write("Arguments = "+directory+" "+output.replace(".root","_"+str(i)+".root")+" "+maxEvents+" "+reportEvery+" "+fileRange+"\n"
                     + "Queue\n")


#Move into .output/ and run newly made condor_submit file
chdir(".output/")
system("condor_submit condor_"+label)
    
    