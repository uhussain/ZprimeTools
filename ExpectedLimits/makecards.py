from sys import *
from ROOT import *
import shutil

cat=""
mchi=""

datamcEvents = {}
SampleOrder = ["DM","ZvvJets","WJets","DiBoson","GJets","TTJets","DYJets","QCD"]
temptext=[]
signal_mx=["10","50","100"]
signal_mv=["100","200","1000","1500","1800","2000","2500","3500"]
category =["08","085","09"]
for cat in category:
    for mx in signal_mx:
        for mv in signal_mv:
            with open("cardtemplate_shape.txt","r") as template:
                temptext=template.readlines()
                with open("zprimeMx"+mx+"_Mv"+mv+"_"+cat+"_shape.txt","w") as newshape:
                    for line2 in temptext:
                        word2 = line2.split()
                        if len(word2)>0:
                            if word2[0] == "shapes":
                                word2[3]="Systematics_"+cat+".root"
                                arg=""
                                for i in range(len(word2)):
                                    arg+=word2[i]+" "
                                newshape.write(arg+"\n")

                            elif word2[0] == "process":
                                if word2[1] == "DM":
                                    word2[1]="    Mx"+mx+"_Mv"+mv
                                    s=20-len(word2[1])
                                    word2[1]+=" "*s
                                    arg=""
                                    for i in range(len(word2)):
                                        s=13-len(word2[i])
                                        if (i == 1):arg+=word2[i]
                                        else:arg+=word2[i]+" "*s
                                    newshape.write(arg+"\n")
                                else:
                                    newshape.write(line2)
                            else:
                                newshape.write(line2)
       
