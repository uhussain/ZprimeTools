#!/bin/bash
category=Pt123Fraction085
rm saveplot
rm ../../../PFConsSyst/Systematics_${category}.root
./rootcom saveplot saveplot
./saveplot TrackerPtFraction_10 TrackerPtFraction_14 TrackerPtFraction_18 EcalPtFraction_14 EcalPtFraction_18 HcalPtFraction_14 HcalPtFraction_18
cd ../../../CMSSW_8_0_26_patch2/src/${category}/
root -l -b savepdf.C
