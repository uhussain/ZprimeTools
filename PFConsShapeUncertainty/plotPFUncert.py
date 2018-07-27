import ROOT
import glob
import os

def plotPFUncert(SystFile,Sample):
    path='/nfs_scratch/uhussain/MonoZprimeJet_postanalyzer_jobsubmission/Uncertainties/AllSystUncertFiles/'
    selection=SystFile.strip(path)
    print SystFile
    print selection
    #selection="Basic." # Basic corresponds to histos; 8(Norm),12(Up),16(Down)
    #selection="Pt123Frac08." # Pt123Fraction08Cut corresponds to histos; 9(Norm),13(Up),17(Down)
    #selection="Pt123Frac085." # Pt123Fraction08Cut corresponds to histos; 10(Norm),14(Up),18(Down)
    #selection="Pt123Frac09." # Pt123Fraction08Cut corresponds to histos; 11(Norm),15(Up),19(Down)
    print selection
    print Sample
    f = ROOT.TFile(SystFile)
    hNorm=ROOT.TH1F(f.Get(Sample+""))
    hUp=ROOT.TH1F(f.Get(Sample+"_TrkUncUp"))
    hDown=ROOT.TH1F(f.Get(Sample+"_TrkUncDown"))
    
    print "TrkUncUpYield: ",hUp.Integral()
    print "NormYield: ",hNorm.Integral()
    print "TrkUncDoYield: ",hDown.Integral()
    
    c1=ROOT.TCanvas()

    nbins = hNorm.GetNbinsX();  
    #Plotting Normalized histos such that Norm = 1, JESUp,JESDown are relative to Norm
    hNorm_Plot = ROOT.TH1F("hNorm_Plot","hNorm_Plot",50,0.0,1.0);
    hUp_Plot = ROOT.TH1F("hUp_Plot","hUp_Plot",50,0.0,1.0);
    hDown_Plot = ROOT.TH1F("hDown_Plot","hDown_Plot",50,0.0,1.0);
    #hNorm.SetDirectory(0)
    #hUp.SetDirectory(0)
    #hDown.SetDirectory(0)
    for i in range(nbins):
        NormContent=hNorm.GetBinContent(i)
        UpContent=hUp.GetBinContent(i)
        DownContent=hDown.GetBinContent(i)
        #Error
        if NormContent!=0:
            NormRatio=(NormContent/NormContent)
            print "NormRatio: ",NormRatio
            hNorm_Plot.SetBinContent(i,NormRatio)
        else: 
            hNorm_Plot.SetBinContent(i,0.0)
        if NormContent!=0:
            UpRatio=(UpContent/NormContent)
            hUp_Plot.SetBinContent(i,UpRatio) 
            print "UpRatio: ",UpRatio
        else: 
            hUp_Plot.SetBinContent(i,UpContent)
        if NormContent!=0:
            DownRatio=(DownContent/NormContent)
            hDown_Plot.SetBinContent(i,DownRatio)
            print "DownRatio: ",DownRatio
        else: 
            hDown_Plot.SetBinContent(i,DownContent)
    for j in range(hNorm_Plot.GetNbinsX()):
        print "bin",j,":",hUp_Plot.GetBinContent(j),",",hNorm_Plot.GetBinContent(j),",",hDown_Plot.GetBinContent(j)
    #Rebin
    #hNorm_Plot.Rebin()#merges 2 bins in one in hNorm_Plot; previous contents of hNorm are lost
    #hUp_Plot.Rebin()
    #hDown_Plot.Rebin()
    hNorm_Plot.SetLineColor(ROOT.kBlue)
    hUp_Plot.SetLineColor(ROOT.kRed)
    hDown_Plot.SetLineColor(ROOT.kBlack)
    hNorm_Plot.GetYaxis().SetTitle("TrkUnc")
    hUp_Plot.GetYaxis().SetTitle("TrkUnc")
    hDown_Plot.GetYaxis().SetTitle("TrkUnc")
    RangeList=[]
    RangeList.append(hUp_Plot.GetBinContent(hUp_Plot.GetMaximumBin()))
    RangeList.append(hNorm_Plot.GetBinContent(hNorm_Plot.GetMaximumBin()))
    RangeList.append(hDown_Plot.GetBinContent(hDown_Plot.GetMaximumBin()))
    print RangeList
    ymax=max(RangeList)
    ymax=ymax*1.05
    hUp_Plot.GetYaxis().SetRangeUser(0.95,ymax) 
    hUp_Plot.Draw("HIST")
    hNorm_Plot.Draw("HIST same")
    hDown_Plot.Draw("HIST same")
    legend=ROOT.TLegend(.65,.70,.85,.90)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.02)
    hUpIntegral=str(round(hUp.Integral(),1))
    hNormIntegral=str(round(hNorm.Integral(),1))
    hDownIntegral=str(round(hDown.Integral(),1))
    legend.AddEntry(hUp_Plot,"TrkUncUp: "+hUpIntegral)
    legend.AddEntry(hNorm_Plot,"Norm: "+hNormIntegral)
    legend.AddEntry(hDown_Plot,"TrkUncDown: "+hDownIntegral)
    legend.Draw()
    selection=selection[:-1]
    file_path="/afs/hep.wisc.edu/home/uhussain/public_html/PTFracUncertPlots/PFConsUncert/TrkUnc/"
    #print file_path
    directory=os.path.join(os.path.dirname(file_path),selection)
    if not os.path.exists(directory):
        os.mkdir(directory,0755)
        print directory
    c1.SaveAs(directory+"/"+Sample+".pdf")
    f.Close()

#JES Samples
#Samples=["data_obs","DM_1GeV","DM_5GeV","DM_10GeV","DM_20GeV","DM_50GeV","DM_100GeV","WJets","ZJets","GJets","DYJets","TTJets","DiBoson","QCD"]
#PDF Samples
#Samples=["DM_Mx10_Mv1000"]
Samples=["data_obs","DM_Mx100_Mv1000","WJets","ZJets","GJets","DYJets","TTJets","QCD"]
#Samples=["PF123PtFraction","h_j1PFCons1PtFraction","h_j1PFCons2PtFraction","h_j1PFCons3PtFraction"]
#Samples=["PF123PtFraction"]
SystFiles = glob.glob("/nfs_scratch/uhussain/MonoZprimeJet_postanalyzer_jobsubmission/Uncertainties/AllSystUncertFiles/*.root")
for f in SystFiles:
   for s in Samples:
        plotPFUncert(f,s)
