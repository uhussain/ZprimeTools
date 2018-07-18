from sys import *
from ROOT import TH2D
from ROOT import gStyle
from ROOT import TLatex
from ROOT import TGaxis
from ROOT import gROOT
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import TColor
from ROOT import TGraph
import ROOT
from array import array
import operator
gROOT.SetBatch(True)

if len(argv) == 1:
    print "usage: python CLplotter.py limits_shape_mchi*085.txt"
    print "replace 085 with and cut you wish to plot"
    exit()


#tdrstyle.setTDRStyle()
gStyle.SetTitleAlign(23)

alpha = 1.0;                                                                                                                                                      
stops = array('d',[0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000])                                                                             
red =   array('d',[  0./255.,  50./255.,  130./255.,  180./255., 200./255.,  215./255., 230./255., 240./255., 255./255.])                                             
green = array('d',[  0./255.,  50./255.,  130./255.,  180./255., 200./255.,  215./255., 230./255., 240./255., 255./255.])                                                
blue =  array('d',[ 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255., 255./255.])                                                  
Idx = TColor.CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);  

#gStyle.SetPalette(Idx)

H_ref = 600;
W_ref = 600;
W = W_ref
H  = H_ref
iPeriod = 4
# references for T, B, L, R                                                                       
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref
B_ratio = 0.1*H_ref
T_ratio = 0.1*H_ref
B_ratio_label = 0.3*H_ref

x1_l = 0.93
y1_l = 0.90
dx_l = 0.20
dy_l = 0.20
x0_l = x1_l-dx_l
y0_l = y1_l-dy_l

signal_mx=[10,50,100]
signal_mv=[100,200,1000,1500,1800,2000,2500,3500]

order_mx = []
order_mv = []

def getScale(masses):
    minnes = {}
    maxes = {}
    maxvalmass=""
    minvalmass=""
    for keys in masses:
        maxes[keys]=max(masses[keys])
        minnes[keys]=min(masses[keys])
    for keys in maxes:
        maxvalmass = max(maxes.iteritems(), key=operator.itemgetter(1))[0]
        minvalmass = min(minnes.iteritems(), key=operator.itemgetter(1))[0]
    relmax=[]
    scale=0.8
    for keys in masses:
        if float(keys) <= 40*0.45:
            if max(masses[keys]) >= maxes[maxvalmass]*scale*0.6:
                relmax.append(max(masses[keys]))
    if len(relmax)>0:
        scale=max(relmax)/(maxes[maxvalmass]*0.6)

    return scale, maxes[maxvalmass], minnes[minvalmass]

class Limit():
    def __init__(self):
        self.mx = array('d')
        self.mv = array('d')
        self.obslim = array('d')
        self.explim = [array('d'),array('d'),array('d'),array('d'),array('d')]

def GetData(filename,ExpectedLimit):
    percent = ["2.5%:","16.0%:","50.0%:","84.0%:","97.5%:"]
    percent2= ['"exp-2":','"exp-1":','"exp0":','"exp+1":','"exp+2":']
    if "Mx" in filename:
        with open(filename,"r") as file:
            text = file.readlines()
            for line in text:
                word = line.split()
                if len(word) > 0:
                    if word[0] == "Observed":
                        ExpectedLimit.obslim.append(float(word[len(word)-1]))
                    elif word[0] == "Expected":
                        for i in range(len(percent)):
                            if word[1] == percent[i]:
                                ExpectedLimit.explim[i].append(float(word[len(word)-1]))
    elif "limits_shape_mchi_" in filename:
        with open(filename,"r") as file:
            text = file.readlines()
            mvBool=False
            for line in text:
                try:
                    if line.split(":")[1] == " {\n":
                        mv = int(float(line.split(":")[0].replace('"','')))
                        if mv == int(ExpectedLimit.mv[0]):
                            if not mvBool:mvBool=True
                            elif mvBool:mvBool=False
                except IndexError:
                    pass
                if mvBool:
                    data=line.split()
                    if "obs" in data[0]:
                        ExpectedLimit.obslim.append(float(data[1].replace(",","")))
                    elif "exp" in data[0]:
                        for i in range(len(percent2)):
                            if percent2[i] in data[0]:
                                ExpectedLimit.explim[i].append(float(data[1].replace(",","")))

def Make2DPlot(Points):
    canvas =TCanvas("c","",800,800)
    canvas.SetMargin(0.15,0.15,0.15,0.08)

    xBins=len(signal_mv)
    yBins=len(signal_mx)
    
    xmax=order_mv[-1]
    xmin=order_mv[0]
    ymax=order_mx[-1]
    ymin=order_mx[0]

    mx_bin,mv_bin={},{}
    for i in range(len(signal_mx)):mx_bin[order_mx[i]]=i+1
    for i in range(len(signal_mv)):mv_bin[order_mv[i]]=i+1
        
    
    graph = TH2D("Expected Limits","",xBins,0,xBins,yBins,0,yBins)
    for mv in signal_mv:
        for mx in signal_mx:
            key = "Mx"+str(mx)+"_Mv"+str(mv)
            #print Points[key].mv[0],Points[key].mx[0],Points[key].explim[2][0]
            graph.SetBinContent(mv_bin[Points[key].mv[0]],mx_bin[Points[key].mx[0]],Points[key].explim[2][0])

    graph.Draw("COLZ")
    graph.SetStats(0)
    
    graph.GetXaxis().SetTitleOffset(999)
    graph.GetXaxis().SetLabelOffset(999)
    graph.GetXaxis().SetTickLength(0)

    graph.GetYaxis().SetTitleOffset(999)
    graph.GetYaxis().SetLabelOffset(999)
    graph.GetYaxis().SetTickLength(0)

    graph.GetZaxis().SetTitle("95% CL limit on #sigma/#sigma_{theory}")
    graph.GetZaxis().SetTitleOffset(1.2)
        
    xaxis = TGaxis(0,0,xBins,0,0,xBins,xBins)
    xaxis.SetTitle("m_{med} [GeV]")
    xaxis.SetLabelFont(42);
    xaxis.SetLabelSize(0);
    xaxis.SetTitleFont(42);
    xaxis.SetTitleSize(0.05);
    xaxis.SetTitleOffset(0.9);
    xaxis.Draw("SAME")

    label=TLatex()
    label.SetTextSize(0.04);
    label.SetTextFont(42)
    
    for i in range(xBins):
        size=float(len(str(order_mv[i])))
        label.DrawLatex(i+0.5/size,-0.15,str(order_mv[i]))

    yaxis = TGaxis(0,0,0,yBins,0,yBins,yBins)
    yaxis.SetTitle("m_{#chi} [GeV]")
    yaxis.SetLabelFont(42);
    yaxis.SetLabelSize(0);
    yaxis.SetTitleFont(42);
    yaxis.SetTitleSize(0.05);
    yaxis.SetTitleOffset(1.2);
    yaxis.Draw("SAME")

    for i in range(yBins):
        size=len(str(order_mx[i]))
        label.DrawLatex(-0.25*size,i+0.5,str(order_mx[i]))


    label.SetTextColor(ROOT.kOrange+10)
    label.SetTextSize(0.028)
    for i in range(xBins):
        for j in range(yBins):
            key="Mx"+str(order_mx[j])+"_Mv"+str(order_mv[i])
            label.DrawLatex(i+0.1,j+0.5,"#bf{"+str(round(Points[key].explim[2][0],3))+"}")
    

    label.SetTextColor(1)
    label.SetTextSize(0.03)
    label.DrawLatex(3,yBins+0.01,"#sqrt{s} = 13 TeV, 35.9 fb^{-1}")
    label.DrawLatex(0.1,yBins+0.01,"#bf{CMS} : #it{Preliminary}")
            
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs("expectedevents2D.png")
    canvas.SaveAs("expectedevents2D.pdf")
    
    
def Make1DPlot(Points):
    color = [ROOT.kBlue,ROOT.kRed,ROOT.kMagenta,ROOT.kCyan,ROOT.kGreen,ROOT.kBlack,ROOT.kYellow-6,ROOT.kRed-8]
    canvas =TCanvas("c","Expected Limits",50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.SetRightMargin(R/W)
    canvas.SetLeftMargin(L/W)
    
    allgraph = []
    legend = TLegend(0.2,0.6,0.45,0.787173,"")
    legend.SetTextSize(0.02)
    legend.SetFillColor(0)
    ##legend.SetNColumns(2)
    i=0
    masses={}
    maxX=0
    minX=0
    shift=0
    for keys in Points:
        maxX=Points[keys].mv[len(Points[keys].mv)-1]
        minX=Points[keys].mv[0]
        graph =TGraph(len(Points[keys].mv),Points[keys].mv,Points[keys].explim[2])
        for j in range(len(Points[keys].mv)):
            if not str(Points[keys].mv[j]) in masses:
                masses[str(Points[keys].mv[j])]=[]
            masses[str(Points[keys].mv[j])].append(Points[keys].explim[2][j])
        graph.SetTitle("")
        graph.SetName(keys)
        try:
            graph.SetLineColor(color[i+shift])
        except IndexError:
            if shift == 0:shift=5
            elif shift == 5:shift=-5
        graph.SetLineWidth(2)
        
        legend.AddEntry(graph,keys,"l")
        allgraph.append(graph)
        i+=1
    scale, maxvalue, minvalue = getScale(masses)
    for i in range(len(allgraph)):
        if i == 0:
            allgraph[i].Draw()
            allgraph[i].GetXaxis().SetRangeUser(minX,maxX)
            allgraph[i].GetYaxis().SetRangeUser(minvalue*(10**-1),maxvalue*(10**scale))
            allgraph[i].GetXaxis().SetTitle("m_{#chi} GeV")
            allgraph[i].GetYaxis().SetTitle("95% CL limit on #sigma/#sigma_{theor}")
            allgraph[i].GetXaxis().SetTitleSize(0.05)
            allgraph[i].GetYaxis().SetTitleSize(0.05)
            allgraph[i].GetXaxis().SetTitleOffset(0.92)
            allgraph[i].GetYaxis().SetTitleOffset(0.92)
            allgraph[i].GetXaxis().SetLabelSize(0.03)
            allgraph[i].GetYaxis().SetLabelSize(0.03)
        else:
            allgraph[i].Draw("SAME")
    canvas.SetLogy()

    texS = TLatex(0.20,0.837173,"#sqrt{s} = 13 TeV, 1.885 fb^{-1}");
    texS.SetNDC();
    texS.SetTextFont(42);
    texS.SetTextSize(0.040);
    texS.Draw();
    texS1 = TLatex(0.12092,0.907173,"#bf{CMS} : #it{Preliminary}");
    texS1.SetNDC();
    texS1.SetTextFont(42);
    texS1.SetTextSize(0.040);
    texS1.Draw();

    legend.Draw()
    #CMS_lumi.CMS_lumi(canvas,iPeriod,iPos)
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs("expectedevents.png")

def main(filenames):

    Points = {}
    if "Mx" in filenames[0]:
        for name in filenames:
            temp= name.replace(".txt","")
            temp= temp.replace("Mx","")
            temp= temp.replace("Mv","")
            temp= temp.replace(".txt","")
            temp= temp.split("_")
            key="Mx"+temp[0]+"_Mv"+temp[1]
            if not key in Points:
                Points[key]=Limit()
            Points[key].mx.append(float(temp[0]))
            Points[key].mv.append(float(temp[1]))
            GetData(name,Points[key])
    elif "limits_shape_mchi_" in filenames[0]:
        for name in filenames:
            temp=name.split("_")
            temp=[temp[3]]
            with open(name,"r") as file:
                text=file.readlines()
                for line in text:
                    try:
                        if line.split(":")[1] == " {\n":
                            mv = str(int(float(line.split(":")[0].replace('"',''))))
                            key="Mx"+temp[0]+"_Mv"+mv
                            if not key in Points:
                                Points[key]=Limit()
                            Points[key].mx.append(float(temp[0]))
                            Points[key].mv.append(float(mv))
                            GetData(name,Points[key])
                    except IndexError:
                        pass
                
        
   
    Make2DPlot(Points)


    
filenames = []
orginize = {}
if "Mx" in argv[1]:
    for i in range(1,len(argv)):
        temp = argv[i].split("_")
        if len(temp) > 1:
            try:
                float(temp[0].replace("Mx",""))
                if not temp[2] in orginize:
                    orginize[temp[2]] = {}
                if not temp[0].replace("Mx","") in orginize[temp[2]]:
                    orginize[temp[2]][temp[0].replace("Mx","")]={}
                orginize[temp[2]][temp[0].replace("Mx","")][temp[1].replace("Mv","")]=argv[i]
            except ValueError:
                pass
    for catType in orginize:
        order_mx = []
        order_mv = []
        for mchi in orginize[catType]:
            order_mx.append(int(mchi))
            if len(order_mv) == 0:
                for mv in orginize[catType][mchi]:
                    order_mv.append(int(mv))
            order_mx.sort()
            order_mv.sort()
        for mchi in order_mx:
            for mv in order_mv:
                filenames.append(orginize[catType][str(mchi)][str(mv)])

elif "limits_shape_mchi" in argv[1]:
    for i in range(1,len(argv)):
        temp = argv[i].split("_")
        if not temp[4] in orginize:
            orginize[temp[4]]={}
        if not temp[3] in orginize[temp[4]]:
            orginize[temp[4]][temp[3]]={}
        with open(argv[i],"r") as file:
            text=file.readlines()
            for line in text:
                try:
                    if line.split(":")[1] == " {\n":
                        mv = str(int(float(line.split(":")[0].replace('"',''))))
                        orginize[temp[4]][temp[3]][mv]=argv[i]
                except IndexError:
                    pass
    
    for catType in orginize:
        order_mx = []
        order_mv = []
        for mchi in orginize[catType]:
            order_mx.append(int(mchi))
            if len(order_mv) == 0:
                for mv in orginize[catType][mchi]:
                    order_mv.append(int(float(mv)))
            order_mx.sort()
            order_mv.sort()
        for mchi in order_mx:
            for mv in order_mv:
                filenames.append(orginize[catType][str(mchi)][str(mv)])
                break

main(filenames)
    
