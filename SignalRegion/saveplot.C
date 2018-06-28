void plotter(const char * variable,std::string name,std::string varname,std::string cat,std::string mchi)
{
  double lumi_1 = 3723.664;
  double lumi_2 = 35900.;
  double lumi_3 = 1885.122;

  std::vector<TH1F*> hs_histo = {};

  std::cout << variable <<std::endl;
  std::ifstream file("postMETdata_final.root");
  if (file)
    {
      file.close();
    }
  else
    {
      system("hadd -f postMETdata_final.root postMETdata_{0..2}.root");
    }

  TCanvas *c = new TCanvas("c", "canvas",800,800);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  //c->SetLeftMargin(0.15);
  //c->SetLogy();
  //c->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.25,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetLogy();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  //opening the data file and adding "h_dileptonM_8" histogram
  TFile *f_datafile_0 = new TFile("postMETdata_final.root");
  TH1F *histo_j1EtaWidth_data_0 = (TH1F*)f_datafile_0->Get(variable);
  histo_j1EtaWidth_data_0->SetName(((std::string("data_obs")+varname)).c_str());

  double integraldata = histo_j1EtaWidth_data_0->Integral();
  std::cout<<"integral of Data here:"<<integraldata<<std::endl;
  
  if (varname.find("jes") == std::string::npos)
    {
      hs_histo.push_back(histo_j1EtaWidth_data_0);
    }

  //opening background WJets Sample file
  std::vector<const char *> WJets_FileNames = {"postWJets_MLM_0.root","postW100to200_0.root","postW200to400_0.root","postW400to600_0.root","postW600to800_0.root","postW800to1200_0.root","postW1200to2500_0.root","postW2500toInf_0.root"};
  std::vector<TFile*> WJets_Files;
  for (int i = 0; i < WJets_FileNames.size(); i++){WJets_Files.push_back(new TFile(WJets_FileNames[i]));}			    
  std::vector<float> WJets_Total = GetTotal(WJets_Files);

  TH1F *histo_j1EtaWidth_WJets_0 = (TH1F*)WJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_W1Jets = (TH1F*)WJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_W2Jets = (TH1F*)WJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_W3Jets = (TH1F*)WJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_W4Jets = (TH1F*)WJets_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_W5Jets = (TH1F*)WJets_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_W6Jets = (TH1F*)WJets_Files[6]->Get(variable);
  TH1F *histo_j1EtaWidth_W7Jets = (TH1F*)WJets_Files[7]->Get(variable);

  double rawWJets = (histo_j1EtaWidth_WJets_0->Integral())+ (histo_j1EtaWidth_W1Jets->Integral())+ (histo_j1EtaWidth_W2Jets->Integral())+ (histo_j1EtaWidth_W3Jets->Integral())+ (histo_j1EtaWidth_W4Jets->Integral())+ (histo_j1EtaWidth_W5Jets->Integral())+ (histo_j1EtaWidth_W6Jets->Integral())+ (histo_j1EtaWidth_W7Jets->Integral());
  std::cout<<"raw WJets events: "<<rawWJets<<std::endl;

   // Scaling = (1/Totalevents)*Luminosity*NNLO-cross-section
  histo_j1EtaWidth_WJets_0->Scale((1.0/WJets_Total[0])*1885.122*50690);
  histo_j1EtaWidth_W1Jets->Scale((1.0/WJets_Total[1])*1885.122*1345);
  histo_j1EtaWidth_W2Jets->Scale((1.0/WJets_Total[2])*1885.122*359.7);
  histo_j1EtaWidth_W3Jets->Scale((1.0/WJets_Total[3])*1885.122*48.91);
  histo_j1EtaWidth_W4Jets->Scale((1.0/WJets_Total[4])*1885.122*12.05);
  histo_j1EtaWidth_W5Jets->Scale((1.0/WJets_Total[5])*1885.122*5.501);
  histo_j1EtaWidth_W6Jets->Scale((1.0/WJets_Total[6])*1885.122*1.329);
  histo_j1EtaWidth_W7Jets->Scale((1.0/WJets_Total[7])*1885.122*0.03216);
  
  
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W1Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W2Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W3Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W4Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W5Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W6Jets);
  histo_j1EtaWidth_WJets_0->Add(histo_j1EtaWidth_W7Jets);

  gPad->Update();

  double integralWJets = histo_j1EtaWidth_WJets_0->Integral();
  std::cout<<"integral of WJets bkg here:"<<integralWJets<<std::endl;

  histo_j1EtaWidth_WJets_0->SetName(((std::string("WJets")+varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_WJets_0);
  std::vector<const char *> Zvv_FileNames = {"postZ100to200_0.root","postZ200to400_0.root","postZ400to600_0.root","postZ600to800_0.root","postZ800to1200_0.root","postZ1200to2500_0.root","postZ2500toInf_0.root"};
  std::vector<TFile *> Zvv_Files;
  for (int i = 0; i < Zvv_FileNames.size(); i++) {Zvv_Files.push_back(new TFile(Zvv_FileNames[i]));}
  std::vector<float> Zvv_Total = GetTotal(Zvv_Files);

  TH1F *histo_j1EtaWidth_100to200 = (TH1F*)Zvv_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_200to400 = (TH1F*)Zvv_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_400to600 = (TH1F*)Zvv_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_600to800 = (TH1F*)Zvv_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_800to1200 = (TH1F*)Zvv_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_1200to2500 = (TH1F*)Zvv_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_2500toInf = (TH1F*)Zvv_Files[6]->Get(variable);

  double rawZvvJets = (histo_j1EtaWidth_100to200->Integral())+ (histo_j1EtaWidth_200to400->Integral())+ (histo_j1EtaWidth_400to600->Integral())+ (histo_j1EtaWidth_600to800->Integral())+ (histo_j1EtaWidth_800to1200->Integral())+ (histo_j1EtaWidth_1200to2500->Integral())+ (histo_j1EtaWidth_2500toInf->Integral()); 
  std::cout<<"raw ZvvJets bkg:"<<rawZvvJets<<std::endl;
  
  // Scaling = (1/Totalevents)*Luminosity*LO-cross-section
  histo_j1EtaWidth_100to200->Scale((1.0/Zvv_Total[0])*1885.122*280.35);
  histo_j1EtaWidth_200to400->Scale((1.0/Zvv_Total[1])*1885.122*77.67);
  histo_j1EtaWidth_400to600->Scale((1.0/Zvv_Total[2])*1885.122*10.73);
  histo_j1EtaWidth_600to800->Scale((1.0/Zvv_Total[3])*1885.122*2.559);
  histo_j1EtaWidth_800to1200->Scale((1.0/Zvv_Total[4])*1885.122*1.1796);
  histo_j1EtaWidth_1200to2500->Scale((1.0/Zvv_Total[5])*1885.122*0.28833);
  histo_j1EtaWidth_2500toInf->Scale((1.0/Zvv_Total[6])*1885.122*0.006945);

  //Add the ZJetsToNuNu histograms to the first one
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_200to400);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_400to600);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_600to800);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_800to1200);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_1200to2500);
  histo_j1EtaWidth_100to200->Add(histo_j1EtaWidth_2500toInf);

  gPad->Update();

  histo_j1EtaWidth_100to200->SetName(((std::string("ZJets")+varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_100to200);

  double integralZvvJets = histo_j1EtaWidth_100to200->Integral();
  std::cout<<"integral of ZvvJets bkg here:"<<integralZvvJets<<std::endl;

  //opening background samples Gamma+jets
  std::vector<const char *> GJets_FileNames = {"postGJets40to100.root","postGJets100to200.root","postGJets200to400.root","postGJets400to600.root","postGJets600toInf.root"};
  std::vector<TFile *> GJets_Files;
  for (int i = 0; i < GJets_FileNames.size(); i++) {GJets_Files.push_back(new TFile(GJets_FileNames[i]));}
  std::vector<float> GJets_Total = GetTotal(GJets_Files);
  
  TH1F *histo_j1EtaWidth_G1Jets = (TH1F*)GJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_G2Jets = (TH1F*)GJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_G3Jets = (TH1F*)GJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_G4Jets = (TH1F*)GJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_G5Jets = (TH1F*)GJets_Files[4]->Get(variable);

  double rawGJets = (histo_j1EtaWidth_G1Jets->Integral())+ (histo_j1EtaWidth_G2Jets->Integral())+ (histo_j1EtaWidth_G3Jets->Integral())+ (histo_j1EtaWidth_G4Jets->Integral())+ (histo_j1EtaWidth_G5Jets->Integral()); 
  std::cout<<"raw GJets bkg:"<<rawGJets<<std::endl;
  //Scaling
  histo_j1EtaWidth_G1Jets->Scale((1.0/GJets_Total[0])*1885.122*17420);
  histo_j1EtaWidth_G2Jets->Scale((1.0/GJets_Total[1])*1885.122*5391);
  histo_j1EtaWidth_G3Jets->Scale((1.0/GJets_Total[2])*1885.122*1168);
  histo_j1EtaWidth_G4Jets->Scale((1.0/GJets_Total[3])*1885.122*132.5);
  histo_j1EtaWidth_G5Jets->Scale((1.0/GJets_Total[4])*1885.122*44.05);
  
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G2Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G3Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G4Jets);
  histo_j1EtaWidth_G1Jets->Add(histo_j1EtaWidth_G5Jets);

  gPad->Update();

  histo_j1EtaWidth_G1Jets->SetName(((std::string("GJets")+varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_G1Jets);
  double integralGJets = histo_j1EtaWidth_G1Jets->Integral();
  std::cout<<"integral of GJets bkg here:"<<integralGJets<<std::endl;

//opening DYJetsToLL backgrounds
  std::vector<const char *> DYJets_FileNames = {"postDY_MLM_0.root","postDY100to200.root","postDY200to400.root","postDY400to600.root","postDY600to800.root","postDY800to1200.root","postDY1200to2500.root","postDY2500toInf.root"};
  std::vector<TFile *> DYJets_Files;
  for (int i = 0; i < DYJets_FileNames.size(); i++) {DYJets_Files.push_back(new TFile(DYJets_FileNames[i]));}
  std::vector<float> DYJets_Total = GetTotal(DYJets_Files);

  TH1F *histo_j1EtaWidth_DY1Jets = (TH1F*)DYJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_DY2Jets = (TH1F*)DYJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_DY3Jets = (TH1F*)DYJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_DY4Jets = (TH1F*)DYJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_DY5Jets = (TH1F*)DYJets_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_DY6Jets = (TH1F*)DYJets_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_DY7Jets = (TH1F*)DYJets_Files[6]->Get(variable);
  TH1F *histo_j1EtaWidth_DY8Jets = (TH1F*)DYJets_Files[7]->Get(variable);

  double rawDY = (histo_j1EtaWidth_DY1Jets->Integral()) + (histo_j1EtaWidth_DY2Jets->Integral())+ (histo_j1EtaWidth_DY3Jets->Integral())+ (histo_j1EtaWidth_DY4Jets->Integral())+ (histo_j1EtaWidth_DY5Jets->Integral()) + (histo_j1EtaWidth_DY6Jets->Integral()) + (histo_j1EtaWidth_DY7Jets->Integral()) + (histo_j1EtaWidth_DY8Jets->Integral());  
  std::cout<<"raw DYJets bkg:"<<rawDY<<std::endl;

  histo_j1EtaWidth_DY1Jets->Scale((1.0/DYJets_Total[0])*1885.122*4895);
  histo_j1EtaWidth_DY2Jets->Scale((1.0/DYJets_Total[1])*1885.122*148);
  histo_j1EtaWidth_DY3Jets->Scale((1.0/DYJets_Total[2])*1885.122*40.94);
  histo_j1EtaWidth_DY4Jets->Scale((1.0/DYJets_Total[3])*1885.122*5.497);
  histo_j1EtaWidth_DY5Jets->Scale((1.0/DYJets_Total[4])*1885.122*1.354);
  histo_j1EtaWidth_DY6Jets->Scale((1.0/DYJets_Total[5])*1885.122*0.6250);
  histo_j1EtaWidth_DY7Jets->Scale((1.0/DYJets_Total[6])*1885.122*0.1511);
  histo_j1EtaWidth_DY8Jets->Scale((1.0/DYJets_Total[7])*1885.122*0.003647);

  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY2Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY3Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY4Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY5Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY6Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY7Jets);
  histo_j1EtaWidth_DY1Jets->Add(histo_j1EtaWidth_DY8Jets);

  gPad->Update();
  
  histo_j1EtaWidth_DY1Jets->SetName(((std::string("DYJets")+varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_DY1Jets);
  double integralDY = histo_j1EtaWidth_DY1Jets->Integral();
  std::cout<<"integral of DYJets bkg here:"<<integralDY<<std::endl;

  //opening background TTJets
  std::vector<const char *> TTJets_FileNames = {"postTTJets_MLM.root"};
  std::vector<TFile *> TTJets_Files;
  for(int i = 0; i < TTJets_FileNames.size(); i++) {TTJets_Files.push_back(new TFile(TTJets_FileNames[i]));}
  std::vector<float> TTJets_Total = GetTotal(TTJets_Files);

  TH1F *histo_j1EtaWidth_TTJets = (TH1F*)TTJets_Files[0]->Get(variable);

  double rawTTJets = histo_j1EtaWidth_TTJets->Integral();
  std::cout<<"raw TTJets here:"<<rawTTJets<<std::endl;

  histo_j1EtaWidth_TTJets->Scale((1.0/TTJets_Total[0])*1885.122*502.2);

  histo_j1EtaWidth_TTJets->SetName(((std::string("TTJets") + varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_TTJets);
  double integralTTJets = histo_j1EtaWidth_TTJets->Integral();
  std::cout<<"integral of integralTTJets bkg here:"<<integralTTJets<<std::endl;

 //addding some backgrounds like WW, WZ, ZZ	

  std::vector<const char *> WW_FileNames = {"postWW.root"};
  std::vector<TFile *> WW_Files;
  for (int i = 0; i < WW_FileNames.size(); i++) {WW_Files.push_back(new TFile(WW_FileNames[i]));}
  std::vector<float> WW_Total = GetTotal(WW_Files);

  TH1F *histo_j1EtaWidth_WW = (TH1F*)WW_Files[0]->Get(variable);
  histo_j1EtaWidth_WW->Scale((1.0/WW_Total[0])*1885.122*118.7);
  histo_j1EtaWidth_WW->SetName(((std::string("WW")+varname)).c_str());
  double integralWW = histo_j1EtaWidth_WW->Integral();
  std::cout<<"integral of WW bkg here:"<<integralWW<<std::endl;
 
  std::vector<const char *> WZ_FileNames = {"postWZ.root"};
  std::vector<TFile *> WZ_Files;
  for (int i = 0; i < WZ_FileNames.size(); i++) {WZ_Files.push_back(new TFile(WZ_FileNames[i]));}
  std::vector<float> WZ_Total = GetTotal(WW_Files);
  
  TH1F *histo_j1EtaWidth_WZ = (TH1F*)WZ_Files[0]->Get(variable);
  histo_j1EtaWidth_WZ->Scale((1.0/WZ_Total[0])*1885.122*47.2);

  histo_j1EtaWidth_WZ->SetName(((std::string("WZ")+varname)).c_str());
  double integralWZ = histo_j1EtaWidth_WZ->Integral();
  std::cout<<"integral of WZ bkg here:"<<integralWZ<<std::endl;

  std::vector<const char *> ZZ_FileNames = {"postZZ.root"};
  std::vector<TFile *> ZZ_Files;
  for (int i = 0; i < ZZ_FileNames.size(); i++) {ZZ_Files.push_back(new TFile(ZZ_FileNames[i]));}
  std::vector<float> ZZ_Total = GetTotal(ZZ_Files);

  TH1F *histo_j1EtaWidth_ZZ = (TH1F*)ZZ_Files[0]->Get(variable);
  histo_j1EtaWidth_ZZ->Scale((1.0/ZZ_Total[0])*1885.122*16.6);

  histo_j1EtaWidth_ZZ->SetName(((std::string("ZZ")+varname)).c_str());
  double integralZZ = histo_j1EtaWidth_ZZ->Integral();
  std::cout<<"integral of ZZ bkg here:"<<integralZZ<<std::endl;
 
  //Ading WW, WZ, ZZ together
  histo_j1EtaWidth_WW->Add(histo_j1EtaWidth_WZ);
  histo_j1EtaWidth_WW->Add(histo_j1EtaWidth_ZZ);
  
  histo_j1EtaWidth_WW->SetName(((std::string("DiBoson")+varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_WW);
  double TotintegralDiBoson = histo_j1EtaWidth_WW->Integral();
  std::cout<<"integral of WW/WZ/ZZ bkg here:"<<TotintegralDiBoson<<std::endl;

 //opening QCD background files (HT-binned samples)
  std::vector<const char *> QJets_FileNames = {"postQCD100to200_0.root","postQCD200to300_0.root","postQCD300to500_0.root","postQCD500to700_0.root","postQCD700to1000_0.root","postQCD1000to1500_0.root","postQCD1500to2000_0.root","postQCD2000toInf_0.root"};
  std::vector<TFile *> QJets_Files;
  for (int i = 0; i < QJets_FileNames.size(); i++) {QJets_Files.push_back(new TFile(QJets_FileNames[i]));}
  std::vector<float> QJets_Total = GetTotal(QJets_Files);

  TH1F *histo_j1EtaWidth_Q1Jets = (TH1F*)QJets_Files[0]->Get(variable);
  TH1F *histo_j1EtaWidth_Q2Jets = (TH1F*)QJets_Files[1]->Get(variable);
  TH1F *histo_j1EtaWidth_Q3Jets = (TH1F*)QJets_Files[2]->Get(variable);
  TH1F *histo_j1EtaWidth_Q4Jets = (TH1F*)QJets_Files[3]->Get(variable);
  TH1F *histo_j1EtaWidth_Q5Jets = (TH1F*)QJets_Files[4]->Get(variable);
  TH1F *histo_j1EtaWidth_Q6Jets = (TH1F*)QJets_Files[5]->Get(variable);
  TH1F *histo_j1EtaWidth_Q7Jets = (TH1F*)QJets_Files[6]->Get(variable);
  TH1F *histo_j1EtaWidth_Q8Jets = (TH1F*)QJets_Files[7]->Get(variable);

  histo_j1EtaWidth_Q1Jets->Scale((1.0/QJets_Total[0])*1885.122*27500000);
  histo_j1EtaWidth_Q2Jets->Scale((1.0/QJets_Total[1])*1885.122*1735000);
  histo_j1EtaWidth_Q3Jets->Scale((1.0/QJets_Total[2])*1885.122*367000);
  histo_j1EtaWidth_Q4Jets->Scale((1.0/QJets_Total[3])*1885.122*29370);
  histo_j1EtaWidth_Q5Jets->Scale((1.0/QJets_Total[4])*1885.122*6524);
  histo_j1EtaWidth_Q6Jets->Scale((1.0/QJets_Total[5])*1885.122*1064);
  histo_j1EtaWidth_Q7Jets->Scale((1.0/QJets_Total[6])*1885.122*121.5);
  histo_j1EtaWidth_Q8Jets->Scale((1.0/QJets_Total[7])*1885.122*25.42);

  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q2Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q3Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q4Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q5Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q6Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q7Jets);
  histo_j1EtaWidth_Q1Jets->Add(histo_j1EtaWidth_Q8Jets);

  gPad->Update();
  
  double integral = histo_j1EtaWidth_Q1Jets->Integral();
  std::cout<<"integral of QCD bkg here:"<<integral<<std::endl;

  histo_j1EtaWidth_Q1Jets->SetName(((std::string("QCD")+varname)).c_str());
  hs_histo.push_back(histo_j1EtaWidth_Q1Jets);
  
  if (mchi == "1GeV")
    {
      TFile *f_signal_1GeVfile = new TFile("postSignal_mchi1GeV.root");
      TH1F *histo_signal_1GeV = (TH1F*)f_signal_1GeVfile ->Get(variable);
      histo_signal_1GeV->Scale((1.0/629)*1885.122*0.056);
      histo_signal_1GeV->SetLineColor(kRed);
      histo_signal_1GeV->SetLineWidth(2);
      histo_signal_1GeV->Draw("HIST SAME");
      
      double integralsignal_1GeV = histo_signal_1GeV->Integral(); 
      std::cout<<"integral of signal_1GeV here:"<<integralsignal_1GeV<<std::endl;
      histo_signal_1GeV->SetName((std::string("DM")+varname).c_str());
      hs_histo.push_back(histo_signal_1GeV);

      hs_save(mchi,cat,variable,hs_histo);
    }
  if (mchi == "5GeV")
    {
      TFile *f_signal_5GeVfile = new TFile("postSignal.root");
      TH1F *histo_signal_5GeV = (TH1F*)f_signal_5GeVfile ->Get(variable);
      std::vector<TFile *> SignalFile = {f_signal_5GeVfile};
      std::vector<float> Signal_Totals = GetTotal(SignalFile);
      histo_signal_5GeV->Scale((1.0/Signal_Totals[0])*1885.122*0.047);
      histo_signal_5GeV->SetLineColor(kBlue);
      histo_signal_5GeV->SetLineWidth(2);
      histo_signal_5GeV->Draw("HIST SAME");
      
      double integralsignal_5GeV = histo_signal_5GeV->Integral(); 
      std::cout<<"integral of signal_5GeV here:"<<integralsignal_5GeV<<std::endl;
      histo_signal_5GeV->SetName((std::string("DM")+varname).c_str());
      hs_histo.push_back(histo_signal_5GeV);
      
      hs_save(mchi,cat,variable,hs_histo);
    }
  if (mchi == "10GeV")
    {
      TFile *f_signal_10GeVfile = new TFile("postSignal_mchi10GeV.root");
      TH1F *histo_signal_10GeV = (TH1F*)f_signal_10GeVfile ->Get(variable);
      histo_signal_10GeV->Scale((1.0/4052)*1885.122*0.04);
      histo_signal_10GeV->SetLineColor(kViolet+1);
      histo_signal_10GeV->SetLineWidth(2);
      histo_signal_10GeV->Draw("HIST SAME");
      
      double integralsignal_10GeV = histo_signal_10GeV->Integral(); 
      std::cout<<"integral of signal_10GeV here:"<<integralsignal_10GeV<<std::endl;
      histo_signal_10GeV->SetName((std::string("DM")+varname).c_str());
      hs_histo.push_back(histo_signal_10GeV);

      hs_save(mchi,cat,variable,hs_histo);
    }
  if (mchi == "20GeV")
    {
      TFile *f_signal_20GeVfile = new TFile("postSignal_mchi20GeV.root");
      TH1F *histo_signal_20GeV = (TH1F*)f_signal_20GeVfile ->Get(variable);
      histo_signal_20GeV->Scale((1.0/9998)*1885.122*0.034);
      histo_signal_20GeV->SetLineColor(kMagenta);
      histo_signal_20GeV->SetLineWidth(2);
      histo_signal_20GeV->Draw("HIST SAME");
      
      double integralsignal_20GeV = histo_signal_20GeV->Integral(); 
      std::cout<<"integral of signal_20GeV here:"<<integralsignal_20GeV<<std::endl;
      histo_signal_20GeV->SetName((std::string("DM")+varname).c_str());
      hs_histo.push_back(histo_signal_20GeV);

      hs_save(mchi,cat,variable,hs_histo);
    }
  if (mchi == "50GeV")
    {
      TFile *f_signal_50GeVfile = new TFile("postSignal_mchi50GeV.root");
      TH1F *histo_signal_50GeV = (TH1F*)f_signal_50GeVfile ->Get(variable);
      histo_signal_50GeV->Scale((1.0/9999)*1885.122*0.025);
      histo_signal_50GeV->SetLineColor(kSpring-1);
      histo_signal_50GeV->SetLineWidth(2);
      histo_signal_50GeV->Draw("HIST SAME");
      
      double integralsignal_50GeV = histo_signal_50GeV->Integral(); 
      std::cout<<"integral of signal_50GeV here:"<<integralsignal_50GeV<<std::endl;
      histo_signal_50GeV->SetName((std::string("DM")+varname).c_str());
      hs_histo.push_back(histo_signal_50GeV);

      hs_save(mchi,cat,variable,hs_histo);
    }
  if (mchi == "100GeV")
    {
      TFile *f_signal_100GeVfile = new TFile("postSignal_mchi100GeV.root");
      TH1F *histo_signal_100GeV = (TH1F*)f_signal_100GeVfile ->Get(variable);
      histo_signal_100GeV->Scale((1.0/9994)*1885.122*0.019);
      histo_signal_100GeV->SetLineColor(kAzure+1);
      histo_signal_100GeV->SetLineWidth(2);
      histo_signal_100GeV->Draw("HIST SAME");
      
      double integralsignal_100GeV = histo_signal_100GeV->Integral(); 
      std::cout<<"integral of signal_100GeV here:"<<integralsignal_100GeV<<std::endl;
      histo_signal_100GeV->SetName((std::string("DM")+varname).c_str());
      hs_histo.push_back(histo_signal_100GeV);

      hs_save(mchi,cat,variable,hs_histo);
    }
}

int main(int argc, const char *argv[])
{
  std::string mchi = std::string(argv[1]);
  std::vector<const char*> variable;
  for (int i = 2; i < argc; i++)
    {
      variable.push_back(argv[i]);
    }
  for (int i = 0; i < variable.size(); i++)
    {
      std::string name = SampleName(variable[i]);
      std::vector<std::string> label = GetName(variable[i]);
      std::string varname = label[0];
      std::string cat = label[1];
      plotter(variable[i],name,varname,cat,mchi);
    } 
  return 0;
}
