{

   TFile sigFile("tmpRootDir/ZPrimeSSM_M-1250_LimitStudy_NoLifetime.root");
   TFile dyFile("tmpRootDir/DYJetsToLeptonLepton_LimitStudy_NoLifetime.root");
   TFile ttjetsFile("tmpRootDir/TTJets_MADGRAPH_LimitStudy_NoLifetime.root");
   TFile ttjets700to1000File("tmpRootDir/TTJets_MADGRAPH_Mtt-700to1000_LimitStudy_NoLifetime.root");
   TFile ttjets1000File("tmpRootDir/TTJets_MADGRAPH_Mtt-1000toInf_LimitStudy_NoLifetime.root");
   TFile wjetsFile("tmpRootDir/WJets_LimitStudy_NoLifetime.root");
   TFile wjets50to70File("tmpRootDir/WJets_Pt-50to70_LimitStudy_NoLifetime.root");
   TFile wjets70to100File("tmpRootDir/WJets_Pt-70to100_LimitStudy_NoLifetime.root");
   TFile wjets100File("tmpRootDir/WJets_Pt-100_LimitStudy_NoLifetime.root");
   TFile wwFile("tmpRootDir/WW_LimitStudy_NoLifetime.root");
   TFile wzFile("tmpRootDir/WZ_LimitStudy_NoLifetime.root");
   TFile zzFile("tmpRootDir/ZZ_LimitStudy_NoLifetime.root");
   TFile qcdFile("tmpRootDir/QCD_LimitStudy_NoLifetime.root");
   TFile dataFile("tmpRootDir/Data_LimitStudy_NoLifetime.root");
   
   TFile outFile("stack_noLifetime_separateBG_1250norm.root", "new");
   

   
   
   //TH1F *dyhist = dyFile->Get("DYJetsToLeptonLeptonHistosDirectory/_elecMuonMetMassNMinus1RebinHisto");
   //dyhist->Scale(6900.28/dyhist->Integral());
   
   TLegend *theLegend = new TLegend(0.65, 0.525, 0.9, 0.825);
   theLegend.SetFillStyle(0);
   theLegend.SetShadowColor(0);
   
   
   
   
   Float_t xbin[] = {0,50,100,150,200,250,300,400,600,900,1500};
   
   TH1F *sigHist = sigFile->Get("ZPrimeSSM_M-1250HistosDirectory/_elecMuonMetMassNMinus1RebinHisto");
   sigHist->SetName("sig");
   sigHist->Scale(21.43/sigHist->Integral());
   sigHist->SetFillColor(kWhite);
   
   theLegend->AddEntry(sigHist, "Z'#rightarrow #tau#tau 1250 GeV", "f");
      
   
   TH1F *dyhist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   dyhist->SetBinContent(1,0.);
   dyhist->SetBinContent(2,3350.);
   dyhist->SetBinContent(3,3590.);
   dyhist->SetBinContent(4,440.);
   dyhist->SetBinContent(5,66.);
   dyhist->SetBinContent(6,15.);
   dyhist->SetBinContent(7,35.);
   dyhist->SetBinContent(8,0.);
   dyhist->SetBinContent(9,2.5);
   dyhist->SetBinContent(10,0.);
   
   
   dyhist->Scale(19.446/19.7);
   
   dyhist->SetFillColor(kRed);

   theLegend->AddEntry(dyhist, "DY->LL", "f");
   
   TH1F *ttjetshist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   ttjetshist->SetBinContent(1, 0.);
   ttjetshist->SetBinContent(2,2.7);
   ttjetshist->SetBinContent(3,41.);
   ttjetshist->SetBinContent(4,70.);
   ttjetshist->SetBinContent(5,83.);
   ttjetshist->SetBinContent(6,53.);
   ttjetshist->SetBinContent(7,57.);
   ttjetshist->SetBinContent(8,28.);
   ttjetshist->SetBinContent(9,4.5);
   ttjetshist->SetBinContent(10,0.2);
   


   ttjetshist->Scale(19.446/19.7);
   
   ttjetshist->SetFillColor(kBlue);
   
   
   theLegend->AddEntry(ttjetshist, "ttjets", "f");

   
   TH1F *wjetshist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   wjetshist->SetBinContent(1, 0.);
   wjetshist->SetBinContent(2,0.7);
   wjetshist->SetBinContent(3,61.);
   wjetshist->SetBinContent(4,146.);
   wjetshist->SetBinContent(5,17.);
   wjetshist->SetBinContent(6,10.);
   wjetshist->SetBinContent(7,8.);
   wjetshist->SetBinContent(8,5.4);
   wjetshist->SetBinContent(9,1.7);
   wjetshist->SetBinContent(10,0.9);
   
   
   wjetshist->Scale(19.446/19.7);
   
   wjetshist->SetFillColor(kGreen);
   
   theLegend->AddEntry(wjetshist, "wjets", "f");

   
   TH1F *wwhist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   wwhist->SetBinContent(1,0.);
   wwhist->SetBinContent(2,27.2);
   wwhist->SetBinContent(3,206.);
   wwhist->SetBinContent(4,215.);
   wwhist->SetBinContent(5,142.);
   wwhist->SetBinContent(6,78.);
   wwhist->SetBinContent(7,75.9);
   wwhist->SetBinContent(8,34.1);
   wwhist->SetBinContent(9,6.2);
   wwhist->SetBinContent(10,1.9);
   
   
   wwhist->Scale(19.446/19.7);
   
   wwhist->SetFillColor(kMagenta);
   theLegend->AddEntry(wwhist, "WW", "f");
   
   TH1F *wzhist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   wzhist->SetBinContent(1,0.);
   wzhist->SetBinContent(2,2.8);
   wzhist->SetBinContent(3,19.9);
   wzhist->SetBinContent(4,27.2);
   wzhist->SetBinContent(5,17.);
   wzhist->SetBinContent(6,10.5);
   wzhist->SetBinContent(7,9.8);
   wzhist->SetBinContent(8,5.1);
   wzhist->SetBinContent(9,1.3);
   wzhist->SetBinContent(10,0.34);
   
   
   wzhist->Scale(19.446/19.7);
   
   wzhist->SetFillColor(kYellow);
   theLegend->AddEntry(wzhist, "WZ", "f");
  
  
   TH1F *zzhist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   zzhist->SetBinContent(1,0.);
   zzhist->SetBinContent(2,1.2);
   zzhist->SetBinContent(3,8.1);
   zzhist->SetBinContent(4,7.3);
   zzhist->SetBinContent(5,5.0);
   zzhist->SetBinContent(6,1.82);
   zzhist->SetBinContent(7,2.5);
   zzhist->SetBinContent(8,1.24);
   zzhist->SetBinContent(9,0.3);
   zzhist->SetBinContent(10,0.0);
   
   
   zzhist->Scale(19.446/19.7);  
   
   zzhist->SetFillColor(kCyan);
   theLegend->AddEntry(zzhist, "ZZ", "f");
   
   //TH1F *qcdhist = qcdFile->Get("QCDHistosDirectory/_elecMuonMetMassNMinus1RebinHisto");
   //qcdhist->Scale(350.96/qcdhist->Integral());
   
   TH1F *qcdhist = new TH1F("", "Elec Muon MEt Mass", 10, xbin);
   
   qcdhist->SetBinContent(1,0.);
   qcdhist->SetBinContent(2,94.);
   qcdhist->SetBinContent(3,313.);
   qcdhist->SetBinContent(4,108.);
   qcdhist->SetBinContent(5,44.);
   qcdhist->SetBinContent(6,42.);
   qcdhist->SetBinContent(7,18.);
   qcdhist->SetBinContent(8,6.);
   qcdhist->SetBinContent(9,4.7.);
   qcdhist->SetBinContent(10,0.);
   
  
   qcdhist->Scale(19.446/19.7);
  
   qcdhist->SetFillColor(kOrange);
   theLegend->AddEntry(qcdhist, "QCD", "f");
  



   THStack *stack = new THStack("stack", "");
   
   stack->Add(zzhist);
   stack->Add(wzhist);
   stack->Add(wwhist);
   stack->Add(ttjetshist);
   
   stack->Add(qcdhist);
   stack->Add(wjetshist);
   
   stack->Add(dyhist);
   
   //stack->Add(sigHist);
   
   stack->SetTitle("massVis (e, #mu, ME_{T})");
   
   stack->Draw();
   
   sigHist->Draw("same");
   
   


   
   TH1F *data_obs = new TH1F("data_obs", "Elec Muon MEt Mass", 10, xbin);
   
   data_obs->SetBinContent(1,0.);
   data_obs->SetBinContent(2,3430.);
   data_obs->SetBinContent(3,4300.);
   data_obs->SetBinContent(4,1000.);
   data_obs->SetBinContent(5,388.);
   data_obs->SetBinContent(6,198.);
   data_obs->SetBinContent(7,176.);
   data_obs->SetBinContent(8,98.);
   data_obs->SetBinContent(9,20.);
   data_obs->SetBinContent(10,3.0);
   
   data_obs->SetBinError(1,0.);
   data_obs->SetBinError(2,60.);
   data_obs->SetBinError(3,70.);
   data_obs->SetBinError(4,30.);
   data_obs->SetBinError(5,20.);
   data_obs->SetBinError(6,14.);
   data_obs->SetBinError(7,13.);
   data_obs->SetBinError(8,10.);
   data_obs->SetBinError(9,4.);
   data_obs->SetBinError(10,1.7);
   
   data_obs->Scale(19.446/19.7);
   
   data_obs->Draw("ep2 same");
   
   theLegend->AddEntry(data_obs, "Data", "ep2");
   
   theLegend->Draw("same");
   
   outFile->Write();
   
   
}


