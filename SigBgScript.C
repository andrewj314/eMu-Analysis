{
   gROOT->Reset();
   c1 = new TCanvas("c1","Elec Muon DCA/Error Sig/BG",200,10,700,500);
   TPad *pad = new TPad("pad","",0,0,1,1);
   pad->SetFillColor(10);
   pad->SetGrid();
   pad->Draw();
   pad->cd();

      // draw a frame to define the range
   TH1F *hr = c1->DrawFrame(0,0,5,100);
   //hr->SetTitle("#frac{#sqrt{IP_{e_{1}}^{2} + IP_{e_{2}}^{2}}}{#sqrt{#sigma_{IP}_{e_{1}}^{2} + #sigma_{IP}_{e_{2}}^{2}}}");
   //hr->SetTitle("DCA (e trk, #mu trk) / #sigma_{DCA}");
   //hr->SetTitle("#Delta p_{T} (e, #mu)");
   //hr->SetTitle("cos#Delta#phi(e, #mu)");
   //hr->SetTitle("cos#Delta#phi(l_{hi-p_{T}}, #Delta p_{T})");
   //hr->SetTitle("cos#Delta#phi(l_{hi-p_{T}}, ME_{T})");
   //hr->SetTitle("Lifetime (e_{1}, e_{2})");   


   //hr->SetTitle("Muon Pt");
    //hr->SetTitle("#frac{#sqrt{(IP_{e} + IP_{#mu})^{2}}}{#sqrt{#sigma_{IP}_{e}^{2} + #sigma_{IP}_{#mu}^{2}}}");
   //hr->SetTitle("#tau Lifetime (DCA OR Impact Parameter)");
   hr->SetTitle("Lifetime");
   //hr->SetXTitle("Min #Delta p_{T} (e, #mu) Threshold");
   //hr->SetXTitle("Max cos#Delta#phi(e, #mu) Threshold");
   //hr->SetXTitle("Min cos#Delta#phi(l_{hi-p_{T}}, ME_{T}) Threshold");
   //hr->SetXTitle("Min DCA / IP Threshold");
   hr->SetXTitle("Min Lifetime Threshold");
   
   hr->SetYTitle("Efficiency");
   //hr->SetAxisColor(kRed);
   //hr->SetLabelColor(kBlack);
   pad->GetFrame()->SetFillColor(10);
   pad->GetFrame()->SetBorderSize(12);

      // create first graph


/*
   gr2 = new TGraphErrors("DataEff.txt", "%lg %lg %lg");
   gr2->SetMarkerColor(kBlack);
   gr2->SetMarkerStyle(8);
   gr2->SetLineColor(kRed);
   gr2->SetName("gr2");
   gr2->Draw("LP"); 
*/   
   gr3 = new TGraphErrors("TTJetsEff.txt", "%lg %lg %lg");
   gr3->SetMarkerColor(kBlue);
   gr3->SetMarkerStyle(8);
   gr3->SetLineColor(kBlue);
   gr3->SetName("gr3");
   gr3->Draw("LP"); 
   
   /*
   gr4 = new TGraphErrors("WJetsEff.txt", "%lg %lg %lg");
   gr4->SetMarkerColor(kMagenta);
   gr4->SetMarkerStyle(8);
   gr4->SetName("gr4");
   gr4->Draw("LP"); 
   
   gr5 = new TGraphErrors("QCDEff.txt", "%lg %lg %lg");
   gr5->SetMarkerColor(kMagenta);
   gr5->SetMarkerStyle(8);
   gr5->SetName("gr5");   
   //gr5->Draw("LP"); 
   
   
   gr6 = new TGraphErrors("DYJetsToLeptonLeptonEff.txt", "%lg %lg %lg");
   gr6->SetMarkerColor(kRed-7);
   gr6->SetMarkerStyle(8);
   gr6->SetName("gr6");   
   //gr6->Draw("LP"); 
   
   gr7 = new TGraphErrors("DYToEEEff.txt", "%lg %lg %lg");
   gr7->SetMarkerColor(kBlue-9);
   gr7->SetName("gr7");
   gr7->SetMarkerStyle(8);
   //gr7->Draw("LP"); 
   
   */
   gr8 = new TGraphErrors("WWEff.txt", "%lg %lg %lg");
   gr8->SetMarkerColor(kGreen);
   gr8->SetMarkerStyle(8);
   gr8->SetName("gr8");
   gr8->Draw("LP"); 
   
   /*
   gr9 = new TGraphErrors("WZEff.txt", "%lg %lg %lg");
   gr9->SetMarkerColor(kCyan);
   gr9->SetMarkerStyle(8);
   gr9->SetName("gr9");
   gr9->Draw("LP"); 
   
   gr10 = new TGraphErrors("ZZEff.txt", "%lg %lg %lg");
   gr10->SetMarkerColor(kPink+10);
   gr10->SetMarkerStyle(8);
   gr10->SetName("gr10");
   gr10->Draw("LP"); 
   */
   
   gr11 = new TGraphErrors("SigEff.txt", "%lg %lg %lg");
   gr11->SetMarkerColor(kRed);
   gr11->SetMarkerStyle(8);
   gr11->SetLineColor(kRed);
   gr11->SetName("gr11");
   gr11->Draw("LP");    
   
   //TF1 *g1 = new TF1("g1", "100*(1-sqrt(1-exp(-x*x/2)))", 0, 8);
   
   TF1 *g1 = new TF1("g1", "100*exp(-x*x/2)", 0, 5);
   
   g1->SetLineColor(kMagenta);
   //g1->Draw("same");
   TLegend *legend = new TLegend(0.7, 0.15, 0.85, 0.25);
   
   
   //legend->AddEntry("gr2", "Data Eff", "LP");
   //legend->AddEntry("g1", "Gaussian Distribution", "L");
   legend->AddEntry("gr3", "TTJets MC Eff", "LP");
   //legend->AddEntry("gr4", "W+Jets Eff", "LP");
   //legend->AddEntry("gr5", "QCD Eff", "LP");
   //legend->AddEntry("gr6", "DYJets->LL Eff", "LP");
   //legend->AddEntry("gr7", "Z->ee Eff", "LP");
   legend->AddEntry("gr8", "WW MC Eff", "LP");
   //legend->AddEntry("gr9", "WZ Eff", "LP");
   //legend->AddEntry("gr10", "ZZ Eff", "lp"); 
   legend->AddEntry("gr11", "Signal MC Eff", "LP");  
   //legend->Draw();
   //create a transparent pad drawn on top of the main pad
  
   c1->cd();
   TPad *overlay = new TPad("overlay","",0,0,1,1);
   overlay->SetFillStyle(4000);
   overlay->SetFillColor(0);
   overlay->SetFrameFillStyle(4000);
   overlay->Draw();
   overlay->cd();
   
   gr1 = new TGraphErrors("SigBgRatio.txt", "%lg %lg %lg");
   gr1->SetMarkerColor(kBlack);
   gr1->SetMarkerStyle(22);
   legend->AddEntry("gr1", "s/sqrt(s+b)", "lp");


   
   
   Double_t xmin = pad->GetUxmin();
   Double_t ymin = 0.;
   Double_t xmax = pad->GetUxmax();
   Double_t ymax = 4;
   TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
   hframe->GetXaxis()->SetLabelOffset(99);
   hframe->GetYaxis()->SetLabelOffset(99);
   gr1->Draw("P");
   legend->Draw();
      
   //Draw an axis on the right side
   TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
   axis->SetLineColor(kBlack);
   axis->SetLabelColor(kBlack);
   axis->SetTitle("#frac{s}{#sqrt{s+b}}");
   axis->Draw();
   
   

   

}
