#define EMuAnalyzer_cxx
#include "EMuAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>

using namespace std;

void EMuAnalyzer::Loop(){
  if (fChain == 0) return;

  Long64_t nentries = _nEvents;
  if(_nEvents < 0) nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
    
  nEvents = nentries;
  
  pdfWeightedEvents = 0;
  pdfWeightedEventsAcc = 0;
  pdfWeightedEventsEId = 0;
  pdfWeightedEventsMuonId = 0;
  pdfWeightedEventsTopology = 0;
  
  isrGammaWeightedEvents = 0;
  isrGammaWeightedEventsAcc = 0;
  isrGammaWeightedEventsEId = 0;
  isrGammaWeightedEventsMuonId = 0;
  isrGammaWeightedEventsTopology = 0;
  
  isrGluonWeightedEvents = 0;
  isrGluonWeightedEventsAcc = 0;
  isrGluonWeightedEventsEId = 0;
  isrGluonWeightedEventsMuonId = 0;
  isrGluonWeightedEventsTopology = 0;
  
  fsrWeightedEvents = 0;
  fsrWeightedEventsAcc = 0;
  fsrWeightedEventsEId = 0;
  fsrWeightedEventsMuonId = 0;
  fsrWeightedEventsTopology = 0;
  
  // Book histos and initialize out tree
  BookHistos();
  
  filterEff = 1;
  branchingFraction = 1;
  _doFactorization = 0;
  _bgEstimate = 0.;
  _bgEstimateError = 0.;
  
  std::cout << " **** Skimming efficiencies for 2011 **** \n";
  
  if (_label == "zprime"){
     branchingFraction = 0.0620;
   //theXSection = 1.9142e-12;	// 500 GeV
    _theXSection = make_pair<float, float>(signalXSection * 1e-12, signalXSection * 1e-12*0.02); // 2% error on xsection
    skimmedEvents = make_pair<int, int>((int)nEvents, 50232);
    theTreeOutName = "ZprimeTree";
  }
  else if(_label == "ztautau"){
    _theXSection = make_pair<float, float>(1.666e-9, 1.666e-9*0.054);  // NNLO for M > 20 with 5.4% error
    skimmedEvents = make_pair<int, int>((int)nEvents, 2032536);
//    _theXSection = 2.289e-9;  // NNLO Z->ll  M > 50
//    skimmedEvents = make_pair<int, int>((int)nEvents, 2661949);
//    _doFactorization = 1;
    theTreeOutName = "ZtautauTree";
    _bgEstimate = _zttEstimate;
    _bgEstimateError = _zttEstimateError;
 }
  else if(_label == "wjets"){
    _theXSection = make_pair<float, float>(31.314e-9, 31.314e-9*0.05);   // NNLO with 5% error
    skimmedEvents = make_pair<int, int>((int)nEvents, 81352581);
    theTreeOutName = "WjetsTree";
    if(_setRegion && _whichRegion.find("wjets") != string::npos) _doFactorization = 0;
    else _doFactorization = 1;
    _bgEstimate = _wjetsEstimate;
    _bgEstimateError = _wjetsEstimateError;
  }
  else if(_label == "ww"){
    //_theXSection = make_pair<float, float>(55.3e-12, 55.3e-12*0.151); // from EWK-10-009 with 15.1% error; update to EWK-11-010
    _theXSection = make_pair<float, float>(43.0e-12, 1.5e-12); // NLO
 skimmedEvents = make_pair<int, int>((int)nEvents, 4225916);
    theTreeOutName = "WWTree";
    _doFactorization = 1;
  }
// for 2011 data
  else if(_label == "wz"){
    //_theXSection = make_pair<float, float>(17.0e-12, 17.0e-12*0.166);  //with 16.6% error; update to EWK-11-010
    _theXSection = make_pair<float, float>(18.2e-12, 0.7e-12);  //NLO
    skimmedEvents = make_pair<int, int>((int)nEvents, 4265243);
    theTreeOutName = "WZTree";
    //_doFactorization = 1;
  }
//  else if(_label == "zccbar"){
//    _theXSection = 2.6e-9 * 0.696;
//    skimmedEvents = make_pair<int, int>((int)nEvents, 729769);
//    theTreeOutName = "ZCCbarTree";
//    //_doFactorization = 1;
//  }
//
  else if(_label == "ttbar"){
    _theXSection = make_pair<float, float>(165e-12, 165e-12*0.06);  // NNLL resummations with 6% error
    skimmedEvents = make_pair<int, int>((int)nEvents, 1089625);
    theTreeOutName = "TTBarTree";
    //_doFactorization = 1;
    _bgEstimate = _ttbarEstimate;
    _bgEstimateError = _ttbarEstimateError;
  }
  else if(_label == "qcd"){
    _theXSection = make_pair<float, float>(0.2966e-3 * 2.855e-4, 0.0);//41X qcd20MuPt15
    skimmedEvents = make_pair<int, int>((int)nEvents, 25080241); // 42X
    theTreeOutName = "IncMuPt15";
    //_theXSection = 2.966e-4 * 0.00118;//41X qcd20MuPt10
    //skimmedEvents = make_pair<int, int>((int)nEvents, 8797418); //42X qcd20MuPt10
    //theTreeOutName = "IncMuPt20";
    _doFactorization = 1;
    if(_setRegion && _whichRegion.find("qcd") != string::npos) _doFactorization = 0;
    _bgEstimate = _qcdEstimate;
    _bgEstimateError = _qcdEstimateError;
  }
  
  else if(_label == "data"){
    _theXSection = make_pair<float, float>(1.0, 0.0);
    theTreeOutName = "data";
  }
  else{
    cout << "label not recogized\n";
    cout << "labels are <zprime> <ztautau> <zee> <wjets> <ttbar> <qcd20to30> <qcd30to80> <qcd80to170>\n";
    exit(1);
  }
 /*
  //get first entry and initialize pdf weights
  GetEntry(1);
  for(unsigned int pdfs = 0; pdfs < PDFWweights->size(); ++pdfs){
    pdfWeightsSum.push_back(0);
    pdfWeightsSumAcc.push_back(0);
    pdfWeightsSumEId.push_back(0);
    pdfWeightsSumMuonId.push_back(0);
    pdfWeightsSumTopology.push_back(0);
  }
*/  
  //eventsAtLumi = make_pair<double, double>(_theXSection.first * _theLumi*1e12, _theXSection.second * _theLumi*1e12);
  eventsAtLumi = make_pair<double, double>(_theXSection.first * _theLumi*1e12, _theXSection.second * _theLumi*1e12);
  skimEff = GetEfficiency(skimmedEvents.first, skimmedEvents.second);
  nSkimmedAtLumi = GetSurvivingEvents(eventsAtLumi, skimEff);
  
  // set variables for histos...
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;
/*    
    for(unsigned int vIt = 0; vIt < PDFWweights->size(); ++vIt){
      if(PDFWweights->at(vIt) < 10) pdfWeightsSum.at(vIt) += PDFWweights->at(vIt);
      else pdfWeightsSum.at(vIt) += 1;
    }
*/    
    //pdfWeightedEvents += GetCentralPDF();
    isrGammaWeightedEvents += ISRGammaWeight;
    isrGluonWeightedEvents += ISRGluonWeight;
    fsrWeightedEvents += FSRWeight;
    
    _eventPUWeight = 1.;
    
    if(_label != "data"){
      if(InTimePU < _thePUWeights.size()) _eventPUWeight = _thePUWeights.at(InTimePU);
      if(InTimePU >= _thePUWeights.size()) _eventPUWeight = _thePUWeights.at(_thePUWeights.size() -1);
    }

    puWeightedEvents += _eventPUWeight; 
    _pileUPDistro->Fill(InTimePU);
    _pileUPReweighted->Fill(InTimePU,_eventPUWeight);
    
    ResetCandCounters();
    if(!passedTrigger())continue;
    nTriggers++;   
    
    float maxElecMuDeltaR = 0;
    float maxElecPt = -10;
    float minElecEta = 3;
    float maxMuPt = -10;
    float minMuonEta = 3;
    
    float minElecHOverEm = 10;
    float minElecDEtaIn = 5;
    float minElecDphiIn = 5;
    float minElecSigmaIEta = 5;
    float maxElec2x5Ratio = -10;
    float maxElec1x5Ratio = -10;
    unsigned int minElecNHits = 25;
    float minElecEcalIsoDr03 = 50;
    float minElecTrkIsoDr03 = 50;
    float minElecEcalIsoDr04 = 50;
    float minElecTrkIsoDr04 = 50;
    float minElecIsoDr03 = 50;
    float minElecIsoDr04 = 50;
    
    float minMuonChi2 = 100;
    unsigned int maxMuonHits = 0;
    float minMuonIp = 10;
    float minMuonEcalIsoDr03 = 50;
    float minMuonTrkIsoDr03 = 50;
    float minMuonEcalIsoDr04 = 50;
    float minMuonTrkIsoDr04 = 50;
    float minMuonEcalIsoDr05 = 50;
    float minMuonTrkIsoDr05 = 50;
    float minMuonIsoDr03 = 50;
    float minMuonIsoDr04 = 50;
    float minMuonIsoDr05 = 50;
    
    float minCosDPhi = 10;
    float minCharge = 5;
    
    unsigned int minNBtagsHiEffTrkCnt = 10;
    unsigned int minNBtagsHiEffTrkCntControl = 10;
    unsigned int minNBtagsHiPurityTrkCnt = 10;
    unsigned int minNBtagsHiEffSimpleSecVtx = 10;
    unsigned int minNBtagsHiPuritySimpleSecVtx = 10;
    unsigned int minNBtagsCombSecVtx = 10;
    
    float minJetSumEt = 10000;
    float minJetMetSumEt = 10000;
    
    float minLeptonMetCosDPhi = 2.;
    
    unsigned int minNBtagsHiEffTrkCntFor2DPlot = 10;
    float minJetSumEtFor2DPlot = 10000;
    float minLeptonMetCosDPhiForPlotVsNBtags = 2.;
    float minLeptonMetCosDPhiForPlotVsEtSum = 2.;
    
    float maxMet = -10;
    
    float max2DZetaCut = -1000.;
    float minZetaVis = 1000.;
    float maxZeta = -1000.;
    float maxMass = -10.;
    float maxMetMass = -10.;
    float maxEMetMass = -10.;
    float maxMuMetMass = -10.;
    float maxTemplateMass = -10.;
    float maxEMuMetMassDataQCD = -10.;
    float maxEMuMetMassDataTTbar = -10.;
    float maxEMuMetMassDataWJets = -10.;
    
    float maxElecMuDeltaRNMinus1 = 0;
    float maxElecPtNMinus1 = -10;
    float minElecEtaNMinus1 = 10;
    float maxMuPtNMinus1 = -10;
    float minMuonEtaNMinus1 = 10;
    
    float minElecHOverEmNMinus1 = 10;
    float minElecDEtaInNMinus1 = 5;
    float minElecDphiInNMinus1 = 5;
    float minElecSigmaIEtaNMinus1 = 5;
    float maxElec2x5RatioNMinus1 = -10;
    float maxElec1x5RatioNMinus1 = -10;
    unsigned int minElecNHitsNMinus1 = 25;
    float minElecEcalIsoDr03NMinus1 = 50;
    float minElecTrkIsoDr03NMinus1 = 50;
    float minElecEcalIsoDr04NMinus1 = 50;
    float minElecTrkIsoDr04NMinus1 = 50;
    float minElecIsoDr03NMinus1 = 50;
    float minElecIsoDr04NMinus1 = 50;
    
    float minMuonChi2NMinus1 = 100;
    unsigned int maxMuonHitsNMinus1 = 0;
    float minMuonIpNMinus1 = 10;
    float minMuonEcalIsoDr03NMinus1 = 50;
    float minMuonTrkIsoDr03NMinus1 = 50;
    float minMuonEcalIsoDr04NMinus1 = 50;
    float minMuonTrkIsoDr04NMinus1 = 50;
    float minMuonEcalIsoDr05NMinus1 = 50;
    float minMuonTrkIsoDr05NMinus1 = 50;
    float minMuonIsoDr03NMinus1 = 50;
    float minMuonIsoDr04NMinus1 = 50;
    float minMuonIsoDr05NMinus1 = 50;
    
    float minCosDPhiNMinus1 = 10;
    float minChargeNMinus1 = 5;
    
    unsigned int minNBtagsHiEffTrkCntNMinus1 = 10;
    unsigned int minNBtagsHiPurityTrkCntNMinus1 = 10;
    unsigned int minNBtagsHiEffSimpleSecVtxNMinus1 = 10;
    unsigned int minNBtagsHiPuritySimpleSecVtxNMinus1 = 10;
    unsigned int minNBtagsCombSecVtxNMinus1 = 10;
    
    float minJetSumEtNMinus1 = 10000;
    float minLeptonMetCosDPhiNMinus1 = 2.;
    float maxMetNMinus1 = -10;
    float max2DZetaCutNMinus1 = -1000.;
    float minZetaVisNMinus1 = 1000.;
    float maxZetaNMinus1 = -1000.;
    unsigned int maxJetsNMinus1 = 100;

    unsigned int itBegin = 0;
    unsigned int itEnd = ePt->size();
    if(_label == "zprime"){
      for(unsigned int mIt = 0; mIt < muonMatched->size(); ++mIt){
        if(eIsMatched(mIt)) elecMatchCounter++;
        if(eIsMatched(mIt) && muonIsMatched(mIt)) muonMatchCounter++ ;
      }
      int theCand = getMatchedCand();
      if(theCand >= 0){
	itBegin = theCand;
	itEnd = itBegin + 1;
      }
      else continue; 							     
    }
    for(unsigned int vIt = itBegin; vIt < itEnd; ++vIt){	     
      bool electronId = passedEId(vIt);
      bool muonId = passedMuId(vIt);
      bool diTau = passedDiTau(vIt);
      //bool elecIso = passedElecIso(vIt);
      bool elecEcalIso = passedElecEcalIso(vIt);
      bool elecTrkIso = passedElecTrkIso(vIt);
      //bool muonIso = passedMuIso(vIt);
      bool muEcalIso = passedMuEcalIso(vIt);
      bool muTrkIso = passedMuTrkIso(vIt);
      bool diTauCharge = passedCharge(vIt);
      
      //define some control regions here...
      bool qcdControl = passedAcceptance(vIt)
	   && passedLooseEId(vIt) && passedMinElecTrkIso(vIt) && passedElecEcalIso(vIt)
           && passedLooseMuId(vIt) && passedMuEcalIso(vIt)
	   && passedCosDPhi(vIt) && passedBTagVeto(vIt) && passedJetSum(vIt) && passedMass(vIt);
	   //&& (eCharge->at(vIt)*muonCharge->at(vIt)) > 0;
      
      bool qcdAltControl = passedAcceptance(vIt)
	   && passedLooseEId(vIt) && passedElecEcalIso(vIt) 
           && passedLooseMuId(vIt) && passedMuEcalIso(vIt) && passedMinMuTrkIso(vIt)
	   && passedCosDPhi(vIt) && passedBTagVeto(vIt) && passedJetSum(vIt) && passedMass(vIt);
      
      bool ttbarControl = passedAcceptance(vIt) && passedEId(vIt) && passedMuId(vIt) 
	   && diTauCharge && passedMet(vIt) && passed2DZetaCut(vIt) && passedMass(vIt)
	   && nBtagsHiEffTrkCnt->at(vIt) >= _btagDiscCut;
      
/*      
      bool ttbarAltControl = passedAcceptance(vIt) && passedEId(vIt) && passedMuId(vIt) 
	   && diTauCharge && passedMet(vIt) && passed2DZetaCut(vIt)  && passedMass(vIt)
	   && jetSumEt->at(vIt) > _jetSumCut;
*/      
      bool ttbarAltControl = passedAcceptance(vIt) && passedEId(vIt) && passedMuId(vIt) 
	   && diTauCharge && passedMass(vIt) && passed2DZetaCut(vIt) 
	   && nJets->at(vIt) >= 2 && mEt->at(vIt) > 40;
      bool muMetRegion = muMEtMt->at(vIt) >= 50. && muMEtMt->at(vIt) <= 100.;
      
      bool wjetsControl = passedAcceptance(vIt)
           && passedLooseEIdForWJets(vIt) && passedElecEcalIso(vIt) && passedMinElecTrkIso(vIt) 
	   //&& passedEId(vIt)
	   && passedMuId(vIt) 
           && passedBTagVeto(vIt) && passedJetSum(vIt) && passedLeptonMETCosDphi(vIt) 
	   && passedMass(vIt) && passedMet(vIt) && muMetRegion;
      
      bool wjetsAltControl = passedAcceptance(vIt)
           && passedLooseEIdForWJets(vIt) && passedElecEcalIso(vIt) 
	   //&& passedEId(vIt)
	   && passedMuId(vIt) 
           && passedBTagVeto(vIt) && passedJetSum(vIt) && passedLeptonMETCosDphi(vIt) && passedMet(vIt) && passedMass(vIt) 
	   && (eCharge->at(vIt)*muonCharge->at(vIt)) > 0 && eMuCosDPhi->at(vIt) > _deltaPhiCut
	   && eMuPZeta->at(vIt) - _zeta2DSlope*eMuPZetaVis->at(vIt) < _zeta2DCut;
	    	   
      //bool ztautauControl = passedDiTauDeltaR(vIt) && passedElecEta(vIt) && passedMuonEta(vIt) && ePt->at(vIt) > 15. 
      //     && muonIsGlobalMuon->at(vIt) > 0 && muonPt->at(vIt) > 15.
      //     && passedHeepId(vIt) && passedMuId(vIt) && diTauCharge && passedBTagVeto(vIt) && passedJetSum(vIt)
      //     && passedCosDPhi(vIt) && passedLeptonMETCosDphi(vIt) && passed2DZetaCut(vIt);
	         
      if(_doFactorization){
        electronId = passedLooseEId(vIt);
	muonId = passedLooseMuId(vIt);
	diTau = passedLooseDiTau(vIt);
	//elecIso = passedLooseElecIso(vIt);
	elecEcalIso = passedLooseElecEcalIso(vIt);
	elecTrkIso = passedLooseElecTrkIso(vIt);
	//muonIso = passedLooseMuIso(vIt);
        muEcalIso = passedLooseMuEcalIso(vIt);
        muTrkIso = passedLooseMuTrkIso(vIt);
	diTauCharge = (eCharge->at(vIt)*muonCharge->at(vIt) < 0 || eCharge->at(vIt)*muonCharge->at(vIt) > 0);
      }
      
      // get ttbar control region taking into account factorization for qcd and WJets
      bool ttbarAltControlFact = passedAcceptance(vIt) && electronId && muonId
	   && diTauCharge && passedMet(vIt) && passed2DZetaCut(vIt)  && passedMass(vIt)
	   && jetSumEt->at(vIt) > _jetSumCut;
      
      
      bool jetCountVar = true;
      if(_setRegion && _whichRegion == "ttbarJetCount") jetCountVar = !passedJetCount(vIt);
      bool wjetsVar = true;
      if(_setRegion && _whichRegion.find("wjets") != string::npos) wjetsVar = muMetRegion;
      //if(_setRegion && _whichRegion == "wjets") wjetsVar = muMetRegion;

            // do N-1 plots first
      if(passedElecPt(vIt) && passedElecEta(vIt) && passedMuonPt(vIt) && passedMuonEta(vIt) && electronId 
        && muonId && diTau && jetCountVar && wjetsVar){
	nMinus1DeltaRCounter++;
        if(eMuDeltaR->at(vIt) > maxElecMuDeltaRNMinus1) maxElecMuDeltaRNMinus1 = eMuDeltaR->at(vIt);
      }
      
      if(passedDiTauDeltaR(vIt) && passedElecEta(vIt) && passedMuonPt(vIt) && passedMuonEta(vIt) && electronId 
        && muonId && diTau && jetCountVar && wjetsVar){
        nMinus1ElecPtCounter++;
        if(ePt->at(vIt) > maxElecPtNMinus1) maxElecPtNMinus1 = ePt->at(vIt);
      }            
      
      if(passedDiTauDeltaR(vIt) && passedElecPt(vIt) && passedMuonPt(vIt) && passedMuonEta(vIt) && electronId 
        && muonId && diTau && jetCountVar && wjetsVar){
	nMinus1ElecEtaCounter++;
	if(fabs(eEta->at(vIt)) < minElecEtaNMinus1) minElecEtaNMinus1 = eEta->at(vIt);
      }
      
      if(passedDiTauDeltaR(vIt) && passedElecPt(vIt) && passedElecEta(vIt) && passedMuonEta(vIt) && electronId 
        && muonId && diTau && jetCountVar && wjetsVar){
	nMinus1MuonPtCounter++;
	if(muonPt->at(vIt) > maxMuPtNMinus1) maxMuPtNMinus1 = muonPt->at(vIt);
      }
      
      if(passedDiTauDeltaR(vIt) && passedElecPt(vIt) && passedElecEta(vIt) && passedMuonPt(vIt) && electronId 
        && muonId && diTau && jetCountVar && wjetsVar){
         nMinus1MuonEtaCounter++;
	 if(fabs(muonEta->at(vIt)) < minMuonEtaNMinus1) minMuonEtaNMinus1 = muonEta->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepDEtaIn(vIt) && passedHeepDPhiIn(vIt) 
        && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) && elecEcalIso && elecTrkIso 
	&& muonId && diTau && jetCountVar && wjetsVar){
        nMinus1ElecEmFracCounter++;
	if(eHOverEm->at(vIt) < minElecHOverEmNMinus1) minElecHOverEmNMinus1 = eHOverEm->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDPhiIn(vIt) 
        && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) && elecEcalIso && elecTrkIso 
	&& muonId && diTau && jetCountVar && wjetsVar){
	nMinus1ElecDEtaInCounter++;
	if(fabs(eDeltaEtaIn->at(vIt)) < minElecDEtaInNMinus1) minElecDEtaInNMinus1 = eDeltaEtaIn->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) 
        && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) && elecEcalIso && elecTrkIso
	&& muonId && diTau && jetCountVar && wjetsVar){
	nMinus1ElecDPhiInCounter++;
	if(fabs(eDeltaPhiIn->at(vIt)) < minElecDphiInNMinus1) minElecDphiInNMinus1 = eDeltaPhiIn->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) 
        && passedHeepDPhiIn(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) && elecEcalIso && elecTrkIso 
	&& muonId && diTau && jetCountVar && wjetsVar){
	nMinus1ElecSigmaIEtaCounter++;
        if(eSigmaIEtaIEta->at(vIt) < minElecSigmaIEtaNMinus1) minElecSigmaIEtaNMinus1 = eSigmaIEtaIEta->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) 
        && passedHeepDPhiIn(vIt) && passedHeepSigmaIEta(vIt) && passedElecMissingHits(vIt) && elecEcalIso && elecTrkIso 
	&& muonId && diTau && jetCountVar && wjetsVar){
        nMinus1ElecClusterRatioCounter++;
        if(eSCE2x5->at(vIt)/eSCE5x5->at(vIt) > maxElec2x5RatioNMinus1) maxElec2x5RatioNMinus1 = eSCE2x5->at(vIt)/eSCE5x5->at(vIt);
        if(eSCE1x5->at(vIt)/eSCE5x5->at(vIt) > maxElec1x5RatioNMinus1) maxElec1x5RatioNMinus1 = eSCE1x5->at(vIt)/eSCE5x5->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) && passedHeepDPhiIn(vIt) 
        && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && elecEcalIso && elecTrkIso 
	&& muonId && diTau && jetCountVar && wjetsVar){
	nMinus1ElecMissingHitsCounter++;
	if(eMissingHits->at(vIt) < minElecNHitsNMinus1) minElecNHitsNMinus1 = eMissingHits->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) 
        && passedHeepDPhiIn(vIt) && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) 
	&& muonId && diTau && passedLooseElecIso(vIt) && jetCountVar && wjetsVar){
	nMinus1ElecIsoCounter++;
	if(heepEcalIso->at(vIt) + heepTrkIso->at(vIt) < minElecIsoDr03NMinus1) minElecIsoDr03NMinus1 = heepEcalIso->at(vIt) + heepTrkIso->at(vIt);
	if(eEcalIsoPat->at(vIt) + eTrkIsoPat->at(vIt) < minElecIsoDr04NMinus1) minElecIsoDr04NMinus1 = eEcalIsoPat->at(vIt) + eTrkIsoPat->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) 
        && passedHeepDPhiIn(vIt) && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) && elecTrkIso
	&& muonId && diTau && passedLooseElecIso(vIt) && jetCountVar && wjetsVar){
	nMinus1ElecEcalIsoCounter++;
        if(heepEcalIso->at(vIt) < minElecEcalIsoDr03NMinus1) minElecEcalIsoDr03NMinus1 = heepEcalIso->at(vIt);
        if(eEcalIsoPat->at(vIt) < minElecEcalIsoDr04NMinus1) minElecEcalIsoDr04NMinus1 = eEcalIsoPat->at(vIt);
      }
      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepEcalDriven(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) 
        && passedHeepDPhiIn(vIt) && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) && elecEcalIso
	&& muonId && diTau && passedLooseElecTrkIso(vIt) && jetCountVar && wjetsVar){
	nMinus1ElecTrkIsoCounter++;
      	if(heepTrkIso->at(vIt) < minElecTrkIsoDr03NMinus1) minElecTrkIsoDr03NMinus1 = heepTrkIso->at(vIt);
      	if(eTrkIsoPat->at(vIt) < minElecTrkIsoDr04NMinus1) minElecTrkIsoDr04NMinus1 = eTrkIsoPat->at(vIt);
      }
      
      
/*      
      if(passedAcceptance(vIt) && passedHeepEt(vIt) && passedHeepHadFrac(vIt) && passedHeepDEtaIn(vIt) && passedHeepDPhiIn(vIt) 
        && passedHeepSigmaIEta(vIt) && passedHeepClusterShape(vIt) && passedElecMissingHits(vIt) 
	&& passedElecEcalIso(vIt) && passedMuId(vIt) && passedDiTau(vIt)){      
        nMinus1ElecTrkIsoCounter++;
      	if(heepTrkIso->at(vIt) < minElecTrkIsoDr03NMinus1) minElecTrkIsoDr03NMinus1 = heepTrkIso->at(vIt);
      	if(eTrkIsoPat->at(vIt) < minElecTrkIsoDr04NMinus1) minElecTrkIsoDr04NMinus1 = eTrkIsoPat->at(vIt);
      }
*/      
      if(passedAcceptance(vIt) && electronId && passedMuNHits(vIt) && passedMuIp(vIt) && passedMuPionVeto(vIt)
        && muEcalIso && muTrkIso && diTau && jetCountVar && wjetsVar){
	nMinus1MuChi2Counter++;
	if(muonNormChiSqrd->at(vIt) < minMuonChi2NMinus1) minMuonChi2NMinus1 = muonNormChiSqrd->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && passedMuChi2(vIt) && passedMuIp(vIt) && passedMuPionVeto(vIt)
        && muEcalIso && muTrkIso && diTau && jetCountVar && wjetsVar){
        nMinus1MuNHitsCounter++;
	if(muonValidHits->at(vIt) > maxMuonHitsNMinus1) maxMuonHitsNMinus1 = muonValidHits->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && passedMuChi2(vIt) && passedMuNHits(vIt) && passedMuPionVeto(vIt)
        && muEcalIso && muTrkIso && diTau && jetCountVar && wjetsVar){
        nMinus1MuIpCounter++;
	if(fabs(muonIpVtx->at(vIt)) < minMuonIpNMinus1) minMuonIpNMinus1 = muonIpVtx->at(vIt); 
      }
      
      if(passedAcceptance(vIt) && electronId && passedMuChi2(vIt) && passedMuNHits(vIt) && passedMuIp(vIt)
        && passedMuPionVeto(vIt) && diTau && passedLooseMuIso(vIt) && jetCountVar && wjetsVar){
	nMinus1MuIsoCounter++;
	if(muonTrkIso->at(vIt) + muonEcalIso->at(vIt) < minMuonIsoDr03NMinus1) minMuonIsoDr03NMinus1 = muonEcalIso->at(vIt) + muonTrkIso->at(vIt);
	if(muonTrkIsoDr04->at(vIt) + muonEcalIsoDr04->at(vIt) < minMuonIsoDr04NMinus1) minMuonIsoDr04NMinus1 = muonTrkIsoDr04->at(vIt) + muonEcalIsoDr04->at(vIt);
	if(muonTrkIsoDr05->at(vIt) + muonEcalIsoDr05->at(vIt) < minMuonIsoDr05NMinus1) minMuonIsoDr05NMinus1 = muonTrkIsoDr05->at(vIt) + muonEcalIsoDr05->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && passedMuChi2(vIt) && passedMuNHits(vIt) && passedMuIp(vIt)
        && passedMuPionVeto(vIt) && muTrkIso && passedLooseMuEcalIso(vIt) && diTau && jetCountVar && wjetsVar){
	nMinus1MuEcalIsoCounter++;
        if(muonEcalIso->at(vIt) < minMuonEcalIsoDr03NMinus1) minMuonEcalIsoDr03NMinus1 = muonEcalIso->at(vIt);
        if(muonEcalIsoDr04->at(vIt) < minMuonEcalIsoDr04NMinus1) minMuonEcalIsoDr04NMinus1 = muonEcalIsoDr04->at(vIt);
        if(muonEcalIsoDr05->at(vIt) < minMuonEcalIsoDr05NMinus1) minMuonEcalIsoDr05NMinus1 = muonEcalIsoDr05->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && passedMuChi2(vIt) && passedMuNHits(vIt) && passedMuIp(vIt)
        && passedMuPionVeto(vIt) && muEcalIso && diTau && jetCountVar && wjetsVar){
	nMinus1MuTrkIsoCounter++;
      	if(muonTrkIso->at(vIt) < minMuonTrkIsoDr03NMinus1) minMuonTrkIsoDr03NMinus1 = muonTrkIso->at(vIt);
      	if(muonTrkIsoDr04->at(vIt) < minMuonTrkIsoDr04NMinus1) minMuonTrkIsoDr04NMinus1 = muonTrkIsoDr04->at(vIt);
      	if(muonTrkIsoDr05->at(vIt) < minMuonTrkIsoDr05NMinus1) minMuonTrkIsoDr05NMinus1 = muonTrkIsoDr05->at(vIt);
      }
/*
      if(passedAcceptance(vIt) && electronId && passedMuChi2(vIt) && passedMuNHits(vIt) && passedMuIp(vIt)
        && passedLooseMuEcalIso(vIt)){
	nMinus1MuTrkIsoCounter++;
      	if(muonTrkIso->at(vIt) < minMuonTrkIsoDr03NMinus1) minMuonTrkIsoDr03NMinus1 = muonTrkIso->at(vIt);
      	if(muonTrkIsoDr04->at(vIt) < minMuonTrkIsoDr04NMinus1) minMuonTrkIsoDr04NMinus1 = muonTrkIsoDr04->at(vIt);
      	if(muonTrkIsoDr05->at(vIt) < minMuonTrkIsoDr05NMinus1) minMuonTrkIsoDr05NMinus1 = muonTrkIsoDr05->at(vIt);
      }
      if(passedAcceptance(vIt) && passedHeepId(vIt) && passedMuChi2(vIt) && passedMuNHits(vIt) && passedMuIp(vIt)
        && passedMuCalo2dComp(vIt) && passedMuEcalIso(vIt) && passedDiTau(vIt)){
        nMinus1MuTrkIsoCounter++;
      	if(muonTrkIso->at(vIt) < minMuonTrkIsoDr03NMinus1) minMuonTrkIsoDr03NMinus1 = muonTrkIso->at(vIt);
      	if(muonTrkIsoDr04->at(vIt) < minMuonTrkIsoDr04NMinus1) minMuonTrkIsoDr04NMinus1 = muonTrkIsoDr04->at(vIt);
      	if(muonTrkIsoDr05->at(vIt) < minMuonTrkIsoDr05NMinus1) minMuonTrkIsoDr05NMinus1 = muonTrkIsoDr05->at(vIt);
      }
*/      
      if(passedAcceptance(vIt) && electronId && muonId && diTauCharge && passedBTagVeto(vIt)
        && passedJetSum(vIt) && passedLeptonMETCosDphi(vIt) && passedMet(vIt) && passedMass(vIt) 
	&& passed2DZetaCut(vIt) && jetCountVar && wjetsVar){
        nMinus1ElecMuDPhiCounter++;
	if(eMuCosDPhi->at(vIt) < minCosDPhiNMinus1) minCosDPhiNMinus1 = eMuCosDPhi->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && passedBTagVeto(vIt)
        && passedJetSum(vIt) && passedLeptonMETCosDphi(vIt) && passedMet(vIt) && passedMass(vIt) 
	&& passed2DZetaCut(vIt) && jetCountVar && wjetsVar){
        nMinus1ElecMuChargeCounter++;
	if(eCharge->at(vIt)*muonCharge->at(vIt) < minChargeNMinus1) minChargeNMinus1 = eCharge->at(vIt)*muonCharge->at(vIt); 
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && diTauCharge
       && passedJetSum(vIt) && passedLeptonMETCosDphi(vIt) && passedMet(vIt) && passedMass(vIt) 
       && passed2DZetaCut(vIt) && jetCountVar && wjetsVar){
       nMinus1BTagCounter++;
       if(nBtagsHiEffTrkCnt->at(vIt) < minNBtagsHiEffTrkCntNMinus1) minNBtagsHiEffTrkCntNMinus1 = nBtagsHiEffTrkCnt->at(vIt);
       if(nBtagsHiPurityTrkCnt->at(vIt) < minNBtagsHiPurityTrkCntNMinus1) minNBtagsHiPurityTrkCntNMinus1 = nBtagsHiPurityTrkCnt->at(vIt);
       if(nBTagsHiEffSimpleSecVtx->at(vIt) < minNBtagsHiEffSimpleSecVtxNMinus1) minNBtagsHiEffSimpleSecVtxNMinus1 = nBTagsHiEffSimpleSecVtx->at(vIt);
       if(nBTagsHiPuritySimpleSecVtx->at(vIt) < minNBtagsHiPuritySimpleSecVtxNMinus1) minNBtagsHiPuritySimpleSecVtxNMinus1 = nBTagsHiPuritySimpleSecVtx->at(vIt);
       if(nBTagsCombSecVtx->at(vIt) < minNBtagsCombSecVtxNMinus1) minNBtagsCombSecVtxNMinus1 = nBTagsCombSecVtx->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && diTauCharge
        && passedBTagVeto(vIt) && passedLeptonMETCosDphi(vIt) && passedMet(vIt) && passedMass(vIt) 
	&& passed2DZetaCut(vIt) && jetCountVar && wjetsVar){
        nMinus1JetSumCounter++;
	if(jetSumEt->at(vIt) < minJetSumEtNMinus1) minJetSumEtNMinus1 = jetSumEt->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && diTauCharge
        && passedBTagVeto(vIt) && passedJetSum(vIt) && passedMet(vIt) && passedMass(vIt) 
	&& passed2DZetaCut(vIt) && jetCountVar && wjetsVar){
        nMinus1LeptonMETCosDPhiCounter++;
	if(leptonMetCosDphi->at(vIt) < minLeptonMetCosDPhiNMinus1) minLeptonMetCosDPhiNMinus1 = leptonMetCosDphi->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && diTauCharge
        && passedBTagVeto(vIt) && passedJetSum(vIt) && passedLeptonMETCosDphi(vIt) && passedMass(vIt) 
	&& passed2DZetaCut(vIt) && jetCountVar && wjetsVar){
        nMinus1METCounter++;
	if(mEt->at(vIt) > maxMetNMinus1) maxMetNMinus1 = mEt->at(vIt);
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && diTauCharge
        && passedBTagVeto(vIt) && passedJetSum(vIt) && passedMet(vIt) && passedLeptonMETCosDphi(vIt) 
	&& passedMass(vIt) && jetCountVar && wjetsVar){
	if(eMuPZeta->at(vIt) - _zeta2DSlope*eMuPZetaVis->at(vIt) > max2DZetaCutNMinus1){
          minZetaVisNMinus1 = eMuPZetaVis->at(vIt);
	  maxZetaNMinus1 = eMuPZeta->at(vIt);
	}
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && passedCosDPhi(vIt) && diTauCharge
        && passedBTagVeto(vIt) && passedJetSum(vIt) && passedMet(vIt) && passedLeptonMETCosDphi(vIt) 
	&& passed2DZetaCut(vIt) && passedMass(vIt) && wjetsVar){
	if(nJets->at(vIt) < maxJetsNMinus1) maxJetsNMinus1 = nJets->at(vIt);
	//cout << "MaxNJets \t" <<  maxJetsNMinus1 << "NJets \t" << nJets << "\n";
      }
      
      if(passedAcceptance(vIt) && electronId && muonId && diTau){
	if(muMEtMt->at(vIt) > maxMuMetMass) maxMuMetMass = muMEtMt->at(vIt);
	if(eMEtMt->at(vIt) > maxEMetMass && jetCountVar) maxEMetMass = eMEtMt->at(vIt);
      }
            
      // Do std Relative Eff Here...
      if(_doFactorization){
        if(passedLooseCuts(vIt)) factorizationCounter++;
      }
     
      //make btag plot for correction
      if(ttbarAltControlFact){
	nBTagControlCounter++;
        if(nBtagsHiEffTrkCnt->at(vIt) < minNBtagsHiEffTrkCntControl){
         minNBtagsHiEffTrkCntControl = nBtagsHiEffTrkCnt->at(vIt);
        }
      }
      
      // Get mass shapes from control regions when possible
      if(_label == "qcd" && (qcdControl || qcdAltControl) && passedElecPt(vIt) && passedMet(vIt) 
         && eMuMetMass->at(vIt) > maxTemplateMass) maxTemplateMass = eMuMetMass->at(vIt);
      
      if(_label == "ttbar" &&  (ttbarControl || ttbarAltControl)  && passedLeptonMETCosDphi(vIt)
         && eMuMetMass->at(vIt) > maxTemplateMass) maxTemplateMass = eMuMetMass->at(vIt);

      if(_label == "wjets" && wjetsControl && passedHeepEt(vIt) && passedElecPt(vIt)
         && eMuMetMass->at(vIt) > maxTemplateMass) maxTemplateMass = eMuMetMass->at(vIt);

      //first get background contribtion from data
      // replace LooseEId for LooseEIdForWJets to remove electron Et requirement 
      
      if(_label == "data"){
        if(qcdControl)dataQCDCounter++;
	if(qcdControl && passedMuTrkIso(vIt)) dataQCDMuIsoCounter++;
	if(qcdControl && passedMinMuTrkIso(vIt) && passedLooseMuTrkIso(vIt)) dataQCDAntiMuIsoCounter++;
	
	if(qcdAltControl)dataQCDAltCounter++;
	if(qcdAltControl && passedCharge(vIt))dataQCDAltOSCounter++;
	if(qcdAltControl && passedElecPt(vIt))dataQCDAltElecPtCounter++;
	if(qcdAltControl && passedElecTrkIso(vIt))dataQCDAltElecIsoCounter++;
	if(qcdAltControl && passedMet(vIt) && passedLeptonMETCosDphi(vIt) && passed2DZetaCut(vIt))dataQCDAltMetCounter++;
	//get qcd mass in control region with met cut included
	if((qcdAltControl || qcdAltControl) && passedElecPt(vIt) && passedMet(vIt) && passedLeptonMETCosDphi(vIt) 
	    && passed2DZetaCut(vIt) && eMuMetMass->at(vIt) > maxEMuMetMassDataQCD) maxEMuMetMassDataQCD = eMuMetMass->at(vIt);
	
	// get ttbar contribution
	if(ttbarControl) dataTTbarCounter++;
	if(ttbarControl && passedJetSum(vIt)) dataTTbarJetSumCounter++;
	if(ttbarControl && passedLeptonMETCosDphi(vIt) && passedCosDPhi(vIt)) dataTTbarAllMETCounter++;
	//get ttbar mass in control region
	if((ttbarControl || ttbarAltControl) && passedLeptonMETCosDphi(vIt)
	   && eMuMetMass->at(vIt) > maxEMuMetMassDataTTbar) maxEMuMetMassDataTTbar = eMuMetMass->at(vIt);
	
	
	if(ttbarAltControl) dataTTbarAltCounter++;
	if(ttbarAltControl && passedBTagVeto(vIt)) dataTTbarAltBtagVetoCounter++;	
        
	//get Wjets contribution via W->munu transverse mass
	if(wjetsControl) dataWJetsCounter++;
	if(wjetsControl && passedHeepEt(vIt) && passedElecPt(vIt)) dataWJetsElecPtCounter++;
	if(wjetsControl && passedCharge(vIt)) dataWJetsOSCounter++;
	if(wjetsControl && passedCosDPhi(vIt) && passed2DZetaCut(vIt)) dataWJetsMETCounter++;
	//get wjets mass in control region
	if(wjetsControl && passedHeepEt(vIt) && passedElecPt(vIt) && passedCosDPhi(vIt) && passed2DZetaCut(vIt)
	   && eMuMetMass->at(vIt) > maxEMuMetMassDataWJets) maxEMuMetMassDataWJets = eMuMetMass->at(vIt);
	
        
	if(wjetsAltControl) dataWJetsAltCounter++;
	if(wjetsAltControl && passedElecTrkIso(vIt)) dataWJetsAltElecIsoCounter++;
	if(wjetsAltControl && passedMinElecTrkIso(vIt) && passedLooseElecTrkIso(vIt)) dataWJetsAltElecAntiIsoCounter++;
	if(wjetsAltControl && muMetRegion) dataWJetsAltMTCounter++;
/*      
	//get Ztautau
        if(ztautauControl) dataZTauTauCounter++;
	if(ztautauControl && passedElecPt(vIt) dataZTauTauElecPtCounter++;
	if(ztautauControl && passedMuonPt(vIt) dataZTauTauMuonPtCounter++;
	if(ztautauControl && passedMet(vIt) && passedLeptonMETCosDphi(vIt) 
	  && passed2DZetaCut(vIt)) dataZTauTauMETCounter++
*/	
      }
      
      if(eMuDeltaR->at(vIt) > maxElecMuDeltaR) maxElecMuDeltaR = eMuDeltaR->at(vIt);
      if(!passedDiTauDeltaR(vIt)) continue;
      deltaRCounter++;
      
      if(ePt->at(vIt) > maxElecPt) maxElecPt = ePt->at(vIt);
      if(!passedElecPt(vIt)) continue;
      elecPtCounter++;
      
      if(fabs(eEta->at(vIt)) < minElecEta) minElecEta = eEta->at(vIt);
      if(!passedElecEta(vIt)) continue;

      elecEtaCounter++;
      
      if(muonPt->at(vIt) > maxMuPt) maxMuPt = muonPt->at(vIt);
      if(!passedMuonPt(vIt)) continue; 
      muonPtCounter++;
      
      if(fabs(muonEta->at(vIt)) < minMuonEta) minMuonEta = muonEta->at(vIt);
      if(!passedMuonEta(vIt)) continue;
      muonEtaCounter++;
  
      if(passedAcceptance(vIt)) accCounter++;
/*      
      // EId HEEP ... 
      if(!passedHeepEt(vIt))continue;
      heepEtCounter++;
      
      if(!passedHeepEcalDriven(vIt)) continue;
      heepEcalDrivenCounter++;  
      
      if(eHOverEm->at(vIt) < minElecHOverEm) minElecHOverEm = eHOverEm->at(vIt);
      if(!passedHeepHadFrac(vIt)) continue;
      heepHadCounter++;  
      
      if(fabs(eDeltaEtaIn->at(vIt)) < minElecDEtaIn) minElecDEtaIn = eDeltaEtaIn->at(vIt);
      if(!passedHeepDEtaIn(vIt)) continue;
      heepDEtaInCounter++;
      
      if(fabs(eDeltaPhiIn->at(vIt)) < minElecDphiIn) minElecDphiIn = eDeltaPhiIn->at(vIt);
      if(!passedHeepDPhiIn(vIt)) continue;
      heepDPhiInCounter++;
      
      if(eSigmaIEtaIEta->at(vIt) < minElecSigmaIEta) minElecSigmaIEta = eSigmaIEtaIEta->at(vIt);
      if(!passedHeepSigmaIEta(vIt)) continue;
      heepSigmaIEtaCounter++;   
      
      if(eSCE2x5->at(vIt)/eSCE5x5->at(vIt) > maxElec2x5Ratio) maxElec2x5Ratio = eSCE2x5->at(vIt)/eSCE5x5->at(vIt);
      if(eSCE1x5->at(vIt)/eSCE5x5->at(vIt) > maxElec1x5Ratio) maxElec1x5Ratio = eSCE1x5->at(vIt)/eSCE5x5->at(vIt);
      if(!passedHeepClusterShape(vIt)) continue; 
      heepClusterRatioCounter++;
      
      if(eMissingHits->at(vIt) < minElecNHits) minElecNHits = eMissingHits->at(vIt);
      if(!passedElecMissingHits(vIt)) continue;
      heepEMissingHitsCounter++;
*/      
/*   No Isolations from heep
      if(!passedHeepEcalHadIso(vIt)) continue;
      heepEcalHadIsoCounter++;	 
      if(!passedHeepCaloIsoDepth2(vIt))continue; 
      heepHadIsoDepth2++;  
      if(!passedHeepTrkIso(vIt)) continue;
      heepTrkIsoCounter++;
*/  

// use heep defaults  
      // EId ...
      if(!passedHeepEcalDriven(vIt)) continue;
      eEcalDrivenCounter++;  
      
      if(eHOverEm->at(vIt) < minElecHOverEm) minElecHOverEm = eHOverEm->at(vIt);
      if(!passedElecEmFrac(vIt)) continue;
      elecEmFracCounter++; 
      
      if(fabs(eDeltaEtaIn->at(vIt)) < minElecDEtaIn) minElecDEtaIn = eDeltaEtaIn->at(vIt);
      if(!passedDEtaIn(vIt)) continue;
      eDEtaInCounter++;
      
      if(fabs(eDeltaPhiIn->at(vIt)) < minElecDphiIn) minElecDphiIn = eDeltaPhiIn->at(vIt);
      if(!passedDPhiIn(vIt)) continue;
      eDPhiInCounter++;
      
      if(eSigmaIEtaIEta->at(vIt) < minElecSigmaIEta) minElecSigmaIEta = eSigmaIEtaIEta->at(vIt);
      if(!passedSigmaIEtaIEta(vIt)) continue;
      eSigmaIEtaCounter++;
      
      if(eSCE2x5->at(vIt)/eSCE5x5->at(vIt) > maxElec2x5Ratio) maxElec2x5Ratio = eSCE2x5->at(vIt)/eSCE5x5->at(vIt);
      if(eSCE1x5->at(vIt)/eSCE5x5->at(vIt) > maxElec1x5Ratio) maxElec1x5Ratio = eSCE1x5->at(vIt)/eSCE5x5->at(vIt);
      if(!passedClusterShape(vIt)) continue;
      eClusterRatioCounter++;
      
      if(eMissingHits->at(vIt) < minElecNHits) minElecNHits = eMissingHits->at(vIt);
      if(!passedElecMissingHits(vIt)) continue;
      eMissingHitsCounter++;
            
      if(heepTrkIso->at(vIt) < minElecTrkIsoDr03) minElecTrkIsoDr03 = heepTrkIso->at(vIt);
      if(eTrkIsoPat->at(vIt) < minElecTrkIsoDr04) minElecTrkIsoDr04 = eTrkIsoPat->at(vIt);
      if(passedElecTrkIso(vIt)) elecTrkIsoTightCounter++;
      if(!elecTrkIso) continue;
      eTrkIsoCounter++;
      
      if(heepEcalIso->at(vIt) < minElecEcalIsoDr03) minElecEcalIsoDr03 = heepEcalIso->at(vIt);
      if(eEcalIsoPat->at(vIt) < minElecEcalIsoDr04) minElecEcalIsoDr04 = eEcalIsoPat->at(vIt);
      if(passedElecEcalIso(vIt)) elecEcalIsoTightCounter++;
      if(!elecEcalIso)continue;
      eEcalIsoCounter++;
      
      //if(heepTrkIso->at(vIt)+heepEcalIso->at(vIt) < minElecIsoDr03) minElecIsoDr03 = heepTrkIso->at(vIt)+heepEcalIso->at(vIt);
      //if(eTrkIsoPat->at(vIt)+eEcalIsoPat->at(vIt) < minElecIsoDr04) minElecIsoDr04 = eTrkIsoPat->at(vIt)+eEcalIsoPat->at(vIt);
      //if(!elecIso) continue;
      //eIsoCounter++;
      
      //get tight counter for factorization
      //if(passedElecIso(vIt)) elecIsoTightCounter++;
      //if(passedAcceptance(vIt) && passedHeepId(vIt)) eHeepIdCounter++;	     
        
      if(passedAcceptance(vIt) && passedEId(vIt)) eIdCounter++;  
  
      // Muon Id ...
      
      if(muonNormChiSqrd->at(vIt) < minMuonChi2) minMuonChi2 = muonNormChiSqrd->at(vIt);
      if(!passedMuChi2(vIt)) continue;
      muChi2Counter++;
      
      if(fabs(muonIpVtx->at(vIt)) < minMuonIp) minMuonIp = muonIpVtx->at(vIt);
      if(!passedMuIp(vIt)) continue;
      muIpCounter++;

            
      if(muonValidHits->at(vIt) > maxMuonHits) maxMuonHits = muonValidHits->at(vIt);
      if(!passedMuNHits(vIt)) continue;
      muNHitsCounter++;
      
      //if(!passedMuCalo2dComp(vIt)) continue;
      if(!passedMuPionVeto(vIt)) continue;
      mu2DCompCounter++;
      
      if(muonTrkIso->at(vIt) < minMuonTrkIsoDr03) minMuonTrkIsoDr03 = muonTrkIso->at(vIt);
      if(muonTrkIsoDr04->at(vIt) < minMuonTrkIsoDr04) minMuonTrkIsoDr04 = muonTrkIsoDr04->at(vIt);
      if(muonTrkIsoDr05->at(vIt) < minMuonTrkIsoDr05) minMuonTrkIsoDr05 = muonTrkIsoDr05->at(vIt);
      if(passedMuTrkIso(vIt)) muTrkIsoTightCounter++;
      if(!muTrkIso) continue;
      muTrkIsoCounter++;
      
      if(muonEcalIso->at(vIt) < minMuonEcalIsoDr03) minMuonEcalIsoDr03 = muonEcalIso->at(vIt);
      if(muonEcalIsoDr04->at(vIt) < minMuonEcalIsoDr04) minMuonEcalIsoDr04 = muonEcalIsoDr04->at(vIt);
      if(muonEcalIsoDr05->at(vIt) < minMuonEcalIsoDr05) minMuonEcalIsoDr05 = muonEcalIsoDr05->at(vIt);
      if(passedMuEcalIso(vIt)) muEcalIsoTightCounter++;
      if(!muEcalIso) continue;
      muEcalIsoCounter++;
      
      //if(muonEcalIso->at(vIt) + muonTrkIso->at(vIt) < minMuonIsoDr03) minMuonIsoDr03 = muonEcalIso->at(vIt)+muonTrkIso->at(vIt);
      //if(muonEcalIsoDr04->at(vIt) + muonTrkIsoDr04->at(vIt)< minMuonIsoDr04) minMuonIsoDr04 = muonEcalIsoDr04->at(vIt)+ muonTrkIsoDr04->at(vIt);
      //if(muonEcalIsoDr05->at(vIt) + muonTrkIsoDr05->at(vIt) < minMuonIsoDr05) minMuonIsoDr05 = muonEcalIsoDr05->at(vIt) + muonTrkIsoDr05->at(vIt) < minMuonIsoDr05;
      //if(!muonIso) continue;
      //muIsoCounter++;
      //get tight counter for factorization
      //if(passedMuIso(vIt)) muIsoTightCounter++;
      
      //if(passedAcceptance(vIt) && passedHeepId(vIt) && passedMuId(vIt)) muIdCounter++;
      if(passedAcceptance(vIt) && passedEId(vIt) && passedMuId(vIt)) muIdCounter++;
        
      // DiTau ...
      //make 2d plots/cand before met cuts
      
      if(nBtagsHiEffTrkCnt->at(vIt) < minNBtagsHiEffTrkCntFor2DPlot){   
        minNBtagsHiEffTrkCntFor2DPlot = nBtagsHiEffTrkCnt->at(vIt);     
        minLeptonMetCosDPhiForPlotVsNBtags = leptonMetCosDphi->at(vIt);
      } 							        
      								        
      if(jetSumEt->at(vIt) < minJetSumEtFor2DPlot){		        
        minJetSumEtFor2DPlot = jetSumEt->at(vIt);		        
        minLeptonMetCosDPhiForPlotVsEtSum = leptonMetCosDphi->at(vIt);
      } 							        
      								        
      //if(passed2DZetaCut(vIt)) zeta2DCounter++;
      if(eMuPZeta->at(vIt) - _zeta2DSlope*eMuPZetaVis->at(vIt) > max2DZetaCut){
        minZetaVis = eMuPZetaVis->at(vIt);
	maxZeta = eMuPZeta->at(vIt);
      }

      if(_setRegion && _whichRegion == "ttbarJetCount" && passedJetCount(vIt)) continue;
      if(_setRegion && _whichRegion.find("wjets") != string::npos && !muMetRegion) continue;
      if(_setRegion && _whichRegion == "ztt" && muMEtMt->at(vIt) > 40.) continue;
      
      if(eMuCosDPhi->at(vIt) < minCosDPhi) minCosDPhi = eMuCosDPhi->at(vIt);
      if(!passedCosDPhi(vIt)) continue;
      eMuDPhiCounter++;	    
      
      if(eCharge->at(vIt)*muonCharge->at(vIt) < minCharge) minCharge = eCharge->at(vIt)*muonCharge->at(vIt); 
      if(!diTauCharge) continue;
      eMuChargeCounter++;
      if(passedCharge(vIt)) chargeTightCounter++;
      
      if(nBtagsHiEffTrkCnt->at(vIt) < minNBtagsHiEffTrkCnt) minNBtagsHiEffTrkCnt = nBtagsHiEffTrkCnt->at(vIt);
      if(nBtagsHiPurityTrkCnt->at(vIt) < minNBtagsHiPurityTrkCnt) minNBtagsHiPurityTrkCnt = nBtagsHiPurityTrkCnt->at(vIt);
      if(nBTagsHiEffSimpleSecVtx->at(vIt) < minNBtagsHiEffSimpleSecVtx) minNBtagsHiEffSimpleSecVtx = nBTagsHiEffSimpleSecVtx->at(vIt);
      if(nBTagsHiPuritySimpleSecVtx->at(vIt) < minNBtagsHiPuritySimpleSecVtx) minNBtagsHiPuritySimpleSecVtx = nBTagsHiPuritySimpleSecVtx->at(vIt);
      if(nBTagsCombSecVtx->at(vIt) < minNBtagsCombSecVtx) minNBtagsCombSecVtx = nBTagsCombSecVtx->at(vIt);
      if(!passedBTagVeto(vIt)) continue;
      btagCounter++;
      
      if(jetSumEt->at(vIt) < minJetSumEt) minJetSumEt = jetSumEt->at(vIt);
      if(jetMETSumEt->at(vIt) < minJetMetSumEt) minJetMetSumEt = jetMETSumEt->at(vIt);
      if(!passedJetSum(vIt)) continue;
      jetSumCounter++;
      
      if(mEt->at(vIt) > maxMet) maxMet = mEt->at(vIt);
      if(!passedMet(vIt)) continue;
      metCounter++;
      
      if(leptonMetCosDphi->at(vIt) < minLeptonMetCosDPhi) minLeptonMetCosDPhi = leptonMetCosDphi->at(vIt);
      if(!passedLeptonMETCosDphi(vIt)) continue;
      leptonMETCosDPhiCounter++;
      
      if(!passed2DZetaCut(vIt))continue;
      zeta2DCounter++;
      
      if(eMuMetMass->at(vIt) > maxMetMass)  maxMetMass = eMuMetMass->at(vIt);     
      if(eMuMass->at(vIt) > maxMass) maxMass = eMuMass->at(vIt);
      if(!passedMass(vIt)) continue;
      massCounter++;
/* 
      if(massCounter > 1) cout << "More than one final pair\t" << _label <<"\t" << massCounter << endl;
      //get surviving event info
      if(eMuMetMass->at(vIt) > 400. && _label == "data"){
        theEventNumberVector.push_back(eventNumber);
        theRunNumberVector.push_back(runNumber);
        cout << "******************************************"
	     << "\n Mass > 400 Info:\n"
	     << "\n Event Number:\t" << eventNumber 
	     << "\n RunNumber:\t" << runNumber
	     << "\n ePt:\t\t" << ePt->at(vIt)
	     << "\n muPt:\t\t" << muonPt->at(vIt)
	     << "\n MET:\t\t" << mEt->at(vIt)
	     << "\n EMuMass:\t" << eMuMass->at(vIt)
	     << "\n EMuMetMass:\t" << eMuMetMass->at(vIt)
	     << "\n******************************************\n";
      }
*/      
  
      //if(passedAcceptance(vIt) && passedHeepId(vIt) && passedMuId(vIt) && passedDiTau(vIt)){
      if(passedAcceptance(vIt) && passedEId(vIt) && passedMuId(vIt) && passedDiTau(vIt)){
        diTauCounter++;
	if(_label == "zprime" || _label == "ztautau" || _label == "ww" || _label == "wz" || _label == "data"){
	  if(eMuMetMass->at(vIt) > maxTemplateMass) maxTemplateMass = eMuMetMass->at(vIt);
        //if(!_doFactorization) FillTemplate(vIt);
          //Template->Fill(eMuMetMass->at(vIt));
	}
      }
    }
    GetEventCounters();
    //GetReWeightedEvents();
    
    //Fill histos ...
    _elecMuDeltaR->Fill(maxElecMuDeltaR, _eventPUWeight);
    _elecPt->Fill(maxElecPt, _eventPUWeight);
    _elecEta->Fill(minElecEta, _eventPUWeight);
    _muonPt->Fill(maxMuPt, _eventPUWeight);
    _muonEta->Fill(minMuonEta, _eventPUWeight);

/*    
    _elecMuDeltaRReWeighted->Fill(maxElecMuDeltaR,PDFWweights->at(0)*ISRGammaWeight*ISRGluonWeight*FSRWeight*eventPUWeight);	
    _elecPtReWeighted->Fill(maxElecPt,PDFWweights->at(0)*ISRGammaWeight*ISRGluonWeight*FSRWeight*eventPUWeight);  	
    _elecEtaReWeighted->Fill(minElecEta,PDFWweights->at(0)*ISRGammaWeight*ISRGluonWeight*FSRWeight*eventPUWeight); 	
    _muonPtReWeighted->Fill(maxMuPt,PDFWweights->at(0)*ISRGammaWeight*ISRGluonWeight*FSRWeight*eventPUWeight);  	
    _muonEtaReWeighted->Fill(minMuonEta,PDFWweights->at(0)*ISRGammaWeight*ISRGluonWeight*FSRWeight*eventPUWeight);
*/  
    _elecHOverE->Fill(minElecHOverEm, _eventPUWeight);       
    _elecDEtaIn->Fill(minElecDEtaIn, _eventPUWeight);        
    _elecDPhiIn->Fill(minElecDphiIn, _eventPUWeight);	
    _elecSigmaIEtaIEta->Fill(minElecSigmaIEta, _eventPUWeight); 
    _elec2x5Over5x5->Fill(maxElec2x5Ratio, _eventPUWeight);	
    _elec1x5Over5x5->Fill(maxElec1x5Ratio, _eventPUWeight);	
    _elecMissingHits->Fill(minElecNHits, _eventPUWeight);	
    _elecEcalIsoDr03->Fill(minElecEcalIsoDr03, _eventPUWeight);	
    _elecTrkIsoDr03->Fill(minElecTrkIsoDr03, _eventPUWeight);	
    _elecEcalIsoDr04->Fill(minElecEcalIsoDr04, _eventPUWeight);	
    _elecTrkIsoDr04->Fill(minElecTrkIsoDr04, _eventPUWeight);	
    _elecIsoDr03->Fill(minElecIsoDr03, _eventPUWeight);   
    _elecIsoDr04->Fill(minElecIsoDr04, _eventPUWeight);   
  
    _muonChi2NdF->Fill(minMuonChi2, _eventPUWeight);	
    _muonIp->Fill(minMuonIp, _eventPUWeight);		
    _muonHits->Fill(maxMuonHits, _eventPUWeight);  	
    _muonEcalIsoDr03->Fill(minMuonEcalIsoDr03, _eventPUWeight);	
    _muonTrkIsoDr03->Fill(minMuonTrkIsoDr03, _eventPUWeight);	
    _muonEcalIsoDr04->Fill(minMuonEcalIsoDr04, _eventPUWeight);	
    _muonTrkIsoDr04->Fill(minMuonTrkIsoDr04, _eventPUWeight);	
    _muonEcalIsoDr05->Fill(minMuonEcalIsoDr05, _eventPUWeight);	
    _muonTrkIsoDr05->Fill(minMuonTrkIsoDr05, _eventPUWeight);	
    _muonIsoDr03->Fill(minMuonIsoDr03, _eventPUWeight);   
    _muonIsoDr04->Fill(minMuonIsoDr04, _eventPUWeight);   
    _muonIsoDr05->Fill(minMuonIsoDr05, _eventPUWeight);   
  
    _elecMuCosDPhi->Fill(minCosDPhi, _eventPUWeight);	
    _elecMuCharge->Fill(minCharge, _eventPUWeight);	
    
    _btagCountHiEffTrkCnt->Fill(minNBtagsHiEffTrkCnt, _eventPUWeight); 	
    _btagCountHiPurityTrkCnt->Fill(minNBtagsHiPurityTrkCnt, _eventPUWeight); 	
    _btagCountSimpleSecVtx->Fill(minNBtagsHiEffSimpleSecVtx, _eventPUWeight); 	
    _btagCountHiPuritySimpleSecVtx->Fill(minNBtagsHiPuritySimpleSecVtx, _eventPUWeight); 	
    _btagCountCombSecVtx->Fill(minNBtagsCombSecVtx, _eventPUWeight); 	
    
    _btagJetTrkCntHiEff->Fill(jetBtagHiDiscByTrkCntHiEff, _eventPUWeight);
    _btagJetTrkCntHiPurity->Fill(jetBtagHiDiscByTrkCntHiPurity, _eventPUWeight);
    _btagJetSimpleSecVtxHiEff->Fill(jetBtagHiDiscBySimpleSecVtxHiEff, _eventPUWeight);
    _btagJetSimpleSecVtxHiPurity->Fill(jetBtagHiDiscBySimpleSecVtxHiPurity, _eventPUWeight);
    _btagJetCombSecVtx->Fill(jetBtagHiDiscByCombSecVtx, _eventPUWeight);

        
    _jetSumEt->Fill(minJetSumEt, _eventPUWeight);  	
    _jetMetSumEt->Fill(minJetMetSumEt, _eventPUWeight);  	
    
    _leptonMETCosDPhi->Fill(minLeptonMetCosDPhi, _eventPUWeight);
    
    _pzeta2D->Fill(minZetaVis, maxZeta, _eventPUWeight);
    _pzeta->Fill(maxZeta, _eventPUWeight);
    _pzetaVis->Fill(minZetaVis, _eventPUWeight);
    
    _leptonMetDPhiVsJetSum->Fill(minJetSumEtFor2DPlot, minLeptonMetCosDPhiForPlotVsEtSum, _eventPUWeight);
    _leptonMetDPhiVsNBtag->Fill(minNBtagsHiEffTrkCntFor2DPlot, minLeptonMetCosDPhiForPlotVsNBtags, _eventPUWeight);
    
    _met->Fill(maxMet, _eventPUWeight);		
    _elecMuMass->Fill(maxMass, _eventPUWeight);
    _elecMuMetMass->Fill(maxMetMass, _eventPUWeight);
    _elecMetMass->Fill(maxEMetMass, _eventPUWeight);
    _muonMetMass->Fill(maxMuMetMass, _eventPUWeight);
    _nJetsElecPt->Fill(countNJets(), _eventPUWeight);
    
    _elecMuMetMassDataQCD->Fill(maxEMuMetMassDataQCD, _eventPUWeight);
    _elecMuMetMassDataTTbar->Fill(maxEMuMetMassDataTTbar, _eventPUWeight);
    _elecMuMetMassDataWJets->Fill(maxEMuMetMassDataWJets, _eventPUWeight);
   
    Template->Fill(maxTemplateMass, _eventPUWeight);
    _btagCountHiEffTrkCntControl->Fill(minNBtagsHiEffTrkCntControl, _eventPUWeight);
    //fill N-1 plots  
    _elecMuDeltaRNMinus1->Fill(maxElecMuDeltaRNMinus1, _eventPUWeight);       
    _elecPtNMinus1->Fill(maxElecPtNMinus1, _eventPUWeight);	
    _elecEtaNMinus1->Fill(minElecEtaNMinus1, _eventPUWeight);		
    _muonPtNMinus1->Fill(maxMuPtNMinus1, _eventPUWeight);     
    _muonEtaNMinus1->Fill(minMuonEtaNMinus1, _eventPUWeight);		
  
    _elecHOverENMinus1->Fill(minElecHOverEmNMinus1, _eventPUWeight);	     
    _elecDEtaInNMinus1->Fill(minElecDEtaInNMinus1, _eventPUWeight);	    
    _elecDPhiInNMinus1->Fill(minElecDphiInNMinus1, _eventPUWeight);	    
    _elecSigmaIEtaIEtaNMinus1->Fill(minElecSigmaIEtaNMinus1, _eventPUWeight); 
    _elec2x5Over5x5NMinus1->Fill(maxElec2x5RatioNMinus1, _eventPUWeight);    
    _elec1x5Over5x5NMinus1->Fill(maxElec1x5RatioNMinus1, _eventPUWeight);    
    _elecMissingHitsNMinus1->Fill(minElecNHitsNMinus1, _eventPUWeight);   
    _elecEcalIsoDr04NMinus1->Fill(minElecEcalIsoDr04NMinus1, _eventPUWeight);		
    _elecTrkIsoDr04NMinus1->Fill(minElecTrkIsoDr04NMinus1, _eventPUWeight);	
    _elecEcalIsoDr03NMinus1->Fill(minElecEcalIsoDr03NMinus1, _eventPUWeight);		
    _elecTrkIsoDr03NMinus1->Fill(minElecTrkIsoDr03NMinus1, _eventPUWeight);	
    _elecIsoDr03NMinus1->Fill(minElecIsoDr03NMinus1, _eventPUWeight);
    _elecIsoDr04NMinus1->Fill(minElecIsoDr04NMinus1, _eventPUWeight);
 
    _muonChi2NdFNMinus1->Fill(minMuonChi2NMinus1, _eventPUWeight);	  
    _muonIpNMinus1->Fill(minMuonIpNMinus1, _eventPUWeight);	
    _muonHitsNMinus1->Fill(maxMuonHitsNMinus1, _eventPUWeight); 	
    _muonEcalIsoDr03NMinus1->Fill(minMuonEcalIsoDr03NMinus1, _eventPUWeight);		
    _muonTrkIsoDr03NMinus1->Fill(minMuonTrkIsoDr03NMinus1, _eventPUWeight);	
    _muonEcalIsoDr04NMinus1->Fill(minMuonEcalIsoDr04NMinus1, _eventPUWeight);		
    _muonTrkIsoDr04NMinus1->Fill(minMuonTrkIsoDr04NMinus1, _eventPUWeight);	
    _muonEcalIsoDr05NMinus1->Fill(minMuonEcalIsoDr05NMinus1, _eventPUWeight);		
    _muonTrkIsoDr05NMinus1->Fill(minMuonTrkIsoDr05NMinus1, _eventPUWeight);	
    _muonIsoDr03NMinus1->Fill(minMuonIsoDr03NMinus1, _eventPUWeight);	  
    _muonIsoDr04NMinus1->Fill(minMuonIsoDr04NMinus1, _eventPUWeight);	  
    _muonIsoDr05NMinus1->Fill(minMuonIsoDr05NMinus1, _eventPUWeight);	  
  
    _elecMuCosDPhiNMinus1->Fill(minCosDPhiNMinus1, _eventPUWeight);	
    _elecMuChargeNMinus1->Fill(minChargeNMinus1, _eventPUWeight);	
    _jetSumEtNMinus1->Fill(minJetSumEtNMinus1, _eventPUWeight); 	
    _metNMinus1->Fill(maxMetNMinus1, _eventPUWeight);		
    _metStdNMinus1->Fill(maxMetNMinus1);		
    _leptonMETCosDPhiNMinus1->Fill(minLeptonMetCosDPhiNMinus1, _eventPUWeight);
  
    _btagCountHiEffTrkCntNMinus1->Fill(minNBtagsHiEffTrkCntNMinus1, _eventPUWeight);
    _pzetaNMinus1->Fill(maxZetaNMinus1, _eventPUWeight);
    _pzetaVisNMinus1->Fill(minZetaVisNMinus1, _eventPUWeight);
    _pzetaPZetaVisNMinus1->Fill(maxZetaNMinus1 - 1.25*minZetaVisNMinus1,_eventPUWeight);
    _nJetsNMinus1->Fill(maxJetsNMinus1, _eventPUWeight);
    
 }
  nEvents = (int)puWeightedEvents;
  GetReport();
  WriteOutFile();

}

void EMuAnalyzer::GetReport(){
  
  std::ofstream outLogFile;
  outLogFile.open(_outLogFile.c_str());
  
  pair<double,pair<double,double> > efficiency;
  //pair<double,pair<double,double> > deltaREff;
  pair<double,pair<double,double> > accEff;
  
  survivingEvents = nMass;
  efficiency = GetEfficiency(survivingEvents, nEvents);
  //deltaREff = GetEfficiency(nDeltaR, nEvents);
  
  accEff = GetEfficiency(nAcc, nEvents);
  pair<double, double> survivingEventsAtLumi = GetSurvivingEvents(nSkimmedAtLumi, efficiency);
  
  // get normalized events in ttbar btag control region
  pair<double,pair<double,double> > efficiencyTTbarControl = GetEfficiency(nBTagControl, nEvents);
  pair<double, double> survivingEventsAtLumiBTagControl = GetSurvivingEvents(nSkimmedAtLumi, efficiencyTTbarControl);
  //
  
  outLogFile.setf(ios::fixed,ios::floatfield);
  outLogFile << " Efficiencies for " << _label << "\n";
  pair<double, double> nFactorizedEvents;
  if(_doFactorization){
    pair<double,double> elecTrkIsoLooseEff = make_pair<double, double>(GetEfficiency(nElecTrkIso,nEMissHits).first,
       max(GetEfficiency(nElecTrkIso,nEMissHits).second.first,
       GetEfficiency(nElecIso,nEMissHits).second.second));
    pair<double, double> elecTrkIsoTightEff = make_pair<double, double>(GetEfficiency(nTightElecTrkIso,nEMissHits).first,
       max(GetEfficiency(nTightElecTrkIso,nEMissHits).second.first,
       GetEfficiency(nTightElecTrkIso,nEMissHits).second.second));
    pair<double,double> elecEcalIsoLooseEff = make_pair<double, double>(GetEfficiency(nElecEcalIso, nElecTrkIso).first,
       max(GetEfficiency(nElecEcalIso, nElecTrkIso).second.first,
       GetEfficiency(nElecEcalIso, nElecTrkIso).second.second));
    pair<double, double> elecEcalIsoTightEff = make_pair<double, double>(GetEfficiency(nTightElecEcalIso,nElecTrkIso).first,
       max(GetEfficiency(nTightElecEcalIso,nElecTrkIso).second.first,
       GetEfficiency(nTightElecEcalIso,nElecTrkIso).second.second));
    pair<double, double> muTrkIsoLooseEff = make_pair<double, double>(GetEfficiency(nMuTrkIso,nMuComp).first,
       max(GetEfficiency(nMuTrkIso,nMuComp).second.first,
       GetEfficiency(nMuTrkIso,nMuComp).second.second));
    pair<double, double> muIsoTrkTightEff = make_pair<double, double>(GetEfficiency(nTightMuonTrkIso,nMuComp).first,
       max(GetEfficiency(nTightMuonTrkIso,nMuComp).second.first,
       GetEfficiency(nTightMuonTrkIso,nMuComp).second.second));
    
    pair<double, double> muEcalIsoLooseEff = make_pair<double, double>(GetEfficiency(nMuCaloIso,nMuTrkIso).first,
       max(GetEfficiency(nMuCaloIso,nMuTrkIso).second.first,
       GetEfficiency(nMuCaloIso,nMuTrkIso).second.second));
    pair<double, double> muIsoEcalTightEff = make_pair<double, double>(GetEfficiency(nTightMuonEcalIso,nMuTrkIso).first,
       max(GetEfficiency(nTightMuonEcalIso,nMuTrkIso).second.first,
       GetEfficiency(nTightMuonEcalIso,nMuTrkIso).second.second));
    
    pair<double, double> chargeLooseEff = make_pair<double, double>(GetEfficiency(nDiTauCharge,nCosDPhi).first,
       max(GetEfficiency(nDiTauCharge,nCosDPhi).second.first,
       GetEfficiency(nDiTauCharge,nCosDPhi).second.second));
    pair<double, double> chargeTightEff = make_pair<double, double>(GetEfficiency(nTightCharge,nCosDPhi).first,
       max(GetEfficiency(nTightCharge,nCosDPhi).second.first,
       GetEfficiency(nTightCharge,nCosDPhi).second.second));
       
    pair<double, double> nAfterElecTrkIso = GetFactorizedEvents(survivingEventsAtLumi,elecTrkIsoLooseEff, elecTrkIsoTightEff);
    pair<double, double> nAfterElecEcalIso = GetFactorizedEvents(nAfterElecTrkIso,elecEcalIsoLooseEff, elecEcalIsoTightEff);
    pair<double, double> nAfterMuTrkIso = GetFactorizedEvents(nAfterElecEcalIso,muTrkIsoLooseEff,muIsoTrkTightEff);
    pair<double, double> nAfterMuCaloIso = GetFactorizedEvents(nAfterMuTrkIso,muEcalIsoLooseEff,muIsoEcalTightEff);
    pair<double, double> nAfterCharge = GetFactorizedEvents(nAfterMuCaloIso,chargeLooseEff,chargeTightEff);
    survivingEventsAtLumi = nAfterCharge;
       
    pair<double, double> nAfterElecTrkIsoBtagControl = GetFactorizedEvents(survivingEventsAtLumiBTagControl,elecTrkIsoLooseEff, elecTrkIsoTightEff);
    pair<double, double> nAfterElecEcalIsoBtagControl = GetFactorizedEvents(nAfterElecTrkIsoBtagControl,elecEcalIsoLooseEff, elecEcalIsoTightEff);
    pair<double, double> nAfterMuTrkIsoBtagControl = GetFactorizedEvents(nAfterElecEcalIsoBtagControl,muTrkIsoLooseEff,muIsoTrkTightEff);
    pair<double, double> nAfterMuCaloIsoBtagControl = GetFactorizedEvents(nAfterMuTrkIsoBtagControl,muEcalIsoLooseEff,muIsoEcalTightEff);
    pair<double, double> nAfterChargeBtagControl = GetFactorizedEvents(nAfterMuCaloIsoBtagControl,chargeLooseEff,chargeTightEff);
    survivingEventsAtLumiBTagControl = nAfterChargeBtagControl;
  }
  if(_label == "data"){
    survivingEventsAtLumi = make_pair<double, double>(survivingEvents, sqrt(survivingEvents));
    survivingEventsAtLumiBTagControl = make_pair<double, double>(nBTagControl, sqrt(nBTagControl));
  }
  //if(_label == "data") survivingEventsAtLumi = make_pair<double, double>(0., 0); // hide the value
  
  unsigned int deltaRDeno = nEvents;
  if(_label == "zprime"){
    outLogFile 
       << "\n EMatched: \t \t" 		<< PrintEfficiency(elecIsMatchedInEvent, nEvents)
       << "\n muonMatched: \t \t"  	<< PrintEfficiency(muonIsMatchedInEvent, elecIsMatchedInEvent)
       << "\n";
    deltaRDeno = muonIsMatchedInEvent;
  }
  //pdfWeightedEvents = GetCentralPDF(pdfWeightsSum);
  outLogFile
       << "\n Trigger: \t \t"		<< PrintEfficiency(nTriggers, nEvents) 
       << "\n deltaR: \t \t" 		<< PrintEfficiency(nDeltaR, nTriggers)
       //<< "\n deltaR: \t \t" 		<< PrintEfficiency(nDeltaR, deltaRDeno)
       << "\n ePt: \t \t \t"  		<< PrintEfficiency(nElePt,nDeltaR)
       << "\n eEta: \t \t \t" 		<< PrintEfficiency(nElecEta,nElePt)
       << "\n muonPt: \t \t"  		<< PrintEfficiency(nMuonPt,nElecEta)
       << "\n muonEta: \t \t" 		<< PrintEfficiency(nMuonEta,nMuonPt)
       << "\n eEcalDrvn: \t \t" 	<< PrintEfficiency(nEcalDrvn,nAcc)
       << "\n eH/E: \t \t \t" 		<< PrintEfficiency(nHadFrac, nEcalDrvn)
       << "\n e|DEtaIn|: \t \t" 	<< PrintEfficiency(nDEtaIn, nHadFrac)
       << "\n e|DPhiIn|: \t \t" 	<< PrintEfficiency(nDPhiIn, nDEtaIn)
       << "\n eSigmaIEta: \t \t" 	<< PrintEfficiency(nSigmaIEta, nDPhiIn)
       << "\n 2x5/5x5: \t \t" 		<< PrintEfficiency(nClusterRatio, nSigmaIEta)
       << "\n eMissHits: \t \t"		<< PrintEfficiency(nEMissHits, nClusterRatio)
       << "\n eTrkIso: \t \t"		<< PrintEfficiency(nTightElecTrkIso, nEMissHits)
       << "\n eEcalIso: \t \t" 		<< PrintEfficiency(nTightElecEcalIso, nElecTrkIso)
/*        
       << "\n HEEP Et: \t \t" 	 	<< PrintEfficiency(nHeepEt, nMuonEta)
       << "\n HEEP EcalDrvn: \t" 	<< PrintEfficiency(nHeepEcalDrvn, nHeepEt)
       << "\n HEEP H/E: \t \t" 		<< PrintEfficiency(nHeepHadFrac, nHeepEcalDrvn)
       << "\n HEEP |DEtaIn|: \t"	<< PrintEfficiency(nHeepDEtaIn, nHeepHadFrac)
       << "\n HEEP |DPhiIn|: \t"	<< PrintEfficiency(nHeepDPhiIn, nHeepDEtaIn)
       << "\n HEEP SigmaIEta: \t"	<< PrintEfficiency(nHeepSigmaIEta, nHeepDPhiIn)
       << "\n HEEP 2x5/5x5:\t \t"	<< PrintEfficiency(nHeepClusterRatio, nHeepSigmaIEta)
       << "\n HEEP MissHits: \t"	<< PrintEfficiency(nHeepMissHits, nHeepClusterRatio)
       << "\n eIso: \t \t \t" 		<< PrintEfficiency(nTightElecIso, nHeepMissHits)
       //<< "\n eTrkIso: \t \t"		<< PrintEfficiency(nElecTrkIso, nElecEcalIso)
*/      

       //<< "\n HEEP CaloIso: \t \t"	<< PrintEfficiency(nHeepCaloIso, nHeepMissHits)
       //<< "\n HEEP HadIso: \t \t"	<< PrintEfficiency(nHeepHadIso, nHeepCaloIso)
       //<< "\n HEEP TrkIso: \t \t"	<< PrintEfficiency(nHeepTrkIso, nHeepHadIso)
       //<< "\n MuChi2: \t \t"		<< PrintEfficiency(nMuChi2, nElecId)
      
       << "\n MuChi2: \t \t"		<< PrintEfficiency(nMuChi2, nElecEcalIso)
       << "\n MuIp: \t \t \t"		<< PrintEfficiency(nMuIp, nMuChi2)
       << "\n MuNHits: \t \t"		<< PrintEfficiency(nMuonHits, nMuIp)
       //<< "\n Mu2DComp: \t \t"		<< PrintEfficiency(nMuComp, nMuonHits)
       << "\n MuPiVeto: \t \t"		<< PrintEfficiency(nMuComp, nMuonHits)
       //<< "\n MuIso: \t \t"		<< PrintEfficiency(nTightMuonIso, nMuComp)
       << "\n MuTrkIso: \t \t"		<< PrintEfficiency(nTightMuonTrkIso, nMuComp)
       << "\n MuEcalIso: \t \t"		<< PrintEfficiency(nTightMuonEcalIso, nMuTrkIso)

       << "\n cosDPhi: \t \t"		<< PrintEfficiency(nCosDPhi, nMuCaloIso)
       << "\n diTauCharge: \t \t"	<< PrintEfficiency(nTightCharge, nCosDPhi)
       << "\n BTag: \t \t \t"		<< PrintEfficiency(nBTag, nDiTauCharge)
       << "\n JetEtSum: \t \t"		<< PrintEfficiency(nJetEtSum, nBTag)
       << "\n MET: \t \t \t"		<< PrintEfficiency(nMet,nJetEtSum)
       << "\n LeptonMETCosDPhi \t"	<< PrintEfficiency(nLeptonMETCosDPhi,nMet)
       << "\n Zeta2DCut: \t \t"		<< PrintEfficiency(nZeta2D,nLeptonMETCosDPhi)
       << "\n Mass: \t \t \t"		<< PrintEfficiency(nMass,nZeta2D)
       
       << "\n \n Trigger: \t \t"	<< PrintEfficiency(nTriggers, nEvents)
       << "\n Acc: \t \t \t"		<< PrintEfficiency(nAcc, nTriggers)
       << "\n eID: \t \t \t"		<< PrintEfficiency(nElecId,nAcc)
       //<< "\n HEEP eID: \t \t"		<< PrintEfficiency(nHeepElecId,nAcc)
       << "\n muonID: \t \t"		<< PrintEfficiency(nMuId,nElecId)
       << "\n diTau: \t \t"		<< PrintEfficiency(nDiTau,nMuId)
              
       << "\n \n Total Events: \t \t" 	<< nEvents 
       
       << "\n Final Nevents Tight: \t" 	<< (int)nDiTau
       
       << "\n Final Nevents Loose: \t" 	<< (int)nLoose
      
       << "\n Cumulative Eff: \t" 	<< setprecision(7) << PrintEfficiency(survivingEvents,nEvents)
       
       << "\n Total Events @" << (int)_theLumi << "pb-1:\t" <<  survivingEventsAtLumi.first 
       << "\t +- \t" <<  survivingEventsAtLumi.second << endl;
       
       if(_label != "data"){
         outLogFile << PrintNMinus1Eff();
	 outLogFile << "NEvetnts after PU ReWeighting " << (int)puWeightedEvents << "\n";
         //outLogFile << PrintReWeightedEvents();
       }
       if(_label == "data") outLogFile << PrintBgFromDataInfo();
       
       outLogFile << "Btag Control Region 2: Needed for Correction of btag eff" << "\n"
       		  << "NBtag CR: \t \t \t" << survivingEventsAtLumiBTagControl.first 
                  << "\t +- \t" <<  survivingEventsAtLumiBTagControl.second 
		  << "\n";
  outLogFile.close();
  //for(unsigned int evt = 0; evt < theEventNumberVector.size(); ++evt){
  //  cout << "RunNumber:EventNumber\t" << theRunNumberVector.at(evt) << "\t" << theEventNumberVector.at(evt) << "\n";
  //}
}

void EMuAnalyzer::BookHistos(){
  outRootFile = new TFile(_outRootFile.c_str(), "Recreate"); 
  
  TDirectory* theDir; 
  std::string histosDir = _label + "HistosDirectory"; 
  if(outRootFile->FindObjectAny(histosDir.c_str()) == 0) theDir = outRootFile->mkdir(histosDir.c_str());
  else theDir = (TDirectory*)outRootFile->FindObjectAny(histosDir.c_str());
  theDir->cd();
  
  _elecMuDeltaR = new TH1F("_elecMuDeltaR", "EMu #Delta R", 50, 0, 2*3.1416);
  _elecPt = new TH1F("_elecPt", "elec Pt", 50, 0, 300);
  _elecEta = new TH1F("_elecEta", "elec Eta", 50, -3., 3.);
  _muonPt = new TH1F("_muonPt", "muon Pt", 50, 0, 300);
  _muonEta = new TH1F("_muonEta", "muon Eta", 50, -3, 3);
  
  _elecHOverE = new TH1F("_elecHOverE", "elec H/E", 25, 0, 0.1);
  _elecDEtaIn = new TH1F("_elecDEtaIn", "elec #Delta#eta", 25, -0.05, 0.05);
  _elecDPhiIn = new TH1F("_elecDPhiIn", "elec #Delta#phi", 25, -0.05, 0.05);
  _elecSigmaIEtaIEta = new TH1F("_elecSigmaIEtaIEta", "elec #sigma_{i#eta i#eta}IEta", 25, 0, 0.05);
  _elec2x5Over5x5 = new TH1F("_elec2x5Over5x5", "elec 2x5 Cluster Ratio", 25, 0, 1);
  _elec1x5Over5x5 = new TH1F("_elec1x5Over5x5", "elec 1x5 Cluster Ratio", 25, 0, 1);
  _elecMissingHits = new TH1F("_elecMissingHits", "elec Missing Hits", 10, 0, 10);
  
  _elecEcalIsoDr03 = new TH1F("_elecEcalIsoDr03", "elec Ecal Iso #Delta R 0.3", 15, 0, 15);
  _elecTrkIsoDr03 = new TH1F("_elecTrkIsoDr03", "elec Trk Iso #Delta R 0.3", 15, 0, 15);
  _elecEcalIsoDr04 = new TH1F("_elecEcalIsoDr04", "elec Ecal Iso #Delta R 0.4", 15, 0, 15);
  _elecTrkIsoDr04 = new TH1F("_elecTrkIsoDr04", "elec Trk Iso #Delta R 0.4", 15, 0, 15);
  _elecIsoDr03 = new TH1F("_elecIsoDr03", "elec Iso (Trk+Ecal) #Delta R 0.3", 15, 0, 15);
  _elecIsoDr04 = new TH1F("_elecIsoDr04", "elec Iso (Trk+Ecal) #Delta R 0.4", 15, 0, 15);
  
  _muonChi2NdF = new TH1F("_muonChi2NdF", "muon #chi/NDF", 25, 0, 25);
  _muonIp = new TH1F("_muonIp", "muon IP", 50, -0.25, 0.25);
  _muonHits = new TH1F("_muonHits", "muon Hits", 100, 0, 100);
  
  _muonEcalIsoDr03 = new TH1F("_muonEcalIsoDr03", "muon Ecal Iso# Delta R 0.3", 15, 0, 15);
  _muonTrkIsoDr03 = new TH1F("_muonTrkIsoDr03", "muon Trk Iso #Delta R 0.3", 15, 0, 15);
  _muonEcalIsoDr04 = new TH1F("_muonEcalIsoDr04", "muon Ecal Iso #Delta R 0.4", 15, 0, 15);
  _muonTrkIsoDr04 = new TH1F("_muonTrkIsoDr04", "muon Trk Iso #Delta R 0.4", 15, 0, 15);
  _muonEcalIsoDr05 = new TH1F("_muonEcalIsoDr05", "muon Ecal Iso #Delta R 0.5", 15, 0, 15);
  _muonTrkIsoDr05 = new TH1F("_muonTrkIsoDr05", "muon Trk Iso #Delta R 0.5", 15, 0, 15);
  _muonIsoDr03 = new TH1F("_muonIsoDr03", "muon Iso (Trk+Ecal) #Delta R 0.3", 15, 0, 15);
  _muonIsoDr04 = new TH1F("_muonIsoDr04", "muon Iso (Trk+Ecal) #Delta R 0.4", 15, 0, 15);
  _muonIsoDr05 = new TH1F("_muonIsoDr05", "muon Iso (Trk+Ecal) #Delta R 0.5", 15, 0, 15);
  
  
  _elecMuCosDPhi = new TH1F("_elecMuCosDPhi", "EMu cos#Delta#phi", 40, -1, 1);
  _elecMuCharge = new TH1F("_elecMuCharge", "EMu charge", 4, -2, 2);
  
  _btagCountHiEffTrkCnt = new TH1F("_btagCountHiEffTrkCnt", "NBtagHiEffTrkCnt", 5, 0, 5);
  _btagCountHiPurityTrkCnt = new TH1F("_btagCountHiPurityTrkCnt", "NBtagHiPurityTrkCnt", 5, 0, 5);
  _btagCountSimpleSecVtx = new TH1F("_btagCountSimpleSecVtx", "NBtagHiEffSimpleSecVtx", 5, 0, 5);
  _btagCountHiPuritySimpleSecVtx = new TH1F("_btagCountHiPuritySimpleSecVtx", "NBtagHiPuritySimpleSecVtx", 5, 0, 5);
  _btagCountCombSecVtx = new TH1F("_btagCountCombSecVtx", "NBtag", 5, 0, 5);
    
  _btagJetTrkCntHiEff = new TH1F("_btagJetTrkCntHiEff", "Btag TrkCnt HiEff", 100, -50, 50);
  _btagJetTrkCntHiPurity = new TH1F("_btagJetTrkCntHiPurity", "Btag TrkCnt HiPurity", 100, -50, 50);
  _btagJetSimpleSecVtxHiEff = new TH1F("_btagJetSimpleSecVtxHiEff", "Btag Simple Sec Vtx HiEff", 50, 0, 10);
  _btagJetSimpleSecVtxHiPurity = new TH1F("_btagJetSimpleSecVtxHiPurity", "Btag Simple Sec Vtx HiPurity", 50, 0, 10);
  _btagJetCombSecVtx = new TH1F("_btagJetCombSecVtx", "Btag Comb Sec Vtx", 100, 0, 1);
  
  _jetSumEt = new TH1F("_jetSumEt", "JetSumEt", 100, 0, 1000);
  _jetMetSumEt = new TH1F("_jetMetSumEt", "JetMETSumEt", 100, 0, 1000);
  _leptonMETCosDPhi = new TH1F("_leptonMETCosDPhi", "LeptonMET cos#Delta#phi", 50, -1, 1);
  
  _met = new TH1F("_met", "MET", 25, 0, 250);
  _elecMuMass = new TH1F("_elecMuMass", "EMu Mass", (int)(_massCut/5.), 0, _massCut);
  //_elecMuMass = new TH1F("_elecMuMass", "EMu Mass", 100, 0, 2000.);
  _elecMuMetMass = new TH1F("_elecMuMetMass", "EMuMet Mass", (int)(_massCut/5.), 0, _massCut);
  _elecMetMass = new TH1F("_elecMetMass", "EMet Mass", 30, 0, 300);
  _muonMetMass = new TH1F("_muonMetMass", "MuonMet Mass", 30, 0, 300);

  _elecMuMetMassDataQCD = new TH1F("_elecMuMetMassDataQCD", "EMuMet Mass QCD CR", (int)(_massCut/5.), 0, _massCut);
  _elecMuMetMassDataTTbar = new TH1F("_elecMuMetMassDataTTbar", "EMuMet Mass TTbar CR", (int)(_massCut/5.), 0, _massCut);
  _elecMuMetMassDataWJets = new TH1F("_elecMuMetMassDataWJets", "EMuMet Mass WJets CR", (int)(_massCut/5.), 0, _massCut);
  
  _pzeta2D = new TH2F("_pzeta2D", "CDF Zeta 2D", 50, 0, 100, 50, -50, 50);
  _pzeta = new TH1F("_pzeta", "CDF Zeta", 50, -50, 50);
  _pzetaVis = new TH1F("_pzetaVis", "CDF ZetaVis", 50, 0, 100);
  
  _leptonMetDPhiVsJetSum = new TH2F("_leptonMetDPhiVsJetSum", "leptonMetDPhiVsJetSum", 50, 0, 500, 20, -1, 1);
  _leptonMetDPhiVsNBtag = new TH2F("_leptonMetDPhiVsNBtag","leptonMetDPhiVsNBtag", 5, 0, 5, 20, -1, 1);
  
  _nJetsElecPt = new TH1F("_nJetsElecPt", "NJetsEvent", 15, 0, 15);
  
  Template = new TH1F("Template", "Mass", 100, 0, 2000); 
  _btagCountHiEffTrkCntControl = new TH1F("_btagCountHiEffTrkCntControl", "NBtagHiEffTrkCntControl", 5, 0, 5);
  
  _elecMuDeltaRReWeighted = new TH1F("_elecMuDeltaRReWeighted", "EMu #Delta R", 50, 0, 2*3.1416);
  _elecPtReWeighted = new TH1F("_elecPtReWeighted", "elec Pt", 100, 0, 300);
  _elecEtaReWeighted = new TH1F("_elecEtaReWeighted", "elec Eta", 50, -3., 3.);
  _muonPtReWeighted = new TH1F("_muonPtReWeighted", "muon Pt", 100, 0, 300);
  _muonEtaReWeighted = new TH1F("_muonEtaReWeighted", "muon Eta", 50, -3, 3);
  
  //_pileUPDistro = new TH1F("_pileUPDistro", "InTimePileUP", _pileUpNBins, _pileUpLowEdge, _pileUpNBins*_pileUpWidth + _pileUpLowEdge);
  //_pileUPReweighted = new TH1F("_pileUPReweighted", "InTimePileUP Reweighed",_pileUpNBins, _pileUpLowEdge,  _pileUpNBins*_pileUpWidth + _pileUpLowEdge); 
  _pileUPDistro = (TH1F*)dataPUHisto->Clone("_pileUPDistro");
  _pileUPReweighted = (TH1F*)dataPUHisto->Clone("_pileUPReweighted");
  
  // N-1 Plots
  
  _elecMuDeltaRNMinus1 = new TH1F("_elecMuDeltaRNMinus1", "EMu #Delta R", 50, 0, 2*3.1416);
  _elecPtNMinus1 = new TH1F("_elecPtNMinus1", "elec Pt", 250, 10, 260);
  _elecEtaNMinus1 = new TH1F("_elecEtaNMinus1", "elec Eta", 50, -3., 3.);
  _muonPtNMinus1 = new TH1F("_muonPtNMinus1", "muon Pt", 250, 10, 260);
  _muonEtaNMinus1 = new TH1F("_muonEtaNMinus1", "muon Eta", 50, -3, 3);
  
  _elecHOverENMinus1 = new TH1F("_elecHOverENMinus1", "elec H/E", 50, 0, 0.25);
  _elecDEtaInNMinus1 = new TH1F("_elecDEtaInNMinus1", "elec #Delta#eta", 20, -0.05, 0.05);
  _elecDPhiInNMinus1 = new TH1F("_elecDPhiInNMinus1", "elec #Delta#phi", 100, -0.25, 0.25);
  _elecSigmaIEtaIEtaNMinus1 = new TH1F("_elecSigmaIEtaIEtaNMinus1", "elec #sigma_{i#eta i#eta}", 100, 0, 0.1);
  _elec2x5Over5x5NMinus1 = new TH1F("_elec2x5Over5x5NMinus1", "elec 2x5 Cluster Ratio", 100, 0, 1);
  _elec1x5Over5x5NMinus1 = new TH1F("_elec1x5Over5x5NMinus1", "elec 1x5 Cluster Ratio", 100, 0, 1);
  _elecMissingHitsNMinus1 = new TH1F("_elecMissingHitsNMinus1", "elec Missing Hits", 10, 0, 10);
  
  _elecEcalIsoDr03NMinus1 = new TH1F("_elecEcalIsoDr03NMinus1", "elec Ecal Iso #Delta R 0.3", 15, 0, 15);
  _elecTrkIsoDr03NMinus1 = new TH1F("_elecTrkIsoDr03NMinus1", "elec Trk Iso #Delta R 0.3", 15, 0, 15);
  _elecEcalIsoDr04NMinus1 = new TH1F("_elecEcalIsoDr04NMinus1", "elec Ecal Iso #Delta R 0.4", 15, 0, 15);
  _elecTrkIsoDr04NMinus1 = new TH1F("_elecTrkIsoDr04NMinus1", "elec Trk Iso #Delta R 0.4", 15, 0, 15);
  _elecIsoDr03NMinus1 = new TH1F("_elecIsoDr03NMinus1", "elec Iso (Trk+Ecal) #Delta R 0.3", 15, 0, 15);
  _elecIsoDr04NMinus1 = new TH1F("_elecIsoDr04NMinus1", "elec Iso (Trk+Ecal) #Delta R 0.4", 15, 0, 15);
  
  _muonChi2NdFNMinus1 = new TH1F("_muonChi2NdFNMinus1", "muon #chi/NDF", 25, 0, 25);
  _muonIpNMinus1 = new TH1F("_muonIpNMinus1", "muon IP", 100, -0.05, 0.05);
  _muonHitsNMinus1 = new TH1F("_muonHitsNMinus1", "muon Hits", 50, 0, 50);
  
  _muonEcalIsoDr03NMinus1 = new TH1F("_muonEcalIsoDr03NMinus1", "muon Ecal Iso# Delta R 0.3", 15, 0, 15);
  _muonTrkIsoDr03NMinus1 = new TH1F("_muonTrkIsoDr03NMinus1", "muon Trk Iso #Delta R 0.3", 15, 0, 15);
  _muonEcalIsoDr04NMinus1 = new TH1F("_muonEcalIsoDr04NMinus1", "muon Ecal Iso #Delta R 0.4", 15, 0, 15);
  _muonTrkIsoDr04NMinus1 = new TH1F("_muonTrkIsoDr04NMinus1", "muon Trk Iso #Delta R 0.4", 15, 0, 15);
  _muonEcalIsoDr05NMinus1 = new TH1F("_muonEcalIsoDr05NMinus1", "muon Ecal Iso #Delta R 0.5", 15, 0, 15);
  _muonTrkIsoDr05NMinus1 = new TH1F("_muonTrkIsoDr05NMinus1", "muon Trk Iso #Delta R 0.5", 15, 0, 15);
  _muonIsoDr03NMinus1 = new TH1F("_muonIsoDr03NMinus1", "muon Iso (Trk+Ecal) #Delta R 0.3", 15, 0, 15);
  _muonIsoDr04NMinus1 = new TH1F("_muonIsoDr04NMinus1", "muon Iso (Trk+Ecal) #Delta R 0.4", 15, 0, 15);
  _muonIsoDr05NMinus1 = new TH1F("_muonIsoDr05NMinus1", "muon Iso (Trk+Ecal) #Delta R 0.5", 15, 0, 15);
  
  _elecMuCosDPhiNMinus1 = new TH1F("_elecMuCosDPhiNMinus1", "EMu cos#Delta#phi", 40, -1, 1);
  _elecMuChargeNMinus1 = new TH1F("_elecMuChargeNMinus1", "EMu charge", 4, -2, 2);
  
  _btagCountHiEffTrkCntNMinus1 = new TH1F("_btagCountHiEffTrkCntNMinus1", "NBtagHiEffTrkCnt", 5, 0, 5);
  _btagCountHiPurityTrkCntNMinus1 = new TH1F("_btagCountHiPurityTrkCntNMinus1", "NBtagHiPurityTrkCnt", 5, 0, 5);
  _btagCountSimpleSecVtxNMinus1 = new TH1F("_btagCountSimpleSecVtxNMinus1", "NBtagHiEffSimpleSecVtx", 5, 0, 5);
  _btagCountHiPuritySimpleSecVtxNMinus1 = new TH1F("_btagCountHiPuritySimpleSecVtxNMinus1", "NBtagHiPuritySimpleSecVtx", 5, 0, 5);
  _btagCountCombSecVtxNMinus1 = new TH1F("_btagCountCombSecVtxNMinus1", "NBtag", 5, 0, 5);
  
  _jetSumEtNMinus1 = new TH1F("_jetSumEtNMinus1", "JetSumEt", 100, 0, 1000);
  _jetMetSumEtNMinus1 = new TH1F("_jetMetSumEtNMinus1", "JetMETSumEt", 100, 0, 1000);
  _leptonMETCosDPhiNMinus1 = new TH1F("_leptonMETCosDPhiNMinus1", "LeptonMET cos#Delta#phi", 50, -1, 1);
  
  _metNMinus1 = new TH1F("_metNMinus1", "MET", 50, 0, 250);
  _metStdNMinus1 = new TH1F("_metStdNMinus1", "MET No PUReWeighted", 50, 0, 250);
  _pzetaNMinus1 = new TH1F("_pzetaNMinus1", "#zeta", 40, -100, 100);
  _pzetaVisNMinus1 = new TH1F("_pzetaVisNMinus1", "#zeta_{vis}", 20, 0, 100);
  _pzetaPZetaVisNMinus1 = new TH1F("_pzetaPZetaVisNMinus1", "#zeta - 1.25 #zeta_{vis}", 40, -100, 100);
  _nJetsNMinus1 = new TH1F("_nJetsNMinus1", "NJets", 10, 0, 10);
}

void EMuAnalyzer::WriteOutFile(){
  WriteTree();
  outRootFile->Write();
  outRootFile->Close();
} 


void EMuAnalyzer::WriteTree(){
  
  //Set Tree for fitter...
  outRootFile->cd();
  theOutTree = new TTree(theTreeOutName.c_str(), "TheEMuTree");
  
  theOutTree->Branch("branchingFraction",&branchingFraction, "branchingFraction/F");
  theOutTree->Branch("FilterEfficiency",&filterEff,"FilterEfficiency/F");
  theOutTree->Branch("xsection", &_theXSection.first, "xsection/F");
  theOutTree->Branch("acceptanceEfficiency",&acceptanceEfficiency,"acceptanceEfficiency/F");
  theOutTree->Branch("electronIDEfficiency",&electronIDEfficiency,"electronIDEfficiency/F");
  theOutTree->Branch("muonIDEfficiency",&muonIDEfficiency,"muonIDEfficiency/F");
  theOutTree->Branch("tauIDEfficiency",&tauIDEfficiency,"tauIDEfficiency/F");
  theOutTree->Branch("topologyEfficiency",&topologyEfficiency,"topologyEfficiency/F");
  
  if(_label == "zprime") acceptanceEfficiency = GetEfficiency(nDiTau, nEvents).first;
  else 
    if(_label == "ztautau" ||_label == "ww" || _label == "wz" || _label == "zz")acceptanceEfficiency = 
    GetEfficiency(nDiTau, nEvents).first * GetEfficiency(skimmedEvents.first, skimmedEvents.second).first;
  else acceptanceEfficiency = _bgEstimate/(_theLumi*_theXSection.first*1e12);
 // else acceptanceEfficiency = GetEfficiency(nDiTau, nEvents).first * GetEfficiency(skimmedEvents.first, skimmedEvents.second).first;
  cout << "Acceptance Eff " << "\t" << _label << "\t" << acceptanceEfficiency << endl;
    
  // these efficiencies are 1; only accepetance is needed
  electronIDEfficiency = 1.;
  tauIDEfficiency = 1.;
  topologyEfficiency = 1.;
  muonIDEfficiency = 1.;														  
  
  theOutTree->Fill();
  
  TDirectory* theDir; 
  std::string dirName = theTreeOutName + "TemplateDirectory"; 
  if(outRootFile->FindObjectAny(dirName.c_str()) == 0) theDir = outRootFile->mkdir(dirName.c_str());
  else theDir = (TDirectory*)outRootFile->FindObjectAny(dirName.c_str());
  
  theDir->cd();
  Template->Write();
}

pair<double, double> EMuAnalyzer::GetFactorizedEvents(pair<double, double> nLoose, pair<double,double> effLoose, 
                                                      pair<double, double> effTight){
  float normFactor = effTight.first/effLoose.first;
  double normFactorSqrd = normFactor*normFactor;
  double sigmaEffLooseSqrd = effLoose.second*effLoose.second;
  double sigmaEffTightSqrd = effTight.second*effTight.second;
  
  double nLooseSqrd = nLoose.first*nLoose.first;
  double sigmaNLooseSqrd = nLoose.second*nLoose.second;
  
  double sigmaNormFactor = sqrt(sigmaEffTightSqrd + (normFactorSqrd*sigmaEffLooseSqrd)/effLoose.first);
  double normNEvents = nLoose.first * normFactor;
  
  double sigmaNormEvents = sqrt((normFactorSqrd*sigmaNLooseSqrd) +(nLooseSqrd*sigmaNormFactor*sigmaNormFactor));
  pair<double, double> theFactorizedNEvents = make_pair<double, double>(normNEvents, sigmaNormEvents);
  
  //std::cout << "Norm NEvents " << normNEvents << "\t" << sigmaNormEvents << "\n";

  return theFactorizedNEvents;

}
template <typename T1, typename T2>
std::pair<double, std::pair<double, double> > EMuAnalyzer::GetEfficiency(T1 passed, T2 total){
  double theEff = double(passed)/double(total);
  double theEffError = sqrt(theEff*(1.-theEff)/double(total));
  pair<double, double> theEffErrors = make_pair<double, double>(theEffError, theEffError);  
  if((theEff + theEffError) >= 1. || (theEff - theEffError <= 0.)){
    //do one bin histo && extract bayes eff
    TH1F* theNumHisto = new TH1F("theNumHisto","theNumHisto",1,0,1);
    theNumHisto->SetBinContent(1,passed);
    theNumHisto->Sumw2();
    TH1F* theDenHisto = new TH1F("theDenHisto","",1,0,1);
    theDenHisto->SetBinContent(1,total);
    theDenHisto->Sumw2();
    TGraphAsymmErrors* bayesEff = new TGraphAsymmErrors();
    bayesEff->BayesDivide(theNumHisto,theDenHisto,"b");
    double effErrorHigh = bayesEff->GetErrorYhigh(0);
    double effErrorLow = bayesEff->GetErrorYlow(0);
    theEffErrors = make_pair<double, double>(effErrorLow, effErrorHigh); 
 
    delete theNumHisto;
    delete theDenHisto;
    delete bayesEff;
  
  }
  pair<double,pair<double,double> > theEffAndError = make_pair<double,pair<double,double> >(theEff,theEffErrors);
  return theEffAndError;
}
template <typename T1, typename T2>
std::string EMuAnalyzer::PrintEfficiency(T1 passed, T2 total){
  std::pair<double, std::pair<double, double> > theEff = GetEfficiency(passed, total);
  
  std::stringstream theEffValue;
  theEffValue.setf(ios::fixed);
  theEffValue.setf(ios::showpoint);
  theEffValue.precision(4);
  theEffValue << theEff.first << "\t +- \t"  << max(theEff.second.first, theEff.second.second);
  return theEffValue.str();
}

std::string EMuAnalyzer::PrintReWeightedEvents(){
  std::stringstream theOutputStream;
  theOutputStream.setf(ios::fixed);
  theOutputStream.setf(ios::showpoint);
  theOutputStream.precision(4);
  theOutputStream 
       << "\n Central PDF\t"<< GetPDFUncertaity(pdfWeightsSum).second/GetCentralPDF(pdfWeightsSum)
       << "\n PDF ReWeighed Events: \t"	<< (float)abs(GetCentralPDF(pdfWeightsSum)-nEvents)/nEvents << "\t" 
       <<  GetPDFUncertaity(pdfWeightsSum).first/GetCentralPDF(pdfWeightsSum) << "\t"  
       << "\t" << GetPDFUncertaity(pdfWeightsSum).second/GetCentralPDF(pdfWeightsSum)
       << "\n Syst: \t\t\t PDF \t\tISRGluon \tISRGamma \tFSR" 
       
       << "\n Xsection:\t \t"
       << fabs(GetCentralPDF(pdfWeightsSum)-nEvents)/(float)nEvents
       << "\t\t" << fabs(isrGluonWeightedEvents-nEvents)/(float)nEvents
       << "\t\t" << fabs(isrGammaWeightedEvents-nEvents)/(float)nEvents
       << "\t\t" << fabs(fsrWeightedEvents-nEvents)/(float)nEvents
       
       << "\n ReWeighed Acc: \t" 
       << GetPerCentError(nAcc, nEvents, (int)pdfWeightedEventsAcc, (int)pdfWeightedEvents)
       << "\t\t" <<   GetPerCentError(nAcc, nEvents, (int)isrGluonWeightedEventsAcc, (int)isrGluonWeightedEvents)
       << "\t\t" <<   GetPerCentError(nAcc, nEvents, (int)isrGammaWeightedEventsAcc, (int)isrGammaWeightedEvents)
       << "\t\t" <<   GetPerCentError(nAcc, nEvents, (int)fsrWeightedEventsAcc, (int)fsrWeightedEvents)
       
       << "\n ReWeighed eId: \t"  
       <<  GetPerCentError(nElecId, nAcc, (int)pdfWeightedEventsEId, (int)pdfWeightedEventsAcc)
       << "\t\t" << GetPerCentError(nElecId, nAcc, (int)isrGluonWeightedEventsEId, (int)isrGluonWeightedEventsAcc)
       << "\t\t" << GetPerCentError(nElecId, nAcc, (int)isrGammaWeightedEventsEId, (int)isrGammaWeightedEventsAcc)
       << "\t\t" << GetPerCentError(nElecId, nAcc, (int)fsrWeightedEventsEId, (int)fsrWeightedEventsAcc)
       
       << "\n ReWeighed muonId: \t"
       << GetPerCentError(nMuId, nElecId, (int)pdfWeightedEventsMuonId, (int)pdfWeightedEventsEId)
       << "\t\t" << GetPerCentError(nMuId, nElecId, (int)isrGluonWeightedEventsMuonId, (int)isrGluonWeightedEventsEId)
       << "\t\t" << GetPerCentError(nMuId, nElecId, (int)isrGammaWeightedEventsMuonId, (int)isrGammaWeightedEventsEId)
       << "\t\t" << GetPerCentError(nMuId, nElecId, (int)fsrWeightedEventsMuonId, (int)fsrWeightedEventsEId)
       
       << "\n ReWeighed Topo: \t"
       << GetPerCentError(nDiTau, nMuId, (int)pdfWeightedEventsTopology, (int)pdfWeightedEventsMuonId)
       << "\t\t" << GetPerCentError(nDiTau, nMuId, (int)isrGluonWeightedEventsTopology, (int)isrGluonWeightedEventsMuonId)
       << "\t\t" << GetPerCentError(nDiTau, nMuId, (int)isrGammaWeightedEventsTopology, (int)isrGammaWeightedEventsMuonId)
       << "\t\t" << GetPerCentError(nDiTau, nMuId, (int)fsrWeightedEventsTopology, (int)fsrWeightedEventsMuonId)
       
       << "\n ReWeighed All: \t"	 
       << GetPerCentError(nDiTau, nEvents, (int)pdfWeightedEventsTopology, (int)pdfWeightedEvents)
       << "\t\t" << GetPerCentError(nDiTau, nEvents, (int)isrGluonWeightedEventsTopology, (int)isrGluonWeightedEvents)
       << "\t\t" << GetPerCentError(nDiTau, nEvents, (int)isrGammaWeightedEventsTopology, (int)isrGammaWeightedEvents)
       << "\t\t" << GetPerCentError(nDiTau, nEvents, (int)fsrWeightedEventsTopology, (int)fsrWeightedEvents)
       << "\t\n";
  return theOutputStream.str();
}

std::string EMuAnalyzer::PrintBgFromDataInfo(){

  pair<double, double> NPureQCD = make_pair<double, double>(nDataQCDAlt, sqrt(nDataQCDAlt));
  pair<double, double> NPureQCDElecPt = GetSurvivingEvents(NPureQCD, GetEfficiency(nDataQCDAltElecPt, nDataQCDAlt));
  pair<double, double> NPureQCDElecIso = GetSurvivingEvents(NPureQCDElecPt, GetEfficiency(nDataQCDAltElecIso, nDataQCDAlt));
  //pair<double, double> NPureQCDElecPt = GetSurvivingEvents(NPureQCDElecIso, GetEfficiency(nDataQCDAltElecEt, nDataQCDAlt));
  //pair<double, double> NPureQCDOS = GetSurvivingEvents(NPureQCDElecPt, GetEfficiency(nDataQCDAltOS, nDataQCDAlt));
  pair<double, double> NPureQCDOS = GetSurvivingEvents(NPureQCDElecIso, GetEfficiency(nDataQCDAltOS, nDataQCDAlt));
  pair<double, double> NPureQCDMet =  GetSurvivingEvents(NPureQCDOS, GetEfficiency(nDataQCDAltMet, nDataQCDAlt));
  
  double muIsoEff = GetEfficiency(nDataQCDMuIso, nDataQCD).first;
  double muAntiIsoEff = GetEfficiency(nDataQCDAntiMuIso, nDataQCD).first;
  
  double muIsoAntiIsoRatio = muIsoEff/muAntiIsoEff;
  double muIsoEffError = max(GetEfficiency(nDataQCDMuIso, nDataQCD).second.first,
                          GetEfficiency(nDataQCDMuIso, nDataQCD).second.second);
  double ratioError = ((1./muAntiIsoEff + (muIsoEff)/(muAntiIsoEff*muAntiIsoEff))*muIsoEffError);
  //std::cout << "Ratio \t" << muIsoAntiIsoRatio << "\t +- \t" << ratioError << endl;
  
  pair<double,double> ratioErrorPair = make_pair<double,double>(ratioError,ratioError);
  pair<double,pair<double,double> >theRatioPair = make_pair<double,pair<double,double> >(muIsoAntiIsoRatio, ratioErrorPair);
  pair<double,double> finalQCDEvents = GetSurvivingEvents(NPureQCDMet, theRatioPair);
  
  //get ttbar estimation; All met correlated variables are done in one step.
  
  pair<double, double> NTTbar = make_pair<double, double>(nDataTTbar, sqrt(nDataTTbar));
  pair<double, double> NTTbarAllMET = GetSurvivingEvents(NTTbar, GetEfficiency(nDataTTbarAllMET, nDataTTbar));
  pair<double, double> NTTbarJeSum = GetSurvivingEvents(NTTbarAllMET, GetEfficiency(nDataTTbarJetSum, nDataTTbar));
  
  //pair<double, double> NTTbarBTagVeto = GetSurvivingEvents(NTTbarZeta2D, GetEfficiency(nDataTTbarAltBtagVeto, nDataTTbarAlt));
  double btagEff = GetEfficiency(nDataTTbarAltBtagVeto, nDataTTbarAlt).first;
  double btagAntiEff = 1. - btagEff;
  double btagRatio = btagEff/btagAntiEff;
  double btagEffError = max(GetEfficiency(nDataTTbarAltBtagVeto, nDataTTbarAlt).second.first,
                          GetEfficiency(nDataTTbarAltBtagVeto, nDataTTbarAlt).second.second);
  double btagRatioError = ((1./btagAntiEff + (btagEff)/(btagAntiEff*btagAntiEff))*btagEffError);
  pair<double,double> bTagRatioErrorPair = make_pair<double,double>(btagRatioError,btagRatioError);
  pair<double,pair<double,double> >btagRatioPair = make_pair<double,pair<double,double> >(btagRatio, bTagRatioErrorPair);
  pair<double,double> finalTTbarEvents = GetSurvivingEvents(NTTbarJeSum, btagRatioPair);
  
  
  //get WJets estimation
  pair<double, double> NWJets = make_pair<double, double>(nDataWJets, sqrt(nDataWJets));
  pair<double, double> NWJetsElecPt = GetSurvivingEvents(NWJets, GetEfficiency(nDataWJetsElecPt, nDataWJets));
  //pair<double, double> NWJetsElecIso = GetSurvivingEvents(NWJetsElecPt, GetEfficiency(nDataWJetsElecIso, nDataWJets));
  
  pair<double, double> NWJetsOS =  GetSurvivingEvents(NWJetsElecPt, GetEfficiency(nDataWJetsOS, nDataWJets));
  pair<double, double> NWJetsMET =  GetSurvivingEvents(NWJetsOS, GetEfficiency(nDataWJetsMET, nDataWJets));
  
  pair<double,double> NWJetsAlt = make_pair<double,double>(nDataWJetsAlt, sqrt(nDataWJetsAlt));
  
  double theMTEff = GetEfficiency(nDataWJetsAltMT, nDataWJetsAlt).first;
  double theMTEffError = max((GetEfficiency(nDataWJetsAltMT, nDataWJetsAlt)).second.first, 
                             (GetEfficiency(nDataWJetsAltMT, nDataWJetsAlt)).second.second);
  double theMTEffInvError = theMTEffError/(theMTEff*theMTEff);
  pair<double, double> theMTEffInvErrorPair = make_pair(theMTEffInvError, theMTEffInvError);
  pair<double, pair<double, double> > theMTEffInverse = make_pair<double, pair<double,double> >
                                                        (1./theMTEff, theMTEffInvErrorPair);
  pair<double,double> NWJetsAltMT = GetSurvivingEvents(NWJetsMET,theMTEffInverse);
  //pair<double,double> NWJetsAltMT = GetSurvivingEvents(NWJetsMET, GetEfficiency(nDataWJetsAltMT, nDataWJetsAlt));
  
  double elecIsoEff = GetEfficiency(nDataWJetsAltElecIso, nDataWJetsAlt).first;
  double elecAntiIsoEff = GetEfficiency(nDataWJetsAltElecAntiIso, nDataWJetsAlt).first;
  
  double elecIsoAntiIsoRatio = elecIsoEff/elecAntiIsoEff;
  double elecIsoEffError = max(GetEfficiency(nDataWJetsAltElecIso, nDataWJetsAlt).second.first,
                          GetEfficiency(nDataWJetsAltElecIso, nDataWJetsAlt).second.second);
  double wjetsElecIsoRatioError = ((1./elecAntiIsoEff + (elecIsoEff)/(elecAntiIsoEff*elecAntiIsoEff))*elecIsoEffError);
  pair<double,double> wjetsElecIsoRatioErrorPair = make_pair<double,double>(wjetsElecIsoRatioError,wjetsElecIsoRatioError);
  pair<double,pair<double,double> >wjetsElecIsoRatioPair = make_pair<double,pair<double,double> >(elecIsoAntiIsoRatio, wjetsElecIsoRatioErrorPair);  
  pair<double, double> finalWJets = GetSurvivingEvents(NWJetsAltMT, wjetsElecIsoRatioPair);
  
  std::stringstream theOutputStream;
  theOutputStream.setf(ios::fixed);
  theOutputStream.setf(ios::showpoint);
  theOutputStream.precision(4);
  theOutputStream 
       << "\n N QCD from Alt Data: \t" 	<< nDataQCDAlt
       << "\n ElecPt: \t \t"		<< PrintEfficiency(nDataQCDAltElecPt, nDataQCDAlt)	        
       << "\n ElecTrkIso: \t \t"	<< PrintEfficiency(nDataQCDAltElecIso, nDataQCDAlt)	        
       << "\n OS: \t \t \t"		<< PrintEfficiency(nDataQCDAltOS, nDataQCDAlt)		        
       //<< "\n ElecEt: \t \t"		<< PrintEfficiency(nDataQCDAltElecEt, nDataQCDAlt)	        
       << "\n Met & Zeta: \t \t"	<< PrintEfficiency(nDataQCDAltMet, nDataQCDAlt) 	        
       << "\n N QCD from Data: \t" 	<< nDataQCD
       << "\n MuTrkIso: \t \t"		<< PrintEfficiency(nDataQCDMuIso, nDataQCD)		        
       << "\n MuAntiTrkIso: \t \t"	<< PrintEfficiency(nDataQCDAntiMuIso, nDataQCD) 	        
       << "\n muTrIso/AntiTrkIso: \t" 	<< muIsoAntiIsoRatio << "\t +- \t" << ratioError	        
       << "\n Final QCD: \t \t"		<< finalQCDEvents.first  << "\t +- \t" << finalQCDEvents.second 
       << "\n"
       << "\n N TTbar from Data: \t" 	<< nDataTTbar
       << "\n JetSum: \t \t"		<< PrintEfficiency(nDataTTbarJetSum, nDataTTbar)
       << "\n All MET: \t \t"		<< PrintEfficiency(nDataTTbarAllMET, nDataTTbar)
       << "\n N after all MET: \t"	<< nDataTTbarAllMET
       << "\n NTTbar 2+Jets: \t"	<< nDataTTbarAlt
       << "\n BTagVeto: \t \t"		<< PrintEfficiency(nDataTTbarAltBtagVeto, nDataTTbarAlt)
       << "\n Btag/AntiBTag: \t \t"	<< btagRatio << "\t +- \t" << btagRatioError
       << "\n Final TTbar: \t \t"	<< finalTTbarEvents.first  << "\t +- \t" << finalTTbarEvents.second
       << "\n" 
       << "\n N WJets from Data: \t" 	<< nDataWJets
       << "\n ElecPt: \t \t"		<< PrintEfficiency(nDataWJetsElecPt, nDataWJets)
       //<< "\n ElecIso: \t \t"		<< PrintEfficiency(nDataWJetsElecIso, nDataWJets)
       << "\n OS: \t \t \t"		<< PrintEfficiency(nDataWJetsOS, nDataWJets)
       << "\n CosDphi & Zeta:\t"	<< PrintEfficiency(nDataWJetsMET, nDataWJets)
       << "\n N WJets from Alt Data:\t" << nDataWJetsAlt
       << "\n ElecTrkIso: \t \t"	<< PrintEfficiency(nDataWJetsAltElecIso, nDataWJetsAlt)
       << "\n ElecAntiTrkIso: \t"	<< PrintEfficiency(nDataWJetsAltElecAntiIso, nDataWJetsAlt)
       << "\n TrkIso/AntiTrkIso: \t"	<< elecIsoAntiIsoRatio << "\t +- \t" << wjetsElecIsoRatioError
       << "\n MT: \t \t \t"		<< PrintEfficiency(nDataWJetsAltMT, nDataWJetsAlt)
       << "\n Final WJets: \t \t"	<< finalWJets.first << "\t +- \t" << finalWJets.second
       << "\n";
  return theOutputStream.str();
}

std::string EMuAnalyzer::PrintNMinus1Eff(){
  std::stringstream theOutputStream;
  theOutputStream.setf(ios::fixed);
  theOutputStream.setf(ios::showpoint);
  theOutputStream.precision(4);
  theOutputStream 			<< "\n N-1 Effs: \n"
       << "\n deltaR: \t \t" 		<< PrintEfficiency(nDiTau, nDeltaRNMinus1)
       << "\n ePt: \t \t \t"  		<< PrintEfficiency(nDiTau, nElePtNMinus1)
       << "\n eEta: \t \t \t" 		<< PrintEfficiency(nDiTau, nElecEtaNMinus1)
       << "\n muonPt: \t \t"  		<< PrintEfficiency(nDiTau, nMuonPtNMinus1)
       << "\n muonEta: \t \t" 		<< PrintEfficiency(nDiTau, nMuonEtaNMinus1)
       
       << "\n eH/E: \t \t \t" 		<< PrintEfficiency(nDiTau, nHadFracNMinus1)
       << "\n e|DEtaIn|: \t \t" 	<< PrintEfficiency(nDiTau, nDEtaInNMinus1)
       << "\n e|DPhiIn|: \t \t" 	<< PrintEfficiency(nDiTau, nDPhiInNMinus1)
       << "\n eSigmaIEta: \t \t" 	<< PrintEfficiency(nDiTau, nSigmaIEtaNMinus1)
       << "\n 2x5/5x5: \t \t" 		<< PrintEfficiency(nDiTau, nClusterRatioNMinus1)
       << "\n eMissHits: \t \t"		<< PrintEfficiency(nDiTau, nEMissHitsNMinus1)
       //<< "\n eIso: \t \t \t" 		<< PrintEfficiency(nDiTau, nElecIsoNMinus1)
       << "\n eTrkIso: \t \t"		<< PrintEfficiency(nDiTau, nElecTrkIsoNMinus1)
       << "\n eEcalIso: \t \t"		<< PrintEfficiency(nDiTau, nElecEcalIsoNMinus1)
      
       << "\n MuChi2: \t \t"		<< PrintEfficiency(nDiTau, nMuChi2NMinus1)
       << "\n MuIp: \t \t \t"		<< PrintEfficiency(nDiTau, nMuIpNMinus1)
       << "\n MuNHits: \t \t"		<< PrintEfficiency(nDiTau, nMuonHitsNMinus1)
       //<< "\n MuIso: \t \t"		<< PrintEfficiency(nDiTau, nMuIsoNMinus1)
       << "\n MuTrkIso: \t \t"		<< PrintEfficiency(nDiTau, nMuTrkIsoNMinus1)
       << "\n MuEcalIso: \t \t"		<< PrintEfficiency(nDiTau, nMuCaloIsoNMinus1)

       << "\n cosDPhi: \t \t"		<< PrintEfficiency(nDiTau, nCosDPhiNMinus1)
       << "\n diTauCharge: \t \t"	<< PrintEfficiency(nDiTau, nDiTauChargeNMinus1)
       << "\n BTag: \t \t \t"		<< PrintEfficiency(nDiTau, nBTagNMinus1)
       << "\n JetEtSum: \t \t"		<< PrintEfficiency(nDiTau, nJetEtSumNMinus1)
       << "\n LeptonMETCosDPhi \t"	<< PrintEfficiency(nDiTau, nLeptonMETCosDPhiNMinus1)
       << "\n MET: \t \t \t"		<< PrintEfficiency(nDiTau, nMetNMinus1)
       << "\n";
       
  return theOutputStream.str();
}

pair<double, double> EMuAnalyzer::GetSurvivingEvents(pair<double, double> nEvents, pair<double, pair<double, double> >efficiency){
  pair<double, double> survivingPair;
  double eff = efficiency.first;
  pair<double, double> effError = efficiency.second;
  double nSurviving = eff*nEvents.first;
  double errorFirstTerm = nEvents.first*effError.first;
  if(effError.first==0.) errorFirstTerm = nEvents.first*effError.second;
  double errorSecondTerm = eff*nEvents.second;
  double nSurvivingError = sqrt((errorFirstTerm*errorFirstTerm)+(errorSecondTerm*errorSecondTerm));
  survivingPair = make_pair<double, double>(nSurviving, nSurvivingError);  
  return survivingPair;
}

double EMuAnalyzer::GetCentralPDF(vector<double> theWeightsVector){
  return theWeightsVector.at(0);
}

std::pair<double, double> EMuAnalyzer::GetPDFUncertaity(vector<double> theWeightsVector){
  std::pair<double, double> thePDFUncertainty;
  double pdf0 = GetCentralPDF(theWeightsVector);
  double deltaXPlusMax = 0.;
  double deltaXMinusMax = 0.;
  //std::vector<double> pdfWeightedEvents;
  for(unsigned int vIt = 1; vIt < theWeightsVector.size();){
    double xplus = theWeightsVector.at(vIt);
    double xminus = theWeightsVector.at(vIt+1);
    
    double deltaPlusCen = xplus - pdf0;
    double deltaMinusCen = xminus - pdf0;
    double theMaxPlusMinusCen = std::max(std::max(deltaPlusCen,deltaMinusCen),0.);
    deltaXPlusMax += theMaxPlusMinusCen*theMaxPlusMinusCen;
    
    double deltaCenPlus = pdf0 - xplus;
    double deltaCenMinus = pdf0 - xminus;
    double theMaxCenPlusMinus = std::max(std::max(deltaCenPlus,deltaCenMinus),0.);
    deltaXMinusMax += theMaxCenPlusMinus*theMaxCenPlusMinus;
    
    vIt += 2;
  }
  deltaXPlusMax = sqrt(deltaXPlusMax);
  deltaXMinusMax = sqrt(deltaXMinusMax);
  thePDFUncertainty = make_pair<double,double>(deltaXPlusMax, deltaXMinusMax);
  return thePDFUncertainty;
}

double EMuAnalyzer::GetPerCentError(unsigned int nomPassed, unsigned int nomTotal, unsigned int systPassed,
                                    unsigned int systTotal){
  std::pair<double, std::pair<double, double> > theNomEff = GetEfficiency(nomPassed, nomTotal);
  std::pair<double, std::pair<double, double> > theSystEff = GetEfficiency(systPassed, systTotal);
  
  return fabs(theNomEff.first - theSystEff.first)/theNomEff.first;
}

void EMuAnalyzer::GetPUWeights(){
  //fisrt get the data distribution
  TFile* dataPUFile = new TFile("/home/eluiggi/work/cms/hiMassTauTau/eMuAnalysis/rootFiles/pileUpCombination.root");
  dataPUHisto = (TH1F*)dataPUFile->Get("pileup");
  dataPUHisto->Scale(1./dataPUHisto->Integral());
  
  //get MC distribution; copy histo from data
  TH1F* mcPUHisto = (TH1F*)dataPUHisto->Clone("mcPUHisto");
  fChain->Project("mcPUHisto", "InTimePU");
  mcPUHisto->Scale(1./mcPUHisto->Integral());
  TH1F* theWeights = (TH1F*)dataPUHisto->Clone("theWeights");
  theWeights->SetName("PUWeights");
  theWeights->Divide(mcPUHisto);
  for(int iBin = 1; iBin < theWeights->GetNbinsX() + 1; iBin++){
     _thePUWeights.push_back(theWeights->GetBinContent(iBin));
  }
}
