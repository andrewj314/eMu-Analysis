#define EMuAnalysis_cxx
#include "EMuAnalysis.h"
#include "pdfUtils.h"
#include "efficiencies.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <utility>
#include <stdio.h>
#include <vector>
#include <math.h>
using namespace std;

void EMuAnalysis::Loop()
{

   bookHistos();
   if (fChain == 0) return;
   
   if(_nEvents < 0 || _nEvents > fChain->GetEntries()) _nEvents = fChain->GetEntries();
   
   cout << "nEvents = " << _nEvents << "\n";
   cout << "fChain->GetEntries = " <<fChain->GetEntries() << "\n";
   counterEvents = 0;
   
   
   //_nEvents = 10;
   
   
   double eventsFraction = (double)_nEvents/(double)fChain->GetEntries();
   
   cout << "eventsFraction = " << eventsFraction << "\n";
   
   _eventPUWeight = 1.0;
   _puWeightedEvents = 0.;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<_nEvents;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      
      if(_source != "Data" && _source != "ZPrimeSSM_M-500" && _source != "ZPrimeSSM_M-750" && _source != "ZPrimeSSM_M-1000" && _source != "ZPrimeSSM_M-1250" && _source != "ZPrimeSSM_M-1500" && _source != "ZPrimeSSM_M-1750"
       && _source != "ZPrimeSSM_M-2000" && _source != "ZPrimeSSM_M-2250" && _source != "ZPrimeSSM_M-2500"){										       
        if(InTimePU < (int)_thePUWeights.size()){
		_eventPUWeight = _thePUWeights.at(InTimePU);
	}
        if(InTimePU >= (int)_thePUWeights.size() && ((int)_thePUWeights.size() != 0)){
		_eventPUWeight = _thePUWeights.at(_thePUWeights.size() -1);
	}    
      }	
           
      
      
      _puWeightedEvents += _eventPUWeight;  
      _pileUPDistro->Fill(InTimePU);									       
      _pileUPReweighted->Fill(InTimePU,_eventPUWeight);	      
   	            
      resetCandCounters();
      
      
      
      


      
      //cross sections: given in pb
      
      if (_source == "ZPrimeSSM_M-500"){
         _theXSection = make_pair<float, float>(2.505, 2.505*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50181);  
      }
      
      else if (_source == "ZPrimeSSM_M-750"){
         _theXSection = make_pair<float, float>(0.472, 0.472*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50116);  
      }      
      
      else if (_source == "ZPrimeSSM_M-1000"){
         _theXSection = make_pair<float, float>(0.129, 0.129*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50400);  
      }   
      
      else if (_source == "ZPrimeSSM_M-1250"){
         _theXSection = make_pair<float, float>(0.041, 0.041*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50170);  
      }               
      
      else if (_source == "ZPrimeSSM_M-1500"){
         _theXSection = make_pair<float, float>(_signalXSection, _signalXSection*0.02); // 2% error on xsection
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50312); 
      }
      
      else if (_source == "ZPrimeSSM_M-1750"){
         _theXSection = make_pair<float, float>(0.0063, 0.0063*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50400);
      }          
      
      else if (_source == "ZPrimeSSM_M-2000"){
         _theXSection = make_pair<float, float>(0.0027, 0.0027*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50038);
      }          
      
      else if (_source == "ZPrimeSSM_M-2250"){
         _theXSection = make_pair<float, float>(0.0012, 0.0012*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50038);
      }                
      
      else if (_source == "ZPrimeSSM_M-2500"){
         _theXSection = make_pair<float, float>(0.0006, 0.0006*0.02);
	 skimmedEvents = make_pair<int, int>((int)_nEvents, 50038);
      }                
      
      else if (_source == "DYToTauTau"){										       
       _theXSection = make_pair<float, float>(1510, 0); // NNLO Z -> ll (M>20) 			       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(3119494.*eventsFraction));		       
      }													       
      else if (_source == "DYToEE"){										       
       _theXSection = make_pair<float, float>(1871, 0); // NNLO Z -> ll (M>20) 			       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(2659625.*eventsFraction));		       
      }	
      else if (_source == "DYJetsToLeptonLepton"){
       _theXSection = make_pair<float, float>(3504, 0); // NNLO Z -> ll (M>20) 			       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(30459503.*eventsFraction));      
      }												       
      else if (_source == "TTJets"){									       
       _theXSection = make_pair<float, float>(245.8, 0); //MCFM NLO 					       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(6923750.*eventsFraction));		       
      }
      else if (_source == "TTJets_MADGRAPH"){									       
       _theXSection = make_pair<float, float>(234, 0); //MCFM NLO 					       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(6474753.*eventsFraction));		       
      }
      else if (_source == "TTJets_MADGRAPH_Mtt-700to1000"){									       
       _theXSection = make_pair<float, float>(17.3, 0); //MCFM NLO 					       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(3082812.*eventsFraction));		       
      }
      else if (_source == "TTJets_MADGRAPH_Mtt-1000toInf"){									       
       _theXSection = make_pair<float, float>(3.28, 0); //MCFM NLO 					       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(1249111.*eventsFraction));		       
      }																			       
      else if(_source == "WJets"){									       
       _theXSection = make_pair<float, float>(36257.2, 0);  //						       
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(18393090.*eventsFraction)); 	       
      }
      else if(_source == "WJets_Pt-50to70"){
       _theXSection = make_pair<float, float>(967.5, 0); //
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(48426609.*eventsFraction));
      }
      else if(_source == "WJets_Pt-70to100"){
       _theXSection = make_pair<float, float>(511.5, 0); //
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(22447541.*eventsFraction));
      }      		
      else if(_source == "WJets_Pt-100"){
       _theXSection = make_pair<float, float>(273.0, 0); //
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(12742382.*eventsFraction));
      }											       
      else if(_source == "WW"){										       
       _theXSection = make_pair<float, float>(54.838, 0); // MCFM NLO inclusive				       
       skimmedEvents = make_pair<int, int>((int)_nEvents, (int)(10000431.*eventsFraction));		       
      }													       
      else if(_source == "WZ"){										       
       _theXSection = make_pair<float, float>(33.2, 0);  //MCFM NLO M_ll > 40 inclusive 		       
       //_theXSection = make_pair<float, float>(0.7192, 0.0277);  //WZ->3leptons+nu with 3.8% error	       
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(10000283.*eventsFraction));		       
      }													       
      else if(_source == "ZZ"){										       
       _theXSection = make_pair<float, float>(17.7, 0);  // NLO inclusive				       
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(9799908.*eventsFraction));		       
      }		      
      else if(_source == "QCD"){
       _theXSection = make_pair<float, float>(3.64E8*3.7E-4, 0); 
       skimmedEvents = make_pair<double, double>((int)_nEvents, (int)(21484602.*eventsFraction));
      }
      
     //get events normilized to lumi
     eventsAtLumi = make_pair<double, double>(_theXSection.first * _theLumi, _theXSection.second * _theLumi);    									          
     skimEff = GetEfficiency(skimmedEvents.first, skimmedEvents.second); 				       
     nSkimmedAtLumi = GetSurvivingEvents(eventsAtLumi, skimEff); 	      
            
      
      resetHistosDefault();
   
     
   
   
      cout << "\n-------Event number: " << jentry << "-------\n\n";
      getCandCounters();
      //cout << "GetEventCounters\n";
      getEventCounters();
      
      fillHistos();      	 
   }
   
   getReport();
   writeOutFile();
   
}

void EMuAnalysis::getCandCounters(){
  unsigned int vItBegin = 0;
  unsigned int vItEnd = elecMatched->size();
  unsigned int tauItBegin = 0;
  unsigned int tauItEnd = min(tauEnergy->size(), muonPx->size());
  unsigned int tauVtxBegin = 0;
  unsigned int tauVtxEnd = genTauVtx_x->size();
  unsigned int muonVtxBegin = 0;
  unsigned int muonVtxEnd = muonVertex_x0->size();
  unsigned int elecVtxBegin = 0;
  unsigned int elecVtxEnd1 = min(elecVertex_x0->size(), elecVertex_y0->size());
  unsigned int elecVtxEnd = min(elecVtxEnd1, (unsigned int)elecVertex_z0->size());
  //cout << "ElecVtxEnd = " << elecVtxEnd << "\n";
  unsigned int zPrimeItBegin = 0;
  unsigned int zPrimeVtx1 = min(zPrimeVtx_x->size(), zPrimeVtx_y->size());
  unsigned int zPrimeItEnd1 = min(zPrimeVtx1, (unsigned int)zPrimeVtx_z->size());
  unsigned int zPrimeItEnd = min(zPrimeItEnd1, (unsigned int)BSx0->size());
  
  //cout << "vItBegin = " << vItBegin << " vItEnd = " << vItEnd << "\n";
  //cout << "tauItEnd = " << tauItEnd << "\n";  
  //counterEvents++;
  
  if(_source == "ZPrimeSSM_M-500" || _source == "ZPrimeSSM_M-750" || _source == "ZPrimeSSM_M-1000" || _source == "ZPrimeSSM_M-1250" || _source == "ZPrimeSSM_M-1500" || _source == "ZPrimeSSM_M-1750" || _source == "ZPrimeSSM_M-2000" || 
  _source == "ZPrimeSSM_M-2250" || _source == "ZPrimeSSM_M-2500"){
      //cout << "ZPRIME\n";
      for(unsigned int mIt = 0; mIt < elecPt->size(); ++mIt){
        if(elecMatched->at((mIt))) _elecMatchedCounter++;
        if(!(elecMatched->at(mIt) && muonMatched->at(mIt))) continue;
	
	_muonMatchedCounter++;
	_matchedCounter++;
        
      }  
      
      int theCand = getMatchedCand();
      if(theCand >= 0){
	vItBegin = theCand;
	vItEnd = vItBegin + 1;
      } 
      							     
  }      
    
  for(unsigned int vIt = vItBegin; vIt < vItEnd; vIt++){ 
	
    getNMinus1CandCounters(vIt);	
    float IsoCorrection = muonPFIsoDR04SumChargedHadronPt->at(vIt) + muonPFIsoDR04SumNeutralHadronPt->at(vIt) + muonPFIsoDR04SumPhotonPt->at(vIt) - 0.5*muonPFIsoDR04SumPUPt->at(vIt);
    //if(IsoCorrection < 0) IsoCorrection = 0;
    float MuonIso = (IsoCorrection)/muonPt->at(vIt); 
    float ElecRelIso = (elecTrkIso->at(vIt) + elecEcalIso->at(vIt) + elecHcalIso->at(vIt))/elecPt->at(vIt);       
    //cout << "Got NMinus1CandCounters\n";
    
    //INCREMENT CONTROL REGION COUNTERS
    /*
    if(passedQCDCR1(vIt)) _QCDCR1Counter++;
    if(passedQCDCR1(vIt) && passedCharge(vIt)) _QCDOSCounter++;
    //if(passedQCDCR1(vIt) && passedElecRelIso(vIt)) _QCDElectronIsoCounter++;
    if(passedQCDCR1(vIt) && passedElecMuonPZeta(vIt) && passedMET(vIt) && passedLeptonMetCosDPhi(vIt)) _QCDTopologyCounter++;
    
    if(passedQCDCR2(vIt)) _QCDCR2Counter++;
    if(passedQCDCR2(vIt) && passedMuonIsoDR04DB(vIt)) _QCDMuonIsoCounter++;
    if(passedQCDCR2(vIt) && !passedMuonIsoDR04DB(vIt)) _QCDMuonAntiIsoCounter++;
    
    
    if(passedTTJetsCR1(vIt)) _TTJetsCR1Counter++;
    if(passedTTJetsCR1(vIt) && passedJetSumEt(vIt)) _TTJetsJetSumEtCounter++;
    if(passedTTJetsCR1(vIt) && passedElecMuonCosDPhi(vIt) && passedLeptonMetCosDPhi(vIt)) _TTJetsTopologyCounter++;
    if(passedTTJetsCR2(vIt) && passedNbTagsCombSecVtx(vIt)) _TTJetsNoBTagsCounter++;
    if(passedTTJetsCR2(vIt) && !passedNbTagsCombSecVtx(vIt)) _TTJetsAtLeastOneBTagCounter++;

    if(passedTTJetsCR2(vIt)) {
      _TTJetsCR2Counter++;
      if(nBTagsCombSecVtx->at(vIt) < _minNBtagsCombSecVtxCR2) _minNBtagsCombSecVtxCR2 = nBTagsCombSecVtx->at(vIt);
    }
    
    if(passedWJetsCR1(vIt)) _WJetsCR1Counter++;
    if(passedWJetsCR1(vIt) && passedCharge(vIt)) _WJetsOSCounter++;
    if(passedWJetsCR1(vIt) && passedElecMuonPZeta(vIt) && passedElecMuonCosDPhi(vIt)) _WJetsTopologyCounter++;
    if(passedWJetsCR2(vIt)) _WJetsCR2Counter++;
    //if(passedWJetsCR2(vIt) && ElecRelIso < 0.1) _WJetsElectronIsoCounter++;
    //if(passedWJetsCR2(vIt) && ElecRelIso > 0.1 && ElecRelIso < 0.2) _WJetsElectronAntiIsoCounter++;
    if(passedWJetsCR2(vIt) && muonMEtMt->at(vIt) > 70 && muonMEtMt->at(vIt) < 120) _WJetsMtCounter++;
    
    if(passedDYToTauTauCR1(vIt)) _DYToTauTauCR1Counter++;
    if(passedDYToTauTauCR2(vIt)) _DYToTauTauCR2Counter++;
    
    */
    /*
    if (_source == "ZPrimeSSM_M-500" || _source == "ZPrimeSSM_M-750" || _source == "ZPrimeSSM_M-1000" || _source == "ZPrimeSSM_M-1250" || _source == "ZPrimeSSM_M-1500" || _source == "ZPrimeSSM_M-1750" || _source == "ZPrimeSSM_M-2000" || 
  _source == "ZPrimeSSM_M-2250" || _source == "ZPrimeSSM_M-2500") {
       if(!(elecMatched->at(vIt) && muonMatched->at(vIt))) continue;
    }
    */
   if(elecMuonMetMass->at(vIt) > _maxElecMuonMetMassBeginning) _maxElecMuonMetMassBeginning = elecMuonMetMass->at(vIt); 
    
   if(!passedElecMuonMetMass(vIt)) continue;
    
   //if(!passedMatched(vIt)) continue; 
    
   /*
   cout << "Muon pt = " << muonPt->at(vIt) << "\n";
   cout << "Muon px = " << muonPx->at(vIt) << "\n";
   cout << "Muon py = " << muonPy->at(vIt) << "\n";
   cout << "Muon pz = " << muonPz->at(vIt) << "\n";
   
   cout << "Elec pt = " << elecPt->at(vIt) << "\n";
   cout << "Elec px = " << elecPx->at(vIt) << "\n";
   cout << "Elec py = " << elecPy->at(vIt) << "\n";
   cout << "Elec pz = " << elecPz->at(vIt) << "\n";
   
   cout << "Met = " << mEt << "\n";
   cout << "MetPx = " << mEtPx << "\n";
   cout << "MetPy = " << mEtPy << "\n";
   cout << "MetPz = " << mEtPz << "\n";
   */
   //cout << "DeltaPtx = " << deltaPtx->at(vIt) << "\n";
   //cout << "DeltaPty = " << deltaPty->at(vIt) << "\n";
   //cout << "DeltaPtz = " << deltaPtz->at(vIt) << "\n";
   
   float deltaPtx = -(muonPx->at(vIt) + elecPx->at(vIt));
   float deltaPty = -(muonPy->at(vIt) + elecPy->at(vIt));
   float deltaPt = sqrt(deltaPtx*deltaPtx + deltaPty*deltaPty);
   
 //  cout << "MEt = " << mEt << " deltaPt = " << deltaPt << " difference = " << mEt - deltaPt << "\n";
   
   if(fabs(mEt - deltaPt) > _maxMetMinusDeltaPt) _maxMetMinusDeltaPt = fabs(mEt - deltaPt);
   if(muonPt->at(vIt) > elecPt->at(vIt)){
     _maxDeltaPtx = muonPhi->at(vIt) - atan2(deltaPty,deltaPtx);
     _maxDeltaPty = deltaPty/deltaPtx;
     _maxMetPx = muonPhi->at(vIt) - atan2(mEtPy,mEtPx);
     _maxMetPy = mEtPy/mEtPx;
  }
  
  else{
     _maxDeltaPtx = elecPhi->at(vIt) - atan2(deltaPty,deltaPtx);
     _maxDeltaPty = deltaPty/deltaPtx;
     _maxMetPx = elecPhi->at(vIt) - atan2(mEtPy,mEtPx);
     _maxMetPy = mEtPy/mEtPx;  
  
  }
  
        //Muons

    
    
    if(muonPt->at(vIt) > _maxMuonPt) _maxMuonPt = muonPt->at(vIt);
    if(!passedMuonPt(vIt)) {
    	continue;
    }
    _muonPtCounter++;
    
 
    if(fabs(muonEta->at(vIt)) < _minMuonEta) _minMuonEta = muonEta->at(vIt);
    if(!passedMuonEta(vIt)){
    	continue;
    }
    _muonEtaCounter++;
    
    
        

	
    if(passedMuonAcc(vIt)) _muonAccCounter++;
    
    if(muonIsGlobalMuon->at(vIt) > _maxIsGlobalMuon) _maxIsGlobalMuon = muonIsGlobalMuon->at(vIt);
    if(!passedMuonIsGlobalMuon(vIt)) continue;
    _muonIsGlobalMuonCounter++;	     
    
    if(elecPt->at(vIt) > _maxElecPt) _maxElecPt = elecPt->at(vIt);
    if(elecEt->at(vIt) > _maxElecEt) _maxElecEt = elecEt->at(vIt);
    if(!passedElecEt(vIt)) continue;
    _elecEtCounter++;	
    
    if(fabs(elecSCEta->at(vIt)) < _minElecEta) _minElecEta = elecSCEta->at(vIt);
    if(!passedElecMaxEta(vIt)){
    	continue;
    }
    _elecEtaCounter++;   	
    
    if(passedElecAcc(vIt)) _elecAccCounter++;

    if(heepPassedEcalDriven->at(vIt) > _maxHeepPassedEcalDriven) _maxHeepPassedEcalDriven = heepPassedEcalDriven->at(vIt);
    
    if(!passedHeepPassedEcalDriven(vIt)){
    	continue;
    }
    _heepPassedEcalDrivenCounter++;	
    
    if(elecMuonDeltaR->at(vIt) > _maxElecMuonDeltaR) _maxElecMuonDeltaR = elecMuonDeltaR->at(vIt);
    if(!passedElecMuonDeltaR(vIt)){
    	continue;
    }
    _elecMuonDeltaRCounter++; 
    
   
    
    if(muonDxyVtx->at(vIt) > _maxMuonDxyVtx) _maxMuonDxyVtx = muonDxyVtx->at(vIt);
    if(!passedMuonDxyVtx(vIt)) continue;
    _muonDxyVtxCounter++;   

    if(muonDzVtx->at(vIt) > _maxMuonDzVtx) _maxMuonDzVtx = muonDzVtx->at(vIt);
    if(!passedMuonDzVtx(vIt)) continue;
    _muonDzVtxCounter++; 
    
    if(muonValidHits->at(vIt) > _maxMuonValidHits) _maxMuonValidHits = muonValidHits->at(vIt);
    if(!passedMuonValidHits(vIt)) continue;
    _muonValidHitsCounter++;
    
    if(muonPixelHits->at(vIt) > _maxMuonPixelHits) _maxMuonPixelHits = muonPixelHits->at(vIt);
    if(!passedMuonPixelHits(vIt)) continue;
    _muonPixelHitsCounter++; 
       
    if(muonMatchedStations->at(vIt) > _maxMuonMatchedStations) _maxMuonMatchedStations = muonMatchedStations->at(vIt);
    if(!passedMuonMatchedStations(vIt)) continue;
    _muonMatchedStationsCounter++;      
    
    if(muonTrkLayersWithHits->at(vIt) > _maxMuonTrkLayersWithHits) _maxMuonTrkLayersWithHits = muonTrkLayersWithHits->at(vIt);
    if(!passedMuonTrkLayersWithHits(vIt)) continue;
    _muonTrkLayersWithHitsCounter++;    
    
    if(!passedMuonDptOverPt(vIt)) continue;
    _muonDptOverPtCounter++;
    
    if(muonNormChiSqrd->at(vIt) > _maxMuonNormChiSqrd) _maxMuonNormChiSqrd = muonNormChiSqrd->at(vIt);
    if(passedMuonNormChiSqrd(vIt)) _muonNormChiSqrdCounter++;
    

    

    

    

    
    
    /*
    if(!passedMuonIsHighPtMuon(vIt)){
    	continue; 
    }
    _muonIsHighPtMuonCounter++;    
    */
    
    if(MuonIso > _maxMuonIsoDR04DB) _maxMuonIsoDR04DB = MuonIso;
    if(!passedMuonIsoDR04DB(vIt)){
        continue;
    }
    _muonIsoDR04DBCounter++;
 
 
    if(passedAllMuonIDCuts(vIt)) _muonIdCounter++;    
    
    
    
    
    //Electrons
    
    //Check to see if there is an electron
    //if(!elecExists(vIt)) continue;
    //_elecExistsCounter++;
    
    
    //if(elecPt->at(vIt) > _maxElecPt) _maxElecPt = elecPt->at(vIt);						     
    //if(!passedElecMinPt(vIt)){
    //	continue;
    //}
    //_elecPtCounter++;
    
    
  
    
    
    //if(passedElecAcc(vIt)) _elecAccCounter++;

  
    
    //HEEP

   

    if(heepPassedDEtaIn->at(vIt) > _maxHeepPassedDEtaIn) _maxHeepPassedDEtaIn = heepPassedDEtaIn->at(vIt);  
    if(!passedHeepPassedDEtaIn(vIt)){
    	continue;
    }
    _heepPassedDEtaInCounter++;
 
   
    if(heepPassedDPhiIn->at(vIt) > _maxHeepPassedDPhiIn) _maxHeepPassedDPhiIn = heepPassedDPhiIn->at(vIt);  
    if(!passedHeepPassedDPhiIn(vIt)){
    	continue;
    }
    _heepPassedDPhiInCounter++;
    
        
    if(heepPassedHadem->at(vIt) > _maxHeepPassedHadem) _maxHeepPassedHadem = heepPassedHadem->at(vIt);
    if(!passedHeepPassedHadem(vIt)){
    	continue;
    }
    _heepPassedHademCounter++;
    
    
    if(heepPassedSigmaIEtaIEta->at(vIt) > _maxHeepPassedSigmaIEtaIEta) _maxHeepPassedSigmaIEtaIEta = heepPassedSigmaIEtaIEta->at(vIt);
    if(!passedHeepPassedSigmaIEtaIEta(vIt)){
    	continue;
    }
    _heepPassedSigmaIEtaIEtaCounter++;
    
    
    if(heepPassed2by5Over5By5->at(vIt) > _maxHeepPassed2by5Over5By5) _maxHeepPassed2by5Over5By5 = heepPassed2by5Over5By5->at(vIt);
    if(!passedHeepPassed2by5Over5By5(vIt)){
    	//cout << "Did not pass HeepPassed2by5Over5By5 cut, continuing\n";
    	continue;
    }
    _heepPassed2by5Over5By5Counter++;
       
    
    if(!passedElecEcalHad1Iso(vIt)) continue;
    _heepPassedEcalHad1IsoCounter++;
         
	 
	/* 
    if(!passedHeepPassedTrkIso(vIt)) continue;
    _heepPassedTrkIsoCounter++; 
         */
    
    if(!passedElecTrkIso(vIt)){
        continue;
    }
    _elecTrkIsoCounter++;    
    
	 
	 
    if(!passedHeepPassedMissingHits(vIt)) continue;
    _heepPassedMissingHitsCounter++;
        
    if(!passedElecDxy(vIt)) continue;
    _elecDxyCounter++;    
    
    
    

    /*
    if(ElecRelIso < _minElecRelIso) _minElecRelIso = ElecRelIso;
    if(!passedElecRelIso(vIt)){
        continue;
    }
    _elecRelIsoCounter++;    
    
    if(!passedElecMissingHits(vIt)) continue;
    _elecMissHitsCounter++; 
   */
    
    if(passedHeepId(vIt)) _heepIdCounter++;
  
    

    
    //Topology
 
       
       
    if(!passedCharge(vIt)){
        continue;
    }
    _chargeCounter++;   


     
   if(elecMuonCosDPhi->at(vIt) < _minElecMuonCosDPhi) _minElecMuonCosDPhi = elecMuonCosDPhi->at(vIt);
    if(!passedElecMuonCosDPhi(vIt)){
    	continue;
    }    
    _elecMuonCosDPhiCounter++;
   
   
    if(mEt > _maxMEt) _maxMEt = mEt;	
    
    
 
    
    			  
    if(!passedMET(vIt)){
    	continue;
    } 		  
    _metCounter++;  
      
    if(!passedNbTagsCombSecVtx(vIt)){
    	continue;
    }
    _nbTagsCombSecVtxCounter++;      
    

    
    if(elecMuonPZeta->at(vIt) > _maxElecMuonPZeta) _maxElecMuonPZeta = elecMuonPZeta->at(vIt);
    if(!passedElecMuonPZeta(vIt)){
    	continue;
    }
    _elecMuonPZetaCounter++;    
      
    
 
    
    
    if(leptonMetCosDphi->at(vIt) < _minLeptonMetCosDPhi) _minLeptonMetCosDPhi = leptonMetCosDphi->at(vIt);
    if(!passedLeptonMetCosDPhi(vIt)){
    	continue;
    }
    _leptonMetCosDPhiCounter++;    
    
    if(leptonDeltaPtCosDphi->at(vIt) < _minLeptonDeltaPtCosDPhi) _minLeptonDeltaPtCosDPhi = leptonDeltaPtCosDphi->at(vIt);
    if(!passedLeptonDeltaPtCosDPhi(vIt)) continue;
    _leptonDeltaPtCosDPhiCounter++;
    
    
    if(!passedThreeLepton(vIt)) continue;
    _threeLeptonCounter++;
    
    
     
    if(!passedJetSumEt(vIt)){
        continue;
    }
    _jetSumEtCounter++;
    
    
    if(!passedNJets(vIt)){
        continue;
    }
    _nJetsCounter++;    
    
    
    if(!passedMuonMetMt(vIt)){
    	continue;
    }
    _muonMetMtCounter++;    
    
    
    if(!passedElecMuonMass(vIt)){
    	continue;
    }
    _elecMuonMassCounter++;
     
    if(!passedElecMetMt(vIt)){
    	continue;
    }
    _elecMetMtCounter++;     
  
    
    if(elecMuonMetMass->at(vIt) > _maxElecMuonMetMass) _maxElecMuonMetMass = elecMuonMetMass->at(vIt);

    
    if(!passedElecMuonMetMass(vIt)){
    	continue;
    }
    _elecMuonMetMassCounter++;
    
    if(elecMuonDeltaPtMass->at(vIt) > _maxElecMuonDeltaPtMass) _maxElecMuonDeltaPtMass = elecMuonDeltaPtMass->at(vIt);
    if(elecE->at(vIt) > _maxElecPtVsMass) _maxElecPtVsMass = elecE->at(vIt);
    if(muonE->at(vIt) > _maxMuonPtVsMass) _maxMuonPtVsMass = muonE->at(vIt); 
    
    if(!passedElecMuonDeltaPtMass(vIt)){
    	continue;
    }
    _elecMuonDeltaPtMassCounter++;      
      
    
    if(!passed2xPLargerLepton(vIt)) continue;
    _2xPLargerLeptonCounter++;
    
    //Values for histos for MEt vs DeltaPt(e, mu)
    
    if(mEt > _maxMEtVsDeltaPt) _maxMEtVsDeltaPt = mEt;
    if(elecPt->at(vIt) > _maxElecPtVsDeltaPt) _maxElecPtVsDeltaPt = elecPt->at(vIt);
    if(muonPt->at(vIt) > _maxMuonPtVsDeltaPt) _maxMuonPtVsDeltaPt = muonPt->at(vIt);        
    
    if(!passedElecMuonDeltaPt(vIt)) continue;
    _elecMuonDeltaPtCounter++;
    
    //cout << "Muon Mother Id = " << muonMotherId->at(vIt) << "\n";
    //cout << "Elec Mother Id = " << eMotherId->at(vIt) << "\n";
    if(eMotherId->at(vIt) > _maxElecMotherId) _maxElecMotherId = eMotherId->at(vIt);
    if(muonMotherId->at(vIt) > _maxMuonMotherId) _maxMuonMotherId = muonMotherId->at(vIt);
    
    
    //Tau veto - if either the muon or the electron comes from a tau, veto the candidate. Used for diagnostics of lifetime cuts
   /*
    if(_source != "ZPrime" && _source.find("SSM") == std::string::npos){
      if(passedMuonMotherId(vIt)) continue;
      if(passedElecMotherId(vIt)) continue;
      //if((passedElecMotherId(vIt) && passedMuonMotherId(vIt)) || (!passedElecMotherId(vIt) && !passedMuonMotherId(vIt))) continue;  
    }
    */
    //Impact Parameter Study
    float ElecIpOverElecIpErr = fabs(elecGSFBSIp->at(vIt)/elecGSFBSIpError->at(vIt));
    float MuonIpOverMuonIpErr = fabs(muonDxyBS->at(vIt)/muonDxyBSError->at(vIt)); 
    float ElecIpError = elecGSFBSIpError->at(vIt);
    float MuonIpError = muonDxyBSError->at(vIt);
    float ElecIp = fabs(elecGSFBSIp->at(vIt));
    float MuonIp = fabs(muonDxyBS->at(vIt));
    float CombinedIpOverIpErr = sqrt((ElecIp*ElecIp + MuonIp*MuonIp + 2*ElecIp*MuonIp)/(ElecIpError*ElecIpError + MuonIpError*MuonIpError));
    //float CombinedIpOverIpErr = ElecIpOverElecIpErr*MuonIpOverMuonIpErr;
    

    
    
    
    if(elecGSFBSIpError->at(vIt) > _maxElecIpError) _maxElecIpError = elecGSFBSIpError->at(vIt);
    if(muonDxyBSError->at(vIt) > _maxMuonIpError) _maxMuonIpError = muonDxyBSError->at(vIt);
    if(elecMissingHits->at(vIt) > _maxElecMissingHits) _maxElecMissingHits = elecMissingHits->at(vIt);
    
    
    if(CombinedIpOverIpErr > _maxCombinedIpOverIpError) _maxCombinedIpOverIpError = CombinedIpOverIpErr;
    if(!passedCombinedIpOverIpError(vIt)){
        continue;
    }
    _combinedIpOverIpErrorCounter++;    
    

    if(ElecIpOverElecIpErr > _maxElecIpOverElecIpError) _maxElecIpOverElecIpError = ElecIpOverElecIpErr;
    if(!passedElecIpOverElecIpError(vIt)){
        continue;
    }
    _elecIpOverElecIpErrorCounter++;
    
    if(MuonIpOverMuonIpErr > _maxMuonIpOverMuonIpError) _maxMuonIpOverMuonIpError = MuonIpOverMuonIpErr;    
    if(!passedMuonIpOverMuonIpError(vIt)){
        continue;
    }
    _muonIpOverMuonIpErrorCounter++;
    

    
    
    float MuonTrkChiSqrd = fabs(muonTrkChiSqrd->at(vIt));
    float ElecTrkChiSqrd = fabs(elecTrkChiSqrd->at(vIt));
    
    if(ElecIp > _maxElecIp) _maxElecIp = ElecIp;
    if(MuonIp > _maxMuonIp) _maxMuonIp = MuonIp;
    
    if(ElecTrkChiSqrd < _minElecTrkChiSqrd) _minElecTrkChiSqrd = ElecTrkChiSqrd;
    if(MuonTrkChiSqrd < _minMuonTrkChiSqrd) _minMuonTrkChiSqrd = MuonTrkChiSqrd;
    
    if(fabs(muonTrkNdof->at(vIt)) > _maxMuonTrkNdof) _maxMuonTrkNdof = fabs(muonTrkNdof->at(vIt));
    if(fabs(elecTrkNdof->at(vIt)) > _maxElecTrkNdof) _maxElecTrkNdof = fabs(elecTrkNdof->at(vIt));
    

    
    
    
    if(elecMuonDistanceOfClosestApproach->at(vIt) > _maxElecMuonDistanceOfClosestApproach){

       _maxElecMuonDistanceOfClosestApproach = elecMuonDistanceOfClosestApproach->at(vIt);
       _maxElecMuonDCACrossingPoint = elecMuonDCACrossingPoint->at(vIt);   


    }    
    
    
    if(elecMuonTwoTrkMinDist->at(vIt) > _maxElecMuonTwoTrkMinDist){
       
       _maxElecMuonTwoTrkMinDist = elecMuonTwoTrkMinDist->at(vIt);
       _maxElecMuonTwoTrkCrossingPoint = elecMuonTwoTrkCrossingPoint->at(vIt);
    
    }
    
    if(eMuDis3D->at(vIt) > _maxElecMuonDis3D) _maxElecMuonDis3D = eMuDis3D->at(vIt);
    if(eMuDis2D->at(vIt) > _maxElecMuonDis2D) _maxElecMuonDis2D = eMuDis2D->at(vIt);
    
    if(elecMuonDCA3DError->at(vIt) > _maxElecMuonDCA3DError) _maxElecMuonDCA3DError = elecMuonDCA3DError->at(vIt);
    if(elecMuonDCA2DError->at(vIt) > _maxElecMuonDCA2DError) _maxElecMuonDCA2DError = elecMuonDCA2DError->at(vIt);
    
    float ElecMuonDis3DOverError = elecMuonTwoTrkMinDist->at(vIt)/elecMuonDCA3DError->at(vIt);
    float ElecMuonDis2DOverError = elecMuonTwoTrkMinDist->at(vIt)/elecMuonDCA2DError->at(vIt);
    
    if(ElecMuonDis3DOverError > _maxElecMuonDis3DOverError) _maxElecMuonDis3DOverError = ElecMuonDis3DOverError;
    
    if(!passedElecMuonDCAOverError3D(vIt)) continue;
    _elecMuonDCAOverError3DCounter++;
    
    if(ElecMuonDis2DOverError > _maxElecMuonDis2DOverError) _maxElecMuonDis2DOverError = ElecMuonDis2DOverError;
    if(!passedElecMuonDCAOverError2D(vIt)) continue;
    _elecMuonDCAOverError2DCounter++;
    
     if(!passedLifetime(vIt)) continue;
    _lifetimeCounter++;
        
    for(unsigned int tauIt = tauItBegin; tauIt < tauItEnd; tauIt++){ 
       if(tauEnergy->at(tauIt) > _maxTauEnergy) _maxTauEnergy = tauEnergy->at(tauIt); 
       
       float MuonP = sqrt(muonPx->at(tauIt)*muonPx->at(tauIt) + muonPy->at(tauIt)*muonPy->at(tauIt) + muonPz->at(tauIt)*muonPz->at(tauIt));
       
       if(MuonP/tauEnergy->at(tauIt) > _maxMuonPOverTauE) _maxMuonPOverTauE = MuonP/tauEnergy->at(tauIt);
       
       //if(genTauVtx_x->at(tauIt) > _maxGenTauVtx_x) _maxGenTauVtx_x = genTauVtx_x->at(tauIt);
       //if(genTauVtx_y->at(tauIt) > _maxGenTauVtx_y) _maxGenTauVtx_y = genTauVtx_y->at(tauIt);
       //if(genTauVtx_z->at(tauIt) > _maxGenTauVtx_z) _maxGenTauVtx_z = genTauVtx_z->at(tauIt);
       
    } 
    
    for(unsigned int tauVtxIt = tauVtxBegin; tauVtxIt < tauVtxEnd; tauVtxIt++){
       if(genTauVtx_x->at(tauVtxIt) > _maxGenTauVtx_x) _maxGenTauVtx_x = genTauVtx_x->at(tauVtxIt);
       if(genTauVtx_y->at(tauVtxIt) > _maxGenTauVtx_y) _maxGenTauVtx_y = genTauVtx_y->at(tauVtxIt);
       if(genTauVtx_z->at(tauVtxIt) > _maxGenTauVtx_z) _maxGenTauVtx_z = genTauVtx_z->at(tauVtxIt);       
       cout << "Tau Vtx x: " << _maxGenTauVtx_x << "\n";
       cout << "Tau Vtx y: " << _maxGenTauVtx_y << "\n";
       cout << "Tau Vtx z: " << _maxGenTauVtx_z << "\n";    
       cout << "Tau Px: " << genTauPx->at(tauVtxIt) << "\n";
       cout << "Tau Py: " << genTauPy->at(tauVtxIt) << "\n";
       cout << "Tau Pz: " << genTauPz->at(tauVtxIt) << "\n";
       cout << "Tau Charge: " << genTauCharge->at(tauVtxIt) << "\n";

       if(otherGenTauVtx_x->at(tauVtxIt) > _maxOtherGenTauVtx_x) _maxOtherGenTauVtx_x = otherGenTauVtx_x->at(tauVtxIt);
       if(otherGenTauVtx_y->at(tauVtxIt) > _maxOtherGenTauVtx_y) _maxOtherGenTauVtx_y = otherGenTauVtx_y->at(tauVtxIt);
       if(otherGenTauVtx_z->at(tauVtxIt) > _maxOtherGenTauVtx_z) _maxOtherGenTauVtx_z = otherGenTauVtx_z->at(tauVtxIt);       
       cout << "Tau Vtx x: " << _maxOtherGenTauVtx_x << "\n";
       cout << "Tau Vtx y: " << _maxOtherGenTauVtx_y << "\n";
       cout << "Tau Vtx z: " << _maxOtherGenTauVtx_z << "\n";               
       cout << "Tau Px: " << otherGenTauPx->at(tauVtxIt) << "\n";
       cout << "Tau Py: " << otherGenTauPy->at(tauVtxIt) << "\n";
       cout << "Tau Pz: " << otherGenTauPz->at(tauVtxIt) << "\n";
       cout << "Tau Charge: " << otherGenTauCharge->at(tauVtxIt) << "\n";

      
    }
    
   
    //cout << "tauIteration\n";
    
    for(unsigned int zPrimeIt = zPrimeItBegin; zPrimeIt < zPrimeItEnd; zPrimeIt++){
    
       //if(zPrimeEnergy->at(zPrimeIt) > _maxZPrimeEnergy) _maxZPrimeEnergy = zPrimeEnergy->at(zPrimeIt);
       float vertexDiff_x = zPrimeVtx_x->at(zPrimeIt) - primaryVertex_x;
       float vertexDiff_y = zPrimeVtx_y->at(zPrimeIt) - primaryVertex_y;
       float vertexDiff_z = zPrimeVtx_z->at(zPrimeIt) - primaryVertex_z;
       float vertexDiff = sqrt(vertexDiff_x*vertexDiff_x + vertexDiff_y*vertexDiff_y + vertexDiff_z*vertexDiff_z);
       float vertexDiffxy = sqrt(vertexDiff_x*vertexDiff_x + vertexDiff_y*vertexDiff_y);
       if(vertexDiff > _maxVertexDiff) _maxVertexDiff = vertexDiff;
       if(vertexDiffxy > _maxVertexDiffxy) _maxVertexDiffxy = vertexDiffxy;
       
       cout << "Generated ZPrime Vtx x = " << zPrimeVtx_x->at(zPrimeIt) << "\n";
       cout << "Generated ZPrime Vtx y = " << zPrimeVtx_y->at(zPrimeIt) << "\n";
       cout << "Generated ZPrime Vtx z = " << zPrimeVtx_z->at(zPrimeIt) << "\n";
       
       
       float beamSpotDiff_x = zPrimeVtx_x->at(zPrimeIt) - BSx0->at(zPrimeIt);
       float beamSpotDiff_y = zPrimeVtx_y->at(zPrimeIt) - BSy0->at(zPrimeIt);
       float beamSpotDiff_z = zPrimeVtx_z->at(zPrimeIt) - BSz0->at(zPrimeIt);
       //cout << "beamSpotDiff_x = " << beamSpotDiff_x << "\n";
       //cout << "zPrimeVtx_x = " << zPrimeVtx_x->at(zPrimeIt) << "\n";
       //cout << "beamSpot_x0 = " << beamSpot_x0 << "\n";
       float beamSpotDiff = sqrt(beamSpotDiff_x*beamSpotDiff_x + beamSpotDiff_y*beamSpotDiff_y + beamSpotDiff_z*beamSpotDiff_z);
       float beamSpotDiffxy = sqrt(beamSpotDiff_x*beamSpotDiff_x + beamSpotDiff_y*beamSpotDiff_y);
       

       //cout << "Beam Spot Diff = " << beamSpotDiff << "\n";
       if(beamSpotDiff > _maxBeamSpotDiff) _maxBeamSpotDiff = beamSpotDiff;
       if(beamSpotDiffxy > _maxBeamSpotDiffxy) _maxBeamSpotDiffxy = beamSpotDiffxy;
       float beamSpotPrimaryzDiff = sqrt(beamSpotDiff_x*beamSpotDiff_x + beamSpotDiff_y*beamSpotDiff_y + vertexDiff_z*vertexDiff_z);
       if(beamSpotPrimaryzDiff > _maxBeamSpotPrimaryzDiff) _maxBeamSpotPrimaryzDiff = beamSpotPrimaryzDiff;
              
       
    }
   
   
   // cout << "ZPrime Iteration\n";
    
    if(fabs(elecEOverP->at(vIt)) > _maxElecEOverP_hiIP) _maxElecEOverP_hiIP = fabs(elecEOverP->at(vIt));
    if(elecSCEta->at(vIt) > _maxElecEta_hiIP) _maxElecEta_hiIP = elecSCEta->at(vIt);
    if(elecEt->at(vIt) > _maxElecEt_hiIP) _maxElecEt_hiIP = elecEt->at(vIt);
    if(elecTrkChiSqrd->at(vIt) > _maxElecTrkChiSqrd_hiIP) _maxElecTrkChiSqrd_hiIP = elecTrkChiSqrd->at(vIt);
    if(!passedElecEOverP(vIt)) continue;
    _elecEOverPCounter++;
    float ElecP = sqrt(elecPx->at(vIt)*elecPx->at(vIt) + elecPy->at(vIt)*elecPy->at(vIt) + elecPz->at(vIt)*elecPz->at(vIt));
    float MuonP = sqrt(muonPx->at(vIt)*muonPx->at(vIt) + muonPy->at(vIt)*muonPy->at(vIt) + muonPz->at(vIt)*muonPz->at(vIt));
    if(ElecP > _maxElecP) _maxElecP = ElecP;
    if(MuonP > _maxMuonP) _maxMuonP = MuonP;
    
    if(elecSCE1x5->at(vIt) > _maxElecSCE1x5) _maxElecSCE1x5 = elecSCE1x5->at(vIt);
    if(elecSCE2x5->at(vIt) > _maxElecSCE2x5) _maxElecSCE2x5 = elecSCE2x5->at(vIt);
    if(elecSCE5x5->at(vIt) > _maxElecSCE5x5) _maxElecSCE5x5 = elecSCE5x5->at(vIt);
    if(elecSCE1x5Over5x5->at(vIt) > _maxElecSCE1x5Over5x5) _maxElecSCE1x5Over5x5 = elecSCE1x5Over5x5->at(vIt);
    if(elecSCE2x5MaxOver5x5->at(vIt) > _maxElecSCE2x5MaxOver5x5) _maxElecSCE2x5MaxOver5x5 = elecSCE2x5MaxOver5x5->at(vIt);
    float Elec1x5EOverP = elecSCE1x5->at(vIt)/ElecP;
    float Elec2x5EOverP = elecSCE2x5->at(vIt)/ElecP;
    float Elec5x5EOverP = elecSCE5x5->at(vIt)/ElecP;
    if(Elec1x5EOverP > _maxElecSCE1x5EOverP) _maxElecSCE1x5EOverP = Elec1x5EOverP;
    if(Elec2x5EOverP > _maxElecSCE2x5EOverP) _maxElecSCE2x5EOverP = Elec2x5EOverP;
    if(Elec5x5EOverP > _maxElecSCE5x5EOverP) _maxElecSCE5x5EOverP = Elec5x5EOverP;
    
    //cout << "Elec ECAL info\n";
    
    if(muonE->at(vIt) > _maxMuonE) _maxMuonE = muonE->at(vIt);

    if(BSx0->at(vIt) > _maxBSx0) _maxBSx0 = BSx0->at(vIt);
    if(BSy0->at(vIt) > _maxBSy0) _maxBSy0 = BSy0->at(vIt);
    if(BSz0->at(vIt) > _maxBSz0) _maxBSz0 = BSz0->at(vIt);   
    
    //cout << "BS info\n";
    /*
    if(muonVertex_x0->at(vIt) > _maxMuonVertex_x0) _maxMuonVertex_x0 = muonVertex_x0->at(vIt); 
    if(muonVertex_y0->at(vIt) > _maxMuonVertex_y0) _maxMuonVertex_y0 = muonVertex_y0->at(vIt); 
    if(muonVertex_z0->at(vIt) > _maxMuonVertex_z0) _maxMuonVertex_z0 = muonVertex_z0->at(vIt);
    */
    //cout << "muon vertices\n";
    for(unsigned int vertexIt = 0; vertexIt < muonDecayLength_x->size(); vertexIt++){
       if(muonDecayLength_x->at(vertexIt) > _maxMuonDecayLength_x) _maxMuonDecayLength_x = muonDecayLength_x->at(vertexIt); 
       if(muonDecayLength_y->at(vertexIt) > _maxMuonDecayLength_y) _maxMuonDecayLength_y = muonDecayLength_y->at(vertexIt); 
       if(muonDecayLength_z->at(vertexIt) > _maxMuonDecayLength_z) _maxMuonDecayLength_z = muonDecayLength_z->at(vertexIt);
    }
    //cout << "muonDecayLengths\n";
    /*
    if(elecVertex_x0->at(vIt) > _maxElecVertex_x0) _maxElecVertex_x0 = elecVertex_x0->at(vIt); 
    if(elecVertex_y0->at(vIt) > _maxElecVertex_y0) _maxElecVertex_y0 = elecVertex_y0->at(vIt); 
    if(elecVertex_z0->at(vIt) > _maxElecVertex_z0) _maxElecVertex_z0 = elecVertex_z0->at(vIt);
    if(elecDecayLength_x->at(vIt) > _maxElecDecayLength_x) _maxElecDecayLength_x = elecDecayLength_x->at(vIt); 
    if(elecDecayLength_y->at(vIt) > _maxElecDecayLength_y) _maxElecDecayLength_y = elecDecayLength_y->at(vIt); 
    if(elecDecayLength_z->at(vIt) > _maxElecDecayLength_z) _maxElecDecayLength_z = elecDecayLength_z->at(vIt);    
    if(elecVx->at(vIt) > _maxElecVx && elecVx->at(vIt) < 1.) _maxElecVx = elecVx->at(vIt); 
    if(elecVy->at(vIt) > _maxElecVy && elecVy->at(vIt) < 1.) _maxElecVy = elecVy->at(vIt); 
    if(elecVz->at(vIt) > _maxElecVz && elecVz->at(vIt) < 1.) _maxElecVz = elecVz->at(vIt); 
  */
    if(passedTopology(vIt)) _topologyCounter++;

    

    

    
    cout << "beam spot x: " << BSx0->at(vIt) << "\n";
    cout << "beam spot y: " << BSy0->at(vIt) << "\n";
    cout << "beam spot z: " << BSz0->at(vIt) << "\n";
   
    cout << "Primary Event Vertex x: " << primaryVertex_x << "\n";
    cout << "Primary Event Vertex y: " << primaryVertex_y << "\n";
    cout << "Primary Event Vertex z: " << primaryVertex_z << "\n";
    cout << "Error on Primary Event Vertex x: " << primaryVertex_xError << "\n";
    cout << "Error on Primary Event Vertex y: " << primaryVertex_yError << "\n";
    cout << "Error on Primary Event Vertex z: " << primaryVertex_zError << "\n";
   
    for(unsigned int elecVtxIt = 0; elecVtxIt < elecVtxEnd; elecVtxIt++){
       cout << "\nElec Candidate: " << elecVtxIt << "\n";
       cout << "elec vertex x: " << elecVertex_x0->at(elecVtxIt) << "\n";
       cout << "elec vertex y: " << elecVertex_y0->at(elecVtxIt) << "\n";
       cout << "elec vertex z: " << elecVertex_z0->at(elecVtxIt) << "\n";
       cout << "elec Pt: " << elecEt->at(elecVtxIt) << "\n";
       cout << "elec Px: " << elecPx->at(elecVtxIt) << "\n";
       cout << "elec Py: " << elecPy->at(elecVtxIt) << "\n";
       cout << "elec Pz: " << elecPz->at(elecVtxIt) << "\n";
       cout << "elec Charge: " << elecCharge->at(elecVtxIt) << "\n";    
    }
    
    for(unsigned int muonVtxIt = 0; muonVtxIt < muonVtxEnd; muonVtxIt++){
       cout << "\nMuon Candidate: "  << muonVtxIt << "\n";
       cout << "muon vertex x: " << muonVertex_x0->at(muonVtxIt) << "\n";
       cout << "muon vertex y: " << muonVertex_y0->at(muonVtxIt) << "\n";
       cout << "muon vertex z: " << muonVertex_z0->at(muonVtxIt) << "\n";  
       cout << "muon Pt: " << muonPt->at(muonVtxIt) << "\n";
       cout << "muon Px: " << muonPx->at(muonVtxIt) << "\n";
       cout << "muon Py: " << muonPy->at(muonVtxIt) << "\n";
       cout << "muon Pz: " << muonPz->at(muonVtxIt) << "\n";        
       cout << "muon Charge: " << muonCharge->at(muonVtxIt) << "\n";
    }
   
    
    cout << "\n";
     
   
    //cout << "elec decay length at event " << jentry << ": " << sqrt(elecDecayLength_x->at(vIt)*elecDecayLength_x->at(vIt) + elecDecayLength_y->at(vIt)*elecDecayLength_y->at(vIt)) << "\n";
    //cout << "muon decay length at event " << jentry << ": " << sqrt(muonDecayLength_x->at(vIt)*muonDecayLength_x->at(vIt) + muonDecayLength_y->at(vIt)*muonDecayLength_y->at(vIt)) << "\n";
    
    
 
  // cout << "PASSED ALL CUTS!\n";   	
  } 
}

void EMuAnalysis::getNMinus1CandCounters(unsigned int theIndex){
  if(passedElecMaxEta(theIndex)){
    _nMinus1ElecPtCounter++;
    if(elecPt->at(theIndex) > _maxElecPtNMinus1) _maxElecPtNMinus1 = elecPt->at(theIndex);
  }
  
  if(passedElecMinPt(theIndex)){
    _nMinus1ElecEtaCounter++;
    if(fabs(elecSCEta->at(theIndex)) < _minElecEtaNMinus1) _minElecEtaNMinus1 = elecSCEta->at(theIndex);
  }
  
  if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedTopology(theIndex) && passedHeepPassedEcalDriven(theIndex)
  && passedHeepPassedCrack(theIndex) && passedHeepPassedDEtaIn(theIndex) && passedHeepPassedDPhiIn(theIndex) && passedHeepPassedHadem(theIndex)
  && passedHeepPassedSigmaIEtaIEta(theIndex) && passedHeepPassed2by5Over5By5(theIndex) && passedElecEt(theIndex)){
    _nMinus1ElecRelIsoCounter++;
    float ElecRelIso = (elecTrkIso->at(theIndex) + elecEcalIso->at(theIndex) + elecHcalIso->at(theIndex))/elecPt->at(theIndex);
    if(ElecRelIso < _minElecRelIsoNMinus1) _minElecRelIsoNMinus1 = ElecRelIso;
  }
													        
 if(passedMuonEta(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedTopology(theIndex)){ 
  _nMinus1MuonPtCounter++;
  if(muonPt->at(theIndex) > _maxMuonPtNMinus1) _maxMuonPtNMinus1 = muonPt->at(theIndex);
  
 }

 if(passedMuonPt(theIndex) && passedMuonEta(theIndex)  && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonMatched(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
 && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex) && passedMuonIsHighPtMuon(theIndex)){ 
  _nMinus1MuonECounter++;
  if(muonE->at(theIndex) > _maxMuonENMinus1) _maxMuonENMinus1 = muonE->at(theIndex);
  
 }

 if(passedMuonPt(theIndex) && passedMuonE(theIndex)  && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonMatched(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
 && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex) && passedMuonIsHighPtMuon(theIndex)){ 
  _nMinus1MuonEtaCounter++;
  if(muonEta->at(theIndex) < _minMuonEtaNMinus1) _minMuonEtaNMinus1 = muonEta->at(theIndex);
  
 }
 
 if(passedMuonPt(theIndex) && passedMuonE(theIndex) && passedMuonEta(theIndex) && passedMuonDzVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonMatched(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
 && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex) && passedMuonIsHighPtMuon(theIndex)){ 
  _nMinus1MuonDxyVtxCounter++;
  if(muonDxyVtx->at(theIndex) > _maxMuonDxyVtxNMinus1) _maxMuonDxyVtxNMinus1 = muonDxyVtx->at(theIndex);
  
 } 
 
 if(passedMuonAcc(theIndex) && passedMuonDxyVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonMatched(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
 && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex) && passedElecAcc(theIndex) &&
 passedHeepId(theIndex) && passedTopology(theIndex)){ 
  _nMinus1MuonDzVtxCounter++;
  if(muonDzVtx->at(theIndex) > _maxMuonDzVtxNMinus1) _maxMuonDzVtxNMinus1 = muonDzVtx->at(theIndex);
 }
 
 
 if(passedMuonAcc(theIndex) && passedMuonDxyVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonMatched(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
 && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex) && passedElecAcc(theIndex) &&
 passedHeepId(theIndex) && passedTopology(theIndex)){ 
  _nMinus1MuonDzBSCounter++;
  if(muonDzBS->at(theIndex) > _maxMuonDzBSNMinus1) _maxMuonDzBSNMinus1 = muonDzBS->at(theIndex);
 } 
 
 if(passedMuonAcc(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedTopology(theIndex) && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) 
 && passedMuonDptOverPt(theIndex) && passedMuonIsHighPtMuon(theIndex)){
   _nMinus1MuonIsoCounter++;
   float IsoCorrection = muonPFIsoDR04SumChargedHadronPt->at(theIndex) + muonPFIsoDR04SumNeutralHadronPt->at(theIndex) + muonPFIsoDR04SumPhotonPt->at(theIndex) - 0.5*muonPFIsoDR04SumPUPt->at(theIndex);
   if(IsoCorrection < 0) IsoCorrection = 0;
   float MuonIso = (IsoCorrection)/muonPt->at(theIndex);
   if(MuonIso < _minMuonIsoNMinus1) _minMuonIsoNMinus1 = MuonIso;
 }
 
 if(passedMuonAcc(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedTopology(theIndex) && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) 
 && passedMuonIsoDR04DB(theIndex) && passedMuonIsHighPtMuon(theIndex)){ 
   _nMinus1MuonDptOverPtCounter++;
   if(muonTrkPt->at(theIndex) > 0 && muonTrkPtError->at(theIndex)/muonTrkPt->at(theIndex) < _minMuonDptOverPtNMinus1) _minMuonDptOverPtNMinus1 = muonTrkPtError->at(theIndex)/muonTrkPt->at(theIndex);
 }
 
 if(passedMuonAcc(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedTopology(theIndex) && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
 && passedMuonIsGlobalMuon(theIndex) && passedMuonNormChiSqrd(theIndex) && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex) && passedMuonTrkLayersWithHits(theIndex) 
 && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex)){  
   _nMinus1MuonIsHighPtMuonCounter++;
   if(muonIsHighPtMuon->at(theIndex) > _maxMuonIsHighPtMuonNMinus1) _maxMuonIsHighPtMuonNMinus1 = muonIsHighPtMuon->at(theIndex);
 }
 
 if(passedHeepPassedCrack(theIndex) && 
 passedHeepPassedDEtaIn(theIndex) && passedHeepPassedDPhiIn(theIndex) && passedHeepPassedHadem(theIndex) && passedHeepPassedSigmaIEtaIEta(theIndex) && 
 passedHeepPassed2by5Over5By5(theIndex) && passedHeepPassedEcalHad1Iso(theIndex) && passedHeepPassedTrkIso(theIndex) && passedHeepPassedEcalDriven(theIndex)){
   _nMinus1ElecEtCounter++;
   if(elecEt->at(theIndex) > _maxElecEtNMinus1) _maxElecEtNMinus1 = elecEt->at(theIndex);
 }

 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1BTagCounter++;
    if(nBTagsCombSecVtxLWP->at(theIndex) < _minNBtagsCombSecVtxNMinus1) _minNBtagsCombSecVtxNMinus1 = nBTagsCombSecVtxLWP->at(theIndex);
 }

 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && /*passedLeptonMetCosDPhi(theIndex) && */ passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    //cout << "PASSED CUTS FOR NMINUS1 COSDPHI\n";
    _nMinus1ElecMuonCosDPhiCounter++;
    if(elecMuonCosDPhi->at(theIndex) < _minElecMuonCosDPhiNMinus1) _minElecMuonCosDPhiNMinus1 = elecMuonCosDPhi->at(theIndex);
    if(leptonMetCosDphi->at(theIndex) < _minLeptonMetCosDPhiNMinus1) _minLeptonMetCosDPhiNMinus1 = leptonMetCosDphi->at(theIndex);
 } 

 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1ElecMuonPZetaCounter++;
    _nMinus1ElecMuonPZetaVisCounter++;
    if(elecMuonPZeta->at(theIndex) > _maxElecMuonPZetaNMinus1) _maxElecMuonPZetaNMinus1 = elecMuonPZeta->at(theIndex);
    if(elecMuonPZetaVis->at(theIndex) > _maxElecMuonPZetaVisNMinus1) _maxElecMuonPZetaVisNMinus1 = elecMuonPZetaVis->at(theIndex);
 } 
  
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1JetSumEtCounter++;
    if(jetSumEt->at(theIndex) < _minJetSumEtNMinus1) _minJetSumEtNMinus1 = jetSumEt->at(theIndex);
 }
 
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) &&  passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1LeptonMetCosDPhiCounter++;
    //if(leptonMetCosDphi->at(theIndex) < _minLeptonMetCosDPhiNMinus1) _minLeptonMetCosDPhiNMinus1 = leptonMetCosDphi->at(theIndex);
    if(leptonDeltaPtCosDphi->at(theIndex) < _minLeptonDeltaPtCosDPhiNMinus1) _minLeptonDeltaPtCosDPhiNMinus1 = leptonDeltaPtCosDphi->at(theIndex);
    
 } 

 if(passedMatched(theIndex) && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex)
 && passedCombinedIpOverIpError(theIndex) && passedElecMuonDCAOverError3D(theIndex) && passedLifetime(theIndex) && passedElecMuonDeltaPt(theIndex) && passed2xPLargerLepton(theIndex)){
    _nMinus1ElecMuonMetMassCounter++;
    if(_source == "TTJets_MADGRAPH" && elecMuonMetMass->at(theIndex) > 700){
      cout << "TTJets Mass = " << elecMuonMetMass->at(theIndex) << " elecPt = " << elecEt->at(theIndex) << " muon Pt = " << muonPt->at(theIndex) << " MEt = " << mEt << " event size = " << elecMatched->size() << "\n";
    }      
    if(elecMuonMetMass->at(theIndex) > _minElecMuonMetMassNMinus1) _minElecMuonMetMassNMinus1 = elecMuonMetMass->at(theIndex);
    if(elecMuonDeltaPtMass->at(theIndex) > _maxElecMuonDeltaPtMassNMinus1) _maxElecMuonDeltaPtMassNMinus1 = elecMuonDeltaPtMass->at(theIndex);
    

    
    
 } 
 
 
if(passedMatched(theIndex) && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex)
 && passedCombinedIpOverIpError(theIndex) && passedElecMuonDCAOverError3D(theIndex) && passedLifetime(theIndex) && passedElecMuonDeltaPt(theIndex) && passed2xPLargerLepton(theIndex)){


    _nMinus1ThreeLeptonCounter++;
    int ThreeLeptons = (eMuMuOneVertex || muEEOneVertex);
    //cout << "Three Leptons = " << ThreeLeptons << "\n";
    
    float ThreeLeptonMass = -10.;
    if(theIndex < muonMuonElecMass->size() && theIndex < elecElecMuonMass->size()){
      muonMuonElecMass->at(theIndex) > elecElecMuonMass->at(theIndex) ? ThreeLeptonMass = muonMuonElecMass->at(theIndex) : ThreeLeptonMass = elecElecMuonMass->at(theIndex);
    }
    if(ThreeLeptonMass > _maxThreeLeptonMassNMinus1) _maxThreeLeptonMassNMinus1 = ThreeLeptonMass; 
    
    if(ThreeLeptons >= _maxThreeLeptonNMinus1 && ThreeLeptons < 2) _maxThreeLeptonNMinus1 = ThreeLeptons;
    
    //cout << "maxThreeLeptonNMinus1 = " << _maxThreeLeptonNMinus1 << "\n";
    
 }  
 if(passedMatched(theIndex) && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex)
 && passedCombinedIpOverIpError(theIndex) && passedElecMuonDCAOverError3D(theIndex) && passedLifetime(theIndex) && passedElecMuonDeltaPt(theIndex) && passed2xPLargerLepton(theIndex) && passedElecMuonMetMass(theIndex)){ 
 
  float muonDptOverPt = fabs(muonTrkPtError->at(theIndex)/muonPt->at(theIndex));
  
  //issues with elec dpt - need to find proper way to express elec dpt - current def is elecpx + muonpx, e.g.
  
  //cout << "Delta Ptx = " << deltaPtx->at(theIndex) << " Delta Pty = " << deltaPty->at(theIndex) << "\n";
  float elecDpt = 5.; //sqrt(deltaPtx->at(theIndex)*deltaPtx->at(theIndex) + deltaPty->at(theIndex)*deltaPty->at(theIndex));
  float elecDptOverPt = fabs(elecDpt/elecEt->at(theIndex));
  
  if(muonDptOverPt > _maxMuonDptOverPtNMinus1) _maxMuonDptOverPtNMinus1 = muonDptOverPt;
  if(elecDptOverPt > _maxElecDptOverPtNMinus1) _maxElecDptOverPtNMinus1 = elecDptOverPt;
  
  if(muonPt->at(theIndex) > _maxMuonPtNMinus1) _maxMuonPtNMinus1 = muonPt->at(theIndex);
  if(elecEt->at(theIndex) > _maxElecEtNMinus1) _maxElecEtNMinus1 = elecEt->at(theIndex);
   
}   
 
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1NJetsCounter++;
    if(nJets->at(theIndex) > _maxNJetsNMinus1) _maxNJetsNMinus1 = nJets->at(theIndex);
 } 
  
 
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex)
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1ElecMuonMassCounter++;
    if(elecMuononMass->at(theIndex) < _minElecMuonMassNMinus1) _minElecMuonMassNMinus1 = elecMuononMass->at(theIndex);
 }  
  
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) 
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1MuonMetMtCounter++;
    if(muonMEtMt->at(theIndex) < _minMuonMetMtNMinus1) _minMuonMetMtNMinus1 = muonMEtMt->at(theIndex);
 }    

 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedNJets(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1ElecMuonChargeCounter++;
    float Charge = elecCharge->at(theIndex) * muonCharge->at(theIndex);
    if(Charge < _minElecMuonChargeNMinus1) _minElecMuonChargeNMinus1 = Charge;
 }     
 
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedCharge(theIndex) && passedNJets(theIndex) && passedJetSumEt(theIndex) && passedCombinedIpOverIpError(theIndex)){
    _nMinus1MEtCounter++;
    if(mEt > _maxMEtNMinus1) _maxMEtNMinus1 = mEt;
 }
  
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedCharge(theIndex) && passedNJets(theIndex) && passedJetSumEt(theIndex) && passedMET(theIndex) &&
 passedElecMuonDCAOverError3D(theIndex)){
    _nMinus1ElecIpOverElecIpErrorCounter++;
    _nMinus1MuonIpOverMuonIpErrorCounter++;
    _nMinus1CombinedIpOverIpErrorCounter++;
    float ElecIpOverElecIpErr = fabs(elecGSFBSIp->at(theIndex)/elecGSFBSIpError->at(theIndex));
    float MuonIpOverMuonIpErr = fabs(muonDxyBS->at(theIndex)/muonDxyBSError->at(theIndex)); 
    float ElecIp = fabs(elecGSFBSIp->at(theIndex));
    float MuonIp = fabs(muonDxyBS->at(theIndex));
    float ElecIpErr = fabs(elecGSFBSIpError->at(theIndex));
    float MuonIpErr = fabs(muonDxyBSError->at(theIndex));
    
    //float CombinedIpOverIpErr = sqrt((ElecIp*ElecIp + MuonIp*MuonIp + 2*ElecIp*MuonIp)/(ElecIpErr*ElecIpErr + MuonIpErr*MuonIpErr));   
    float   CombinedIpOverIpErr = ElecIpOverElecIpErr*MuonIpOverMuonIpErr;
    
    if(ElecIpOverElecIpErr > _maxElecIpOverElecIpErrorNMinus1) _maxElecIpOverElecIpErrorNMinus1 = ElecIpOverElecIpErr;
    if(MuonIpOverMuonIpErr > _maxMuonIpOverMuonIpErrorNMinus1) _maxMuonIpOverMuonIpErrorNMinus1 = MuonIpOverMuonIpErr;
    if(CombinedIpOverIpErr > _maxCombinedIpOverIpErrorNMinus1) {
      //cout << "Min combined Ip err/Ip = " << _minCombinedIpErrorOverIpNMinus1 << "\n";
      //cout << "CombinedIpErrOverIp = " << CombinedIpErrOverIp << "\n";
      _maxCombinedIpOverIpErrorNMinus1 = CombinedIpOverIpErr;
      //cout << "Min combined Ip err/Ip = " << _minCombinedIpErrorOverIpNMinus1 << "\n\n";
    }
   
 }
 
 
 if(passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedElecAcc(theIndex) && passedHeepId(theIndex) && passedElecMuonMass(theIndex) 
 && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedElecMuonMetMass(theIndex) && passedElecMetMt(theIndex) && passedMuonMetMt(theIndex)
 && passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedCharge(theIndex) && passedNJets(theIndex) && passedJetSumEt(theIndex) && passedMET(theIndex) &&
 passedCombinedIpOverIpError(theIndex)){
    
    _nMinus1ElecMuonDCAOverError3DCounter++;
    
    float ElecMuonDis3DOverError = elecMuonTwoTrkMinDist->at(theIndex)/elecMuonDCA3DError->at(theIndex);
    
    if(ElecMuonDis3DOverError > _maxElecMuonDCAOverError3DNMinus1) _maxElecMuonDCAOverError3DNMinus1 = ElecMuonDis3DOverError;

 
 }
 
 
 
} 

void EMuAnalysis::getReport(){

  pair<double, pair<double,double> > efficiency = GetEfficiency(_nElecEOverP, _puWeightedEvents);
  pair<double, pair<double,double> > startEfficiency = GetEfficiency(_puWeightedEvents, _puWeightedEvents);
  
  pair<double, double> survivingEventsAtLumi = GetSurvivingEvents(nSkimmedAtLumi, efficiency);
  pair<double, double> survivingEventsAtStart = GetSurvivingEvents(nSkimmedAtLumi, startEfficiency);
  
  double _FirstDenom = _puWeightedEvents;
  
  if (_source == "ZPrimeSSM_M-500" || _source == "ZPrimeSSM_M-750" || _source == "ZPrimeSSM_M-1000" || _source == "ZPrimeSSM_M-1250" || _source == "ZPrimeSSM_M-1500" || _source == "ZPrimeSSM_M-1750" || _source == "ZPrimeSSM_M-2000" || 
  _source == "ZPrimeSSM_M-2250" || _source == "ZPrimeSSM_M-2500"){
  	_FirstDenom = (double)_nMatched;	
  }
  
  if (_source == "Data"){
  	survivingEventsAtLumi.first = _nElecEOverP;
	survivingEventsAtLumi.second = sqrt(_nElecEOverP);
	
	survivingEventsAtStart.first = _puWeightedEvents;
	survivingEventsAtStart.second = _puWeightedEvents;
  }
  
  std::ofstream outLogFile;
  outLogFile.open(_outLogFile.c_str());    
  outLogFile.setf(ios::fixed,ios::floatfield);
  
  outLogFile << "\t Efficiencies \n";
  outLogFile  
  

  << "\n _nEvents: \t \t \t"		 << _nEvents
  << "\n _puWeightedEvents: \t \t" 	 << _puWeightedEvents
  << "\n Counter Events: \t \t"		 << counterEvents
  << "\n _nElecMuonMetMass: \t \t"	 << _nElecMuonMetMass
  

  << "\n _nMatched: \t \t \t" 		 << _nMatched
  << "\n _nMuonPt: \t \t \t"		 << _nMuonPt
  << "\n _nMuonEta: \t \t \t"		 << _nMuonEta
  << "\n _nMuonAcc: \t \t \t"		 << _nMuonAcc
  << "\n _nMuonIsGlobal: \t \t"		 << _nMuonIsGlobalMuon
  << "\n _nElecEt: \t \t \t"		 << _nElecEt
  << "\n _nElecEta: \t \t \t"		 << _nElecEta
  << "\n _nElecAcc: \t \t \t"	 	 << _nElecAcc
  << "\n _nElecIsEcalDriven: \t \t"	 << _nHeepPassedEcalDriven
  << "\n _nEMuDeltaR: \t \t \t"		 << _nElecMuonDeltaR
  << "\n _nMuonDxy: \t \t \t"		 << _nMuonDxyVtx
  << "\n _nMuonDz: \t \t \t"		 << _nMuonDzVtx
  << "\n _nMuonValidHits: \t \t"	 << _nMuonValidHits
  << "\n _nMuonPixelHits: \t \t"	 << _nMuonPixelHits
  << "\n _nMuonMatchedStations: \t"	 << _nMuonMatchedStations
  << "\n _nMuonTrackerLayersWithHits: \t"<< _nMuonTrkLayersWithHits
  << "\n _nMuonDptOverPt: \t \t"	 << _nMuonDptOverPt
  << "\n _nMuonIsoDR04DB: \t \t" 	 << _nMuonIsoDR04DB
  

  << "\n _nMuonIsHighPtMuon: \t \t" 	 << _nMuonIsHighPtMuon
  << "\n _nMuonId: \t \t \t"		 << _nMuonId
  
  << "\n"
  << "\n _nHeepPassedDEtaIn: \t \t"	 << _nHeepPassedDEtaIn
  << "\n _nHeepPassedDPhiIn: \t \t"	 << _nHeepPassedDPhiIn
  << "\n _nHeepPassedHadem: \t \t"	 << _nHeepPassedHadem
  << "\n _nHeepPassedSigmaIEtaIEta: \t"	 << _nHeepPassedSigmaIEtaIEta
  << "\n _nHeepPassed2by5Over5By5: \t"   << _nHeepPassed2by5Over5By5
  << "\n _nHeepPassedEcalHad1Iso: \t"	 << _nHeepPassedEcalHad1Iso
  << "\n _nHeepPassedTrkIso: \t \t"	 << _nHeepPassedTrkIso
  << "\n _nElecTrkIso: \t \t \t:"	 << _nElecTrkIso
  << "\n _nHeepPassedMissingHits: \t"	 << _nHeepPassedMissingHits
  << "\n _nHeepPassedElecDxy: \t \t"	 << _nElecDxy
  << "\n _nHeepId: \t \t \t"		 << _nHeepId
  << "\n"
  << "\n _nCharge: \t \t \t"		 << _nCharge
  << "\n _nElecMuonCosDPhi: \t \t"	 << _nElecMuonCosDPhi
  << "\n _nMEt: \t \t \t"		 << _nMEt
  << "\n _nNbTagsCombSecVtx: \t \t" 	 << _nNbTagsCombSecVtx
  << "\n _nElecMuonPZeta: \t \t"	 << _nElecMuonPZeta
  << "\n _nLeptonMetCosDPhi: \t \t"	 << _nLeptonMetCosDPhi
  << "\n _nLeptonDeltaPtCosDPhi: \t"	 << _nLeptonDeltaPtCosDPhi
  << "\n _nJetSumEt: \t \t \t "	 	 << _nJetSumEt
  << "\n _nNjets: \t \t \t"		 << _nNJets
  << "\n _nMuonMetMt: \t \t \t"	 	 << _nMuonMetMt	
  << "\n _nElecMuonMass: \t \t"	 	 << _nElecMuonMass
  << "\n _nElecMetMt: \t \t \t"	 	 << _nElecMetMt
  << "\n _nElecMuonMetMass: \t \t"	 << _nElecMuonMetMass
  << "\n _nElecMuonDeltaPtMass: \t"	 << _nElecMuonDeltaPtMass
  << "\n _n2xPLargerLepton: \t \t"	 << _n2xPLargerLepton
  << "\n _nElecMuonDeltaPt: \t \t"	 << _nElecMuonDeltaPt
  << "\n _nElecIpOverElecIpError:\t"	 << _nElecIpOverElecIpError
  << "\n _nMuonIpOverMuonIpError:\t"     << _nMuonIpOverMuonIpError
  << "\n _nCombinedIpOverIpError:\t"	 << _nCombinedIpOverIpError  
  << "\n _nElecMuonDCAOverError3D:\t"	 << _nElecMuonDCAOverError3D
  << "\n _nElecMuonDCAOverError2D:\t"	 << _nElecMuonDCAOverError2D
  << "\n _nLifetime: \t \t \t \t"	 << _nLifetime
  << "\n _nTopology: \t \t \t \t"	 << _nTopology
  
  << "\n"

  

  
  
  << "\n"
  << "\n Matching: \t \t \t \t"	   	<< PrintEfficiency(_nMatched, (double)_nEvents)
  << "\n"
   
  << "\n Electron_ip/Error: \t \t \t"   << PrintEfficiency(_nElecIpOverElecIpError, (double)_nEvents)
  << "\n Muon_ip/Error: \t \t \t"	<< PrintEfficiency(_nMuonIpOverMuonIpError, (double)_nEvents)
  << "\n Elec Muon Met Mass: \t \t \t"	<< PrintEfficiency(_nElecMuonMetMass, _FirstDenom)
  << "\n"
  << "\n Muon Pt: \t \t \t \t" 		<< PrintEfficiency(_nMuonPt, _FirstDenom)
  << "\n Muon Eta: \t \t \t \t" 	<< PrintEfficiency(_nMuonEta, _nMuonPt)
  << "\n Muon Is Global Muon: \t \t \t" << PrintEfficiency(_nMuonIsGlobalMuon, _nMuonEta)
  << "\n"
  << "\n ElecEt: \t \t \t \t"		<< PrintEfficiency(_nElecEt, _nMuonIsGlobalMuon)  

  << "\n Elec Eta: \t \t \t \t"		<< PrintEfficiency(_nElecEta, _nElecEt)
  << "\n HeepPassedEcalDriven: \t \t"   << PrintEfficiency(_nHeepPassedEcalDriven, _nElecEta)
  << "\n Elec Muon Delta R: \t \t \t" 	<< PrintEfficiency(_nElecMuonDeltaR, _nHeepPassedEcalDriven)  
  << "\n Muon Dxy Vtx: \t \t \t \t" 	<< PrintEfficiency(_nMuonDxyVtx, _nElecMuonDeltaR)
  << "\n Muon Dz Vtx: \t \t \t \t" 	<< PrintEfficiency(_nMuonDzVtx, _nMuonDxyVtx)
  << "\n Muon Valid Hits: \t \t \t" 	<< PrintEfficiency(_nMuonValidHits, _nMuonDzVtx)
  << "\n Muon Pixel Hits: \t \t \t" 	<< PrintEfficiency(_nMuonPixelHits, _nMuonValidHits)
  << "\n Muon Matched Stations: \t \t " << PrintEfficiency(_nMuonMatchedStations, _nMuonPixelHits)
  << "\n Muon Track Layers with Hits: \t \t" << PrintEfficiency(_nMuonTrkLayersWithHits, _nMuonMatchedStations)
  << "\n Muon Dpt Over Pt: \t \t \t"    << PrintEfficiency(_nMuonDptOverPt, _nMuonTrkLayersWithHits)
  << "\n Muon Iso DR04 DB: \t \t \t" 	<< PrintEfficiency(_nMuonIsoDR04DB, _nMuonDptOverPt)
  
  << "\n"
  << "\n HeepPassedDEtaIn: \t \t \t"	<< PrintEfficiency(_nHeepPassedDEtaIn, _nMuonIsoDR04DB)
  << "\n HeepPassedDPhiIn: \t \t \t"	<< PrintEfficiency(_nHeepPassedDPhiIn, _nHeepPassedDEtaIn)
  << "\n HeepPassedHadem: \t \t \t"	<< PrintEfficiency(_nHeepPassedHadem, _nHeepPassedDPhiIn)
  << "\n HeepPassedSigmaIEtaIEta: \t \t"	<< PrintEfficiency(_nHeepPassedSigmaIEtaIEta, _nHeepPassedHadem)
  << "\n HeepPassed2by5Over5By5: \t \t"	<< PrintEfficiency(_nHeepPassed2by5Over5By5, _nHeepPassedSigmaIEtaIEta)
  << "\n HeepPassedEcalHad1Iso: \t \t"  << PrintEfficiency(_nHeepPassedEcalHad1Iso, _nHeepPassed2by5Over5By5)
  << "\n Elec Trk Iso: \t \t \t"	<< PrintEfficiency(_nElecTrkIso, _nHeepPassedEcalHad1Iso)
  << "\n HeepPassedMissingHits: \t \t"	<< PrintEfficiency(_nHeepPassedMissingHits, _nElecTrkIso)
  << "\n HeepPassedElecDxy: \t \t \t"	<< PrintEfficiency(_nElecDxy, _nHeepPassedMissingHits)
  << "\n"
  << "\n Charge Requirement: \t \t \t" 	<< PrintEfficiency(_nCharge, _nElecDxy)  
  << "\n Elec Muon Cos DPhi: \t \t \t" 	<< PrintEfficiency(_nElecMuonCosDPhi, _nCharge)
  << "\n MEt: \t \t \t \t \t" 		<< PrintEfficiency(_nMEt, _nElecMuonCosDPhi)
  << "\n nBtagsCombSecVtx: \t \t \t" 	<< PrintEfficiency(_nNbTagsCombSecVtx, _nMEt)
  << "\n Elec Muon PZeta: \t \t \t" 	<< PrintEfficiency(_nElecMuonPZeta, _nNbTagsCombSecVtx)
  << "\n Lepton-Met Cos DPhi: \t \t \t" << PrintEfficiency(_nLeptonMetCosDPhi, _nElecMuonPZeta)
  << "\n Lepton-DeltaPt Cos DPhi: \t"	<< PrintEfficiency(_nLeptonDeltaPtCosDPhi, _nLeptonMetCosDPhi)
  << "\n Three Lepton One Vertex: \t"	<< PrintEfficiency(_nThreeLepton, _nLeptonDeltaPtCosDPhi)
  << "\n"
  
  
   

  
  
  << "\n"
  
  << "\n Jet Sum Et: \t \t \t \t" 		<< PrintEfficiency(_nJetSumEt, _nThreeLepton)
  << "\n nJets: \t \t \t \t" 		<< PrintEfficiency(_nNJets, _nJetSumEt)
  << "\n Muon Met Mt: \t \t \t \t" 	<< PrintEfficiency(_nMuonMetMt, _nNJets)
  << "\n Elec Muon Mass: \t \t \t" 	<< PrintEfficiency(_nElecMuonMass, _nMuonMetMt)
  
  << "\n Elec Met Mt: \t \t \t \t" 	<< PrintEfficiency(_nElecMetMt, _nElecMuonMass)  
  << "\n Elec Muon Met Mass: \t \t \t" 	<< PrintEfficiency(_nElecMuonMetMass, _nElecMetMt)
  << "\n Elec Muon DeltaPt Mass: \t \t"	<< PrintEfficiency(_nElecMuonDeltaPtMass, _nElecMuonMetMass)
  << "\n 2x p(larger p lepton): \t \t"  << PrintEfficiency(_n2xPLargerLepton, _nElecMuonDeltaPtMass)
  << "\n Elec Muon Delta P: \t \t \t"	<< PrintEfficiency(_nElecMuonDeltaPt, _n2xPLargerLepton)
  << "\n Combined Ip/Error: \t \t \t"   << PrintEfficiency(_nCombinedIpOverIpError, _nElecMuonDeltaPt)
  << "\n Elec Ip/Error: \t \t \t"	<< PrintEfficiency(_nElecIpOverElecIpError, _nCombinedIpOverIpError)
  << "\n Muon Ip/Error: \t \t \t"	<< PrintEfficiency(_nMuonIpOverMuonIpError, _nElecIpOverElecIpError)
  << "\n Elec Muon DCA/Error 3D: \t \t"	<< PrintEfficiency(_nElecMuonDCAOverError3D, _nMuonIpOverMuonIpError)
  << "\n Elec Muon DCA/Error 2D: \t \t" << PrintEfficiency(_nElecMuonDCAOverError2D, _nElecMuonDCAOverError3D)
  << "\n Tau Lifetime: \t \t \t"	<< PrintEfficiency(_nLifetime, _nElecMuonDCAOverError2D)
  << "\n Elec E over P: \t \t \t"	<< PrintEfficiency(_nElecEOverP,  _nLifetime)
  


  << "\n" 
  << "\n Muon Acc: \t \t \t \t" 	<< PrintEfficiency(_nMuonAcc, _FirstDenom)
  << "\n Elec Acc: \t \t \t \t"		<< PrintEfficiency(_nElecAcc, _nMuonAcc)
  << "\n Muon ID: \t \t \t \t" 		<< PrintEfficiency(_nMuonId, _nElecAcc)  
  << "\n HeepID: \t \t \t \t" 		<< PrintEfficiency(_nHeepId, _nMuonId)	

  << "\n All Topology: \t \t \t \t" 	<< PrintEfficiency(_nTopology, _nHeepId)
  << "\n All Cuts: \t \t \t \t" 		<< PrintEfficiency(_nElecEOverP, _FirstDenom)
  << "\n"
  << "\n Total Events: \t \t \t \t"	<< _nEvents 
  << "\n Events @Start: \t \t \t \t"	<< setprecision(2) << survivingEventsAtStart.first << "\t +- \t" << setprecision(2) << survivingEventsAtStart.second
  << "\n Events @" << (int)_theLumi << "pb-1:\t \t \t" <<  setprecision(2) << survivingEventsAtLumi.first 
  << "\t +- \t" <<  setprecision(2) << survivingEventsAtLumi.second
  << "\n ReWeighted Events: \t \t \t" << _puWeightedEvents 
  << "\n Events Passing: \t \t \t" 	<< _nElecEOverP
  << "\n";

  
  outLogFile.close();
}





void EMuAnalysis::fillHistos(){

  _tauEnergyHisto->Fill(_maxTauEnergy, _eventPUWeight);
  _zPrimeEnergyHisto->Fill(_maxZPrimeEnergy, _eventPUWeight);
  _muonPOverTauEHisto->Fill(_maxMuonPOverTauE, _eventPUWeight);
  _muonPOverTauEVsMuonIpOverErrorHisto->Fill(_maxMuonIpOverMuonIpError, _maxMuonPOverTauE, _eventPUWeight);
  _muonPOverTauEVsCosDPhiHisto->Fill(_minElecMuonCosDPhi, _maxMuonPOverTauE, _eventPUWeight);
  
  float DecayLength = sqrt(_maxMuonDecayLength_x*_maxMuonDecayLength_x + _maxMuonDecayLength_y*_maxMuonDecayLength_y + _maxMuonDecayLength_z*_maxMuonDecayLength_z);
  _muonPOverTauEVsTauDecayLengthHisto->Fill(DecayLength, _maxMuonPOverTauE, _eventPUWeight);
  _tauEnergyVsTauDecayLengthHisto->Fill(DecayLength, _maxTauEnergy, _eventPUWeight);
  _vertexDiffHisto->Fill(_maxVertexDiff, _eventPUWeight);
  _vertexDiffxyHisto->Fill(_maxVertexDiffxy, _eventPUWeight);
  
  _impactParameterVsVertexDiffHisto->Fill(_maxVertexDiff, _maxCombinedIpOverIpError, _eventPUWeight);
  
  _beamSpotDiffHisto->Fill(_maxBeamSpotDiff, _eventPUWeight);
  _beamSpotDiffxyHisto->Fill(_maxBeamSpotDiffxy, _eventPUWeight);
  
  _beamSpotPrimaryzDiffHisto->Fill(_maxBeamSpotPrimaryzDiff, _eventPUWeight);
  
  _tauEnergyVsElecMuonMetMassHisto->Fill(_minElecMuonMetMassNMinus1, _maxTauEnergy, _eventPUWeight);
  _zPrimeEnergyVsElecMuonMetMassHisto->Fill(_minElecMuonMetMassNMinus1, _maxZPrimeEnergy, _eventPUWeight);
  _tauEnergyVsMuonEtaHisto->Fill(_minMuonEta, _maxTauEnergy, _eventPUWeight);
  _zPrimeEnergyVsMuonEtaHisto->Fill(_minMuonEta, _maxZPrimeEnergy, _eventPUWeight);
  _tauEnergyVsElecEtaHisto->Fill(_minElecEta, _maxTauEnergy, _eventPUWeight);
  _zPrimeEnergyVsElecEtaHisto->Fill(_minElecEta, _maxZPrimeEnergy, _eventPUWeight);  
  _zPrimeEnergyVs2xMuonPtHisto->Fill(2*_maxMuonPt, _maxZPrimeEnergy, _eventPUWeight);
  _zPrimeEnergyVs2xElecEtHisto->Fill(2*_maxElecEt, _maxZPrimeEnergy, _eventPUWeight);
  float maxPt = max(_maxMuonPt, _maxElecEt);
  float maxP = max(_maxMuonP, _maxElecP);
  _zPrimeEnergyVs2xMaxLeptonPtHisto->Fill(2*maxPt, _maxZPrimeEnergy, _eventPUWeight);
  _zPrimeEnergyVs2xMaxLeptonPHisto->Fill(2*maxP, _maxZPrimeEnergy, _eventPUWeight);

  _nVtxHisto->Fill(nVtx,_eventPUWeight);				      
  _elecPtHisto->Fill(_maxElecPt,_eventPUWeight); 			      
  _elecEtaHisto->Fill(_minElecEta,_eventPUWeight);			      
  _elecDeltaPhiInBarrelHisto->Fill(_minElecDeltaPhiInBarrel,_eventPUWeight);		      
  _elecDeltaPhiInEndCapHisto->Fill(_minElecDeltaPhiInEndCap,_eventPUWeight);		      
  _elecDeltaEtaInBarrelHisto->Fill(_minElecDeltaEtaInBarrel,_eventPUWeight);		      
  _elecDeltaEtaInEndCapHisto->Fill(_minElecDeltaEtaInEndCap,_eventPUWeight);		      
  _elecSigmaIEtaIEtaBarrelHisto->Fill(_minElecSigmaIEtaBarrel,_eventPUWeight);	      
  _elecSigmaIEtaIEtaEndCapHisto->Fill(_minElecSigmaIEtaEndCap,_eventPUWeight);	      
  _elecHadFracBarrelHisto->Fill(_minElecHadFracBarrel,_eventPUWeight);		      
  _elecHadFracEndCapHisto->Fill(_minElecHadFracEndCap,_eventPUWeight);		      
  _elecEOverPBarrelHisto->Fill(_minElecEOverPBarrel,_eventPUWeight);		      
  _elecEOverPEndCapHisto->Fill(_minElecEOverPEndCap,_eventPUWeight);		      		      
  _elecMissingHitsHisto->Fill(_maxElecMissingHits,_eventPUWeight);		      
  _elecPFDR03ChargedHadIsoHisto->Fill(_minPFIsoChargedHadDr03,_eventPUWeight);	      
  _elecPFDR04ChargedHadIsoHisto->Fill(_minPFIsoChargedHadDr04,_eventPUWeight);	      
  _elecPFDR03NeutralHadIsoHisto->Fill(_minPFIsoNeutralHadDr03,_eventPUWeight);	      
  _elecPFDR04NeutralHadIsoHisto->Fill(_minPFIsoNeutralHadDr04,_eventPUWeight);	      
  _elecPFDR03PhotonIsoHisto->Fill(_minPFIsoPhotonDr03,_eventPUWeight);		      
  _elecPFDR04PhotonIsoHisto->Fill(_minPFIsoPhotonDr04,_eventPUWeight);		      
  _elecPFDR03RelIsoHisto->Fill(_minPFRelIsoDr03,_eventPUWeight);		      
  _elecPFDR04RelIsoHisto->Fill(_minPFRelIsoDr04,_eventPUWeight);
  
  _muonPtHisto->Fill(_maxMuonPt,_eventPUWeight);
  _muonEtaHisto->Fill(_minMuonEta,_eventPUWeight);
  _muonPtVsMuonEtaHisto->Fill(_minMuonEta, _maxMuonPt, _eventPUWeight);
  _muonEVsMuonEtaHisto->Fill(_minMuonEta, _maxMuonE, _eventPUWeight);

  _muonDxyVtxHisto->Fill(_maxMuonDxyVtx,_eventPUWeight);
  _muonDzVtxHisto->Fill(_maxMuonDzVtx,_eventPUWeight);
  _muonIsGlobalMuonHisto->Fill(_maxIsGlobalMuon,_eventPUWeight);
  _muonNormChiSqrdHisto->Fill(_maxMuonNormChiSqrd,_eventPUWeight);
  _muonValidHitsHisto->Fill(_maxMuonValidHits,_eventPUWeight);
  _muonMatchedStationsHisto->Fill(_maxMuonMatchedStations,_eventPUWeight);
  _muonPixelHitsHisto->Fill(_maxMuonPixelHits,_eventPUWeight);
  _muonTrkLayersWithHitsHisto->Fill(_maxMuonTrkLayersWithHits,_eventPUWeight);
  _muonIsoDR04DBHisto->Fill(_maxMuonIsoDR04DB, _eventPUWeight);
  _muonDptOverPtHisto->Fill(_minMuonDptOverPt, _eventPUWeight);
  _muonIsHighPtMuonHisto->Fill(_maxMuonIsHighPtMuon, _eventPUWeight);



  _heepPassedEcalDrivenHisto->Fill(_maxHeepPassedEcalDriven,_eventPUWeight);
  _heepPassedCrackHisto->Fill(_maxHeepPassedCrack,_eventPUWeight); 
  _heepPassedDEtaInHisto->Fill(_maxHeepPassedDEtaIn,_eventPUWeight); 
  _heepPassedDPhiInHisto->Fill(_maxHeepPassedDPhiIn,_eventPUWeight); 
  _heepPassedHademHisto->Fill(_maxHeepPassedHadem,_eventPUWeight); 
  _heepPassedSigmaIEtaIEtaHisto->Fill(_maxHeepPassedSigmaIEtaIEta,_eventPUWeight); 
  _heepPassedEcalHad1IsoHisto->Fill(_maxHeepPassedEcalHad1Iso,_eventPUWeight); 
  _heepPassedTrkIsoHisto->Fill(_maxHeepPassedTrkIso,_eventPUWeight); 
  _elecEtHisto->Fill(_maxElecEt,_eventPUWeight);      
  
  _elecSCE1x5Histo->Fill(_maxElecSCE1x5, _eventPUWeight);
  _elecSCE2x5Histo->Fill(_maxElecSCE2x5, _eventPUWeight);
  _elecSCE5x5Histo->Fill(_maxElecSCE5x5, _eventPUWeight);
  _elecSCE1x5Over5x5Histo->Fill(_maxElecSCE1x5Over5x5, _eventPUWeight);
  _elecSCE2x5MaxOver5x5Histo->Fill(_maxElecSCE2x5MaxOver5x5, _eventPUWeight);
  _elecSCE1x5EOverPHisto->Fill(_maxElecSCE1x5EOverP, _eventPUWeight);
  _elecSCE2x5EOverPHisto->Fill(_maxElecSCE2x5EOverP, _eventPUWeight);
  _elecSCE5x5EOverPHisto->Fill(_maxElecSCE5x5EOverP, _eventPUWeight);
  
  _muonEHisto->Fill(_maxMuonE, _eventPUWeight);
  
  _elecSCE1x5VsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE1x5, _eventPUWeight);
  _elecSCE2x5VsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE2x5, _eventPUWeight);
  _elecSCE5x5VsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE5x5, _eventPUWeight);
  _elecSCE1x5Over5x5VsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE1x5Over5x5, _eventPUWeight);
  _elecSCE2x5MaxOver5x5VsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE2x5MaxOver5x5, _eventPUWeight);
  _elecSCE1x5EOverPVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE1x5EOverP, _eventPUWeight);
  _elecSCE2x5EOverPVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE2x5EOverP, _eventPUWeight);
  _elecSCE5x5EOverPVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecSCE5x5EOverP, _eventPUWeight);
  
  _elecSCE1x5EOverPVsEOverP_hiIPHisto->Fill(_maxElecEOverP_hiIP, _maxElecSCE1x5EOverP, _eventPUWeight);
  _elecSCE2x5EOverPVsEOverP_hiIPHisto->Fill(_maxElecEOverP_hiIP, _maxElecSCE2x5EOverP, _eventPUWeight);
  _elecSCE5x5EOverPVsEOverP_hiIPHisto->Fill(_maxElecEOverP_hiIP, _maxElecSCE5x5EOverP, _eventPUWeight);
  
  	     
  
  _MetMinusDeltaPtHisto->Fill(_maxMetMinusDeltaPt, _eventPUWeight);  	
  _MetxVsDeltaPtxHisto->Fill(_maxDeltaPtx, _maxMetPx, _eventPUWeight);				            
  _MetyVsDeltaPtyHisto->Fill(_maxDeltaPty, _maxMetPy, _eventPUWeight);				            
		      	      
  _elecIpOverElecIpErrorHisto->Fill(_maxElecIpOverElecIpError, _eventPUWeight);
  _muonIpOverMuonIpErrorHisto->Fill(_maxMuonIpOverMuonIpError, _eventPUWeight);
  _combinedIpOverIpErrorHisto->Fill(_maxCombinedIpOverIpError, _eventPUWeight);	
  _elecIpErrorHisto->Fill(_maxElecIpError, _eventPUWeight);
  _muonIpErrorHisto->Fill(_maxMuonIpError, _eventPUWeight);
  _elecIpOverErrVsElecIpHisto->Fill(_maxElecIp, _maxElecIpOverElecIpError, _eventPUWeight);
  _elecIpHisto->Fill(_maxElecIp, _eventPUWeight);
  
  _elecEOverP_hiIPHisto->Fill(_maxElecEOverP_hiIP, _eventPUWeight);
  _elecEta_hiIPHisto->Fill(_maxElecEta_hiIP, _eventPUWeight);
  _elecEOverPVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecEOverP_hiIP, _eventPUWeight);
  _elecEtaVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecEta_hiIP, _eventPUWeight);
  _elecEtVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecEt_hiIP, _eventPUWeight);
  _elecEOverPVsMuonIpOverError_hiIPHisto->Fill(_maxMuonIpOverMuonIpError, _maxElecEOverP_hiIP, _eventPUWeight);
  _elecTrkChiSqrdVsElecIpOverError_hiIPHisto->Fill(_maxElecIpOverElecIpError, _maxElecTrkChiSqrd_hiIP, _eventPUWeight);
  _elecTrkChiSqrdVsElecEOverP_hiIPHisto->Fill(_maxElecEOverP_hiIP, _maxElecTrkChiSqrd_hiIP, _eventPUWeight);
  _elecEOverPVsElecMissingHits_hiIPHisto->Fill(_maxElecMissingHits, _maxElecEOverP_hiIP, _eventPUWeight);
  _elecEOverPVsElecIpError_hiIPHisto->Fill(_maxElecIpError, _maxElecEOverP_hiIP, _eventPUWeight);
  _elecPVsElecIpError_hiIPHisto->Fill(_maxElecIpError, _maxElecP, _eventPUWeight);
  _elecPVsElecTrkChiSqrd_hiIPHisto->Fill(_maxElecTrkChiSqrd_hiIP, _maxElecP, _eventPUWeight);
  _elecEOverPVsElecP_hiIPHisto->Fill(_maxElecP, _maxElecEOverP_hiIP, _eventPUWeight);
  
  
  _elecMuonDistanceOfClosestApproachHisto->Fill(_maxElecMuonDistanceOfClosestApproach, _eventPUWeight);
  _elecMuonDCAVsCrossingPointHisto->Fill(_maxElecMuonDCACrossingPoint, _maxElecMuonDistanceOfClosestApproach, _eventPUWeight);
  _elecMuonTwoTrkMinDistHisto->Fill(_maxElecMuonTwoTrkMinDist, _eventPUWeight);
  _elecMuonTwoTrkMinDistVsCrossingPointHisto->Fill(_maxElecMuonTwoTrkCrossingPoint, _maxElecMuonTwoTrkMinDist, _eventPUWeight);
  _eMuDis3DHisto->Fill(_maxElecMuonDis3D, _eventPUWeight);
  _eMuDis2DHisto->Fill(_maxElecMuonDis2D, _eventPUWeight);
  _elecMuonDCA3DErrorHisto->Fill(_maxElecMuonDCA3DError, _eventPUWeight);
  _elecMuonDCA2DErrorHisto->Fill(_maxElecMuonDCA2DError, _eventPUWeight);
  _eMuDis3DOverErrorHisto->Fill(_maxElecMuonDis3DOverError, _eventPUWeight);
  _eMuDis2DOverErrorHisto->Fill(_maxElecMuonDis2DOverError, _eventPUWeight);
  
  
  _elecTrkChiSqrdHisto->Fill(_minElecTrkChiSqrd, _eventPUWeight);
  _muonTrkChiSqrdHisto->Fill(_minMuonTrkChiSqrd, _eventPUWeight);
  
  _IpVsElecTrkChiSqrdHisto->Fill(_minElecTrkChiSqrd, _maxCombinedIpOverIpError, _eventPUWeight);
  _IpVsMuonTrkChiSqrdHisto->Fill(_minMuonTrkChiSqrd, _maxCombinedIpOverIpError, _eventPUWeight);
  _elecIpVsElecTrkChiSqrdHisto->Fill(_minElecTrkChiSqrd, _maxElecIpOverElecIpError, _eventPUWeight);
  _muonIpVsMuonTrkChiSqrdHisto->Fill(_minMuonTrkChiSqrd, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _elecIpErrVsElecTrkChiSqrdHisto->Fill(_minElecTrkChiSqrd, _maxElecIpError, _eventPUWeight);
  _muonIpErrVsMuonTrkChiSqrdHisto->Fill(_minMuonTrkChiSqrd, _maxMuonIpError, _eventPUWeight);
  _elecIpOverErrVsElecIpErrHisto->Fill(_maxElecIpError, _maxElecIpOverElecIpError, _eventPUWeight);
  _muonIpOverErrVsMuonIpErrHisto->Fill(_maxMuonIpError, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _muonIpOverErrVsElecIpOverErrHisto->Fill(_maxElecIpOverElecIpError, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _muonIpOverErrVsMuonPixelHitsHisto->Fill(_maxMuonPixelHits, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _muonEtaVsMuonPixelHitsHisto->Fill(_maxMuonPixelHits, _minMuonEta, _eventPUWeight);
  _muonIpOverErrVsMuonEtaHisto->Fill(_minMuonEta, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _muonPtVsMuonIpOverErrHisto->Fill(_maxMuonIpOverMuonIpError, _maxMuonPt, _eventPUWeight);
  _elecEtVsElecIpOverErrHisto->Fill(_maxElecIpOverElecIpError, _maxElecEt, _eventPUWeight);
  _MEtVsElecIpOverErrHisto->Fill(_maxElecIpOverElecIpError, _maxMEt, _eventPUWeight);
  _MEtVsMuonIpOverErrHisto->Fill(_maxMuonIpOverMuonIpError, _maxMEt, _eventPUWeight);
  _MEtVsCombinedIpOverErrHisto->Fill(_maxCombinedIpOverIpError, _maxMEt, _eventPUWeight);
  _muonIpOverErrVsMuonIpHisto->Fill(_maxMuonIp, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _muonIpHisto->Fill(_maxMuonIp, _eventPUWeight);
  
  _DCAOverErrorVsElecIpOverErrorHisto->Fill(_maxElecIpOverElecIpError, _maxElecMuonDis3DOverError, _eventPUWeight);
  _DCAOverErrorVsMuonIpOverErrorHisto->Fill(_maxMuonIpOverMuonIpError, _maxElecMuonDis3DOverError, _eventPUWeight);
  _DCAOverErrorVsCombinedIpOverErrorHisto->Fill(_maxCombinedIpOverIpError, _maxElecMuonDis3DOverError, _eventPUWeight);
  
  
  _elecIpOverErrVsElecMissingHitsHisto->Fill(_maxElecMissingHits, _maxElecIpOverElecIpError, _eventPUWeight);
  _muonIpOverErrVsMuonTrkLayersWithHitsHisto->Fill(_maxMuonTrkLayersWithHits, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _muonIpOverErrVsMuonValidHitsHisto->Fill(_maxMuonValidHits, _maxMuonIpOverMuonIpError, _eventPUWeight);		      

  _muonIpOverErrVsMuonTrkNdofHisto->Fill(_maxMuonTrkNdof, _maxMuonIpOverMuonIpError, _eventPUWeight);
  _elecIpOverErrVsElecTrkNdofHisto->Fill(_maxElecTrkNdof, _maxElecIpOverElecIpError, _eventPUWeight);
  
  _elecMotherIdHisto->Fill(_maxElecMotherId, _eventPUWeight);
  _muonMotherIdHisto->Fill(_maxMuonMotherId, _eventPUWeight);
  _combinedMotherIdHisto->Fill(_maxElecMotherId + _maxMuonMotherId, _eventPUWeight);
  
  _elecPtVsMuonPtHisto->Fill(_maxMuonPt, _maxElecPt, _eventPUWeight);
  _elecEtVsMuonPtHisto->Fill(_maxMuonPt, _maxElecEt, _eventPUWeight);
  
			      
  _elecMuonCosDPhiHisto->Fill(_minElecMuonCosDPhi,_eventPUWeight); 
  _elecMuonDeltaRHisto->Fill(_maxElecMuonDeltaR,_eventPUWeight); 
  _elecMuonPZetaHisto->Fill(_maxElecMuonPZeta,_eventPUWeight); 
  _leptonMetCosDPhiHisto->Fill(_minLeptonMetCosDPhi,_eventPUWeight); 
  _leptonDeltaPtCosDPhiHisto->Fill(_minLeptonDeltaPtCosDPhi, _eventPUWeight);
  _MEtHisto->Fill(_maxMEt,_eventPUWeight);   
  _elecMuonMetMassHisto->Fill(_maxElecMuonMetMass, _eventPUWeight);
  _elecMuonMetMassBeginningHisto->Fill(_maxElecMuonMetMassBeginning, _eventPUWeight);
  _elecMuonDeltaPtMassHisto->Fill(_maxElecMuonDeltaPtMass, _eventPUWeight);
  
  _elecMuonMetMassVs2xLeptonPtHisto->Fill(2*max(_maxElecPtVsMass, _maxMuonPtVsMass), _maxElecMuonMetMass, _eventPUWeight); 
  _elecMuonMetMassVsElecMuonDeltaPtMassHisto->Fill(_maxElecMuonDeltaPtMassNMinus1, _minElecMuonMetMassNMinus1, _eventPUWeight); 

  _BSy0VsBSx0Histo->Fill(_maxBSx0, _maxBSy0, _eventPUWeight);
  _muonVertex_x0Histo->Fill(_maxMuonVertex_x0, _eventPUWeight);
  _muonVertex_y0Histo->Fill(_maxMuonVertex_y0, _eventPUWeight);
  _muonVertex_z0Histo->Fill(_maxMuonVertex_z0, _eventPUWeight);
  _muonDecayLength_xHisto->Fill(_maxMuonDecayLength_x, _eventPUWeight);
  _muonDecayLength_yHisto->Fill(_maxMuonDecayLength_y, _eventPUWeight);
  _muonDecayLength_zHisto->Fill(_maxMuonDecayLength_z, _eventPUWeight);
  _elecVertex_x0Histo->Fill(_maxElecVertex_x0, _eventPUWeight);
  _elecVertex_y0Histo->Fill(_maxElecVertex_y0, _eventPUWeight);
  _elecVertex_z0Histo->Fill(_maxElecVertex_z0, _eventPUWeight);
  _elecDecayLength_xHisto->Fill(_maxElecDecayLength_x, _eventPUWeight);
  _elecDecayLength_yHisto->Fill(_maxElecDecayLength_y, _eventPUWeight);
  _elecDecayLength_zHisto->Fill(_maxElecDecayLength_z, _eventPUWeight);  
  _elecVxHisto->Fill(_maxElecVx, _eventPUWeight);
  _elecVyHisto->Fill(_maxElecVy, _eventPUWeight);
  _elecVzHisto->Fill(_maxElecVz, _eventPUWeight);
  
  float MuonDecayLength = sqrt(_maxMuonDecayLength_x*_maxMuonDecayLength_x + _maxMuonDecayLength_y*_maxMuonDecayLength_y);
  float ElecDecayLength = sqrt(_maxElecDecayLength_x*_maxElecDecayLength_x + _maxElecDecayLength_y*_maxElecDecayLength_y);
  
  _muonDecayLengthHisto->Fill(MuonDecayLength, _eventPUWeight);
  _elecDecayLengthHisto->Fill(ElecDecayLength, _eventPUWeight);
  
  _MEtVsLeptonDeltaPHisto->Fill(fabs(_maxMuonPtVsDeltaPt - _maxElecPtVsDeltaPt), _maxMEtVsDeltaPt, _eventPUWeight);
  
  _elecMuonDCAOverErrorVsCosDPhiHisto->Fill(_minElecMuonCosDPhi, _maxElecMuonDis3DOverError, _eventPUWeight);

  //N-1 histos
  _elecPtNMinus1Histo->Fill(_maxElecPtNMinus1,_eventPUWeight);
  _elecEtaNMinus1Histo->Fill(_minElecEtaNMinus1,_eventPUWeight);

  _muonPtNMinus1Histo->Fill(_maxMuonPtNMinus1,_eventPUWeight);
  _muonENMinus1Histo->Fill(_maxMuonENMinus1,_eventPUWeight);
  _muonEtaNMinus1Histo->Fill(_minMuonEtaNMinus1,_eventPUWeight);
  _muonDxyVtxNMinus1Histo->Fill(_maxMuonDxyVtxNMinus1,_eventPUWeight);
  _muonDzVtxNMinus1Histo->Fill(_maxMuonDzVtxNMinus1, _eventPUWeight);
  _muonDzBSNMinus1Histo->Fill(_maxMuonDzBSNMinus1, _eventPUWeight);
  _muonIsoNMinus1Histo->Fill(_minMuonIsoNMinus1, _eventPUWeight);
  _muonDptOverPtNMinus1Histo->Fill(_minMuonDptOverPtNMinus1, _eventPUWeight);
  _muonIsHighPtMuonNMinus1Histo->Fill(_maxMuonIsHighPtMuonNMinus1, _eventPUWeight);

  _elecIpOverElecIpErrorNMinus1Histo->Fill(_maxElecIpOverElecIpErrorNMinus1, _eventPUWeight);
  _muonIpOverMuonIpErrorNMinus1Histo->Fill(_maxMuonIpOverMuonIpErrorNMinus1, _eventPUWeight);
  _combinedIpOverIpErrorNMinus1Histo->Fill(_maxCombinedIpOverIpErrorNMinus1, _eventPUWeight);

  _elecEtNMinus1Histo->Fill(_maxElecEtNMinus1,_eventPUWeight);   
  _nBTagNMinus1Histo->Fill(_minNBtagsCombSecVtxNMinus1, _eventPUWeight); 
  _elecMuonCosDPhiNMinus1Histo->Fill(_minElecMuonCosDPhiNMinus1, _eventPUWeight);
  _elecMuonPZetaNMinus1Histo->Fill(_maxElecMuonPZetaNMinus1, _eventPUWeight);  
  _elecMuonPZetaVisNMinus1Histo->Fill(_maxElecMuonPZetaVisNMinus1, _eventPUWeight);
  _jetSumEtNMinus1Histo->Fill(_minJetSumEtNMinus1, _eventPUWeight);
  _leptonMetCosDPhiNMinus1Histo->Fill(_minLeptonMetCosDPhiNMinus1, _eventPUWeight);
  _leptonDeltaPtCosDPhiNMinus1Histo->Fill(_minLeptonDeltaPtCosDPhiNMinus1, _eventPUWeight);
  _leptonMetCosDPhiVsLeptonDeltaPtCosDPhiNMinus1Histo->Fill(_minLeptonDeltaPtCosDPhiNMinus1, _minLeptonMetCosDPhiNMinus1, _eventPUWeight);
  _elecMuonCosDPhiVsLeptonMetCosDPhiNMinus1Histo->Fill(_minLeptonMetCosDPhiNMinus1, _minElecMuonCosDPhiNMinus1, _eventPUWeight);
  
  _elecMuonMetMassNMinus1Histo->Fill(_minElecMuonMetMassNMinus1, 1);
  _elecMuonMetMassNMinus1RebinHisto->Fill(_minElecMuonMetMassNMinus1, _eventPUWeight);
  _elecMuonDeltaPtMassNMinus1Histo->Fill(_maxElecMuonDeltaPtMassNMinus1, _eventPUWeight);
  _elecMuonDeltaPtMassNMinus1RebinHisto->Fill(_maxElecMuonDeltaPtMassNMinus1, _eventPUWeight);  
  _nJetsNMinus1Histo->Fill(_maxNJetsNMinus1, _eventPUWeight);
  _elecMuonMassNMinus1Histo->Fill(_minElecMuonMassNMinus1, _eventPUWeight);
  _muonMetMtNMinus1Histo->Fill(_minMuonMetMtNMinus1, _eventPUWeight);
  _elecMuonChargeNMinus1Histo->Fill(_minElecMuonChargeNMinus1, _eventPUWeight);
  _MEtNMinus1Histo->Fill(_maxMEtNMinus1, _eventPUWeight);
  
  _threeLeptonNMinus1Histo->Fill(_maxThreeLeptonNMinus1, _eventPUWeight);
  _threeLeptonMassNMinus1Histo->Fill(_maxThreeLeptonMassNMinus1, _eventPUWeight);
  
  _muonDeltaPtOverPtNMinus1Histo->Fill(_maxMuonDptOverPtNMinus1, _eventPUWeight);
  _elecDeltaPtOverPtNMinus1Histo->Fill(_maxElecDptOverPtNMinus1, _eventPUWeight);
  
  _muonDeltaPtOverPtVsPtNMinus1Histo->Fill(_maxMuonPtNMinus1, _maxMuonDptOverPtNMinus1, _eventPUWeight);
  _elecDeltaPtOverPtVsPtNMinus1Histo->Fill(_maxElecEtNMinus1, _maxElecDptOverPtNMinus1, _eventPUWeight);
  
  _PZetaVsPZetaVisHisto->Fill(_maxElecMuonPZetaVisNMinus1, _maxElecMuonPZetaNMinus1, _eventPUWeight);
  _nBTagVsJetSumEtHisto->Fill(_minJetSumEtNMinus1, _minNBtagsCombSecVtxNMinus1, _eventPUWeight);
  _nBTagVsMEtHisto->Fill(_maxMEtNMinus1, _minNBtagsCombSecVtxNMinus1, _eventPUWeight);
  _nBTagVsMuonMetMtHisto->Fill(_minMuonMetMtNMinus1, _minNBtagsCombSecVtxNMinus1, _eventPUWeight);
  _nBTagVsNJetsHisto->Fill(_maxNJetsNMinus1, _minNBtagsCombSecVtxNMinus1, _eventPUWeight);
  
  _nBTagCR2Histo->Fill(_minNBtagsCombSecVtxCR2, _eventPUWeight);

}

void EMuAnalysis::bookHistos(){

  outRootFile = new TFile(_outRootFileName.c_str(), "Recreate");
  TDirectory* theDir; 
  std::string histosDir = _source + "HistosDirectory"; 
  if(outRootFile->FindObjectAny(histosDir.c_str()) == 0) theDir = outRootFile->mkdir(histosDir.c_str());
  else theDir = (TDirectory*)outRootFile->FindObjectAny(histosDir.c_str());
  theDir->cd();
  
  _pileUPDistro = (TH1F*)dataPUHisto->Clone("_pileUPDistro");
  _pileUPReweighted = (TH1F*)dataPUHisto->Clone("_pileUPReweighted");  
  
  _tauEnergyHisto = new TH1F("_tauEnergyHisto", "Tau Energy", 200, 0, 4000);
  _muonPOverTauEHisto = new TH1F("_muonPOverTauEHisto", "p_{#mu}/E_{#tau}", 100, 0, 1.);
  _muonPOverTauEVsMuonIpOverErrorHisto = new TH2F("_muonPOverTauEVsMuonIpOverErrorHisto", "p_{#mu}/E_{#tau} vs IP_{#mu}/#sigma_{IP}", 100, 0, 10, 100, 0, 1.);
  _muonPOverTauEVsCosDPhiHisto = new TH2F("_muonPOverTauEVsCosDPhiHisto", "p_{#mu}/E_{#tau} vs Cos#Delta#phi(e, #mu)", 100, -1, -0.95, 100, 0, 1.);
  _muonPOverTauEVsTauDecayLengthHisto = new TH2F("_muonPOverTauEVsTauDecayLengthHisto", "p_{#mu}/E_{#tau} vs Tau Decay Length", 100, 0, 1, 100, 0, 1);
  _tauEnergyVsTauDecayLengthHisto = new TH2F("_tauEnergyVsTauDecayLengthHisto", "#tau Energy vs #tau Decay Length", 100, 0, 10, 400, 0, 4000);
  
  _zPrimeEnergyHisto = new TH1F("_zPrimeEnergyHisto", "Z' Energy", 200, 0, 4000);
  _vertexDiffHisto = new TH1F("_vertexDiffHisto", "3D Distance Between Generated ZPrime Vertex and Reconstructed Primary", 100, 0, 0.05);
  _vertexDiffxyHisto = new TH1F("_vertexDiffxyHisto", "xy Distance Between Generated ZPrime Vertex and Reconstructed Primary", 100, 0, 0.01);
  _impactParameterVsVertexDiffHisto = new TH2F("_impactParameterVsVertexDiffHisto", "Combined Ip/#sigma vs 3D Distance Between Generated ZPrime Vertex and Reco Vertex", 100, 0, 0.1, 100, 0, 30.);
  _beamSpotDiffHisto = new TH1F("_beamSpotDiffHisto", "3D Distance Between Generated ZPrime Vertex and Beam Spot", 100, 0, 1.);
  _beamSpotDiffxyHisto = new TH1F("_beamSpotDiffxyHisto", "xy Distance Between Generated ZPrime Vertex and Beam Spot", 100, 0, 0.01);
  _beamSpotPrimaryzDiffHisto = new TH1F("_beamSpotPrimaryzDiffHisto", "3D Distance Between Generated ZPrime Vertex and Beam Spot (z from Primary)", 100, 0, 0.05);

  
  _tauEnergyVsElecMuonMetMassHisto = new TH2F("_tauEnergyVsElecMuonMetMassHisto", "Tau Energy vs MassVis(e, #mu, ME_{T})", 200, 0, 2000, 400, 0, 4000);
  _zPrimeEnergyVsElecMuonMetMassHisto = new TH2F("_zPrimeEnergyVsElecMuonMetMassHisto", "Z' Energy vs MassVis(e, #mu, ME_{T})", 200, 0, 2000, 400, 0, 4000);
  _tauEnergyVsMuonEtaHisto = new TH2F("_tauEnergyVsMuonEtaHisto", "Tau Energy vs Muon |#eta|", 100, -2.5, 2.5, 400, 0, 4000);
  _zPrimeEnergyVsMuonEtaHisto = new TH2F("_zPrimeEnergyVsMuonEtaHisto", "Z' Energy vs Muon |#eta|", 100, -2.5, 2.5, 400, 0, 4000);
  _tauEnergyVsElecEtaHisto = new TH2F("_tauEnergyVsElecEtaHisto", "Tau Energy vs Elec |#eta|", 100, -2.5, 2.5, 400, 0, 4000);
  _zPrimeEnergyVsElecEtaHisto = new TH2F("_zPrimeEnergyVsElecEtaHisto", "Z' Energy vs Elec |#eta|", 100, -2.5, 2.5, 400, 0, 4000);  
  _zPrimeEnergyVs2xMuonPtHisto = new TH2F("_zPrimeEnergyVs2xMuonPtHisto", "Z' Energy vs 2 x #mu p_{T}", 100, 0, 1000, 400, 0, 4000);
  _zPrimeEnergyVs2xElecEtHisto = new TH2F("_zPrimeEnergyVs2xElecEtHisto", "Z' Energy vs 2 x Elec e_{T}", 100, 0, 1000, 400, 0, 4000);
  _zPrimeEnergyVs2xMaxLeptonPtHisto = new TH2F("_zPrimeEnergyVs2xMaxLeptonPtHisto", "Z' Energy vs 2 x Max(#mu p_{T} , Elec e_{T})", 200, 0, 2000, 400, 0, 4000);
  _zPrimeEnergyVs2xMaxLeptonPHisto = new TH2F("_zPrimeEnergyVs2xMaxLeptonPHisto", "Z' Energy vs 2 x Max(#mu p_ , Elec p)", 200, 0, 2000, 400, 0, 4000);
  
  
  _nVtxHisto = new TH1I("_nVtxHisto", "NVtx", 35, 0., 35);
  _elecPtHisto = new TH1F("_elecPtHisto", "elec p_{T}", 50, 0., 250.);
  _elecEtaHisto = new TH1F(" _elecEtaHisto", "elec #eta", 50, -2.5, 2.5);
  _elecDeltaPhiInBarrelHisto = new TH1F("_elecDeltaPhiInBarrelHisto", "elec #Delta#phi_{in} Barrel", 200, -0.5, 0.5);
  _elecDeltaPhiInEndCapHisto = new TH1F("_elecDeltaPhiInEndCapHisto", "elec #Delta#phi_{in} EndCap", 200, -0.5, 0.5);
  _elecDeltaEtaInBarrelHisto = new TH1F("_elecDeltaEtaInBarrelHisto", "elec #Delta#eta_{in} Barrel", 100, -0.1, 0.1);
  _elecDeltaEtaInEndCapHisto = new TH1F("_elecDeltaEtaInEndCapHisto", "elec #Delta#eta_{in} EndCap", 100, -0.1, 0.1);
  _elecSigmaIEtaIEtaBarrelHisto = new TH1F("_elecSigmaIEtaIEtaBarrelHisto", "elec #sigma_{i#etai#eta} Barrel", 100, 0, 0.5);
  _elecSigmaIEtaIEtaEndCapHisto = new TH1F("_elecSigmaIEtaIEtaEndCapHisto", "elec #sigma_{i#etai#eta} EndCap", 100, 0, 0.5);
  _elecTrkChiSqrdHisto = new TH1F("_elecTrkChiSqrdHisto", "Elec Trk Chi Sqrd", 100, 0, 20);
  _elecHadFracBarrelHisto = new TH1F("_elecHadFracBarrelHisto", "elec H/E", 100, 0, 10);
  _elecHadFracEndCapHisto = new TH1F("_elecHadFracEndCapHisto", "elec H/E", 100, 0, 10);
  _elecEOverPBarrelHisto = new TH1F("_elecEOverPBarrelHisto", "elec (1/E - 1/p)", 200, -1., 1.);
  _elecEOverPEndCapHisto = new TH1F("_elecEOverPEndCapHisto", "elec (1/E - 1/p)", 200, -1., 1.);
  _elecMissingHitsHisto = new TH1I("_elecMissingHitsHisto", "elec Missing Hits", 10, 0, 10);
  _elecPFDR03ChargedHadIsoHisto = new TH1F("_elecPFDR03ChargedHadIsoHisto", "elec Charged Had PFIso #Delta R = 0.3", 100, 0., 100.);
  _elecPFDR03NeutralHadIsoHisto = new TH1F("_elecPFDR03NeutralHadIsoHisto", "elec Neutral Had PFIso #Delta R = 0.3", 100, 0., 100.);
  _elecPFDR03PhotonIsoHisto = new TH1F("_elecPFDR03PhotonIsoHisto", "elec Photon PFIso #Delta R = 0.3", 100, 0., 100.);
  _elecPFDR03RelIsoHisto = new TH1F("_elecPFDR03RelIsoHisto", "elec Charged + Neutral + Photon RelPFIso #Delta R = 0.3", 100, 0., 10.);
  _elecPFDR04ChargedHadIsoHisto = new TH1F("_elecPFDR04ChargedHadIsoHisto", "elec Charged Had PFIso #Delta R = 0.4", 100, 0., 100.);
  _elecPFDR04NeutralHadIsoHisto = new TH1F("_elecPFDR04NeutralHadIsoHisto", "elec Neutral Had PFIso #Delta R = 0.4", 100, 0., 100.);
  _elecPFDR04PhotonIsoHisto = new TH1F("_elecPFDR04PhotonIsoHisto", "elec Photon PFIso #Delta R = 0.4", 100, 0., 100.);
  _elecPFDR04RelIsoHisto = new TH1F("_elecPFDR04RelIsoHisto", "elec Charged + Neutral + Photon RelPFIso #Delta R = 0.4", 100, 0., 10.);  
  _elecRelIsoHisto = new TH1F("_elecRelIsoHisto", "Electron Iso: (Ecal + Hcal + Trk)/Pt", 500, 0., 5.); 
  _elecTrkIsoHisto = new TH1F("_elecTrkIsoHisto", "Electron Track Isolation", 100, 0., 5.);
  
  _muonPtHisto = new TH1F("_muonPtHisto", "Muon p_{T}", 50, 0., 200.);
  _muonEtaHisto = new TH1F("_muonEtaHisto", "Muon Eta", 50, -2.5, 2.5);
  _muonPtVsMuonEtaHisto = new TH2F("_muonPtVsMuonEtaHisto", "Muon P_{T} vs Muon #eta", 50, -2.5, 2.5, 200, 0, 1000);
  _muonEVsMuonEtaHisto = new TH2F("_muonEVsMuonEtaHisto", "Muon E vs Muon #eta", 50, -2.5, 2.5, 200, 0, 1000);
  
  
  _muonDxyVtxHisto = new TH1F("_muonDxyVtxHisto", "Muon Dxy Vtx", 1000, -0.2, 0.2);
  _muonDzVtxHisto = new TH1F("_muonDzVtxHisto", "Muon Dz Vtx", 1000, -5, 5);
   
  _muonIsGlobalMuonHisto = new TH1I("_muonIsGlobalMuonHisto", "Is Global Muon", 2, 0, 2);
  _muonMatchedHisto = new TH1I("_muonMatchedHisto", "Muon Matched", 2, 0, 2);
  _muonNormChiSqrdHisto = new TH1F("_muonNormChiSqrdHisto", "Muon Norm Chi Sqrd", 50, 0, 15);
  _muonTrkChiSqrdHisto = new TH1F("_muonTrkChiSqrdHisto", "Muon Trk Chi Sqrd", 100, 0, 20);
  _muonValidHitsHisto = new TH1I("_muonValidHitsHisto", "Muon Valid Hits", 50, 0, 10);
  _muonMatchedStationsHisto = new TH1I("_muonMatchedStationsHisto", "Muon Matched Stations", 50, 0, 20);
  _muonPixelHitsHisto = new TH1I("_muonPixelHitsHisto", "Muon Pixel Hits", 50, 0, 10);
  _muonTrkLayersWithHitsHisto = new TH1F("_muonTrkLayersWithHitsHisto", "Muon Track Layers with Hits", 50, 0, 50);
  _muonIsoDR04DBHisto = new TH1F("_muonIsoDR04DBHisto", "Muon Iso: (Charged + Neutral + Photon)/Pt", 500, 0., 5.); 
  _muonDptOverPtHisto = new TH1F("_muonDptOverPtHisto", "Muon TrkPt_Err/TrkPt", 100, 0., 10.);
  _muonIsHighPtMuonHisto = new TH1I("_muonIsHighPtMuonHisto", "Is High Pt Muon", 2, 0, 2);
   
  _elecMatchedHisto = new TH1I("_elecMatchedHisto", "Elec Matched", 2, 0, 2);
  _heepPassedPtHisto = new TH1I("_heepPassedPtHisto", "HEEP Passed Pt", 2, 0, 2);
  _heepPassedDetEtaHisto = new TH1I("_heepPassedDetEtaHisto", "HEEP Passed Delta Eta", 2, 0, 2);
  _heepPassedCrackHisto = new TH1I("_heepPassedCrackHisto", "HEEP Passed Crack", 2, 0, 2);
  _heepPassedDEtaInHisto = new TH1I("_heepPassedDEtaInHisto", "HEEP Passed Delta Eta In", 2, 0, 2);
  _heepPassedDPhiInHisto = new TH1I("_heepPassedDPhiInHisto", "HEEP Passed Delta Phi In", 2, 0, 2);
  _heepPassedHademHisto = new TH1I("_heepPassedHademHisto", "HEEP Passed HadEM", 2, 0, 2);
  _heepPassedSigmaIEtaIEtaHisto = new TH1I("_heepPassedSigmaIEtaIEtaHisto", "HEEP Passed SigmaIEtaIEta", 2, 0, 2);
  _heepPassed2by5Over5By5Histo = new TH1I("_heepPassed2by5Over5By5Histo", "HEEP Passed 2 by 5 Over 5 by 5", 2, 0, 2);
  _heepPassedEcalHad1IsoHisto = new TH1I("_heepPassedEcalHad1IsoHisto", "HEEP Passed Ecal Had 1 Iso", 2, 0, 2);
  //_heepPassedHad2IsoHisto = new TH1F("_heepPassedHad2IsoHisto", "HEEP Passed Had 2 Iso", 2, 0, 2);
  _heepPassedTrkIsoHisto = new TH1I("_heepPassedTrkIsoHisto", "HEEP Passed Track Iso", 2, 0, 2);
  _heepPassedEcalDrivenHisto = new TH1I("_heepPassedEcalDrivenHisto", "HEEP Passed Ecal Driven", 2, 0, 2);  
  /*
  _elecIpOverElecIpErrorHisto = new TH1F("_elecIpOverElecIpErrorHisto", "|Electron Ip/err|", 50, 0, 10);
  _muonIpOverMuonIpErrorHisto = new TH1F("_muonIpOverMuonIpErrorHisto", "|Muon Ip/err|", 50, 0, 10);
  _combinedIpOverIpErrorHisto = new TH1F("_combinedIpOverIpErrorHisto", "|e-Ip/err| + |mu-Ip/err|", 50, 0, 10);
  */
  _elecIpOverElecIpErrorHisto = new TH1F("_elecIpOverElecIpErrorHisto", "|Electron Ip/err|", 100, 0, 5);
  _muonIpOverMuonIpErrorHisto = new TH1F("_muonIpOverMuonIpErrorHisto", "|Muon Ip/err|", 100, 0, 5);
  _combinedIpOverIpErrorHisto = new TH1F("_combinedIpOverIpErrorHisto", "sqrt((e-Ip^{2} + #mu-Ip^{2})/(#sigma_{e-ip}^{2} + #sigma_{#mu-ip}^{2}))", 100, 0, 5);  
  _elecIpErrorHisto = new TH1F("_elecIpErrorHisto", "Elec Ip Error", 100, 0, 0.1);
  _muonIpErrorHisto = new TH1F("_muonIpErrorHisto", "Muon Ip Error", 100, 0, 0.01);
  
  
  _elecSCE1x5Histo = new TH1F("_elecSCE1x5Histo", "Elec SCE 1x5", 100, 0, 300);
  _elecSCE2x5Histo = new TH1F("_elecSCE2x5Histo", "Elec SCE 2x5", 100, 0, 300);
  _elecSCE5x5Histo = new TH1F("_elecSCE5x5Histo", "Elec SCE 5x5", 100, 0, 300);
  _elecSCE1x5Over5x5Histo = new TH1F("_elecSCE1x5Over5x5Histo", "Elec SCE 1x5", 100, 0, 1);
  _elecSCE2x5MaxOver5x5Histo = new TH1F("_elecSCE2x5MaxOver5x5Histo", "Elec SCE 1x5", 100, 0, 1);
  _elecSCE1x5EOverPHisto = new TH1F("_elecSCE1x5EOverPHisto", "Elec SCE/p 1x5", 100, 0, 1);
  _elecSCE2x5EOverPHisto = new TH1F("_elecSCE2x5EOverPHisto", "Elec SCE/p 2x5", 100, 0, 1);
  _elecSCE5x5EOverPHisto = new TH1F("_elecSCE5x5EOverPHisto", "Elec SCE/p 5x5", 100, 0, 1);
  _muonEHisto = new TH1F("_muonEHisto", "Muon E", 100, 0, 1000);
  
  _elecSCE1x5VsElecIpOverError_hiIPHisto = new TH2F("_elecSCE1x5VsElecIpOverError_hiIPHisto", "Elec SCE 1x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 300);
  _elecSCE2x5VsElecIpOverError_hiIPHisto = new TH2F("_elecSCE2x5VsElecIpOverError_hiIPHisto", "Elec SCE 2x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 300);
  _elecSCE5x5VsElecIpOverError_hiIPHisto = new TH2F("_elecSCE5x5VsElecIpOverError_hiIPHisto", "Elec SCE 5x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 300);
  _elecSCE1x5Over5x5VsElecIpOverError_hiIPHisto = new TH2F("_elecSCE5x5VsElecIpOverError_hiIPHisto", "Elec SCE 1x5/5x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 100);
  _elecSCE2x5MaxOver5x5VsElecIpOverError_hiIPHisto = new TH2F("_elecSCE5x5VsElecIpOverError_hiIPHisto", "Elec SCE 2x5/5x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 100);
  _elecSCE1x5EOverPVsElecIpOverError_hiIPHisto = new TH2F("_elecSCE1x5EOverPVsElecIpOverError_hiIPHisto", "Elec SCE/p 1x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 1);
  _elecSCE2x5EOverPVsElecIpOverError_hiIPHisto = new TH2F("_elecSCE2x5EOverPVsElecIpOverError_hiIPHisto", "Elec SCE/p 2x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 1);
  _elecSCE5x5EOverPVsElecIpOverError_hiIPHisto = new TH2F("_elecSCE5x5EOverPVsElecIpOverError_hiIPHisto", "Elec SCE/p 5x5 vs Elec IP/#sigma", 100, 0, 10, 100, 0, 1);
  
  _elecSCE1x5EOverPVsEOverP_hiIPHisto = new TH2F("_elecSCE1x5EOverPVsEOverP_hiIPHisto", "Elec SCE/p 1x5 vs E/p", 100, 0, 10, 100, 0, 1);
  _elecSCE2x5EOverPVsEOverP_hiIPHisto = new TH2F("_elecSCE2x5EOverPVsEOverP_hiIPHisto", "Elec SCE/p 2x5 vs E/p", 100, 0, 10, 100, 0, 1);
  _elecSCE5x5EOverPVsEOverP_hiIPHisto = new TH2F("_elecSCE5x5EOverPVsEOverP_hiIPHisto", "Elec SCE/p 5x5 vs E/p", 100, 0, 10, 100, 0, 1);
  
  
  _elecMuonDistanceOfClosestApproachHisto = new TH1F("_elecMuonDistanceOfClosestApproachHisto", "Distance Of Closest Approach (e, #mu) in RPhi", 100, 0, 0.1);
  _elecMuonDCAVsCrossingPointHisto = new TH2F("_elecMuonDCAVsCrossingPointHisto", "DCA(e, #mu) vs R(track crossing)", 100, 0, 10, 100, 0, 0.1);
  _elecMuonTwoTrkMinDistHisto = new TH1F("_elecMuonTwoTrkMinDistHisto", "DCA (e trk, #mu trk)", 100, 0, 0.05);
  _elecMuonTwoTrkMinDistVsCrossingPointHisto = new TH2F("_elecMuonTwoTrkMinDistVsCrossingPointHisto", "DCA (e trk, #mu trk) vs R(track crossing)", 100, 0, 10, 100, 0, 0.05);
  _eMuDis3DHisto = new TH1F("_eMuDis3DHisto", "Distance of Closest Approach (e, #mu) in 3D", 100, 0, 0.01);
  _eMuDis2DHisto = new TH1F("_eMuDis2DHisto", "Distance of Closest Approach (e, #mu) in 2D", 100, 0, 0.01);
  _elecMuonDCA3DErrorHisto = new TH1F("_elecMuonDCA3DErrorHisto", "Error on 3D DCA (e, #mu)", 100, 0, 0.05);
  _elecMuonDCA2DErrorHisto = new TH1F("_elecMuonDCA2DErrorHisto", "Error on 2D DCA (e, #mu)", 100, 0, 0.05);
  _eMuDis3DOverErrorHisto = new TH1F("_eMuDis3DOverErrorHisto", "DCA/Error (e, #mu) in 3D", 100, 0, 50);
  _eMuDis2DOverErrorHisto = new TH1F("_eMuDis2DOverErrorHisto", "DCA/Error (e, #mu) in 2D", 100, 0, 50);
  
  
  _IpVsElecTrkChiSqrdHisto = new TH2F("_IpVsElecTrkChiSqrdHisto", "IP vs Elec Trk Chi Sqrd", 100, 0, 20, 100, 0, 10);
  _IpVsMuonTrkChiSqrdHisto = new TH2F("_IpVsMuonTrkChiSqrdHisto", "IP vs Muon Trk Chi Sqrd", 100, 0, 20, 100, 0, 10);
 
 
  _DCAOverErrorVsElecIpOverErrorHisto = new TH2F("_DCAOverErrorVsElecIpOverErrorHisto", "DCA/#sigma(e, #mu) vs IP/#sigma(e)", 100, 0, 10, 100, 0, 10);
  _DCAOverErrorVsMuonIpOverErrorHisto = new TH2F("_DCAOverErrorVsMuonIpOverErrorHisto", "DCA/#sigma(e, #mu) vs IP/#sigma(#mu)", 100, 0, 10, 100, 0, 10);
  _DCAOverErrorVsCombinedIpOverErrorHisto = new TH2F("_DCAOverErrorVsCombinedIpOverErrorHisto", "DCA/#sigma(e, #mu) vs IP/#sigma(e + #mu)", 100, 0, 10, 100, 0, 10);
  
  
  
  _elecIpVsElecTrkChiSqrdHisto = new TH2F("_elecIpVsElecTrkChiSqrdHisto", "ElecIp/Err vs Elec Trk Chi Sqrd", 100, 0, 20, 100, 0, 10);
  _muonIpVsMuonTrkChiSqrdHisto = new TH2F("_muonIpVsMuonTrkChiSqrdHisto", "MuonIp/Err vs Muon Trk Chi Sqrd", 100, 0, 20, 100, 0, 10);
  _elecIpErrVsElecTrkChiSqrdHisto = new TH2F("_elecIpErrVsElecTrkChiSqrdHisto", "ElecIpErr vs Elec Trk Chi Sqrd", 100, 0, 20, 100, 0, 0.1);
  _muonIpErrVsMuonTrkChiSqrdHisto = new TH2F("_muonIpErrVsMuonTrkChiSqrdHisto", "MuonIpErr vs Muon Trk Chi Sqrd", 100, 0, 20, 100, 0, 0.01);
  _elecIpOverErrVsElecIpErrHisto = new TH2F("_elecIpOverErrVsElecIpErrHisto", "ElecIpOverErr vs ElecIpErr", 100, 0, 0.1, 100, 0, 10);
  _muonIpOverErrVsMuonIpErrHisto = new TH2F("_muonIpOverErrVsMuonIpErrHisto", "MuonIpOverErr vs MuonIpErr", 100, 0, 0.01, 100, 0, 10);
  _muonIpOverErrVsElecIpOverErrHisto = new TH2F("_muonIpOverErrVsElecIpOverErrHisto", "Muon IP/#sigma vs Elec IP/#sigma", 100, 0, 10, 100, 0, 10);
  _muonIpOverErrVsMuonPixelHitsHisto = new TH2F("_muonIpOverErrVsMuonPixelHitsHisto", "MuonIp/Err vs Muon Pixel Hits", 10, 0, 10, 100, 0, 10);
  _muonEtaVsMuonPixelHitsHisto = new TH2F("_muonEtaVsMuonPixelHitsHisto", "Muon Eta vs Muon Pixel Hits", 10, 0, 10, 50, -2.5, 2.5);
  _muonIpOverErrVsMuonEtaHisto = new TH2F("_muonIpOverErrVsMuonEtaHisto", "MuonIp/Err vs Muon Eta", 50, -2.5, 2.5, 100, 0, 10);
  
  _elecIpOverErrVsElecMissingHitsHisto = new TH2F("_elecIpOverErrVsElecMissingHitsHisto", "ElecIp/Err vs Elec Missing Hits", 10, 0, 10, 100, 0, 10);
  _muonIpOverErrVsMuonTrkLayersWithHitsHisto = new TH2F("_muonIpOverErrVsMuonTrkLayersWithHitsHisto", "MuonIp/Err vs Muon Track Layers with Hits", 50, 0, 50, 100, 0, 10);
  _muonIpOverErrVsMuonValidHitsHisto = new TH2F("_muonIpOverErrVsMuonValidHitsHisto", "MuonIp/Err vs Muon Valid Hits", 10, 0, 10, 100, 0, 10);
  
  _muonPtVsMuonIpOverErrHisto = new TH2F("_muonPtVsMuonIpOverErrHisto", "Muon Pt vs Muon Ip/Err", 100, 0, 10, 100, 0, 100);
  _elecEtVsElecIpOverErrHisto = new TH2F("_elecEtVsElecIpOverErrHisto", "Elec Et vs Elec Ip/Err", 100, 0, 10, 100, 0, 100);
  _MEtVsElecIpOverErrHisto = new TH2F("_MEtVsElecIpOverErrHisto", "MEt vs Elec Ip/Err", 100, 0, 10, 100, 0, 100);
  _MEtVsMuonIpOverErrHisto = new TH2F("_MEtVsMuonIpOverErrHisto", "MEt vs Muon Ip/Err", 100, 0, 10, 100, 0, 100);
  _MEtVsCombinedIpOverErrHisto = new TH2F("_MEtVsCombinedIpOverErrHisto", "MEt vs Combined Ip/Err", 100, 0, 10, 100, 0, 100);
  
  _elecEOverP_hiIPHisto = new TH1F("_elecEOverP_hiIPHisto", "Electron E/p (hi IP)", 100, 0, 10);
  _elecEta_hiIPHisto = new TH1F("_elecEta_hiIPHisto", "Electron SC #eta (hi IP)", 100, -2.5, 2.5);
  _elecEOverPVsElecIpOverError_hiIPHisto = new TH2F("_elecEOverPVsElecIpOverError_hiIPHisto", "Electron E/p vs Electron IP/#sigma", 100, 0, 10, 100, 0, 10);
  _elecEtaVsElecIpOverError_hiIPHisto = new TH2F("_elecEtaVsElecIpOverError_hiIPHisto", "Electron SC #eta vs Electron IP/#sigma", 100, 0, 10, 100, -2.5, 2.5);  
  _elecEtVsElecIpOverError_hiIPHisto = new TH2F("_elecEtVsElecIpOverError_hiIPHisto", "Electron E_{T} vs Electron IP/#sigma", 100, 0, 10, 100, 0, 300);  
  _elecEOverPVsMuonIpOverError_hiIPHisto = new TH2F("_elecEOverPVsMuonIpOverError_hiIPHisto", "Electron E/p vs Muon IP/#sigma", 100, 0, 10, 100, 0, 10);
  _elecTrkChiSqrdVsElecIpOverError_hiIPHisto = new TH2F("_elecTrkChiSqrdVsElecIpOverError_hiIPHisto", "Electron Track #chi^{2} vs Electron IP/#sigma", 100, 0, 10, 100, 0, 10);
  _elecTrkChiSqrdVsElecEOverP_hiIPHisto = new TH2F("_elecTrkChiSqrdVsElecEOverP_hiIPHisto", "Electron Track #chi^{2} vs Electron E/p", 100, 0, 10, 100, 0, 10);
  _elecEOverPVsElecMissingHits_hiIPHisto = new TH2F("_elecEOverPVsElecMissingHits_hiIPHisto", "Electron E/p vs n(electron missing hits)", 5, 0, 5, 100, 0, 10);
  _elecEOverPVsElecIpError_hiIPHisto = new TH2F("_elecEOverPVsElecIpError_hiIPHisto", "Electron E/p vs #sigma_{IP}", 100, 0, 0.1, 100, 0, 10);
  _elecPVsElecIpError_hiIPHisto = new TH2F("_elecPVsElecIpError_hiIPHisto", "Electron p vs #sigma_{IP}", 100, 0, 0.05, 100, 0, 300);
  _elecPVsElecTrkChiSqrd_hiIPHisto = new TH2F("_elecPVsElecTrkChiSqrd_hiIPHisto", "Electron p vs Electron Track #chi^{2}", 100, 0, 10, 100, 0, 300);
  _elecEOverPVsElecP_hiIPHisto = new TH2F("_elecEOverPVsElecP_hiIPHisto", "Electron E/p vs Electron p", 100, 0, 300, 100, 0, 10);
  
  _muonIpOverErrVsMuonIpHisto = new TH2F("_muonIpOverErrVsMuonIpHisto", "Muon Ip/Err vs Muon Ip", 100, 0, 10, 100, 0, 10);
  _elecIpOverErrVsElecIpHisto = new TH2F("_elecIpOverErrVsElecIpHisto", "Elec Ip/Err vs Elec Ip", 100, 0, 10, 100, 0, 10);
  
  _muonIpOverErrVsMuonTrkNdofHisto = new TH2F("_muonIpOverErrVsMuonTrkNdofHisto", "MuonIp/Err vs Muon Trk #DOF", 50, 0, 50, 100, 0, 10);
  _elecIpOverErrVsElecTrkNdofHisto = new TH2F("_elecIpOverErrVsElecTrkNdofHisto", "ElecIp/Err vs Elec Trk #DOF", 50, 0, 50, 100, 0, 10);
  
  
  _muonIpHisto = new TH1F("_muonIpHisto", "Muon IP", 1000, -0.1, 0.1);
  _elecIpHisto = new TH1F("_elecIpHisto", "Elec IP", 1000, -0.1, 0.1);
  
  _BSy0VsBSx0Histo = new TH2F("_BSy0VsBSx0Histo", "Beam Spot Y vs. X", 100, 0, 0.24402, 100, 0.39282, 0.39285);
  
  _muonVertex_x0Histo = new TH1F("_muonVertex_x0Histo", "Muon Vertex x_{0}", 100, 0, 1);
  _muonVertex_y0Histo = new TH1F("_muonVertex_y0Histo", "Muon Vertex y_{0}", 100, 0, 1);
  _muonVertex_z0Histo = new TH1F("_muonVertex_z0Histo", "Muon Vertex z_{0}", 100, 0, 1);
  _muonDecayLength_xHisto = new TH1F("_muonDecayLength_xHisto", "Muon Decay x Length (l_{x})", 100, 0, 1);
  _muonDecayLength_yHisto = new TH1F("_muonDecayLength_yHisto", "Muon Decay y Length (l_{y})", 100, 0, 1);
  _muonDecayLength_zHisto = new TH1F("_muonDecayLength_zHisto", "Muon Decay z Length (l_{z})", 100, 0, 1);
  
  _muonDecayLengthHisto = new TH1F("_muonDecayLengthHisto", "Muon Decay Length", 100, 0, 0.1);
  
  _elecVertex_x0Histo = new TH1F("_elecVertex_x0Histo", "Elec Vertex x_{0}", 100, 0, 1);
  _elecVertex_y0Histo = new TH1F("_elecVertex_y0Histo", "Elec Vertex y_{0}", 100, 0, 1);
  _elecVertex_z0Histo = new TH1F("_elecVertex_z0Histo", "Elec Vertex z_{0}", 100, 0, 1);
  _elecDecayLength_xHisto = new TH1F("_elecDecayLength_xHisto", "Elec Decay x Length (l_{x})", 100, 0, 1);
  _elecDecayLength_yHisto = new TH1F("_elecDecayLength_yHisto", "Elec Decay y Length (l_{y})", 100, 0, 1);
  _elecDecayLength_zHisto = new TH1F("_elecDecayLength_zHisto", "Elec Decay z Length (l_{z})", 100, 0, 1);  
  
  _elecDecayLengthHisto = new TH1F("_elecDecayLengthHisto", "Elec Decay Length", 100, 0, 0.05);
  
  _elecVxHisto = new TH1D("_elecVxHisto", "Elec Vx", 100, 0, 1);
  _elecVyHisto = new TH1D("_elecVyHisto", "Elec Vy", 100, 0, 1);
  _elecVzHisto = new TH1D("_elecVzHisto", "Elec Vz", 100, 0, 1);
  
  _elecMotherIdHisto = new TH1I("_elecMotherIdHisto", "Electron Mother Id", 101, -50, 50);
  _muonMotherIdHisto = new TH1I("_muonMotherIdHisto", "Muon Mother Id", 101, -50, 50);
  _combinedMotherIdHisto = new TH1I("_combinedMotherIdHisto", "Elec Mother Id + Muon Mother Id", 1001, 0, 1000);
  
  _elecPtVsMuonPtHisto = new TH2F("_elecPtVsMuonPtHisto", "Elec Pt Vs Muon Pt", 50, 0, 250, 50, 0, 250);
  _elecEtVsMuonPtHisto = new TH2F("_elecEtVsMuonPtHisto", "Elec Et Vs Muon Pt", 50, 0, 250, 50, 0, 250);
 
  _elecEtHisto = new TH1F("_elecEtHisto", "Electron Et", 5, 0, 100); 
  
  _MetMinusDeltaPtHisto = new TH1F("_MetMinusDeltaPtHisto", "MEt - #Delta p_{T}", 100, 0, 500); 
  _MetxVsDeltaPtxHisto = new TH2F("_MetxVsDeltaPtxHisto", "MEt_{x} vs #Delta p_{T_{x}}", 100, -5, 5, 100, -5, 5);
  _MetyVsDeltaPtyHisto = new TH2F("_MetyVsDeltaPtyHisto", "MEt_{y} vs #Delta p_{T_{y}}", 100, -5, 5, 100, -5, 5);
  
  _elecMuonCosDPhiHisto = new TH1F("_elecMuonCosDPhiHisto", "Elec-Muon Cos #Delta #phi", 50, -2, 2);
  _elecMuonDeltaRHisto = new TH1F("_elecMuonDeltaRHisto", "Elec-Muon #Delta R", 100, -20, 20);
  _elecMuonPZetaHisto = new TH1F("_elecMuonPZetaHisto", "Elec-Muon P #zeta", 100, -20, 20);
  _leptonMetCosDPhiHisto = new TH1F("_leptonMetCosDPhiHisto", "Lepton-MEt Cos #Delta #phi", 50, -2, 2);
  _leptonDeltaPtCosDPhiHisto = new TH1F("_leptonDeltaPtCosDPhiHisto", "Lepton-#Delta p_{T} Cos #Delta #phi", 50, -2, 2);
  _MEtHisto = new TH1F("_MEtHisto", "MEt", 100, 0, 500);     
  _elecMuonMetMassHisto = new TH1F("_elecMuonMetMassHisto", "M_{T} (e, #mu, ME_{T})", 150, 0, 1500);
  _elecMuonMetMassBeginningHisto = new TH1F("_elecMuonMetMassBeginningHisto", "M_{T} (e, #mu, ME_{T})", 150, 0, 1500);
  
  _elecMuonDeltaPtMassHisto = new TH1F("_elecMuonDeltaPtMassHisto", "M_{T} (e, #mu, #Delta p_{T})", 150, 0, 1500);
  
  _elecMuonMetMassVs2xLeptonPtHisto = new TH2F("_elecMuonMetMassVs2xLeptonPtHisto", "M_{T} (e, #mu ME_{T}) vs 2 x p of larger p lepton", 100,
  0, 1000, 100, 0, 1000);
  _elecMuonMetMassVsElecMuonDeltaPtMassHisto = new TH2F("_elecMuonMetMassVsElecMuonDeltaPtHisto", "M_{T} (e, #mu ME_{T}) vs M_{T} (e, #mu, #Delta p_{T})", 150, 0, 1500, 150, 0, 1500);
  
  _MEtVsLeptonDeltaPHisto = new TH2F("_MEtVsLeptonDeltaPHisto", "MEt vs #Delta p (e, #mu)", 50, 0, 500, 100, 0, 500);
  _elecMuonDCAOverErrorVsCosDPhiHisto = new TH2F("_elecMuonDCAOverErrorVsCosDPhiHisto", "DCA/#sigma (e, #mu) vs Cos #Delta#phi(e, #mu)", 100, -1, -0.9, 100, 0, 10);
  // N-1 Histos  
  
  _elecPtNMinus1Histo = new TH1F("_elecPtNMinus1Histo", "elec p_{T}", 25, 17., 267.);
  _elecEtaNMinus1Histo = new TH1F("_elecEtaNMinus1Histo", "elec #eta", 50, -2.5, 2.5);

  _muonPtNMinus1Histo = new TH1F("_muonPtNMinus1Histo", "Muon p_{T}", 40, 0, 200);
  _muonENMinus1Histo = new TH1F("_muonENMinus1Histo", "Muon E", 50, 0, 100);
  _muonEtaNMinus1Histo = new TH1F("_muonEtaNMinus1Histo", "Muon Eta NMinus1", 50, -2.5, 2.5);
  _muonDxyVtxNMinus1Histo = new TH1F("_muonDxyVtxNMinus1Histo", "Muon Dxy BS NMinus1", 100, -0.5, 0.5);
  _muonDzVtxNMinus1Histo = new TH1F("_muonDzVtxNMinus1Histo", "Muon Dz Vtx NMinus1", 1000, -1, 1);
  _muonDzBSNMinus1Histo = new TH1F("_muonDzBSNMinus1Histo", "Muon Dz BS NMinus1", 1000, -100, 100);
  _muonIsoNMinus1Histo = new TH1F("_muonIsoNMinus1Histo", "Muon Iso NMinus1", 100, 0, 20);
  _muonDptOverPtNMinus1Histo = new TH1F("_muonDptOverPtNMinus1Histo", "Muon TrkPt_Err/TrkPt NMinus1", 100, 0., 1.);
  _muonIsHighPtMuonNMinus1Histo = new TH1I("_muonIsHighPtMuonNMinus1Histo", "Is High Pt Muon NMinus1", 2, 0, 2); 
   
  _elecIpOverElecIpErrorNMinus1Histo = new TH1F("_elecIpOverElecIpErrorNMinus1Histo", "|Electron Ip/err| NMinus1", 100, 0, 5);
  _muonIpOverMuonIpErrorNMinus1Histo = new TH1F("_muonIpOverMuonIpErrorNMinus1Histo", "|Muon Ip/err| NMinus1", 100, 0, 5);

  //_elecIpOverElecIpErrorNMinus1Histo = new TH1F("_elecIpOverElecIpErrorNMinus1Histo", "|Electron Ip/err| NMinus1", 50, 0, 10);
  //_muonIpOverMuonIpErrorNMinus1Histo = new TH1F("_muonIpOverMuonIpErrorNMinus1Histo", "|Muon Ip/err| NMinus1", 50, 0, 10);
  _combinedIpOverIpErrorNMinus1Histo = new TH1F("_combinedIpOverIpErrorNMinus1Histo", "sqrt(|e-Ip/err|^2 + |mu-Ip/err|^2) NMinus1", 100, 0, 5);
  
  
  
 
  _elecEtNMinus1Histo = new TH1F("_elecEtNMinus1Histo", "Elec Et NMinus1", 20, 0, 100);   
  
  _nBTagNMinus1Histo = new TH1F("_nBTagNMinus1Histo", "N B-tags NMinus1", 10, 0, 10);   
  _elecMuonCosDPhiNMinus1Histo = new TH1F("_elecMuonCosDPhiNMinus1Histo", "Elec Muon Cos #Delta #phi NMinus1", 200, -1, -0.8);
  _elecMuonPZetaNMinus1Histo = new TH1F("_elecMuonPZetaNMinus1Histo", "Elec Muon PZeta NMinus1", 300, -150, 150);
  _elecMuonPZetaVisNMinus1Histo = new TH1F("_elecMuonPZetaVisNMinus1Histo", "Elec Muon PZeta Vis NMinus1", 100, 0, 100);
  
  _jetSumEtNMinus1Histo = new TH1F("_jetSumEtNMinus1Histo", "Jet Sum Et NMinus1", 100, 0, 1000);
  _leptonMetCosDPhiNMinus1Histo = new TH1F("_leptonMetCosDPhiNMinus1Histo", "Lepton-Met Cos #Delta #phi NMinus1", 100, -1, 1);
  _leptonDeltaPtCosDPhiNMinus1Histo = new TH1F("_leptonDeltaPtCosDPhiNMinus1Histo", "Lepton-#Delta p_{T} Cos #Delta #phi NMinus1", 100, -1, 1);
  
  _leptonMetCosDPhiVsLeptonDeltaPtCosDPhiNMinus1Histo = new TH2F("_leptonMetCosDPhiVsLeptonDeltaPtCosDPhiNMinus1Histo", "Lepton-Met Cos #Delta #phi vs Lepton-#Delta p_{T} Cos #Delta #phi NMinus1", 100, -1, 1, 100, -1, 1);
  _elecMuonCosDPhiVsLeptonMetCosDPhiNMinus1Histo = new TH2F(" _elecMuonCosDPhiVsLeptonMetCosDPhiNMinus1Histo", "Elec Muon Cos#Delta#phi vs Lepton-Met Cos #Delta#phi", 200, -1, 1, 200, -1, 1);
  
  _elecMuonMetMassNMinus1Histo = new TH1F("_elecMuonMetMassNMinus1Histo", "Elec Muon MEt Mass", 100, 0, 1500);
  _elecMuonDeltaPtMassNMinus1Histo = new TH1F("_elecMuonDeltaPtMassNMinus1Histo", "Elec Muon #Delta p_{T} Mass", 100, 0, 1500);
  
  Float_t xbin[] = {0,50,100,150,200,250,300,400,600,900,1500};
  _elecMuonMetMassNMinus1RebinHisto = new TH1F("_elecMuonMetMassNMinus1RebinHisto", "Elec Muon MEt Mass", 10, xbin);
  _elecMuonDeltaPtMassNMinus1RebinHisto = new TH1F("_elecMuonDeltaPtMassNMinus1RebinHisto", "Elec Muon Delta Pt Mass", 10, xbin);
  
  _nJetsNMinus1Histo = new TH1F("_nJetsNMinus1Histo", "NJets", 10, 0, 10);
  
  
  _elecMuonMassNMinus1Histo = new TH1F("_elecMuonMassNMinus1Histo", "Elec Muon Mass NMinus1", 1000, 0, 2000);
  _muonMetMtNMinus1Histo = new TH1F("_muonMetMtNMinus1Histo", "Muon Met NMinus1", 100, 0, 200);
  
  _elecMuonChargeNMinus1Histo = new TH1F("_elecMuonChargeNMinus1Histo", "ElecCharge * MuonCharge", 3, -1, 2);
  _MEtNMinus1Histo = new TH1F("_MEtNMinus1Histo", "MEt NMinus1", 150, 0, 150);
  
  _threeLeptonNMinus1Histo = new TH1I("_threeLeptonNMinus1Histo", "Three Leptons with two from same vertex?", 2, 0, 2);
  _threeLeptonMassNMinus1Histo = new TH1F("_threeLeptonMassNMinus1Histo", "Mass of tri-lepton System where two like leptons come from same vertex", 100, 0, 500);
  
  
  _muonDeltaPtOverPtNMinus1Histo = new TH1F("_muonDeltaPtOverPtNMinus1Histo", "#Delta p_{T}/p_{T} (#mu)", 100, 0., 1.);
  _elecDeltaPtOverPtNMinus1Histo = new TH1F("_elecDeltaPtOverPtNMinus1Histo", "#Delta p_{T}/p_{T} (e)", 100, 0., 1.);
  
  _muonDeltaPtOverPtVsPtNMinus1Histo = new TH2F("_muonDeltaPtOverPtVsPtNMinus1Histo", "#Delta p_{T}/p_{T} (#mu) vs p_{T} (#mu)", 100, 0., 1000., 100, 0., 1.);
  _elecDeltaPtOverPtVsPtNMinus1Histo = new TH2F("_elecDeltaPtOverPtVsPtNMinus1Histo", "#Delta p_{T}/p_{T} (e) vs p_{T} (e)", 100, 0., 500., 100, 0., 1.);
  
  
  
  _nBTagCR2Histo = new TH1F("_nBTagCR2Histo", "N B-tags CR2", 10, 0, 10); 
  
  _PZetaVsPZetaVisHisto = new TH2F("_PZetaVsPZetaVisHisto", "PZeta vs. PZeta-Visible", 60, 0, 300, 120, -200, 400);  
  _nBTagVsJetSumEtHisto = new TH2F("_nBTagVsJetSumEtHisto", "n B-tags vs. Jet Sum Et", 100, 0, 1000, 10, 0, 10);
  _nBTagVsMEtHisto = new TH2F("_nBTagVsMEtHisto", "n B-tags vs. MEt", 40, 0, 200, 10, 0, 10);
  _nBTagVsMuonMetMtHisto = new TH2F("_nBTagVsMuonMetMtHisto", "n B-tags vs. Mt(#mu, MEt)", 40, 0, 200, 10, 0, 10);
  _nBTagVsNJetsHisto = new TH2F("_nBTagVsNJetsHisto", "n B-tags vs. n Jets", 10, 0, 10, 10, 0, 10);

}

void EMuAnalysis::writeOutFile(){
  outRootFile->Write();        
  outRootFile->Close();    
} 

void EMuAnalysis::GetPUWeights(){
  //fisrt get the data distribution
  TFile* dataPUFile = new TFile("/usr/users/ajohnson/CMS/elecTauTauAnalysis/puRootFiles/run2012AB_190456-196531.root");
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
