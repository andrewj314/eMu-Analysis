//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  8 14:30:32 2013 by ROOT version 5.32/00
// from TTree HMTTree/EMuTree
// found on file: eMuAnalysis_10_1_6jf.root
//////////////////////////////////////////////////////////

#ifndef EMuAnalysis_h
#define EMuAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2.h>
#include <TH1.h>
#include <math.h>
#include <map>
#include <sstream>
#include <iostream>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include "puUtils.h"
#include <TAxis.h>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class EMuAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Long64_t        runNumber;
   Long64_t        eventNumber;
   Int_t           InTimePU;
   Int_t           OutTimePlusPU;
   Int_t           OutTimeMinusPU;
   Bool_t          passedHLTMu5Elec17;
   Bool_t          passedHLTMu11Elec8;
   Bool_t          passedHLTMu8Elec17;
   Bool_t          passedHLTMu17Elec8;
   Bool_t          passedHLTIsoMu24;
   std::vector<double>  *PDFWeights;
   Double_t        ISRGluonWeight;
   Double_t        ISRGammaWeight;
   Double_t        FSRWeight;
   Float_t         mEt;
   Float_t	   mEtPx;
   Float_t	   mEtPy;
   Float_t	   mEtPz;
   Float_t	   beamSpot_x0;
   Float_t	   beamSpot_y0;
   Float_t	   beamSpot_z0;
   Double_t	   primaryVertex_x;
   Double_t	   primaryVertex_y;
   Double_t	   primaryVertex_z;
   Double_t	   primaryVertex_xError;
   Double_t	   primaryVertex_yError;
   Double_t	   primaryVertex_zError;
   
   Double_t	   rho;
   UInt_t          nVtx;
   Bool_t	   eMuMuOneVertex;
   Bool_t	   muEEOneVertex;
   std::vector<float> *elecElecMuonMass;
   std::vector<float> *muonMuonElecMass;
   std::vector<float> *muonMuonMass;
   std::vector<float> *elecElecMass;
   std::vector<float> *tauEnergy;
   std::vector<float> *zPrimeEnergy;
   std::vector<float> *zPrimeVtx_x;
   std::vector<float> *zPrimeVtx_y;
   std::vector<float> *zPrimeVtx_z;
   
   std::vector<float> *threeLeptonMuonPt;
   std::vector<float> *threeLeptonOtherMuonPt;
   std::vector<float> *threeLeptonMuonEta;
   std::vector<float> *threeLeptonOtherMuonEta;
   std::vector<int> *threeLeptonIsGlobalMuon;
   std::vector<int> *threeLeptonOtherMuonIsGlobalMuon;
   std::vector<float> *threeLeptonMuonDxyVtx;
   std::vector<float> *threeLeptonOtherMuonDxy;
   std::vector<float> *threeLeptonMuonDzVtx;
   std::vector<float> *threeLeptonOtherMuonDzVtx;
   std::vector<int> *threeLeptonMuonChamberHits;
   std::vector<int> *threeLeptonOtherMuonChamberHits;
   std::vector<int> *threeLeptonMuonPixelHits;
   std::vector<int> *threeLeptonOtherMuonPixelHits;
   std::vector<int> *threeLeptonMuonMatchedStations;
   std::vector<int> *threeLeptonOtherMuonMatchedStations;
   std::vector<int> *threeLeptonMuonTrkLayersWithHits;
   std::vector<int> *threeLeptonOtherMuonTrkLayersWithHits;
   std::vector<float> *threeLeptonMuonTrkPt;
   std::vector<float> *threeLeptonOtherMuonTrkPt;
   std::vector<float> *threeLeptonMuonTrkPtError;
   std::vector<float> *threeLeptonOtherMuonTrkPtError;
   std::vector<float> *threeLeptonMuonPFIsoDR04SumChargedHadronPt;
   std::vector<float> *threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt;
   std::vector<float> *threeLeptonMuonPFIsoDR04SumCNeutralHadronPt;
   std::vector<float> *threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt;	    
   std::vector<float> *threeLeptonMuonPFIsoDR04SumCPhotonPt;
   std::vector<float> *threeLeptonOtherMuonPFIsoDR04SumPhotonPt;	    
   std::vector<float> *threeLeptonMuonPFIsoDR04SumPUPt;
   std::vector<float> *threeLeptonOtherMuonPFIsoDR04SumPUPt;
   std::vector<float> *threeLeptonMuMuDCAOverError;

   std::vector<float> *threeLeptonElecEt;
   std::vector<float> *threeLeptonOtherElecEt;
   std::vector<float> *threeLeptonElecSCEta;
   std::vector<float> *threeLeptonOtherElecSCEta;
   std::vector<int> *threeLeptonElecEcalDriven;
   std::vector<int> *threeLeptonOtherElecEcalDriven;
   std::vector<int> *threeLeptonElecDEtaIn;
   std::vector<int> *threeLeptonOtherElecDEtaIn;
   std::vector<int> *threeLeptonElecDPhiIn;
   std::vector<int> *threeLeptonOtherElecDPhiIn;
   std::vector<int> *threeLeptonElecHadem;
   std::vector<int> *threeLeptonOtherElecHadem;
   std::vector<int> *threeLeptonElecSigmaIEtaIEta;
   std::vector<int> *threeLeptonOtherElecSigmaIEtaIEta;
   std::vector<int> *threeLeptonElec2by5Over5By5;
   std::vector<int> *threeLeptonOtherElec2by5Over5By5;
   std::vector<int> *threeLeptonElecEcalHad1Iso;
   std::vector<int> *threeLeptonOtherElecEcalHad1Iso;
   std::vector<int> *threeLeptonElecTrkIso;
   std::vector<int> *threeLeptonOtherElecTrkIso;
   std::vector<int> *threeLeptonElecMissingHits;
   std::vector<int> *threeLeptonOtherElecMissingHits;
   std::vector<int> *threeLeptonElecDxy;
   std::vector<int> *threeLeptonOtherElecDxy;
   std::vector<float> *elecElecDCAOverError;   
   
   
   
   
   Float_t         jetBtagHiDiscByTrkCntHiEff;
   Float_t         jetBtagHiDiscByTrkCntHiPurity;
   Float_t         jetBtagHiDiscBySimpleSecVtxHiEff;
   Float_t         jetBtagHiDiscBySimpleSecVtxHiPurity;
   Float_t         jetBtagHiDiscByCombSecVtx;
   std::vector<float>   *jetPt;
   std::vector<float>   *jetEt;
   std::vector<float>   *jetE;
   std::vector<float>   *jetEta;
   std::vector<float>   *jetPhi;
   std::vector<unsigned int> *muonMatched;
   std::vector<unsigned int> *muonMotherId;
   std::vector<float>   *muonPt;
   std::vector<float>   *muonPx;
   std::vector<float>   *muonPy;
   std::vector<float>   *muonPz;
   std::vector<float>   *muonE;
   std::vector<float>   *muonEta;
   std::vector<float>   *muonPhi;
   std::vector<float>   *muonCharge;
   std::vector<double>  *muonDxyVtx;
   std::vector<double>  *muonDxyBS;
   std::vector<double>  *muonDxyBSError;
   std::vector<double>  *muonDxyError;
   std::vector<double>  *muonDzVtx;
   std::vector<double>  *muonDzBS;
   std::vector<double>  *muonDzError;
   std::vector<double>  *muonVx;
   std::vector<double>  *muonVy;
   std::vector<double>  *muonVz;
   std::vector<double>  *muonVertex_x0;
   std::vector<double>  *muonVertex_y0;
   std::vector<double>  *muonVertex_z0;
   std::vector<double>  *muonDecayLength_x;
   std::vector<double>  *muonDecayLength_y;
   std::vector<double>  *muonDecayLength_z;
   
   std::vector<unsigned int> *muonValidHits;
   std::vector<double>  *muonNormChiSqrd;
   std::vector<double>  *muonTrkChiSqrd;
   std::vector<int>     *muonTrkNdof;
   std::vector<int>     *muonChambersValidHits;
   std::vector<int>     *muonMatchedStations;
   std::vector<int>     *muonPixelHits;
   std::vector<int>     *muonTrkLayersWithHits;
   std::vector<float>   *muonPFIsoDR03SumChargedHadronPt;
   std::vector<float>   *muonPFIsoDR03SumChargedParticlePt;
   std::vector<float>   *muonPFIsoDR03SumNeutralHadronPt;
   std::vector<float>   *muonPFIsoDR03SumPhotonPt;
   std::vector<float>   *muonPFIsoDR03SumPUPt;
   std::vector<float>   *muonPFIsoDR04SumChargedHadronPt;
   std::vector<float>   *muonPFIsoDR04SumChargedParticlePt;
   std::vector<float>   *muonPFIsoDR04SumNeutralHadronPt;
   std::vector<float>   *muonPFIsoDR04SumPhotonPt;
   std::vector<float>   *muonPFIsoDR04SumPUPt;
   std::vector<float>   *muonCaloComp;
   std::vector<int>     *muonIsGlobalMuon;
   std::vector<int>     *muonisTrackerMuon;
   std::vector<int>     *muonIsGlobalPromptTight;
   std::vector<int>     *muonIsTMLastStationLoose;
   std::vector<int>     *muonIsTMLastStationTight;
   std::vector<int>     *muonIsTM2DCompatibilityLoose;
   std::vector<int>     *muonIsTM2DCompatibilityTight;
   std::vector<float>   *muonPionVeto;
   
   std::vector<int>     *muonIsHighPtMuon;
   std::vector<float>   *muonTrkPt;
   std::vector<float>   *muonTrkPtError;
   
   std::vector<int>     *elecMatched;
   std::vector<unsigned int> *eMotherId;
   std::vector<float>   *elecE;
   std::vector<float>   *elecEt;
   std::vector<float>   *elecPt;
   std::vector<float>   *elecPx;
   std::vector<float>   *elecPy;
   std::vector<float>   *elecPz;
   std::vector<double>   *elecVertex_x0;
   std::vector<double>   *elecVertex_y0;
   std::vector<double>   *elecVertex_z0;
   std::vector<double>   *elecVx;
   std::vector<double>   *elecVy;
   std::vector<double>   *elecVz;
   std::vector<float>   *elecDecayLength_x;
   std::vector<float>   *elecDecayLength_y;
   std::vector<float>   *elecDecayLength_z;
   std::vector<float>   *elecCharge;
   std::vector<float>   *elecEta;
   std::vector<float>   *elecSCEta;
   std::vector<float>   *elecPhi;
   std::vector<float>   *elecSigmaEtaEta;
   std::vector<float>   *elecSigmaIEtaIEta;
   std::vector<float>   *elecEOverP;
   std::vector<float>   *elecHOverEm;
   std::vector<float>   *elecDeltaPhiIn;
   std::vector<float>   *elecDeltaEtaIn;
   std::vector<float>   *elecFBrem;
   std::vector<float>   *elecEOverInvP;
   std::vector<double>  *elecTrkChiSqrd;
   std::vector<int>     *elecTrkNdof;
   std::vector<float>   *elecSCE1x5;
   std::vector<float>   *elecSCE2x5;
   std::vector<float>   *elecSCE5x5;
   std::vector<float>   *elecSCE1x5Over5x5;
   std::vector<float>   *elecSCE2x5MaxOver5x5;
   std::vector<float>   *elecTrkIso;
   std::vector<float>   *elecEcalIso;
   std::vector<float>   *elecHcalIso;
   std::vector<float>   *elecHcalIsoDepth1;
   std::vector<float>   *elecHcalIsoDepth2;
   std::vector<float>   *elecEcalHcalIsoDepth1;
   std::vector<float>   *gsfEcalIso;
   std::vector<float>   *gsfHcalIso;
   std::vector<float>   *gsfEt;
   std::vector<float>   *gsfEnergy;
   std::vector<float>   *gsfCaloEnergy;
   std::vector<float>   *elecPFIsoDR03SumChargedHadronPt;
   std::vector<float>   *elecPFIsoDR03SumNeutralHadronPt;
   std::vector<float>   *elecPFIsoDR03SumPhotonPt;
   std::vector<float>   *elecGSFIp;
   std::vector<float>   *elecGSFIpError;
   std::vector<float>   *elecGSFBSIp;
   std::vector<float>   *elecGSFBSIpError;   
   std::vector<float>   *elecCTFIp;
   std::vector<float>   *elecCTFIpError;
   std::vector<float>   *elecMuonDistanceOfClosestApproach;
   std::vector<float>   *elecMuonDCACrossingPoint;
   std::vector<float>   *elecMuonTwoTrkMinDist;
   std::vector<float>   *elecMuonTwoTrkCrossingPoint;
   std::vector<float>   *eMuDis3D;
   std::vector<float>   *eMuDis2D;
   std::vector<float>   *elecMuonDCA3DError;
   std::vector<float>   *elecMuonDCA2DError;
   std::vector<int>     *elecMissingHits;
   std::vector<int>     *elecClass;
   std::vector<bool>    *elecIsEcalDriven;
   std::vector<bool>    *elecEcalDriven;
   std::vector<bool>    *elecIsTrkDriven;
   std::vector<bool>    *elecInEE;
   std::vector<bool>    *elecInEB;
   std::vector<bool>    *elecInGap;
   std::vector<int>     *heepPassedEt;
   std::vector<int>     *heepPassedPt;
   std::vector<int>     *heepPassedDetEta;
   std::vector<int>     *heepPassedCrack;
   std::vector<int>     *heepPassedDEtaIn;
   std::vector<int>     *heepPassedDPhiIn;
   std::vector<int>     *heepPassedHadem;
   std::vector<int>     *heepPassedSigmaIEtaIEta;
   std::vector<int>     *heepPassed2by5Over5By5;
   std::vector<int>     *heepPassedEcalHad1Iso;
   std::vector<int>     *heepPassedHad2Iso;
   std::vector<int>     *heepPassedTrkIso;
   std::vector<int>     *heepPassedEcalDriven;
   std::vector<int>	*heepPassedMissingHits;
   std::vector<int>	*heepPassedElecDxy;
   std::vector<int>     *heepPassedAllCuts;
   std::vector<float>   *deltaPtx;
   std::vector<float>   *deltaPty;
   std::vector<float>   *deltaPtz;
   
   std::vector<int>     *elecMVAStatus;
   std::vector<float>   *elecMVAOut;
   std::vector<float>   *elecMuononMass;
   std::vector<float>   *elecMuonCosDPhi;
   std::vector<float>   *elecMuonDeltaR;
   std::vector<float>   *elecMuonMetMass;
   std::vector<float>   *elecMuonDeltaPtMass;
   std::vector<float>   *elecMuonCollMass;
   std::vector<float>   *elecMEtMt;
   std::vector<float>   *muonMEtMt;
   std::vector<float>   *elecMuonPZetaVis;
   std::vector<float>   *elecMuonPZeta;
   std::vector<unsigned int> *nBtagsHiEffTrkCnt;
   std::vector<unsigned int> *nBtagsHiPurityTrkCnt;
   std::vector<unsigned int> *nBTagsHiEffSimpleSecVtx;
   std::vector<unsigned int> *nBTagsHiPuritySimpleSecVtx;
   std::vector<unsigned int> *nBTagsCombSecVtxLWP;
   std::vector<float>   *jetSumEt;
   std::vector<float>   *jetMETSumEt;
   std::vector<float>   *extraTrkPtSum;
   std::vector<float>   *leptonMetCosDphi;
   std::vector<float>   *leptonDeltaPtCosDphi;
   std::vector<unsigned int> *nJets;
   
   std::vector<float> *BSx0;
   std::vector<float> *BSy0;
   std::vector<float> *BSz0;

   std::vector<float> *genTauVtx_x;
   std::vector<float> *genTauVtx_y;
   std::vector<float> *genTauVtx_z;
   std::vector<float> *otherGenTauVtx_x;
   std::vector<float> *otherGenTauVtx_y;
   std::vector<float> *otherGenTauVtx_z;   
   std::vector<float> *genTauPx;   
   std::vector<float> *genTauPy;   
   std::vector<float> *genTauPz;   
   std::vector<float> *otherGenTauPx;   
   std::vector<float> *otherGenTauPy;   
   std::vector<float> *otherGenTauPz;   
   std::vector<float> *genTauCharge;   
   std::vector<float> *otherGenTauCharge;   


   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_InTimePU;   //!
   TBranch        *b_OutTimePlusPU;   //!
   TBranch        *b_OutTimeMinusPU;   //!
   TBranch        *b_passedHLTMu5Elec17;   //!
   TBranch        *b_passedHLTMu11Elec8;   //!
   TBranch        *b_passedHLTMu8Elec17;   //!
   TBranch        *b_passedHLTMu17Elec8;   //!
   TBranch        *b_passedHLTIsoMu24;   //!
   TBranch        *b_PDFWeights;   //!
   TBranch        *b_ISRGluonWeight;   //!
   TBranch        *b_ISRGammaWeight;   //!
   TBranch        *b_FSRWeight;   //!
   TBranch        *b_mEt;   //!
   TBranch	  *b_mEtPx;   //!
   TBranch	  *b_mEtPy;   //!
   TBranch	  *b_mEtPz;   //!
   TBranch	  *b_beamSpot_x0;   //!
   TBranch	  *b_beamSpot_y0;   
   TBranch	  *b_beamSpot_z0;
   TBranch	  *b_rho;   //!
   TBranch        *b_nVtx;   //!
   TBranch	  *b_eMuMuOneVertex;   //!
   TBranch	  *b_muEEOneVertex;   //!
   TBranch	  *b_elecElecMuonMass;   //!
   TBranch	  *b_muonMuonElecMass;   //!
   
   TBranch        *b_muonMuonMass;   //!
   TBranch        *b_elecElecMass;   //!
   TBranch	  *b_tauEnergy;   //!
   TBranch	  *b_zPrimeEnergy;   //!
   TBranch	  *b_zPrimeVtx_x;   //!
   TBranch	  *b_zPrimeVtx_y;   //!
   TBranch	  *b_zPrimeVtx_z;   //!
   
   TBranch        *b_threeLeptonMuonPt;   //!
   TBranch        *b_threeLeptonOtherMuonPt;   //!
   TBranch        *b_threeLeptonMuonEta;   //!
   TBranch        *b_threeLeptonOtherMuonEta;   //!
   TBranch        *b_threeLeptonIsGlobalMuon;   //!
   TBranch        *b_threeLeptonOtherMuonIsGlobalMuon;   //!
   TBranch        *b_threeLeptonMuonDxyVtx;   //!
   TBranch        *b_threeLeptonOtherMuonDxy;   //!
   TBranch        *b_threeLeptonMuonDzVtx;   //!
   TBranch        *b_threeLeptonOtherMuonDzVtx;   //!
   TBranch        *b_threeLeptonMuonChamberHits;   //!
   TBranch        *b_threeLeptonOtherMuonChamberHits;   //!
   TBranch        *b_threeLeptonMuonPixelHits;   //!
   TBranch        *b_threeLeptonOtherMuonPixelHits;   //!
   TBranch        *b_threeLeptonMuonMatchedStations;   //!
   TBranch        *b_threeLeptonOtherMuonMatchedStations;   //!
   TBranch        *b_threeLeptonMuonTrkLayersWithHits;   //!
   TBranch        *b_threeLeptonOtherMuonTrkLayersWithHits;   //!
   TBranch        *b_threeLeptonMuonTrkPt;   //!
   TBranch        *b_threeLeptonOtherMuonTrkPt;   //!
   TBranch        *b_threeLeptonMuonTrkPtError;   //!
   TBranch        *b_threeLeptonOtherMuonTrkPtError;   //!
   TBranch        *b_threeLeptonMuonPFIsoDR04SumChargedHadronPt;   //!
   TBranch        *b_threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt;   //!
   TBranch        *b_threeLeptonMuonPFIsoDR04SumCNeutralHadronPt;   //!
   TBranch        *b_threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt;   //!	    
   TBranch        *b_threeLeptonMuonPFIsoDR04SumCPhotonPt;   //!
   TBranch        *b_threeLeptonOtherMuonPFIsoDR04SumPhotonPt;   //!	    
   TBranch        *b_threeLeptonMuonPFIsoDR04SumPUPt;   //!
   TBranch        *b_threeLeptonOtherMuonPFIsoDR04SumPUPt;   //!
   TBranch        *b_threeLeptonMuMuDCAOverError;   //!

   TBranch        *b_threeLeptonElecEt;   //!
   TBranch        *b_threeLeptonOtherElecEt;   //!
   TBranch        *b_threeLeptonElecSCEta;   //!
   TBranch        *b_threeLeptonOtherElecSCEta;   //!
   TBranch        *b_threeLeptonElecEcalDriven;   //!
   TBranch        *b_threeLeptonOtherElecEcalDriven;   //!
   TBranch        *b_threeLeptonElecDEtaIn;   //!
   TBranch        *b_threeLeptonOtherElecDEtaIn;   //!
   TBranch        *b_threeLeptonElecDPhiIn;   //!
   TBranch        *b_threeLeptonOtherElecDPhiIn;   //!
   TBranch        *b_threeLeptonElecHadem;   //!
   TBranch        *b_threeLeptonOtherElecHadem;   //!
   TBranch        *b_threeLeptonElecSigmaIEtaIEta;   //!
   TBranch        *b_threeLeptonOtherElecSigmaIEtaIEta;   //!
   TBranch        *b_threeLeptonElec2by5Over5By5;   //!
   TBranch        *b_threeLeptonOtherElec2by5Over5By5;   //!
   TBranch        *b_threeLeptonElecEcalHad1Iso;   //!
   TBranch        *b_threeLeptonOtherElecEcalHad1Iso;   //!
   TBranch        *b_threeLeptonElecTrkIso;   //!
   TBranch        *b_threeLeptonOtherElecTrkIso;   //!
   TBranch        *b_threeLeptonElecMissingHits;   //!
   TBranch        *b_threeLeptonOtherElecMissingHits;   //!
   TBranch        *b_threeLeptonElecDxy;   //!
   TBranch        *b_threeLeptonOtherElecDxy;   //!
   TBranch        *b_elecElecDCAOverError;   //!      
   
   TBranch        *b_jetBtagHiDiscByTrkCntHiEff;   //!
   TBranch        *b_jetBtagHiDiscByTrkCntHiPurity;   //!
   TBranch        *b_jetBtagHiDiscBySimpleSecVtxHiEff;   //!
   TBranch        *b_jetBtagHiDiscBySimpleSecVtxHiPurity;   //!
   TBranch        *b_jetBtagHiDiscByCombSecVtx;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_muonMatched;   //!
   TBranch        *b_muonMotherId;   //!
   TBranch        *b_muonPt;   //!
   TBranch        *b_muonPx;   //!
   TBranch        *b_muonPy;   //!
   TBranch        *b_muonPz;   //!
   TBranch        *b_muonE;   //!
   TBranch        *b_muonEta;   //!
   TBranch        *b_muonPhi;   //!
   TBranch        *b_muonCharge;   //!
   TBranch        *b_muonDxyVtx;   //!
   TBranch        *b_muonDxyBS;   //!
   TBranch 	  *b_muonDxyBSError; //!
   TBranch        *b_muonDxyError;   //!
   TBranch        *b_muonDzVtx;   //!
   TBranch        *b_muonDzBS;   //!
   TBranch        *b_muonDzError;   //!
   TBranch        *b_muonVx;   //!
   TBranch        *b_muonVy;   //!
   TBranch        *b_muonVz;   //!
   TBranch	  *b_muonVertex_x0;   
   TBranch	  *b_muonVertex_y0;   
   TBranch	  *b_muonVertex_z0;   
   TBranch	  *b_muonDecayLength_x;
   TBranch	  *b_muonDecayLength_y;
   TBranch	  *b_muonDecayLength_z;
   
   TBranch        *b_muonValidHits;   //!
   TBranch        *b_muonNormChiSqrd;   //!
   TBranch	  *b_muonTrkChiSqrd;   //!
   TBranch	  *b_muonTrkNdof;   //!
   TBranch        *b_muonChambersValidHits;   //!
   TBranch        *b_muonMatchedStations;   //!
   TBranch        *b_muonPixelHits;   //!
   TBranch        *b_muonTrkLayersWithHits;   //!
   TBranch        *b_muonPFIsoDR03SumChargedHadronPt;   //!
   TBranch        *b_muonPFIsoDR03SumChargedParticlePt;   //!
   TBranch        *b_muonPFIsoDR03SumNeutralHadronPt;   //!
   TBranch        *b_muonPFIsoDR03SumPhotonPt;   //!
   TBranch        *b_muonPFIsoDR03SumPUPt;   //!
   TBranch        *b_muonPFIsoDR04SumChargedHadronPt;   //!
   TBranch        *b_muonPFIsoDR04SumChargedParticlePt;   //!
   TBranch        *b_muonPFIsoDR04SumNeutralHadronPt;   //!
   TBranch        *b_muonPFIsoDR04SumPhotonPt;   //!
   TBranch        *b_muonPFIsoDR04SumPUPt;   //!
   TBranch        *b_muonCaloComp;   //!
   TBranch        *b_muonIsGlobalMuon;   //!
   TBranch        *b_muonisTrackerMuon;   //!
   TBranch        *b_muonIsGlobalPromptTight;   //!
   TBranch        *b_muonIsTMLastStationLoose;   //!
   TBranch        *b_muonIsTMLastStationTight;   //!
   TBranch        *b_muonIsTM2DCompatibilityLoose;   //!
   TBranch        *b_muonIsTM2DCompatibilityTight;   //!
   TBranch        *b_muonPionVeto;   //!
   
   TBranch	  *b_muonIsHighPtMuon;   //!
   TBranch 	  *b_muonTrkPt;   //!
   TBranch	  *b_muonTrkPtError;   //!
   
   TBranch        *b_elecMatched;   //!
   TBranch        *b_eMotherId;   //!
   TBranch        *b_elecE;   //!
   TBranch        *b_elecEt;   //!
   TBranch        *b_elecPt;   //!
   TBranch        *b_elecPx;   //!
   TBranch        *b_elecPy;   //!
   TBranch        *b_elecPz;   //!
   TBranch	  *b_elecVertex_x0;   //!
   TBranch	  *b_elecVertex_y0;   //!
   TBranch	  *b_elecVertex_z0;   //!
   TBranch	  *b_elecDecayLength_x;   //!
   TBranch	  *b_elecDecayLength_y;   //!
   TBranch	  *b_elecDecayLength_z;   //!
   TBranch	  *b_elecVx;   //!
   TBranch	  *b_elecVy;   //!
   TBranch	  *b_elecVz;   //!
   
   TBranch        *b_elecCharge;   //!
   TBranch        *b_elecEta;   //!
   TBranch	  *b_elecSCEta;   //!
   TBranch        *b_elecPhi;   //!
   TBranch        *b_elecSigmaEtaEta;   //!
   TBranch        *b_elecSigmaIEtaIEta;   //!
   TBranch        *b_elecEOverP;   //!
   TBranch        *b_elecHOverEm;   //!
   TBranch        *b_elecDeltaPhiIn;   //!
   TBranch        *b_elecDeltaEtaIn;   //!
   TBranch        *b_elecFBrem;   //!
   TBranch	  *b_elecTrkChiSqrd;   //!
   TBranch	  *b_elecTrkNdof;   //!
   TBranch        *b_elecEOverInvP;   //!
   TBranch        *b_elecSCE1x5;   //!
   TBranch        *b_elecSCE2x5;   //!
   TBranch        *b_elecSCE5x5;   //!
   TBranch        *b_elecSCE1x5Over5x5;   //!
   TBranch        *b_elecSCE2x5MaxOver5x5;   //!
   TBranch        *b_elecTrkIso;   //!
   TBranch        *b_elecEcalIso;   //!
   TBranch        *b_elecHcalIso;   //!
   TBranch        *b_elecHcalIsoDepth1;   //!
   TBranch        *b_elecHcalIsoDepth2;   //!
   TBranch        *b_elecEcalHcalIsoDepth1;   //!
   TBranch	  *b_gsfEcalIso;   //!
   TBranch	  *b_gsfHcalIso;   //!
   TBranch	  *b_gsfEt;   //!
   TBranch	  *b_gsfEnergy;   //!
   TBranch	  *b_gsfCaloEnergy;   //!
   TBranch        *b_elecPFIsoDR03SumChargedHadronPt;   //!
   TBranch        *b_elecPFIsoDR03SumNeutralHadronPt;   //!
   TBranch        *b_elecPFIsoDR03SumPhotonPt;   //!
   TBranch        *b_elecGSFIp;   //!
   TBranch        *b_elecGSFIpError;   //!
   TBranch        *b_elecGSFBSIp;   //!
   TBranch        *b_elecGSFBSIpError;   //!   
   TBranch        *b_elecCTFIp;   //!
   TBranch        *b_elecCTFIpError;   //!
   TBranch	  *b_elecMuonDistanceOfClosestApproach;   //!
   TBranch	  *b_elecMuonDCACrossingPoint;   //!
   TBranch	  *b_elecMuonTwoTrkMinDist;   //!
   TBranch	  *b_elecMuonTwoTrkCrossingPoint;   //!
   TBranch	  *b_eMuDis3D;   //!
   TBranch	  *b_eMuDis2D;   //!
   TBranch	  *b_elecMuonDCA3DError;   //!
   TBranch	  *b_elecMuonDCA2DError;   //!
   TBranch        *b_elecMissingHits;   //!
   TBranch        *b_elecClass;   //!
   TBranch        *b_elecIsEcalDriven;   //!
   TBranch	  *b_elecEcalDriven;   //!
   TBranch        *b_elecIsTrkDriven;   //!
   TBranch        *b_elecInEE;   //!
   TBranch        *b_elecInEB;   //!
   TBranch        *b_elecInGap;   //!
   TBranch        *b_heepPassedEt;   //!
   TBranch        *b_heepPassedPt;   //!
   TBranch        *b_heepPassedDetEta;   //!
   TBranch        *b_heepPassedCrack;   //!
   TBranch        *b_heepPassedDEtaIn;   //!
   TBranch        *b_heepPassedDPhiIn;   //!
   TBranch        *b_heepPassedHadem;   //!
   TBranch        *b_heepPassedSigmaIEtaIEta;   //!
   TBranch        *b_heepPassed2by5Over5By5;   //!
   TBranch        *b_heepPassedEcalHad1Iso;   //!
   TBranch        *b_heepPassedHad2Iso;   //!
   TBranch        *b_heepPassedTrkIso;   //!
   TBranch        *b_heepPassedEcalDriven;   //!
   TBranch	  *b_heepPassedMissingHits;   //!
   TBranch	  *b_heepPassedElecDxy;   //!
   TBranch	  *b_heepPassedAllCuts;   //!
   TBranch	  *b_deltaPtx;   //!
   TBranch	  *b_deltaPty;   //!
   TBranch	  *b_deltaPtz;   //!
   TBranch        *b_elecMVAStatus;   //!
   TBranch        *b_elecMVAOut;   //!
   TBranch        *b_elecMuononMass;   //!
   TBranch        *b_elecMuonCosDPhi;   //!
   TBranch        *b_elecMuonDeltaR;   //!
   TBranch        *b_elecMuonMetMass;   //!
   TBranch	  *b_elecMuonDeltaPtMass;   //!
   TBranch        *b_elecMuonCollMass;   //!
   TBranch        *b_elecMEtMt;   //!
   TBranch        *b_muonMEtMt;   //!
   TBranch        *b_elecMuonPZetaVis;   //!
   TBranch        *b_elecMuonPZeta;   //!
   TBranch        *b_nBtagsHiEffTrkCnt;   //!
   TBranch        *b_nBtagsHiPurityTrkCnt;   //!
   TBranch        *b_nBTagsHiEffSimpleSecVtx;   //!
   TBranch        *b_nBTagsHiPuritySimpleSecVtx;   //!
   TBranch        *b_nBTagsCombSecVtxLWP;   //!
   TBranch        *b_jetSumEt;   //!
   TBranch        *b_jetMETSumEt;   //!
   TBranch        *b_extraTrkPtSum;   //!
   TBranch        *b_leptonMetCosDphi;   //!
   TBranch	  *b_leptonDeltaPtCosDphi;   //!
   TBranch        *b_nJets;   //!
   
   TBranch	  *b_BSx0;   //!
   TBranch	  *b_BSy0;   //!
   TBranch	  *b_BSz0;   //!
   
   TBranch	  *b_primaryVertex_x;   //!
   TBranch	  *b_primaryVertex_y;   //!
   TBranch	  *b_primaryVertex_z;   //!
   TBranch	  *b_primaryVertex_xError;   //!
   TBranch	  *b_primaryVertex_yError;   //!
   TBranch	  *b_primaryVertex_zError;   //!
   TBranch	  *b_genTauVtx_x;   //!
   TBranch	  *b_genTauVtx_y;   //!
   TBranch	  *b_genTauVtx_z;   //!
   TBranch	  *b_otherGenTauVtx_x;   //!
   TBranch	  *b_otherGenTauVtx_y;   //!
   TBranch	  *b_otherGenTauVtx_z;   //!   
   TBranch	  *b_genTauPx;   //!
   TBranch	  *b_genTauPy;   //!
   TBranch	  *b_genTauPz;   //!
   TBranch	  *b_otherGenTauPx;   //!
   TBranch	  *b_otherGenTauPy;   //!
   TBranch	  *b_otherGenTauPz;   //!
   TBranch	  *b_genTauCharge;   //!
   TBranch	  *b_otherGenTauCharge;   //!


   EMuAnalysis(TTree *tree=0);
   virtual ~EMuAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   void SetSignalXSection(double theXSection);
   void SetLumi(double theLumi);
   void SetNEvents(int nEvents);
   void SetSource(std::string theSource);
   
   void SetOutputLogFileName(std::string theLogFileName);
   void SetOutputRootFileName(std::string theRootFileName);
   
   void SetElecAccCuts(double minElecPtCut, double maxElecEtaCut);
   
   void SetElecIdGlobalCuts(int maxElecMissingHitsCut,
  			    double maxElecPFRelIsoCut);
   void SetElecIdBarrelCuts(double maxElecDeltaPhiInEBCut, double maxElecDeltaEtaInEBCut,
     			    double maxElecSigmaIEtaIEtaEBCut, double maxElecHadFracEBCut, 
			    double maxOneOverEMinusOneOverPEBCut, double maxElecDxyBarrelCut);
   void SetElecIdEndCapCuts(double maxElecDeltaPhiInBBCut, double maxElecDeltaEtaInBBCut,
     			    double maxElecSigmaIEtaIEtaBBCut, double maxElecHadFracBBCut, 
			    double maxOneOverEMinusOneOverPBBCut, double maxElecDxyEndCapCut);  
			    
   void SetMuonMotherIdCut(double muonMotherIdCut);			    
   void SetElecMotherIdCut(double eMotherIdCut);			    
			    
   void SetMuonAccCuts(double minMuonPtCut, double maxMuonEtaCut);			     
   
   void SetMuonIdCuts(double maxMuonDxyVtxCut, double maxMuonDzVtxCut,
			      double minMuonIsGlobalMuonCut, double minMuonMatchedCut, double maxMuonNormChiSqrdCut, double minMuonValidHitsCut,
			      double minMuonMatchedStationsCut, double minMuonPixelHitsCut, double minMuonTrkLayersWithHitsCut, double minMuonIsoDR04DBCut,
			      double maxMuonIsoDR04DBCut, double minMuonIsHighPtMuonCut, double maxMuonDptOverPtCut);
			      
   void SetHeepCuts(double minElecMatchedCut, double minHeepPassedEtCut, double minHeepPassedDetEtaCut, double minHeepPassedCrackCut, double minHeepPassedDEtaInCut, 
   			      double minHeepPassedDPhiInCut, double minHeepPassedHademCut, double minHeepPassedSigmaIEtaIEtaCut, double minHeepPassed2by5Over5By5Cut, 
			      double minHeepPassedEcalDrivenCut, double minHeepPassedMissingHitsCut, double minHeepPassedElecDxyCut, double minElecEtCut, double minElecRelIsoCut, double maxElecRelIsoCut, double minElecTrkIsoCut, double maxElecTrkIsoCut, double minHeepPassedEcalHad1IsoCut,
			      double minHeepPassedTrkIsoCut);      
  
  void SetImpactParameterCuts(double minMuonIpOverMuonIpErrorCut, double minElecIpOverElecIpErrorCut, double minCombinedIpOverIpErrorCut, double maxCombinedIpOverIpErrorCut, double minElecMuonDCAOverError2DCut, 
  			      double maxElecMuonDCAOverError2DCut, double minElecMuonDCAOverError3DCut, double maxElecMuonDCAOverError3DCut, double minLifetimeCut);
  
  void SetTopologyCuts(double elecMuonMassCut, double minElecMuonCosDPhiCut, double maxElecMuonCosDPhiCut, double elecMuonDeltaRCut, double minElecMuonMetMassCut,
  		       double maxElecMuonMetMassCut, double elecMuonDeltaPtMassCut, double elecMetMtCut, double minMuonMetMtCut, double maxMuonMetMtCut, double elecMuonPZetaVisCut,
		       double minElecMuonPZetaCut, double maxElecMuonPZetaCut, double pZetaSlope, double minNBtagsHiEffTrkCntCut, double maxNBtagsHiEffTrkCntCut, double minNBtagsHiPurityTrkCntCut,
		       double maxNBtagsHiPurityTrkCntCut, double minNBTagsHiEffSimpleSecVtxCut, double maxNBTagsHiEffSimpleSecVtxCut, double minNBTagsHiPuritySimpleSecVtxCut, 
		       double maxNBTagsHiPuritySimpleSecVtxCut, double minNBTagsCombSecVtxCut, double maxNBTagsCombSecVtxCut, double leptonMetCosDPhiCut, double maxLeptonDeltaPtCosDPhiCut, double MEtCut, double nJetsCut, 
		       double minChargeCut, double maxChargeCut, double jetSumEtCut, double minElecMuonDeltaPtCut, double min2xPLargerLeptonCut, double minElecEOverPCut, double maxElecEOverPCut, double maxThreeLeptonCut);
  
  
  /*
  void SetTopologyCuts(double elecMuonMassCut, double minElecMuonCosDPhiCut, double maxElecMuonCosDPhiCut, double elecMuonDeltaRCut, 
  		       double elecMuonMetMassCut, double elecMetMtCut, double minMuonMetMtCut, double maxMuonMetMtCut, double pZetaCentroidCut, double pZetaVisCentroidCut,
		       double pZetaSemiMajorAxisCut, double pZetaSemiMinorAxisCut, double pZetaEllipseAngleCut, double minNBtagsHiEffTrkCntCut, double maxNBtagsHiEffTrkCntCut, double minNBtagsHiPurityTrkCntCut,
		       double maxNBtagsHiPurityTrkCntCut, double minNBTagsHiEffSimpleSecVtxCut, double maxNBTagsHiEffSimpleSecVtxCut, double minNBTagsHiPuritySimpleSecVtxCut, 
		       double maxNBTagsHiPuritySimpleSecVtxCut, double minNBTagsCombSecVtxCut, double maxNBTagsCombSecVtxCut, double leptonMetCosDPhiCut, double MEtCut, double nJetsCut, 
		       double minChargeCut, double maxChargeCut, double jetSumEtCut);
  		*/       
  
private:
   TFile* outRootFile;


   bool elecExists(unsigned int theIndex);
   bool passedMatched(unsigned int theIndex);
   bool passedElecMinPt(unsigned int theIndex);
   bool passedElecMaxEta(unsigned int theIndex);
   bool passedElecAcc(unsigned int theIndex);
   bool passedElecDeltaPhiIn(unsigned int theIndex);
   bool passedElecDeltaEtaIn(unsigned int theIndex);
   bool passedElecSigmaIEtaIEta(unsigned int theIndex);
   bool passedElecHOverE(unsigned int theIndex);
   bool passedElecEOverInvP(unsigned int theIndex);
   bool passedElecRelIso(unsigned int theIndex);
   bool passedElecTrkIso(unsigned int theIndex);
   bool passedElecMissingHits(unsigned int theIndex);
   //bool passedElecIsoDr(unsigned int theIndex);
   bool passedElecPFRelIso(unsigned int theIndex);
   bool passedElecGlobalId(unsigned int theIndex);
   bool passedElecId(unsigned int theIndex);
   bool elecInBarrel(unsigned int theIndex);
   bool elecInEndCap(unsigned int theIndex);
   bool elecInGaps(unsigned int theIndex);
   bool passedElecMatched(unsigned int theIndex);
   bool passedElecMotherId(unsigned int theIndex);  
   bool passedElec2x5Over5x5(unsigned int theIndex);
   bool passedElecEcalHad1Iso(unsigned int theIndex);
   bool passedElecDxy(unsigned int theIndex); 

   bool passedElecMuonMass(unsigned int theIndex);
   bool passedElecMuonCosDPhi(unsigned int theIndex);
   bool passedElecMuonDeltaR(unsigned int theIndex);
   bool passedElecMuonMetMass(unsigned int theIndex);
   bool passedElecMuonDeltaPtMass(unsigned int theIndex);
   bool passedElecMetMt(unsigned int theIndex);
   bool passedMuonMetMt(unsigned int theIndex);
   bool passedElecMuonPZetaVis(unsigned int theIndex);
   bool passedElecMuonPZeta(unsigned int theIndex);
   bool passedNbTagsHiEffTrkCnt(unsigned int theIndex);
   bool passedNbTagsHiPurityTrkCnt(unsigned int theIndex);
   bool passedNbTagsHiEffSimpleSecVtx(unsigned int theIndex);
   bool passedNbTagsHiPuritySimpleSecVtx(unsigned int theIndex);
   bool passedNbTagsCombSecVtx(unsigned int theIndex);
   bool passedLeptonMetCosDPhi(unsigned int theIndex);
   bool passedLeptonDeltaPtCosDPhi(unsigned int theIndex);
   bool passedMET(unsigned int theIndex);
   bool passedThreeLepton(unsigned int theIndex);
   bool passedNJets(unsigned int theIndex);
   bool passedCharge(unsigned int theIndex);
   bool passedJetSumEt(unsigned int theIndex);
   bool passedElecMuonDeltaPt(unsigned int theIndex);
   bool passed2xPLargerLepton(unsigned int theIndex);
   bool passedElecEOverP(unsigned int theIndex);
   bool passedTopology(unsigned int theIndex);
   
   
   bool passedMuonPt(unsigned int theIndex);
   bool passedMuonE(unsigned int theIndex);
   bool passedMuonEta(unsigned int theIndex);
   bool passedMuonAcc(unsigned int theIndex);
   bool passedMuonDxyVtx(unsigned int theIndex);
   bool passedMuonDzVtx(unsigned int theIndex);
   bool passedMuonIsGlobalMuon(unsigned int theIndex);
   bool passedMuonMatched(unsigned int theIndex);
   bool passedMuonMotherId(unsigned int theIndex);
   bool passedMuonNormChiSqrd(unsigned int theIndex);
   bool passedMuonValidHits(unsigned int theIndex);
   bool passedMuonMatchedStations(unsigned int theIndex);
   bool passedMuonPixelHits(unsigned int theIndex);
   bool passedMuonTrkLayersWithHits(unsigned int theIndex);
   bool passedMuonIsoDR03NoDB(unsigned int theIndex);
   bool passedMuonIsoDR03DB(unsigned int theIndex);
   bool passedMuonIsoDR04NoDB(unsigned int theIndex);
   bool passedMuonIsoDR04DB(unsigned int theIndex);
   
   bool passedMuonIsHighPtMuon(unsigned int theIndex);
   bool passedMuonDptOverPt(unsigned int theIndex);
   
   bool passedAllMuonIDCuts(unsigned int theIndex);
   
   bool passedHeepPassedEt(unsigned int theIndex);
   bool passedHeepPassedCrack(unsigned int theIndex);
   bool passedHeepPassedDEtaIn(unsigned int theIndex);
   bool passedHeepPassedDPhiIn(unsigned int theIndex);
   bool passedHeepPassedHadem(unsigned int theIndex);
   bool passedHeepPassedSigmaIEtaIEta(unsigned int theIndex);
   bool passedHeepPassed2by5Over5By5(unsigned int theIndex);
   bool passedHeepPassedEcalHad1Iso(unsigned int theIndex);
   //bool passedHeepPassedHad2Iso(unsigned int theIndex);   
   bool passedHeepPassedTrkIso(unsigned int theIndex);
   bool passedHeepPassedEcalDriven(unsigned int theIndex);
   bool passedHeepPassedMissingHits(unsigned int theIndex);
   bool passedHeepPassedElecDxy(unsigned int theIndex);
   bool passedElecEt(unsigned int theIndex);
   bool passedHeepId(unsigned int theIndex);
   
   bool passedElecIpOverElecIpError(unsigned int theIndex);
   bool passedMuonIpOverMuonIpError(unsigned int theIndex);
   bool passedCombinedIpOverIpError(unsigned int theIndex);
   bool passedElecMuonDCAOverError2D(unsigned int theIndex);
   bool passedElecMuonDCAOverError3D(unsigned int theIndex);
   bool passedLifetime(unsigned int theIndex);
   
   //CONTROL REGIONS
   
   bool passedQCDCR1(unsigned int theIndex);
   bool passedQCDCR2(unsigned int theIndex);
   bool passedTTJetsCR1(unsigned int theIndex);
   bool passedTTJetsCR2(unsigned int theIndex);
   bool passedWJetsCR1(unsigned int theIndex);
   bool passedWJetsCR2(unsigned int theIndex);
   bool passedDYToTauTauCR1(unsigned int theIndex);
   bool passedDYToTauTauCR2(unsigned int theIndex);


   float _eventPUWeight;
   std::vector<float> _thePUWeights;
   Long64_t _nEvents;
   int counterEvents;
   std::string _source;
   std::pair<int, int> skimmedEvents;
   std::pair<double, std::pair<double, double> > skimEff;
   double _signalXSection;
   std::pair<float, float> _theXSection; // xsection, error
   double _theLumi;   
   std::pair<double, double> nSkimmedAtLumi;
   std::pair<double, double> eventsAtLumi;
      
   
   std::string _outLogFile;
   std::string _outRootFileName;   
      
   //private functions

   
   void getEventCounters();
   void getCandCounters();
   void getNMinus1CandCounters(unsigned int theIndex);
   void resetCandCounters();
   void resetHistosDefault();
   void bookHistos();
   void fillHistos();
   void getReport();
   void writeOutFile();
   void GetPUWeights();  
   
   void GetEventWeight();
   int getMatchedCand();   
   
   //cuts from config
   double _minElecPtCut;
   double _maxElecEtaCut;  
   double _maxElecDeltaPhiInEBCut;
   double _maxElecDeltaPhiInBBCut;
   double _maxElecDeltaEtaInEBCut;
   double _maxElecDeltaEtaInBBCut;
   double _maxElecSigmaIEtaIEtaEBCut;
   double _maxElecSigmaIEtaIEtaBBCut;
   double _maxElecHadFracEBCut;
   double _maxElecHadFracBBCut;
   double _maxOneOverEMinusOneOverPEBCut;
   double _maxOneOverEMinusOneOverPBBCut;
   double _minElecRelIsoCut;
   double _maxElecRelIsoCut;
   double _minElecTrkIsoCut;
   double _maxElecTrkIsoCut;
   int _maxElecMissingHitsCut;
   //double _elecIsoDr;
   double _maxElecPFRelIsoCut;   
   double _eMotherIdCut;
   double _maxElecDxyBarrelCut;
   double _maxElecDxyEndCapCut;
      
      
   double _minMuonPtCut;
   double _minMuonECut;
   double _maxMuonEtaCut;
   double _maxMuonDxyVtxCut;
   double _maxMuonDzVtxCut;
   double _minMuonIsGlobalMuonCut;
   double _minMuonMatchedCut;
   double _muonMotherIdCut;
   double _maxMuonNormChiSqrdCut;
   double _minMuonValidHitsCut;
   double _minMuonMatchedStationsCut;
   double _minMuonPixelHitsCut;
   double _minMuonTrkLayersWithHitsCut;
   double _maxMuonIsoDR03NoDBCut;
   double _maxMuonIsoDR03DBCut;
   double _maxMuonIsoDR04NoDBCut;
   double _minMuonIsoDR04DBCut;
   double _maxMuonIsoDR04DBCut;
   double _minMuonIsHighPtMuonCut;
   double _maxMuonDptOverPtCut;
   
   
   double _minElecMatchedCut;
   double _isoCut;

   double _minHeepPassedEtCut;
   double _minHeepPassedDetEtaCut;
   double _minHeepPassedCrackCut;
   double _minHeepPassedDEtaInCut;
   double _minHeepPassedDPhiInCut;
   double _minHeepPassedHademCut;
   double _minHeepPassedSigmaIEtaIEtaCut;
   double _minHeepPassed2by5Over5By5Cut;
   double _minHeepPassedEcalHad1IsoCut;
  // double _minHeepPassedHad2IsoCut;
   double _minHeepPassedTrkIsoCut;
   double _minHeepPassedEcalDrivenCut;
   double _minHeepPassedMissingHitsCut;
   double _minHeepPassedElecDxyCut;
   double _minElecEtCut; 
   
   double _minElecIpOverElecIpErrorCut;
   double _minMuonIpOverMuonIpErrorCut;
   double _minCombinedIpOverIpErrorCut;
   double _maxCombinedIpOverIpErrorCut;
   double _minElecMuonDCAOverError2DCut;
   double _maxElecMuonDCAOverError2DCut;
   double _minElecMuonDCAOverError3DCut;
   double _maxElecMuonDCAOverError3DCut;
   double _minLifetimeCut;
   
   double _elecMuonMassCut;
   double _minElecMuonCosDPhiCut;
   double _maxElecMuonCosDPhiCut;
   double _elecMuonDeltaRCut;
   double _minElecMuonMetMassCut;
   double _maxElecMuonMetMassCut;
   double _elecMuonDeltaPtMassCut;
   double _elecMetMtCut;
   double _minMuonMetMtCut;
   double _maxMuonMetMtCut;
   
   double _elecMuonPZetaVisCut;
   double _minElecMuonPZetaCut;
   double _maxElecMuonPZetaCut;
   double _pZetaSlope;
   
   double _pZetaCentroidCut;
   double _pZetaVisCentroidCut;
   double _pZetaSemiMajorAxisCut;
   double _pZetaSemiMinorAxisCut;
   double _pZetaEllipseAngleCut;
   
   unsigned int _minNbTagsHiEffTrkCntCut;
   unsigned int _maxNbTagsHiEffTrkCntCut;
   unsigned int _minNbTagsHiPurityTrkCntCut;
   unsigned int _maxNbTagsHiPurityTrkCntCut;   
   unsigned int _minNbTagsHiEffSimpleSecVtxCut;
   unsigned int _maxNbTagsHiEffSimpleSecVtxCut;
   unsigned int _minNbTagsHiPuritySimpleSecVtxCut;
   unsigned int _maxNbTagsHiPuritySimpleSecVtxCut;   
   unsigned int _minNbTagsCombSecVtxCut;
   unsigned int _maxNbTagsCombSecVtxCut;
   double _leptonMetCosDPhiCut;
   double _maxLeptonDeltaPtCosDPhiCut;
   double _nJetsCut;
   double _jetSumEtCut;
   double _metCut;
   double _maxThreeLeptonCut;
   double _minChargeCut;
   double _maxChargeCut;    
   double _minElecMuonDeltaPtCut;
   double _min2xPLargerLeptonCut;
   double _minElecEOverPCut;
   double _maxElecEOverPCut;
    
    
   //candidate counters
   unsigned int _elecExistsCounter;
   unsigned int _elecPtCounter;
   unsigned int _elecEtaCounter;
   unsigned int _elecAccCounter;  
   unsigned int _elecDeltaPhiInCounter;
   unsigned int _elecDeltaEtaInCounter;
   unsigned int _elecSigmaIEtaIEtaCounter;
   unsigned int _elecTrkChiSqrdCounter;
   unsigned int _elecHadFracCounter;
   unsigned int _elecEOverPCounter;
   unsigned int _elecRelIsoCounter;
   unsigned int _elecTrkIsoCounter;
   unsigned int _elecMissHitsCounter;
   unsigned int _elecPFRelIsoCounter;
   unsigned int _elecGlobalIdCounter;
   unsigned int _elecIdCounter;   
   unsigned int _elecMotherIdCounter;
   unsigned int _elecHOverECounter;
   unsigned int _elec2x5Over5x5Counter;
   unsigned int _elecEcalHad1IsoCounter;
   unsigned int _elecDxyCounter;
   
   unsigned int _muonPtCounter;
   unsigned int _muonECounter;
   unsigned int _muonEtaCounter;
   unsigned int _muonAccCounter;
   unsigned int _muonDxyVtxCounter;
   unsigned int _muonDzVtxCounter;
   unsigned int _muonIsGlobalMuonCounter;
   unsigned int _muonMatchedCounter;
   unsigned int _muonMotherIdCounter;
   unsigned int _muonNormChiSqrdCounter;
   unsigned int _muonTrkChiSqrdCounter;
   unsigned int _muonValidHitsCounter;
   unsigned int _muonMatchedStationsCounter;
   unsigned int _muonPixelHitsCounter;
   unsigned int _muonTrkLayersWithHitsCounter;
   unsigned int _muonIsoDR03NoDBCounter;
   unsigned int _muonIsoDR03DBCounter;
   unsigned int _muonIsoDR04NoDBCounter;
   unsigned int _muonIsoDR04DBCounter;
   unsigned int _muonIsHighPtMuonCounter;
   unsigned int _muonDptOverPtCounter;
   unsigned int _muonIdCounter;
   
   unsigned int _elecMatchedCounter;
   unsigned int _heepPassedEtCounter;
   unsigned int _heepPassedCrackCounter;
   unsigned int _heepPassedDEtaInCounter;
   unsigned int _heepPassedDPhiInCounter;
   unsigned int _heepPassedHademCounter;
   unsigned int _heepPassedSigmaIEtaIEtaCounter;
   unsigned int _heepPassed2by5Over5By5Counter;
   unsigned int _heepPassedEcalHad1IsoCounter;
   //unsigned int _heepPassedHad2IsoCounter;
   unsigned int _heepPassedTrkIsoCounter;
   unsigned int _heepPassedEcalDrivenCounter;
   unsigned int _heepPassedMissingHitsCounter;
   unsigned int _heepPassedElecDxyCounter;
   unsigned int _elecEtCounter;    
   unsigned int _heepIdCounter;

   unsigned int _elecIpOverElecIpErrorCounter;
   unsigned int _muonIpOverMuonIpErrorCounter;
   unsigned int _combinedIpOverIpErrorCounter;
   unsigned int _elecMuonDCAOverError2DCounter;
   unsigned int _elecMuonDCAOverError3DCounter;
   unsigned int _lifetimeCounter;

   unsigned int _matchedCounter;
   
   unsigned int _elecMuonMassCounter;
   unsigned int _elecMuonCosDPhiCounter;
   unsigned int _elecMuonDeltaRCounter;
   unsigned int _elecMuonMetMassCounter;
   unsigned int _elecMuonDeltaPtMassCounter;
   unsigned int _elecMetMtCounter;
   unsigned int _muonMetMtCounter;
   unsigned int _elecMuonPZetaVisCounter;
   unsigned int _elecMuonPZetaCounter;
   unsigned int _nbTagsHiEffTrkCntCounter;
   unsigned int _nbTagsHiPurityTrkCntCounter;
   unsigned int _nbTagsHiEffSimpleSecVtxCounter;
   unsigned int _nbTagsHiPuritySimpleSecVtxCounter;
   unsigned int _nbTagsCombSecVtxCounter;
   unsigned int _leptonMetCosDPhiCounter;
   unsigned int _leptonDeltaPtCosDPhiCounter;
   unsigned int _metCounter;
   unsigned int _threeLeptonCounter;
   unsigned int _nJetsCounter;
   unsigned int _jetSumEtCounter;
   unsigned int _chargeCounter;
   unsigned int _elecMuonDeltaPtCounter;
   unsigned int _2xPLargerLeptonCounter;
   unsigned int _topologyCounter;
   
   unsigned int _lowPtLeptonAntiIsoCounter;
   
   //CONTROL REGION EFFICIENCY COUNTERS
   unsigned int _QCDCR1Counter;
   unsigned int _QCDCR2Counter;
   unsigned int _TTJetsCR1Counter;
   unsigned int _TTJetsCR2Counter;
   unsigned int _WJetsCR1Counter;
   unsigned int _WJetsCR2Counter;
   unsigned int _DYToTauTauCR1Counter;
   unsigned int _DYToTauTauCR2Counter;
   
   //QCD Efficiencies
   unsigned int _QCDOSCounter;
   unsigned int _QCDElectronIsoCounter;
   unsigned int _QCDTopologyCounter;
   unsigned int _QCDMuonIsoCounter;
   unsigned int _QCDMuonAntiIsoCounter;

   //TTJets Efficiencies
   unsigned int _TTJetsJetSumEtCounter;
   unsigned int _TTJetsTopologyCounter;
   unsigned int _TTJetsNoBTagsCounter;
   unsigned int _TTJetsAtLeastOneBTagCounter;
  
   //WJets Efficiencies
   unsigned int _WJetsOSCounter;
   unsigned int _WJetsTopologyCounter;
   unsigned int _WJetsElectronIsoCounter;
   unsigned int _WJetsElectronAntiIsoCounter;
   unsigned int _WJetsMtCounter;

   //DYToTauTau Efficiencies
   unsigned int _DYToTauTauMETCounter;   
   
   //N-1 candidate counters
   unsigned int _nMinus1ElecPtCounter;
   unsigned int _nMinus1ElecEtaCounter;  
   
   unsigned int _nMinus1MuonPtCounter;
   unsigned int _nMinus1MuonECounter;
   unsigned int _nMinus1MuonEtaCounter;
   unsigned int _nMinus1MuonDxyVtxCounter;
   unsigned int _nMinus1MuonDzVtxCounter;
   unsigned int _nMinus1MuonDzBSCounter;
   unsigned int _nMinus1MuonIsGlobalMuonCounter;
   unsigned int _nMinus1MuonMatchedCounter;
   unsigned int _nMinus1MuonNormChiSqrdCounter;
   unsigned int _nMinus1MuonValidHitsCounter;
   unsigned int _nMinus1MuonMatchedStationsCounter;
   unsigned int _nMinus1MuonPixelHitsCounter;
   unsigned int _nMinus1MuonTrkLayersWithHitsCounter;
   unsigned int _nMinus1MuonIsoDR03NoDBCounter;
   unsigned int _nMinus1MuonIsoDR03DBCounter;
   unsigned int _nMinus1MuonIsoDR04NoDBCounter;
   unsigned int _nMinus1MuonIsoDR04DBCounter;   
   unsigned int _nMinus1MuonIsoCounter;
   unsigned int _nMinus1MuonIsHighPtMuonCounter;
   unsigned int _nMinus1MuonDptOverPtCounter;
   
   unsigned int _nMinus1ElecMatchedCounter;
   unsigned int _nMinus1HeepPassedPtCounter;
   unsigned int _nMinus1HeepPassedDetEtaCounter;
   unsigned int _nMinus1HeepPassedCrackCounter;
   unsigned int _nMinus1HeepPassedDEtaInCounter;
   unsigned int _nMinus1HeepPassedDPhiInCounter;
   unsigned int _nMinus1HeepPassedHademCounter;
   unsigned int _nMinus1HeepPassedSigmaIEtaIEtaCounter;
   unsigned int _nMinus1HeepPassed2by5Over5By5Counter;
   unsigned int _nMinus1HeepPassedEcalHad1IsoCounter;
   //unsigned int _nMinus1HeepPassedHad2IsoCounter;
   unsigned int _nMinus1HeepPassedTrkIsoCounter;
   unsigned int _nMinus1HeepPassedEcalDrivenCounter;
   unsigned int _nMinus1ElecEtCounter;       
   unsigned int _nMinus1ElecRelIsoCounter;
   
   unsigned int _nMinus1ElecIpOverElecIpErrorCounter;
   unsigned int _nMinus1MuonIpOverMuonIpErrorCounter;
   unsigned int _nMinus1CombinedIpOverIpErrorCounter;
   unsigned int _nMinus1ElecMuonDCAOverError2DCounter;
   unsigned int _nMinus1ElecMuonDCAOverError3DCounter;
   unsigned int _nMinus1ElecTrkChiSqrdCounter;
   unsigned int _nMinus1MuonTrkChiSqrdCounter;
   
   unsigned int _nMinus1MEtCounter;
   unsigned int _nMinus1BTagCounter;
   unsigned int _nMinus1ElecMuonCosDPhiCounter;
   unsigned int _nMinus1ElecMuonPZetaCounter;
   unsigned int _nMinus1ElecMuonPZetaVisCounter;
   unsigned int _nMinus1JetSumEtCounter;
   unsigned int _nMinus1LeptonMetCosDPhiCounter;
   unsigned int _nMinus1ElecMuonMetMassCounter;
   unsigned int _nMinus1NJetsCounter;
   unsigned int _nMinus1ElecMuonMassCounter;
   unsigned int _nMinus1MuonMetMtCounter;
   unsigned int _nMinus1ElecMuonChargeCounter;   
   
   unsigned int _nMinus1ThreeLeptonCounter;
   
   
   
   //event counter floats due to pileup reweighting
   
   double _puWeightedEvents;
   double _nMatched;
   double _nElecExists;   
   double _nElecPt;
   double _nElecEta;
   double _nElecAcc;
   double _nElecDeltaPhiIn;
   double _nElecDeltaEtaIn;
   double _nElecSigmaIEtaIEta;
   double _nElecHadFrac;
   double _nElecTrkChiSqrd;
   double _nElecEOverP;
   double _nElecRelIso;
   double _nElecTrkIso;
   double _nElecMissHits;
   double _nElecPFRelIso;
   double _nElecGlobalId;
   double _nElecId;   
   double _nElecMotherId;
   double _nElecHOverE;
   double _nElec2x5Over5x5;
   double _nElecEcalHad1Iso;
   double _nElecDxy;
   
   double _nMuonPt;
   double _nMuonE;
   double _nMuonEta;
   double _nMuonAcc;
   double _nMuonDxyVtx;
   double _nMuonDzVtx;
   double _nMuonIsGlobalMuon;
   double _nMuonMatched;
   double _nMuonMotherId;
   double _nMuonNormChiSqrd;
   double _nMuonTrkChiSqrd;
   double _nMuonMatchedStations;
   double _nMuonValidHits;
   double _nMuonPixelHits;
   double _nMuonTrkLayersWithHits;
   double _nMuonIsoDR03NoDB;
   double _nMuonIsoDR03DB;
   double _nMuonIsoDR04NoDB;
   double _nMuonIsoDR04DB;
   double _nMuonIsHighPtMuon;
   double _nMuonDptOverPt;
   double _nMuonId;
  
   
   double _nElecMatched;
   double _nHeepPassedEt;
   double _nHeepPassedCrack;
   double _nHeepPassedDEtaIn;
   double _nHeepPassedDPhiIn;
   double _nHeepPassedHadem;
   double _nHeepPassedSigmaIEtaIEta;
   double _nHeepPassed2by5Over5By5;
   double _nHeepPassedEcalHad1Iso;
   //double _nHeepPassedHad2Iso;
   double _nHeepPassedTrkIso;
   double _nHeepPassedEcalDriven;
   double _nHeepPassedMissingHits;
   double _nHeepPassedElecDxy;
   double _nElecEt;      
   double _nHeepId;
   
   double _nElecIpOverElecIpError;
   double _nMuonIpOverMuonIpError;
   double _nCombinedIpOverIpError;
   double _nElecMuonDCAOverError2D;
   double _nElecMuonDCAOverError3D;
   double _nLifetime;
   
   double _nElecMuonMass;
   double _nElecMuonCosDPhi;
   double _nElecMuonDeltaR;
   double _nElecMuonMetMass;
   double _nElecMuonDeltaPtMass;
   double _nElecMetMt;
   double _nMuonMetMt;
   double _nElecMuonPZetaVis;
   double _nElecMuonPZeta;
   double _nNbTagsHiEffTrkCnt;
   double _nNbTagsHiPurityTrkCnt;
   double _nNbTagsHiEffSimpleSecVtx;
   double _nNbTagsHiPuritySimpleSecVtx;
   double _nNbTagsCombSecVtx;
   double _nLeptonMetCosDPhi;
   double _nLeptonDeltaPtCosDPhi;
   double _nMEt;
   double _nThreeLepton;
   double _nNJets;
   double _nJetSumEt;
   double _nCharge;
   double _nElecMuonDeltaPt;
   double _n2xPLargerLepton;
   double _nTopology;
   
   double _nLowPtLeptonAntiIso;
   
   double _nDummy;

   //Control Region Event Counters
   double _nQCDCR1;
   double _nQCDCR2;
   double _nTTJetsCR1;
   double _nTTJetsCR2;
   double _nWJetsCR1;
   double _nWJetsCR2;
   double _nDYToTauTauCR1;
   double _nDYToTauTauCR2;
   
   double _nQCDOS;
   double _nQCDElectronIso;
   double _nQCDTopology;
   double _nQCDMuonIso;
   double _nQCDMuonAntiIso;
  

   double _nTTJetsJetSumEt;
   double _nTTJetsTopology;
   double _nTTJetsNoBTags;
   double _nTTJetsAtLeastOneBTag;
   
   double _nWJetsOS;
   double _nWJetsTopology;
   double _nWJetsElectronIso;
   double _nWJetsElectronAntiIso;
   double _nWJetsMt;


   //N-1 Event Counters
   double _nElecPtNMinus1;
   double _nElecEtaNMinus1;
 
   double _nMuonPtNMinus1;
   double _nMuonENMinus1;
   double _nMuonEtaNMinus1;
   double _nMuonDxyVtxNMinus1;
   double _nMuonDzVtxNMinus1;
   double _nMuonDzBSNMinus1;
   double _nMuonIsGlobalMuonNMinus1;
   double _nMuonMatchedNMinus1;
   double _nMuonNormChiSqrdNMinus1;
   double _nMuonMatchedStationsNMinus1;
   double _nMuonValidHitsNMinus1;
   double _nMuonPixelHitsNMinus1;
   double _nMuonTrkLayersWithHitsNMinus1;
   double _nMuonIsoDR03NoDBNMinus1;
   double _nMuonIsoDR03DBNMinus1;
   double _nMuonIsoDR04NoDBNMinus1;
   double _nMuonIsoDR04DBNMinus1;
   double _nMuonIsoNMinus1; 
   double _nMuonIsHighPtMuonNMinus1;
   double _nMuonDptOverPtNMinus1;
  
   double _nElecMatchedNMinus1;
   double _nHeepPassedPtNMinus1;
   double _nHeepPassedDetEtaNMinus1;
   double _nHeepPassedCrackNMinus1;
   double _nHeepPassedDEtaInNMinus1;
   double _nHeepPassedDPhiInNMinus1;
   double _nHeepPassedHademNMinus1;
   double _nHeepPassedSigmaIEtaIEtaNMinus1;
   double _nHeepPassed2by5Over5By5NMinus1;
   double _nHeepPassedEcalHad1IsoNMinus1;
   double _nHeepPassedTrkIsoNMinus1;
   //double _nHeepPassedHad2IsoNMinus1;
   double _nHeepPassedEcalDrivenNMinus1;
   double _nElecEtNMinus1;         
   double _nElecRelIsoNMinus1;
    
   double _nElecIpOverElecIpErrorNMinus1;
   double _nMuonIpOverMuonIpErrorNMinus1;   
   double _nCombinedIpOverIpErrorNMinus1; 
   double _nElecMuonDCAOverError2DNMinus1;
   double _nElecMuonDCAOverError3DNMinus1;  
   double _nElecTrkChiSqrdNMinus1;
   double _nMuonTrkChiSqrdNMinus1;
      
   double _nMEtNMinus1;
   double _nBTagsNMinus1;   
   double _nElecMuonCosDPhiNMinus1;
   double _nElecMuonPZetaNMinus1;
   double _nElecMuonPZetaVisNMinus1;
   double _nElecMuonMassNMinus1;
   double _nMuonMetMtNMinus1;
   double _nJetSumEtNMinus1;
   double _nLeptonMetCosDPhiNMinus1;
   double _nElecMuonMetMassNMinus1;
   double _nNJetsNMinus1;
   double _nElecMuonChargeNMinus1;   
   
   
 
   //values for histos
   float _maxElecPt;
   float _maxElecP;
   float _minElecEta;           
   float _minElecDeltaPhiInBarrel;
   float _minElecDeltaEtaInBarrel;
   float _minElecSigmaIEtaBarrel;
   float _minElecHadFracBarrel;
   float _minElecTrkChiSqrd;
   float _maxElecTrkNdof;
   float _minElecEOverPBarrel;
   float _minElecDeltaPhiInEndCap;
   float _minElecDeltaEtaInEndCap;
   float _minElecSigmaIEtaEndCap;
   float _minElecHadFracEndCap;
   float _minElecEOverPEndCap;
   float _minElecRelIso;
   float _minElecTrkIso;
   float _maxElecMissingHits;
   float _minPFIsoChargedHadDr03;
   float _minPFIsoNeutralHadDr03;
   float _minPFIsoPhotonDr03;
   float _minPFRelIsoDr03;
   float _minPFIsoChargedHadDr04;
   float _minPFIsoNeutralHadDr04;
   float _minPFIsoPhotonDr04;
   float _minPFRelIsoDr04;  
   
   float _maxHeepPassedEcalDriven;
   float _maxHeepPassedCrack;
   float _maxHeepPassedDEtaIn;
   float _maxHeepPassedDPhiIn;
   float _maxHeepPassedHadem;
   float _maxHeepPassedSigmaIEtaIEta;
   float _maxHeepPassed2by5Over5By5;
   float _maxHeepPassedEcalHad1Iso;
   float _maxHeepPassedTrkIso;
   float _maxElecEt;
   
   
   float _maxElecSCE1x5;
   float _maxElecSCE2x5;
   float _maxElecSCE5x5;
   float _maxElecSCE1x5Over5x5;
   float _maxElecSCE2x5MaxOver5x5;
   float _maxElecSCE1x5EOverP;
   float _maxElecSCE2x5EOverP;
   float _maxElecSCE5x5EOverP;
   
   
   float _maxMetMinusDeltaPt;
   
   float _maxDeltaPtx;
   float _maxDeltaPty;
   float _maxMetPx;
   float _maxMetPy; 
    
   float _maxTauEnergy; 
   float _maxMuonPOverTauE;
   float _maxZPrimeEnergy;
   float _maxVertexDiff;
   float _maxVertexDiffxy;
   
   double _maxBeamSpotDiff;
   double _maxBeamSpotDiffxy;
   double _maxBeamSpotPrimaryzDiff;
   
    
   float _maxElecIpOverElecIpError;
   float _maxMuonIpOverMuonIpError;
   float _maxCombinedIpOverIpError; 
   float _maxElecIp;
   float _maxMuonIp;
   float _maxElecIpError;
   float _maxMuonIpError;
   float _maxElecMuonDistanceOfClosestApproach;
   float _maxElecMuonDCACrossingPoint;   
   float _maxElecMuonTwoTrkMinDist;
   float _maxElecMuonTwoTrkCrossingPoint;
   float _maxElecMuonDis3D;
   float _maxElecMuonDis2D;
   float _maxElecMuonDCA3DError;
   float _maxElecMuonDCA2DError;
   float _maxElecMuonDis3DOverError;
   float _maxElecMuonDis2DOverError;
   float _maxLifetime;
   float _maxElecEOverP_hiIP;
   float _maxElecEta_hiIP;
   float _maxElecEt_hiIP;
   float _maxElecTrkChiSqrd_hiIP;
 
   float _maxMuonPt;
   float _maxMuonP;
   float _maxMuonE;
   float _minMuonEta;
   float _maxMuonDxyVtx;
   float _maxMuonDzVtx;
   float _maxMuonNormChiSqrd;
   float _minMuonTrkChiSqrd;
   float _maxMuonTrkNdof;
   float _maxIsGlobalMuon;
   float _maxMuonValidHits;
   float _maxMuonMatchedStations;
   float _maxMuonPixelHits;
   float _maxMuonTrkLayersWithHits;
   float _maxMuonIsoDR04DB;
   float _maxMuonIsHighPtMuon;
   float _minMuonDptOverPt;
   
      									      
   int _minCharge;
   float _minElecMuonCosDPhi;
   float _maxElecMuonDeltaR;
   float _maxElecMuonPZeta;
   float _maxElecMuonPZetaVis;
   float _minLeptonMetCosDPhi;	
   float _minLeptonDeltaPtCosDPhi;								      
   float _maxMEt;  
   unsigned int _maxThreeLepton;										      
   unsigned int _minNBtagsHiEffTrkCnt;  
   float _maxElecMuonMetMass;
   float _maxElecMuonDeltaPtMass;
   float _maxElecMuonMetMassBeginning;
   
   float _maxElecPtVsMass;
   float _maxMuonPtVsMass;	   
   
   float _maxBSx0;
   float _maxBSy0;
   float _maxBSz0;
   
   float _maxMuonVertex_x0;
   float _maxMuonVertex_y0;
   float _maxMuonVertex_z0;
   float _maxMuonDecayLength_x;
   float _maxMuonDecayLength_y;
   float _maxMuonDecayLength_z;
   
   float _maxElecVertex_x0;
   float _maxElecVertex_y0;
   float _maxElecVertex_z0;
   float _maxElecDecayLength_x;
   float _maxElecDecayLength_y;
   float _maxElecDecayLength_z; 
   double _maxElecVx;  
   double _maxElecVy;  
   double _maxElecVz;  
   
   float _maxPrimaryVertex_x;
   float _maxPrimaryVertex_y;
   float _maxPrimaryVertex_z;
   float _maxPrimaryVertex_xError;
   float _maxPrimaryVertex_yError;
   float _maxPrimaryVertex_zError;
   float _maxGenTauVtx_x;
   float _maxGenTauVtx_y;
   float _maxGenTauVtx_z;
   float _maxOtherGenTauVtx_x;
   float _maxOtherGenTauVtx_y;
   float _maxOtherGenTauVtx_z;   
   
   float _maxElecMotherId;
   float _maxMuonMotherId;
   
   float _maxMEtVsDeltaPt;
   float _maxElecPtVsDeltaPt;
   float _maxMuonPtVsDeltaPt;
   
   float _maxElecDeltaPtxNMinus1;
   float _maxElecDeltaPtyNMinus1;
   float _maxMuonDeltaPtNMinus1;
   
   float _maxMuonDptOverPtNMinus1;
   float _maxElecDptOverPtNMinus1;
   
   float _maxElecPtNMinus1;
   float _maxElecEtNMinus1;
   float _minElecEtaNMinus1;
   float _minElecRelIsoNMinus1;
   
   float _maxMuonPtNMinus1;
   float _maxMuonENMinus1;
   float _minMuonEtaNMinus1;
   float _maxMuonDxyVtxNMinus1;
   float _maxMuonDzVtxNMinus1; 
   float _maxMuonDzBSNMinus1;
   float _minMuonIsoNMinus1;
   float _maxMuonIsHighPtMuonNMinus1;
   float _minMuonDptOverPtNMinus1;
   
   float _maxElecIpOverElecIpErrorNMinus1;
   float _maxMuonIpOverMuonIpErrorNMinus1;
   float _maxCombinedIpOverIpErrorNMinus1;
   float _maxElecMuonDCAOverError2DNMinus1;
   float _maxElecMuonDCAOverError3DNMinus1;
   
   
   									      
   int _minChargeNMinus1;										      
   float _maxMEtNMinus1;	
   int _maxThreeLeptonNMinus1;	
   float _maxThreeLeptonMassNMinus1;								      
   unsigned int _minNBtagsHiEffTrkCntNMinus1;      
   unsigned int _minNBtagsCombSecVtxNMinus1;
   float _maxElecMuonPZetaNMinus1;
   float _maxElecMuonPZetaVisNMinus1;
   float _minElecMuonCosDPhiNMinus1;
   float _minElecMuonMassNMinus1;
   float _minMuonMetMtNMinus1;
   float _minJetSumEtNMinus1;
   float _minLeptonMetCosDPhiNMinus1;
   float _minLeptonDeltaPtCosDPhiNMinus1;
   float _minElecMuonMetMassNMinus1;
   float _maxElecMuonDeltaPtMassNMinus1;
   float _maxNJetsNMinus1;
   float _minElecMuonChargeNMinus1;
   
   unsigned int _minNBtagsCombSecVtxCR2;
   
      
   //histograms
   
   TH1F* dataPUHisto;
   TH1F* _pileUPDistro;
   TH1F* _pileUPReweighted;   
   
   TH1F* _tauEnergyHisto;
   TH1F* _muonPOverTauEHisto;
   TH2F* _muonPOverTauEVsMuonIpOverErrorHisto;
   TH2F* _muonPOverTauEVsCosDPhiHisto;
   TH2F* _muonPOverTauEVsTauDecayLengthHisto;
   TH1F* _zPrimeEnergyHisto;
   TH1F* _vertexDiffHisto;
   TH1F* _vertexDiffxyHisto;
   TH2F* _impactParameterVsVertexDiffHisto;
   TH2F* _tauEnergyVsTauDecayLengthHisto;
   
   TH1F* _beamSpotDiffHisto;
   TH1F* _beamSpotDiffxyHisto;
   
   TH1F* _beamSpotPrimaryzDiffHisto;
   
   TH2F* _tauEnergyVsElecMuonMetMassHisto;
   TH2F* _zPrimeEnergyVsElecMuonMetMassHisto;
   TH2F* _tauEnergyVsMuonEtaHisto;
   TH2F* _zPrimeEnergyVsMuonEtaHisto;
   TH2F* _tauEnergyVsElecEtaHisto;
   TH2F* _zPrimeEnergyVsElecEtaHisto;
   TH2F* _zPrimeEnergyVs2xMuonPtHisto;
   TH2F* _zPrimeEnergyVs2xElecEtHisto;
   TH2F* _zPrimeEnergyVs2xMaxLeptonPtHisto;
   TH2F* _zPrimeEnergyVs2xMaxLeptonPHisto;
   TH1I* _nVtxHisto;  
   TH1F* _elecPtHisto;
   TH1F* _elecEtaHisto;  
   TH1F* _elecDeltaPhiInBarrelHisto;
   TH1F* _elecDeltaEtaInBarrelHisto;
   TH1F* _elecSigmaIEtaIEtaBarrelHisto;
   TH1F* _elecHadFracBarrelHisto;
   TH1F* _elecEOverPBarrelHisto;
   TH1F* _elecDeltaPhiInEndCapHisto;
   TH1F* _elecDeltaEtaInEndCapHisto;
   TH1F* _elecSigmaIEtaIEtaEndCapHisto;
   TH1F* _elecHadFracEndCapHisto;
   TH1F* _elecEOverPEndCapHisto;
   TH1F* _elecEcalHisto;
   TH1F* _elecHcalHisto;
   TH1I* _elecMissingHitsHisto;
   TH1F* _elecPFDR03ChargedHadIsoHisto;
   TH1F* _elecPFDR03NeutralHadIsoHisto;
   TH1F* _elecPFDR03PhotonIsoHisto;
   TH1F* _elecPFDR03RelIsoHisto;
   TH1F* _elecPFDR04ChargedHadIsoHisto;
   TH1F* _elecPFDR04NeutralHadIsoHisto;
   TH1F* _elecPFDR04PhotonIsoHisto;
   TH1F* _elecPFDR04RelIsoHisto;   
   TH1F* _elecRelIsoHisto;
   TH1F* _elecTrkIsoHisto;
   
   TH1F* _elecSCE1x5Histo;
   TH1F* _elecSCE2x5Histo;
   TH1F* _elecSCE5x5Histo;
   TH1F* _elecSCE1x5Over5x5Histo;
   TH1F* _elecSCE2x5MaxOver5x5Histo;
   TH1F* _elecSCE1x5EOverPHisto;
   TH1F* _elecSCE2x5EOverPHisto;
   TH1F* _elecSCE5x5EOverPHisto;
   
   TH2F* _elecSCE1x5VsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE2x5VsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE5x5VsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE1x5Over5x5VsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE2x5MaxOver5x5VsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE1x5EOverPVsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE2x5EOverPVsElecIpOverError_hiIPHisto;
   TH2F* _elecSCE5x5EOverPVsElecIpOverError_hiIPHisto;
   
   TH2F* _elecSCE1x5EOverPVsEOverP_hiIPHisto;
   TH2F* _elecSCE2x5EOverPVsEOverP_hiIPHisto;
   TH2F* _elecSCE5x5EOverPVsEOverP_hiIPHisto;
   
   TH1F* _muonPtHisto;
   TH1F* _muonEHisto;
   TH1F* _muonEtaHisto;
   TH2F* _muonPtVsMuonEtaHisto;
   TH2F* _muonEVsMuonEtaHisto;
   TH1F* _muonDxyVtxHisto;
   TH1F* _muonDzVtxHisto;
   TH1I* _muonIsGlobalMuonHisto;
   TH1I* _muonMatchedHisto;
   TH1F* _muonNormChiSqrdHisto;
   TH1I* _muonValidHitsHisto;
   TH1I* _muonMatchedStationsHisto;
   TH1I* _muonPixelHitsHisto;
   TH1F* _muonTrkLayersWithHitsHisto;
   TH1F* _muonIsoDR04DBHisto;
   TH1I* _muonIsHighPtMuonHisto;
   TH1F* _muonDptOverPtHisto;
   
   TH1I* _elecMatchedHisto;
   TH1I* _heepPassedPtHisto;
   TH1I* _heepPassedDetEtaHisto;
   TH1I* _heepPassedCrackHisto;
   TH1I* _heepPassedDEtaInHisto;
   TH1I* _heepPassedDPhiInHisto;
   TH1I* _heepPassedHademHisto;
   TH1I* _heepPassedSigmaIEtaIEtaHisto;
   TH1I* _heepPassed2by5Over5By5Histo;
   TH1I* _heepPassedEcalHad1IsoHisto;
   //TH1F* _heepPassedHad2IsoHisto;
   TH1I* _heepPassedTrkIsoHisto;
   TH1I* _heepPassedEcalDrivenHisto;
   TH1F* _elecEtHisto; 
   
   TH1F* _MetMinusDeltaPtHisto;  
   TH2F* _MetxVsDeltaPtxHisto;
   TH2F* _MetyVsDeltaPtyHisto;
   
   TH1F* _elecIpOverElecIpErrorHisto;
   TH1F* _muonIpOverMuonIpErrorHisto; 
   TH1F* _combinedIpOverIpErrorHisto;
   TH1F* _elecTrkChiSqrdHisto;
   TH1F* _muonTrkChiSqrdHisto;
   TH1F* _elecIpErrorHisto;
   TH1F* _muonIpErrorHisto;
   
   TH1F* _elecMuonDistanceOfClosestApproachHisto;
   TH2F* _elecMuonDCAVsCrossingPointHisto;
   
   TH1F* _elecMuonTwoTrkMinDistHisto;
   TH2F* _elecMuonTwoTrkMinDistVsCrossingPointHisto;
   
   TH1F* _eMuDis3DHisto;
   TH1F* _eMuDis2DHisto;
   TH1F* _elecMuonDCA3DErrorHisto;
   TH1F* _elecMuonDCA2DErrorHisto;
   TH1F* _eMuDis3DOverErrorHisto;
   TH1F* _eMuDis2DOverErrorHisto;
   
   TH2F* _elecIpVsElecTrkChiSqrdHisto;
   TH2F* _elecIpOverErrVsElecTrkNdofHisto;
   TH2F* _muonIpVsMuonTrkChiSqrdHisto;
   TH2F* _muonIpOverErrVsMuonTrkNdofHisto;
   TH2F* _elecIpErrVsElecTrkChiSqrdHisto;
   TH2F* _muonIpErrVsMuonTrkChiSqrdHisto;
   TH2F* _elecIpOverErrVsElecIpErrHisto;
   TH2F* _muonIpOverErrVsMuonIpErrHisto;
   TH2F* _muonIpOverErrVsElecIpOverErrHisto;
   TH2F* _IpVsElecTrkChiSqrdHisto;
   TH2F* _IpVsMuonTrkChiSqrdHisto;
   
   TH2F* _elecIpOverErrVsElecMissingHitsHisto;
   TH2F* _muonIpOverErrVsMuonTrkLayersWithHitsHisto;
   TH2F* _muonIpOverErrVsMuonValidHitsHisto;
   
   TH2F* _muonIpOverErrVsMuonPixelHitsHisto;
   TH2F* _muonEtaVsMuonPixelHitsHisto;
   TH2F* _muonIpOverErrVsMuonEtaHisto;
   TH2F* _muonPtVsMuonIpOverErrHisto;
   TH2F* _elecEtVsElecIpOverErrHisto;
   TH2F* _MEtVsElecIpOverErrHisto;
   TH2F* _MEtVsMuonIpOverErrHisto;
   TH2F* _MEtVsCombinedIpOverErrHisto;
   
   TH1F* _elecEOverP_hiIPHisto;
   TH1F* _elecEta_hiIPHisto;
   TH2F* _elecEOverPVsElecIpOverError_hiIPHisto;
   TH2F* _elecEtaVsElecIpOverError_hiIPHisto;
   TH2F* _elecEtVsElecIpOverError_hiIPHisto;
   TH2F* _elecEOverPVsMuonIpOverError_hiIPHisto;
   TH2F* _elecTrkChiSqrdVsElecIpOverError_hiIPHisto;
   TH2F* _elecTrkChiSqrdVsElecEOverP_hiIPHisto;
   TH2F* _elecEOverPVsElecMissingHits_hiIPHisto;
   TH2F* _elecEOverPVsElecIpError_hiIPHisto;
   TH2F* _elecPVsElecIpError_hiIPHisto;
   TH2F* _elecPVsElecTrkChiSqrd_hiIPHisto;
   TH2F* _elecEOverPVsElecP_hiIPHisto;

   TH2F* _muonIpOverErrVsMuonIpHisto;
   TH2F* _elecIpOverErrVsElecIpHisto;

   TH2F* _elecPtVsElecIpOverErrHisto;
   TH2F* _BSy0VsBSx0Histo;
   
   TH1F* _muonVertex_x0Histo;
   TH1F* _muonVertex_y0Histo;
   TH1F* _muonVertex_z0Histo;
   TH1F* _muonDecayLength_xHisto;
   TH1F* _muonDecayLength_yHisto;
   TH1F* _muonDecayLength_zHisto;
   
   TH1F* _elecVertex_x0Histo;
   TH1F* _elecVertex_y0Histo;
   TH1F* _elecVertex_z0Histo;
   TH1F* _elecDecayLength_xHisto;
   TH1F* _elecDecayLength_yHisto;
   TH1F* _elecDecayLength_zHisto; 
   
   TH1D* _elecVxHisto;
   TH1D* _elecVyHisto;
   TH1D* _elecVzHisto;
   
   
   TH1F* _muonDecayLengthHisto;
   TH1F* _elecDecayLengthHisto;  
   
   TH1F* _elecIpHisto;
   TH1F* _muonIpHisto;
   
   TH1F* _lifetimeHisto;
   
   TH2F* _elecPtVsMuonPtHisto;
   TH2F* _elecEtVsMuonPtHisto;
   
   TH1I* _elecMotherIdHisto;
   TH1I* _muonMotherIdHisto;
   TH1I* _combinedMotherIdHisto;
   
   TH1F* _elecMuonCosDPhiHisto;
   TH1F* _elecMuonDeltaRHisto;
   TH1F* _elecMuonPZetaHisto;
   TH1F* _leptonMetCosDPhiHisto;
   TH1F* _leptonDeltaPtCosDPhiHisto;
   TH1F* _MEtHisto;
   TH1I* _threeLeptonHisto;
   TH1I* _nBTagHisto;
   TH1F* _elecMuonMetMassHisto;
   TH1F* _elecMuonDeltaPtMassHisto;
   TH1F* _elecMuonMetMassBeginningHisto;
   
   TH2F* _elecMuonMetMassVs2xLeptonPtHisto;
   TH2F* _elecMuonMetMassVsElecMuonDeltaPtMassHisto;
 
   TH1F* _elecMEtMassHisto;
   
   TH2F* _elecMuonDCAOverErrorVsCosDPhiHisto;
   
   TH2F* _DCAOverErrorVsElecIpOverErrorHisto;
   TH2F* _DCAOverErrorVsMuonIpOverErrorHisto;
   TH2F* _DCAOverErrorVsCombinedIpOverErrorHisto;
   
 
   //N-1 Histos
   TH1F* _elecPtNMinus1Histo;
   TH1F* _elecEtaNMinus1Histo;
   TH1F* _elecTrkChiSqrdNMinus1Histo;
   
   TH1F* _muonPtNMinus1Histo;
   TH1F* _muonENMinus1Histo;
   TH1F* _muonEtaNMinus1Histo;
   TH1F* _muonDxyVtxNMinus1Histo;
   TH1F* _muonDzVtxNMinus1Histo;
   TH1F* _muonDzBSNMinus1Histo;
   TH1F* _muonIsGlobalMuonNMinus1Histo;
   TH1F* _muonMatchedNMinus1Histo;
   TH1F* _muonNormChiSqrdNMinus1Histo;
   TH1F* _muonTrkChiSqrdNMinus1Histo;
   TH1F* _muonValidHitsNMinus1Histo;
   TH1F* _muonMatchedStationsNMinus1Histo;
   TH1F* _muonPixelHitsNMinus1Histo;
   TH1F* _muonTrkLayersWithHitsNMinus1Histo;
   TH1F* _muonIsoNMinus1Histo;
   TH1I* _muonIsHighPtMuonNMinus1Histo;
   TH1F* _muonDptOverPtNMinus1Histo;
   
   
   TH1F* _elecMatchedNMinus1Histo;
   TH1F* _heepPassedPtNMinus1Histo;
   TH1F* _heepPassedDetEtaNMinus1Histo;
   TH1F* _heepPassedCrackNMinus1Histo;
   TH1F* _heepPassedDEtaInNMinus1Histo;
   TH1F* _heepPassedDPhiInNMinus1Histo;
   TH1F* _heepPassedHademNMinus1Histo;
   TH1F* _heepPassedSigmaIEtaIEtaNMinus1Histo;
   TH1F* _heepPassed2by5Over5By5NMinus1Histo;
   TH1F* _heepPassedEcalHad1IsoNMinus1Histo;
   //TH1F* _heepPassedHad2IsoNMinus1Histo;
   TH1F* _heepPassedTrkIsoNMinus1Histo;
   TH1F* _heepPassedEcalDrivenNMinus1Histo;
   TH1F* _elecEtNMinus1Histo;       
   TH1F* _elecRelIsoNMinus1Histo;
   
   TH1F* _elecIpOverElecIpErrorNMinus1Histo;
   TH1F* _muonIpOverMuonIpErrorNMinus1Histo;
   TH1F* _combinedIpOverIpErrorNMinus1Histo;
   TH1F* _eMuDis2DOverErrorNMinus1Histo;
   TH1F* _eMuDis3DOverErrorNMinus1Histo;
   
   TH1F* _metNMinus1Histo;
   TH1F* _nBTagNMinus1Histo;   
   TH1F* _elecMuonCosDPhiNMinus1Histo;
   TH1F* _elecMuonPZetaNMinus1Histo;
   TH1F* _elecMuonPZetaVisNMinus1Histo;
   TH1F* _elecMuonMassNMinus1Histo;
   TH1F* _muonMetMtNMinus1Histo;
   TH1F* _jetSumEtNMinus1Histo;
   TH1F* _leptonMetCosDPhiNMinus1Histo;
   TH1F* _leptonDeltaPtCosDPhiNMinus1Histo;
   TH2F* _leptonMetCosDPhiVsLeptonDeltaPtCosDPhiNMinus1Histo;
   TH2F* _elecMuonCosDPhiVsLeptonMetCosDPhiNMinus1Histo;
   
   TH1F* _elecMuonMetMassNMinus1Histo;
   TH1F* _elecMuonMetMassNMinus1RebinHisto;
   TH1F* _elecMuonDeltaPtMassNMinus1Histo;
   TH1F* _elecMuonDeltaPtMassNMinus1RebinHisto;   
   TH1F* _nJetsNMinus1Histo;
   TH1F* _elecMuonChargeNMinus1Histo;
   TH1F* _MEtNMinus1Histo;
   TH1I* _threeLeptonNMinus1Histo;
   TH1F* _threeLeptonMassNMinus1Histo;
   TH1F* _nBTagCR2Histo;
   
   TH2F* _PZetaVsPZetaVisHisto;
   TH2F* _nBTagVsJetSumEtHisto;
   TH2F* _nBTagVsMEtHisto;
   TH2F* _nBTagVsMuonMetMtHisto;
   TH2F* _nBTagVsNJetsHisto;   
   
   TH2F* _MEtVsLeptonDeltaPHisto;
   
   TH2F* _elecDeltaPtOverPtVsPtNMinus1Histo;
   TH2F* _muonDeltaPtOverPtVsPtNMinus1Histo;
   
   TH1F* _elecDeltaPtOverPtNMinus1Histo;
   TH1F* _muonDeltaPtOverPtNMinus1Histo;
   

   
};

#endif

#ifdef EMuAnalysis_cxx
EMuAnalysis::EMuAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      std::cout << "TREE == 0\n";
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("eMuAnalysis_10_1_6jf.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("eMuAnalysis_10_1_6jf.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("eMuAnalysis_10_1_6jf.root:/analyzeHiMassTau");
      dir->GetObject("HMTTree",tree);

   }
   Init(tree);
   GetPUWeights();
   
   //set float counters to zero
   _nMatched = 0.;
   _nElecExists = 0.;   
   _nElecPt = 0.;
   _nElecEta = 0.;
   _nElecAcc = 0.;   
   _nElecDeltaPhiIn = 0.;
   _nElecDeltaEtaIn = 0.;
   _nElecSigmaIEtaIEta = 0.;
   _nElecHadFrac = 0.;
   _nElecTrkChiSqrd = 0.;
   _nElecEOverP = 0.;
   _nElecRelIso = 0.;
   _nElecTrkIso = 0.;
   _nElecMissHits = 0.;
   _nElecPFRelIso = 0.;
   _nElecGlobalId = 0.;
   _nElecId = 0.;
   _nElecMotherId = 0.;   
   _nElecHOverE = 0.;
   _nElec2x5Over5x5 = 0.;
   _nElecEcalHad1Iso = 0.;
   _nElecDxy = 0.;
      
   _nMuonPt = 0.;
   _nMuonE = 0.;
   _nMuonEta = 0.;
   _nMuonAcc = 0.;
   _nMuonDxyVtx = 0.;
   _nMuonDzVtx = 0.;
   _nMuonIsGlobalMuon = 0.;
   _nMuonMatched = 0.;
   _nMuonMotherId = 0.;
   _nMuonNormChiSqrd = 0.;
   _nMuonTrkChiSqrd = 0.;
   _nMuonValidHits = 0.;
   _nMuonMatchedStations = 0.;
   _nMuonPixelHits = 0.;
   _nMuonTrkLayersWithHits = 0.;
   _nMuonIsoDR03NoDB = 0.;
   _nMuonIsoDR03DB = 0.;
   _nMuonIsoDR04NoDB = 0.;
   _nMuonIsoDR04DB = 0.;
   _nMuonIsHighPtMuon = 0.;
   _nMuonDptOverPt = 0.;
   _nMuonId = 0.;
   
   
   _nElecMatched = 0.;
   _nHeepPassedEt = 0.;
   _nHeepPassedCrack = 0.;
   _nHeepPassedDEtaIn = 0.;
   _nHeepPassedDPhiIn = 0.;
   _nHeepPassedHadem = 0.;
   _nHeepPassedSigmaIEtaIEta = 0.;
   _nHeepPassed2by5Over5By5 = 0.;
   _nHeepPassedEcalHad1Iso = 0.;
   _nHeepPassedEcalDriven = 0.;
   _nHeepPassedTrkIso = 0.;
   _nHeepPassedEcalDriven = 0.;
   _nHeepPassedMissingHits = 0.;
   _nHeepPassedElecDxy = 0.;
   _nElecEt = 0.;      
   _nHeepId = 0.; 
   
   _nElecIpOverElecIpError = 0.;
   _nMuonIpOverMuonIpError = 0.;
   _nCombinedIpOverIpError = 0.;
   _nElecMuonDCAOverError3D = 0.;
   _nElecMuonDCAOverError2D = 0.;
   _nLifetime = 0.;

   _nElecMuonMass = 0.;
   _nElecMuonCosDPhi = 0.;
   _nElecMuonDeltaR = 0.;
   _nElecMuonMetMass = 0.;
   _nElecMuonDeltaPtMass = 0.;
   _nElecMetMt = 0.;
   _nMuonMetMt = 0.;
   _nElecMuonPZetaVis = 0.;
   _nElecMuonPZeta = 0.;
   _nNbTagsHiEffTrkCnt = 0.;
   _nNbTagsHiPurityTrkCnt = 0.;
   _nNbTagsHiEffSimpleSecVtx = 0.;
   _nNbTagsHiPuritySimpleSecVtx = 0.;
   _nNbTagsCombSecVtx = 0.;
   _nLeptonMetCosDPhi = 0.;
   _nLeptonDeltaPtCosDPhi = 0.;
   _nMEt = 0.;
   _nThreeLepton = 0.;
   _nNJets = 0.;
   _nJetSumEt = 0.;
   _nCharge = 0.;
   _nElecMuonDeltaPt = 0.;
   _n2xPLargerLepton = 0.;
   _nTopology = 0.;
   
   _nLowPtLeptonAntiIso = 0.;
   
   _nQCDCR1 = 0.;
   _nQCDCR2 = 0.;
   _nTTJetsCR1 = 0.;
   _nTTJetsCR2 = 0.;
   _nWJetsCR1 = 0.;
   _nWJetsCR2 = 0.;
   _nDYToTauTauCR1 = 0.;
   _nDYToTauTauCR2 = 0.;
   
   _nQCDOS = 0.;
   _nQCDElectronIso = 0.;
   _nQCDTopology = 0.;
   _nQCDMuonIso = 0.;
   _nQCDMuonAntiIso = 0.;
   
   _nTTJetsJetSumEt = 0.;
   _nTTJetsTopology = 0.;
   _nTTJetsNoBTags = 0.;
   _nTTJetsAtLeastOneBTag = 0.;
   
   _nWJetsOS = 0.;
   _nWJetsTopology = 0.;
   _nWJetsElectronIso = 0.;
   _nWJetsElectronAntiIso = 0.;
   _nWJetsMt = 0.;
   
   
   //set N-1 float counters to zero
   
   _nElecPtNMinus1 = 0.;
   _nElecEtaNMinus1 = 0.;   
   _nMEtNMinus1 = 0.;
   _nBTagsNMinus1 = 0.;
   _nElecMuonCosDPhiNMinus1 = 0.;
   _nElecMuonPZetaNMinus1 = 0.;
   _nElecMuonPZetaVisNMinus1 = 0.;
   _nJetSumEtNMinus1 = 0.;
   _nLeptonMetCosDPhiNMinus1 = 0.;
   _nElecMuonMetMassNMinus1 = 0.;
   _nNJetsNMinus1 = 0.;
   _nElecMuonMassNMinus1 = 0.;
   _nMuonMetMtNMinus1 = 0.;   
   _nElecMuonChargeNMinus1 = 0.;
   _nElecTrkChiSqrdNMinus1 = 0.;
   
   _nMuonPtNMinus1 = 0.;
   _nMuonENMinus1 = 0.;
   _nMuonEtaNMinus1 = 0.;
   _nMuonDxyVtxNMinus1 = 0.;
   _nMuonDzVtxNMinus1 = 0.;
   _nMuonDzBSNMinus1 = 0.;
   _nMuonIsGlobalMuonNMinus1 = 0.;
   _nMuonMatchedNMinus1 = 0.;
   _nMuonNormChiSqrdNMinus1 = 0.;
   _nMuonTrkChiSqrdNMinus1 = 0.;
   _nMuonMatchedStationsNMinus1 = 0.;
   _nMuonValidHitsNMinus1 = 0.;
   _nMuonPixelHitsNMinus1 = 0.;
   _nMuonTrkLayersWithHitsNMinus1 = 0.;
   _nMuonTrkLayersWithHits = 0.;
   _nMuonIsoDR03NoDBNMinus1 = 0.;
   _nMuonIsoDR03DBNMinus1 = 0.;
   _nMuonIsoDR04NoDBNMinus1 = 0.;
   _nMuonIsoDR04DBNMinus1 = 0.;   
   _nMuonIsoNMinus1 = 0.;
   _nMuonIsHighPtMuonNMinus1 = 0.;
   _nMuonDptOverPtNMinus1 = 0.;
   
   _nElecMatchedNMinus1 = 0.;
   _nHeepPassedPtNMinus1 = 0.;
   _nHeepPassedDetEtaNMinus1 = 0.;
   _nHeepPassedCrackNMinus1 = 0.;
   _nHeepPassedDEtaInNMinus1 = 0.;
   _nHeepPassedDPhiInNMinus1 = 0.;
   _nHeepPassedHademNMinus1 = 0.;
   _nHeepPassedSigmaIEtaIEtaNMinus1 = 0.;
   _nHeepPassed2by5Over5By5NMinus1 = 0.;
   _nHeepPassedEcalHad1IsoNMinus1 = 0.;
   _nHeepPassedTrkIsoNMinus1 = 0.;
   _nHeepPassedEcalDrivenNMinus1 = 0.;
   _nElecEtNMinus1 = 0.; 
   _nElecRelIsoNMinus1 = 0.;    
   
   _nElecIpOverElecIpErrorNMinus1 = 0.;
   _nMuonIpOverMuonIpErrorNMinus1 = 0.;     
   _nCombinedIpOverIpErrorNMinus1 = 0.;  
   _nElecMuonDCAOverError2DNMinus1 = 0.;
   _nElecMuonDCAOverError3DNMinus1 = 0.;
}

EMuAnalysis::~EMuAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EMuAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EMuAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EMuAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PDFWeights = 0;
   jetPt = 0;
   jetEt = 0;
   jetE = 0;
   jetEta = 0;
   jetPhi = 0;
   muonMatched = 0;
   muonMotherId = 0;
   muonPt = 0;
   muonPx = 0;
   muonPy = 0;
   muonPz = 0;
   muonE = 0;
   muonEta = 0;
   muonPhi = 0;
   muonCharge = 0;
   muonDxyVtx = 0;
   muonDxyBS = 0;
   muonDxyBSError = 0;
   muonDxyError = 0;
   muonDzVtx = 0;
   muonDzBS = 0;
   muonDzError = 0;
   muonVx = 0;
   muonVy = 0;
   muonVz = 0;
   muonValidHits = 0;
   muonNormChiSqrd = 0;
   muonTrkChiSqrd = 0;
   muonTrkNdof = 0;
   muonChambersValidHits = 0;
   muonMatchedStations = 0;
   muonPixelHits = 0;
   muonTrkLayersWithHits = 0;
   muonPFIsoDR03SumChargedHadronPt = 0;
   muonPFIsoDR03SumChargedParticlePt = 0;
   muonPFIsoDR03SumNeutralHadronPt = 0;
   muonPFIsoDR03SumPhotonPt = 0;
   muonPFIsoDR03SumPUPt = 0;
   muonPFIsoDR04SumChargedHadronPt = 0;
   muonPFIsoDR04SumChargedParticlePt = 0;
   muonPFIsoDR04SumNeutralHadronPt = 0;
   muonPFIsoDR04SumPhotonPt = 0;
   muonPFIsoDR04SumPUPt = 0;
   muonCaloComp = 0;
   muonIsGlobalMuon = 0;
   muonisTrackerMuon = 0;
   muonIsGlobalPromptTight = 0;
   muonIsTMLastStationLoose = 0;
   muonIsTMLastStationTight = 0;
   muonIsTM2DCompatibilityLoose = 0;
   muonIsTM2DCompatibilityTight = 0;
   muonPionVeto = 0;
   muonIsHighPtMuon = 0;
   muonTrkPt = 0;
   muonTrkPtError = 0;
   elecMatched = 0;
   eMotherId = 0;
   elecE = 0;
   elecEt = 0;
   elecPt = 0;
   elecPx = 0;
   elecPy = 0;
   elecPz = 0;
   elecCharge = 0;
   elecEta = 0;
   elecSCEta = 0;
   elecPhi = 0;
   elecSigmaEtaEta = 0;
   elecSigmaIEtaIEta = 0;
   elecTrkChiSqrd = 0;
   elecTrkNdof = 0;
   elecEOverP = 0;
   elecHOverEm = 0;
   elecDeltaPhiIn = 0;
   elecDeltaEtaIn = 0;
   elecFBrem = 0;
   elecEOverInvP = 0;
   elecSCE1x5 = 0;
   elecSCE2x5 = 0;
   elecSCE5x5 = 0;
   elecSCE1x5Over5x5 = 0;
   elecSCE2x5MaxOver5x5 = 0;
   elecTrkIso = 0;
   elecEcalIso = 0;
   elecHcalIso = 0;
   elecHcalIsoDepth1 = 0;
   elecHcalIsoDepth2 = 0;
   elecEcalHcalIsoDepth1 = 0;
   gsfEcalIso = 0;
   gsfHcalIso = 0;
   gsfEt = 0;
   gsfEnergy = 0;
   gsfCaloEnergy = 0;
   elecPFIsoDR03SumChargedHadronPt = 0;
   elecPFIsoDR03SumNeutralHadronPt = 0;
   elecPFIsoDR03SumPhotonPt = 0;
   elecGSFIp = 0;
   elecGSFIpError = 0;
   elecGSFBSIp = 0;
   elecGSFBSIpError = 0;   
   elecCTFIp = 0;
   elecCTFIpError = 0;
   elecMuonDistanceOfClosestApproach = 0;
   elecMuonDCACrossingPoint = 0;
   elecMuonTwoTrkMinDist = 0;
   elecMuonTwoTrkCrossingPoint = 0;
   eMuDis3D = 0;
   eMuDis2D = 0;
   elecMuonDCA3DError = 0;
   elecMuonDCA2DError = 0;
   elecElecMuonMass = 0;
   muonMuonElecMass = 0;
   
   muonMuonMass = 0;
   elecElecMass = 0;
   tauEnergy = 0;
   zPrimeEnergy = 0;
   zPrimeVtx_x = 0;
   zPrimeVtx_y = 0;
   zPrimeVtx_z = 0;
   
   threeLeptonMuonPt = 0;
   threeLeptonOtherMuonPt = 0;
   threeLeptonMuonEta = 0;
   threeLeptonOtherMuonEta = 0;
   threeLeptonIsGlobalMuon = 0;
   threeLeptonOtherMuonIsGlobalMuon = 0;
   threeLeptonMuonDxyVtx = 0;
   threeLeptonOtherMuonDxy = 0;
   threeLeptonMuonDzVtx = 0;
   threeLeptonOtherMuonDzVtx = 0;
   threeLeptonMuonChamberHits = 0;
   threeLeptonOtherMuonChamberHits = 0;
   threeLeptonMuonPixelHits = 0;
   threeLeptonOtherMuonPixelHits = 0;
   threeLeptonMuonMatchedStations = 0;
   threeLeptonOtherMuonMatchedStations = 0;
   threeLeptonMuonTrkLayersWithHits = 0;
   threeLeptonOtherMuonTrkLayersWithHits = 0;
   threeLeptonMuonTrkPt = 0;
   threeLeptonOtherMuonTrkPt = 0;
   threeLeptonMuonTrkPtError = 0;
   threeLeptonOtherMuonTrkPtError = 0;
   threeLeptonMuonPFIsoDR04SumChargedHadronPt = 0;
   threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt = 0;
   threeLeptonMuonPFIsoDR04SumCNeutralHadronPt = 0;
   threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt = 0;	    
   threeLeptonMuonPFIsoDR04SumCPhotonPt = 0;
   threeLeptonOtherMuonPFIsoDR04SumPhotonPt = 0;	    
   threeLeptonMuonPFIsoDR04SumPUPt = 0;
   threeLeptonOtherMuonPFIsoDR04SumPUPt = 0;
   threeLeptonMuMuDCAOverError = 0;

   threeLeptonElecEt = 0;
   threeLeptonOtherElecEt = 0;
   threeLeptonElecSCEta = 0;
   threeLeptonOtherElecSCEta = 0;
   threeLeptonElecEcalDriven = 0;
   threeLeptonOtherElecEcalDriven = 0;
   threeLeptonElecDEtaIn = 0;
   threeLeptonOtherElecDEtaIn = 0;
   threeLeptonElecDPhiIn = 0;
   threeLeptonOtherElecDPhiIn = 0;
   threeLeptonElecHadem = 0;
   threeLeptonOtherElecHadem = 0;
   threeLeptonElecSigmaIEtaIEta = 0;
   threeLeptonOtherElecSigmaIEtaIEta = 0;
   threeLeptonElec2by5Over5By5 = 0;
   threeLeptonOtherElec2by5Over5By5 = 0;
   threeLeptonElecEcalHad1Iso = 0;
   threeLeptonOtherElecEcalHad1Iso = 0;
   threeLeptonElecTrkIso = 0;
   threeLeptonOtherElecTrkIso = 0;
   threeLeptonElecMissingHits = 0;
   threeLeptonOtherElecMissingHits = 0;
   threeLeptonElecDxy = 0;
   threeLeptonOtherElecDxy = 0;
   elecElecDCAOverError = 0;      
   
   
   elecMissingHits = 0;
   elecClass = 0;
   elecIsEcalDriven = 0;
   elecEcalDriven = 0;
   elecIsTrkDriven = 0;
   elecInEE = 0;
   elecInEB = 0;
   elecInGap = 0;
   heepPassedEt = 0;
   heepPassedPt = 0;
   heepPassedDetEta = 0;
   heepPassedCrack = 0;
   heepPassedDEtaIn = 0;
   heepPassedDPhiIn = 0;
   heepPassedHadem = 0;
   heepPassedSigmaIEtaIEta = 0;
   heepPassed2by5Over5By5 = 0;
   heepPassedEcalHad1Iso = 0;
   heepPassedHad2Iso = 0;
   heepPassedTrkIso = 0;
   heepPassedEcalDriven = 0;
   heepPassedMissingHits = 0;
   heepPassedElecDxy = 0;
   heepPassedAllCuts = 0;
   deltaPtx = 0;
   deltaPty = 0;
   deltaPtz = 0;
   elecMVAStatus = 0;
   elecMVAOut = 0;
   elecMuononMass = 0;
   elecMuonCosDPhi = 0;
   elecMuonDeltaR = 0;
   elecMuonMetMass = 0;
   elecMuonDeltaPtMass = 0;
   elecMuonCollMass = 0;
   elecMEtMt = 0;
   muonMEtMt = 0;
   elecMuonPZetaVis = 0;
   elecMuonPZeta = 0;
   nBtagsHiEffTrkCnt = 0;
   nBtagsHiPurityTrkCnt = 0;
   nBTagsHiEffSimpleSecVtx = 0;
   nBTagsHiPuritySimpleSecVtx = 0;
   nBTagsCombSecVtxLWP = 0;
   jetSumEt = 0;
   jetMETSumEt = 0;
   extraTrkPtSum = 0;
   leptonMetCosDphi = 0;
   leptonDeltaPtCosDphi = 0;
   nJets = 0;
   BSx0 = 0;
   BSy0 = 0;
   BSz0 = 0;
   beamSpot_x0 = 0;
   beamSpot_y0 = 0;
   beamSpot_z0 = 0;
   muonVertex_x0 = 0;
   muonVertex_y0 = 0;
   muonVertex_z0 = 0;
   muonDecayLength_x = 0;
   muonDecayLength_y = 0;
   muonDecayLength_z = 0;
   elecVertex_x0 = 0;
   elecVertex_y0 = 0;
   elecVertex_z0 = 0;
   elecDecayLength_x = 0;
   elecDecayLength_y = 0;
   elecDecayLength_z = 0;   
   elecVx = 0;
   elecVy = 0;
   elecVz = 0;
   
   primaryVertex_x = 0;
   primaryVertex_y = 0;
   primaryVertex_z = 0;
   primaryVertex_xError = 0;
   primaryVertex_yError = 0;
   primaryVertex_zError = 0;
   genTauVtx_x = 0;
   genTauVtx_y = 0;
   genTauVtx_z = 0;
   otherGenTauVtx_x = 0;
   otherGenTauVtx_y = 0;
   otherGenTauVtx_z = 0;
   genTauPx = 0;  
   genTauPy = 0;  
   genTauPz = 0;  
   otherGenTauPx = 0;  
   otherGenTauPy = 0;  
   otherGenTauPz = 0;  
   genTauCharge = 0;  
   otherGenTauCharge = 0;  

 
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("InTimePU", &InTimePU, &b_InTimePU);
   fChain->SetBranchAddress("OutTimePlusPU", &OutTimePlusPU, &b_OutTimePlusPU);
   fChain->SetBranchAddress("OutTimeMinusPU", &OutTimeMinusPU, &b_OutTimeMinusPU);
   fChain->SetBranchAddress("passedHLTMu5Elec17", &passedHLTMu5Elec17, &b_passedHLTMu5Elec17);
   fChain->SetBranchAddress("passedHLTMu11Elec8", &passedHLTMu11Elec8, &b_passedHLTMu11Elec8);
   fChain->SetBranchAddress("passedHLTMu8Elec17", &passedHLTMu8Elec17, &b_passedHLTMu8Elec17);
   fChain->SetBranchAddress("passedHLTMu17Elec8", &passedHLTMu17Elec8, &b_passedHLTMu17Elec8);
   fChain->SetBranchAddress("passedHLTIsoMu24", &passedHLTIsoMu24, &b_passedHLTIsoMu24);
   fChain->SetBranchAddress("PDFWeights", &PDFWeights, &b_PDFWeights);
   fChain->SetBranchAddress("ISRGluonWeight", &ISRGluonWeight, &b_ISRGluonWeight);
   fChain->SetBranchAddress("ISRGammaWeight", &ISRGammaWeight, &b_ISRGammaWeight);
   fChain->SetBranchAddress("FSRWeight", &FSRWeight, &b_FSRWeight);
   fChain->SetBranchAddress("mEt", &mEt, &b_mEt);
   fChain->SetBranchAddress("mEtPx", &mEtPx, &b_mEtPx);
   fChain->SetBranchAddress("mEtPy", &mEtPy, &b_mEtPy);
   fChain->SetBranchAddress("mEtPz", &mEtPz, &b_mEtPz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("eMuMuOneVertex", &eMuMuOneVertex, &b_eMuMuOneVertex);
   fChain->SetBranchAddress("muEEOneVertex", &muEEOneVertex, &b_muEEOneVertex);
   fChain->SetBranchAddress("elecElecMuonMass", &elecElecMuonMass, &b_elecElecMuonMass);
   fChain->SetBranchAddress("muonMuonElecMass", &muonMuonElecMass, &b_muonMuonElecMass);
   
   fChain->SetBranchAddress("muonMuonMass", &muonMuonMass, &b_muonMuonMass);
   fChain->SetBranchAddress("elecElecMass", &elecElecMass, &b_elecElecMass);
   fChain->SetBranchAddress("tauEnergy", &tauEnergy, &b_tauEnergy);
   fChain->SetBranchAddress("zPrimeEnergy", &zPrimeEnergy, &b_zPrimeEnergy);
   fChain->SetBranchAddress("zPrimeVtx_x", &zPrimeVtx_x, &b_zPrimeVtx_x);
   fChain->SetBranchAddress("zPrimeVtx_y", &zPrimeVtx_y, &b_zPrimeVtx_y);
   fChain->SetBranchAddress("zPrimeVtx_z", &zPrimeVtx_z, &b_zPrimeVtx_z);
   
   fChain->SetBranchAddress("threeLeptonMuonPt", &threeLeptonMuonPt, &b_threeLeptonMuonPt);
   fChain->SetBranchAddress("threeLeptonOtherMuonPt", &threeLeptonOtherMuonPt, &b_threeLeptonOtherMuonPt);
   fChain->SetBranchAddress("threeLeptonMuonEta", &threeLeptonMuonEta, &b_threeLeptonMuonEta);
   fChain->SetBranchAddress("threeLeptonOtherMuonEta", &threeLeptonOtherMuonEta, &b_threeLeptonOtherMuonEta);
   fChain->SetBranchAddress("threeLeptonIsGlobalMuon", &threeLeptonIsGlobalMuon, &b_threeLeptonIsGlobalMuon);
   fChain->SetBranchAddress("threeLeptonOtherMuonIsGlobalMuon", &threeLeptonOtherMuonIsGlobalMuon, &b_threeLeptonOtherMuonIsGlobalMuon);
   fChain->SetBranchAddress("threeLeptonMuonDxyVtx", &threeLeptonMuonDxyVtx, &b_threeLeptonMuonDxyVtx);
   fChain->SetBranchAddress("threeLeptonOtherMuonDxy", &threeLeptonOtherMuonDxy, &b_threeLeptonOtherMuonDxy);
   fChain->SetBranchAddress("threeLeptonMuonDzVtx", &threeLeptonMuonDzVtx, &b_threeLeptonMuonDzVtx);
   fChain->SetBranchAddress("threeLeptonOtherMuonDzVtx", &threeLeptonOtherMuonDzVtx, &b_threeLeptonOtherMuonDzVtx);
   fChain->SetBranchAddress("threeLeptonMuonChamberHits", &threeLeptonMuonChamberHits, &b_threeLeptonMuonChamberHits);
   fChain->SetBranchAddress("threeLeptonOtherMuonChamberHits", &threeLeptonOtherMuonChamberHits, &b_threeLeptonOtherMuonChamberHits);
   fChain->SetBranchAddress("threeLeptonMuonPixelHits", &threeLeptonMuonPixelHits, &b_threeLeptonMuonPixelHits);
   fChain->SetBranchAddress("threeLeptonOtherMuonPixelHits", &threeLeptonOtherMuonPixelHits, &b_threeLeptonOtherMuonPixelHits);
   fChain->SetBranchAddress("threeLeptonMuonMatchedStations", &threeLeptonMuonMatchedStations, &b_threeLeptonMuonMatchedStations);
   fChain->SetBranchAddress("threeLeptonOtherMuonMatchedStations", &threeLeptonOtherMuonMatchedStations, &b_threeLeptonOtherMuonMatchedStations);
   fChain->SetBranchAddress("threeLeptonMuonTrkLayersWithHits", &threeLeptonMuonTrkLayersWithHits, &b_threeLeptonMuonTrkLayersWithHits);
   fChain->SetBranchAddress("threeLeptonOtherMuonTrkLayersWithHits", &threeLeptonOtherMuonTrkLayersWithHits, &b_threeLeptonOtherMuonTrkLayersWithHits);
   fChain->SetBranchAddress("threeLeptonMuonTrkPt", &threeLeptonMuonTrkPt, &b_threeLeptonMuonTrkPt);
   fChain->SetBranchAddress("threeLeptonOtherMuonTrkPt", &threeLeptonOtherMuonTrkPt, &b_threeLeptonOtherMuonTrkPt);
   fChain->SetBranchAddress("threeLeptonMuonTrkPtError", &threeLeptonMuonTrkPtError, &b_threeLeptonMuonTrkPtError);
   fChain->SetBranchAddress("threeLeptonOtherMuonTrkPtError", &threeLeptonOtherMuonTrkPtError, &b_threeLeptonOtherMuonTrkPtError);
   fChain->SetBranchAddress("threeLeptonMuonPFIsoDR04SumChargedHadronPt", &threeLeptonMuonPFIsoDR04SumChargedHadronPt, &b_threeLeptonMuonPFIsoDR04SumChargedHadronPt);
   fChain->SetBranchAddress("threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt", &threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt, &b_threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt);
   fChain->SetBranchAddress("threeLeptonMuonPFIsoDR04SumCNeutralHadronPt", &threeLeptonMuonPFIsoDR04SumCNeutralHadronPt, &b_threeLeptonMuonPFIsoDR04SumCNeutralHadronPt);
   fChain->SetBranchAddress("threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt", &threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt, &b_threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt);	    
   fChain->SetBranchAddress("threeLeptonMuonPFIsoDR04SumCPhotonPt", &threeLeptonMuonPFIsoDR04SumCPhotonPt, &b_threeLeptonMuonPFIsoDR04SumCPhotonPt);
   fChain->SetBranchAddress("threeLeptonOtherMuonPFIsoDR04SumPhotonPt", &threeLeptonOtherMuonPFIsoDR04SumPhotonPt , &b_threeLeptonOtherMuonPFIsoDR04SumPhotonPt);
   fChain->SetBranchAddress("threeLeptonMuonPFIsoDR04SumPUPt", &threeLeptonMuonPFIsoDR04SumPUPt, &b_threeLeptonMuonPFIsoDR04SumPUPt);
   fChain->SetBranchAddress("threeLeptonOtherMuonPFIsoDR04SumPUPt", &threeLeptonOtherMuonPFIsoDR04SumPUPt, &b_threeLeptonOtherMuonPFIsoDR04SumPUPt);
   fChain->SetBranchAddress("threeLeptonMuMuDCAOverError", &threeLeptonMuMuDCAOverError, &b_threeLeptonMuMuDCAOverError);

   fChain->SetBranchAddress("threeLeptonElecEt", &threeLeptonElecEt, &b_threeLeptonElecEt);
   fChain->SetBranchAddress("threeLeptonOtherElecEt", &threeLeptonOtherElecEt, &b_threeLeptonOtherElecEt);
   fChain->SetBranchAddress("threeLeptonElecSCEta", &threeLeptonElecSCEta, &b_threeLeptonElecSCEta);
   fChain->SetBranchAddress("threeLeptonOtherElecSCEta", &threeLeptonOtherElecSCEta, &b_threeLeptonOtherElecSCEta);
   fChain->SetBranchAddress("threeLeptonElecEcalDriven", &threeLeptonElecEcalDriven, &b_threeLeptonElecEcalDriven);
   fChain->SetBranchAddress("threeLeptonOtherElecEcalDriven", &threeLeptonOtherElecEcalDriven, &b_threeLeptonOtherElecEcalDriven);
   fChain->SetBranchAddress("threeLeptonElecDEtaIn", &threeLeptonElecDEtaIn, &b_threeLeptonElecDEtaIn);
   fChain->SetBranchAddress("threeLeptonOtherElecDEtaIn", &threeLeptonOtherElecDEtaIn, &b_threeLeptonOtherElecDEtaIn);
   fChain->SetBranchAddress("threeLeptonElecDPhiIn", &threeLeptonElecDPhiIn, &b_threeLeptonElecDPhiIn);
   fChain->SetBranchAddress("threeLeptonOtherElecDPhiIn", &threeLeptonOtherElecDPhiIn, &b_threeLeptonOtherElecDPhiIn);
   fChain->SetBranchAddress("threeLeptonElecHadem", &threeLeptonElecHadem, &b_threeLeptonElecHadem);
   fChain->SetBranchAddress("threeLeptonOtherElecHadem", &threeLeptonOtherElecHadem, &b_threeLeptonOtherElecHadem);
   fChain->SetBranchAddress("threeLeptonElecSigmaIEtaIEta", &threeLeptonElecSigmaIEtaIEta, &b_threeLeptonElecSigmaIEtaIEta);
   fChain->SetBranchAddress("threeLeptonOtherElecSigmaIEtaIEta", &threeLeptonOtherElecSigmaIEtaIEta, &b_threeLeptonOtherElecSigmaIEtaIEta);
   fChain->SetBranchAddress("threeLeptonElec2by5Over5By5", &threeLeptonElec2by5Over5By5, &b_threeLeptonElec2by5Over5By5);
   fChain->SetBranchAddress("threeLeptonOtherElec2by5Over5By5", &threeLeptonOtherElec2by5Over5By5, &b_threeLeptonOtherElec2by5Over5By5);
   fChain->SetBranchAddress("threeLeptonElecEcalHad1Iso", &threeLeptonElecEcalHad1Iso, &b_threeLeptonElecEcalHad1Iso);
   fChain->SetBranchAddress("threeLeptonOtherElecEcalHad1Iso", &threeLeptonOtherElecEcalHad1Iso, &b_threeLeptonOtherElecEcalHad1Iso);
   fChain->SetBranchAddress("threeLeptonElecTrkIso", &threeLeptonElecTrkIso, &b_threeLeptonElecTrkIso);
   fChain->SetBranchAddress("threeLeptonOtherElecTrkIso", &threeLeptonOtherElecTrkIso, &b_threeLeptonOtherElecTrkIso);
   fChain->SetBranchAddress("threeLeptonElecMissingHits", &threeLeptonElecMissingHits, &b_threeLeptonElecMissingHits);
   fChain->SetBranchAddress("threeLeptonOtherElecMissingHits", &threeLeptonOtherElecMissingHits, &b_threeLeptonOtherElecMissingHits);
   fChain->SetBranchAddress("threeLeptonElecDxy", &threeLeptonElecDxy, &b_threeLeptonElecDxy);
   fChain->SetBranchAddress("threeLeptonOtherElecDxy", &threeLeptonOtherElecDxy, &b_threeLeptonOtherElecDxy);
   fChain->SetBranchAddress("elecElecDCAOverError", & elecElecDCAOverError, &b_elecElecDCAOverError);
   
   
   fChain->SetBranchAddress("jetBtagHiDiscByTrkCntHiEff", &jetBtagHiDiscByTrkCntHiEff, &b_jetBtagHiDiscByTrkCntHiEff);
   fChain->SetBranchAddress("jetBtagHiDiscByTrkCntHiPurity", &jetBtagHiDiscByTrkCntHiPurity, &b_jetBtagHiDiscByTrkCntHiPurity);
   fChain->SetBranchAddress("jetBtagHiDiscBySimpleSecVtxHiEff", &jetBtagHiDiscBySimpleSecVtxHiEff, &b_jetBtagHiDiscBySimpleSecVtxHiEff);
   fChain->SetBranchAddress("jetBtagHiDiscBySimpleSecVtxHiPurity", &jetBtagHiDiscBySimpleSecVtxHiPurity, &b_jetBtagHiDiscBySimpleSecVtxHiPurity);
   fChain->SetBranchAddress("jetBtagHiDiscByCombSecVtx", &jetBtagHiDiscByCombSecVtx, &b_jetBtagHiDiscByCombSecVtx);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEt", &jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("muonMatched", &muonMatched, &b_muonMatched);
   fChain->SetBranchAddress("muonMotherId", &muonMotherId, &b_muonMotherId);
   fChain->SetBranchAddress("muonPt", &muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonPx", &muonPx, &b_muonPx);
   fChain->SetBranchAddress("muonPy", &muonPy, &b_muonPy);
   fChain->SetBranchAddress("muonPz", &muonPz, &b_muonPz);
   fChain->SetBranchAddress("muonE", &muonE, &b_muonE);
   fChain->SetBranchAddress("muonEta", &muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonPhi", &muonPhi, &b_muonPhi);
   fChain->SetBranchAddress("muonCharge", &muonCharge, &b_muonCharge);
   fChain->SetBranchAddress("muonDxyVtx", &muonDxyVtx, &b_muonDxyVtx);
   fChain->SetBranchAddress("muonDxyBS", &muonDxyBS, &b_muonDxyBS);
   fChain->SetBranchAddress("muonDxyBSError", &muonDxyBSError, &b_muonDxyBSError);
   fChain->SetBranchAddress("muonDxyError", &muonDxyError, &b_muonDxyError);
   fChain->SetBranchAddress("muonDzVtx", &muonDzVtx, &b_muonDzVtx);
   fChain->SetBranchAddress("muonDzBS", &muonDzBS, &b_muonDzBS);
   fChain->SetBranchAddress("muonDzError", &muonDzError, &b_muonDzError);
   fChain->SetBranchAddress("muonVx", &muonVx, &b_muonVx);
   fChain->SetBranchAddress("muonVy", &muonVy, &b_muonVy);
   fChain->SetBranchAddress("muonVz", &muonVz, &b_muonVz);
   fChain->SetBranchAddress("muonValidHits", &muonValidHits, &b_muonValidHits);
   fChain->SetBranchAddress("muonNormChiSqrd", &muonNormChiSqrd, &b_muonNormChiSqrd);
   fChain->SetBranchAddress("muonTrkChiSqrd", &muonTrkChiSqrd, &b_muonTrkChiSqrd);
   fChain->SetBranchAddress("muonTrkNdof", &muonTrkNdof, &b_muonTrkNdof);
   fChain->SetBranchAddress("muonChambersValidHits", &muonChambersValidHits, &b_muonChambersValidHits);
   fChain->SetBranchAddress("muonMatchedStations", &muonMatchedStations, &b_muonMatchedStations);
   fChain->SetBranchAddress("muonPixelHits", &muonPixelHits, &b_muonPixelHits);
   fChain->SetBranchAddress("muonTrkLayersWithHits", &muonTrkLayersWithHits, &b_muonTrkLayersWithHits);
   fChain->SetBranchAddress("muonPFIsoDR03SumChargedHadronPt", &muonPFIsoDR03SumChargedHadronPt, &b_muonPFIsoDR03SumChargedHadronPt);
   fChain->SetBranchAddress("muonPFIsoDR03SumChargedParticlePt", &muonPFIsoDR03SumChargedParticlePt, &b_muonPFIsoDR03SumChargedParticlePt);
   fChain->SetBranchAddress("muonPFIsoDR03SumNeutralHadronPt", &muonPFIsoDR03SumNeutralHadronPt, &b_muonPFIsoDR03SumNeutralHadronPt);
   fChain->SetBranchAddress("muonPFIsoDR03SumPhotonPt", &muonPFIsoDR03SumPhotonPt, &b_muonPFIsoDR03SumPhotonPt);
   fChain->SetBranchAddress("muonPFIsoDR03SumPUPt", &muonPFIsoDR03SumPUPt, &b_muonPFIsoDR03SumPUPt);
   fChain->SetBranchAddress("muonPFIsoDR04SumChargedHadronPt", &muonPFIsoDR04SumChargedHadronPt, &b_muonPFIsoDR04SumChargedHadronPt);
   fChain->SetBranchAddress("muonPFIsoDR04SumChargedParticlePt", &muonPFIsoDR04SumChargedParticlePt, &b_muonPFIsoDR04SumChargedParticlePt);
   fChain->SetBranchAddress("muonPFIsoDR04SumNeutralHadronPt", &muonPFIsoDR04SumNeutralHadronPt, &b_muonPFIsoDR04SumNeutralHadronPt);
   fChain->SetBranchAddress("muonPFIsoDR04SumPhotonPt", &muonPFIsoDR04SumPhotonPt, &b_muonPFIsoDR04SumPhotonPt);
   fChain->SetBranchAddress("muonPFIsoDR04SumPUPt", &muonPFIsoDR04SumPUPt, &b_muonPFIsoDR04SumPUPt);
   fChain->SetBranchAddress("muonCaloComp", &muonCaloComp, &b_muonCaloComp);
   fChain->SetBranchAddress("muonIsGlobalMuon", &muonIsGlobalMuon, &b_muonIsGlobalMuon);
   fChain->SetBranchAddress("muonisTrackerMuon", &muonisTrackerMuon, &b_muonisTrackerMuon);
   fChain->SetBranchAddress("muonIsGlobalPromptTight", &muonIsGlobalPromptTight, &b_muonIsGlobalPromptTight);
   fChain->SetBranchAddress("muonIsTMLastStationLoose", &muonIsTMLastStationLoose, &b_muonIsTMLastStationLoose);
   fChain->SetBranchAddress("muonIsTMLastStationTight", &muonIsTMLastStationTight, &b_muonIsTMLastStationTight);
   fChain->SetBranchAddress("muonIsTM2DCompatibilityLoose", &muonIsTM2DCompatibilityLoose, &b_muonIsTM2DCompatibilityLoose);
   fChain->SetBranchAddress("muonIsTM2DCompatibilityTight", &muonIsTM2DCompatibilityTight, &b_muonIsTM2DCompatibilityTight);
   fChain->SetBranchAddress("muonPionVeto", &muonPionVeto, &b_muonPionVeto);
   fChain->SetBranchAddress("muonIsHighPtMuon", &muonIsHighPtMuon, &b_muonIsHighPtMuon);
   fChain->SetBranchAddress("muonTrkPt", &muonTrkPt, &b_muonTrkPt);
   fChain->SetBranchAddress("muonTrkPtError", &muonTrkPtError, &b_muonTrkPtError);
   fChain->SetBranchAddress("elecMatched", &elecMatched, &b_elecMatched);
   fChain->SetBranchAddress("eMotherId", &eMotherId, &b_eMotherId);
   fChain->SetBranchAddress("elecE", &elecE, &b_elecE);
   fChain->SetBranchAddress("elecEt", &elecEt, &b_elecEt);
   fChain->SetBranchAddress("elecPt", &elecPt, &b_elecPt);
   fChain->SetBranchAddress("elecPx", &elecPx, &b_elecPx);
   fChain->SetBranchAddress("elecPy", &elecPy, &b_elecPy);
   fChain->SetBranchAddress("elecPz", &elecPz, &b_elecPz);
   fChain->SetBranchAddress("elecCharge", &elecCharge, &b_elecCharge);
   fChain->SetBranchAddress("elecEta", &elecEta, &b_elecEta);
   fChain->SetBranchAddress("elecSCEta", &elecSCEta, &b_elecSCEta);
   fChain->SetBranchAddress("elecPhi", &elecPhi, &b_elecPhi);
   fChain->SetBranchAddress("elecSigmaEtaEta", &elecSigmaEtaEta, &b_elecSigmaEtaEta);
   fChain->SetBranchAddress("elecSigmaIEtaIEta", &elecSigmaIEtaIEta, &b_elecSigmaIEtaIEta);
   fChain->SetBranchAddress("elecTrkChiSqrd", &elecTrkChiSqrd, &b_elecTrkChiSqrd);
   fChain->SetBranchAddress("elecTrkNdof", &elecTrkNdof, &b_elecTrkNdof);
   fChain->SetBranchAddress("elecEOverP", &elecEOverP, &b_elecEOverP);
   fChain->SetBranchAddress("elecHOverEm", &elecHOverEm, &b_elecHOverEm);
   fChain->SetBranchAddress("elecDeltaPhiIn", &elecDeltaPhiIn, &b_elecDeltaPhiIn);
   fChain->SetBranchAddress("elecDeltaEtaIn", &elecDeltaEtaIn, &b_elecDeltaEtaIn);
   fChain->SetBranchAddress("elecFBrem", &elecFBrem, &b_elecFBrem);
   fChain->SetBranchAddress("elecEOverInvP", &elecEOverInvP, &b_elecEOverInvP);
   fChain->SetBranchAddress("elecSCE1x5", &elecSCE1x5, &b_elecSCE1x5);
   fChain->SetBranchAddress("elecSCE2x5", &elecSCE2x5, &b_elecSCE2x5);
   fChain->SetBranchAddress("elecSCE5x5", &elecSCE5x5, &b_elecSCE5x5);
   fChain->SetBranchAddress("elecSCE1x5Over5x5", &elecSCE1x5Over5x5, &b_elecSCE1x5Over5x5);
   fChain->SetBranchAddress("elecSCE2x5MaxOver5x5", &elecSCE2x5MaxOver5x5, &b_elecSCE2x5MaxOver5x5);
   fChain->SetBranchAddress("elecTrkIso", &elecTrkIso, &b_elecTrkIso);
   fChain->SetBranchAddress("elecEcalIso", &elecEcalIso, &b_elecEcalIso);
   fChain->SetBranchAddress("elecHcalIso", &elecHcalIso, &b_elecHcalIso);
   fChain->SetBranchAddress("elecHcalIsoDepth1", &elecHcalIsoDepth1, &b_elecHcalIsoDepth1);
   fChain->SetBranchAddress("elecHcalIsoDepth2", &elecHcalIsoDepth2, &b_elecHcalIsoDepth2);
   fChain->SetBranchAddress("elecEcalHcalIsoDepth1", &elecEcalHcalIsoDepth1, &b_elecEcalHcalIsoDepth1);
   fChain->SetBranchAddress("gsfEcalIso", &gsfEcalIso, &b_gsfEcalIso);
   fChain->SetBranchAddress("gsfHcalIso", &gsfHcalIso, &b_gsfHcalIso);
   fChain->SetBranchAddress("gsfEt", &gsfEt, &b_gsfEt);
   fChain->SetBranchAddress("gsfEnergy", &gsfEnergy, &b_gsfEnergy);
   fChain->SetBranchAddress("gsfCaloEnergy", &gsfCaloEnergy, &b_gsfCaloEnergy);
   fChain->SetBranchAddress("elecPFIsoDR03SumChargedHadronPt", &elecPFIsoDR03SumChargedHadronPt, &b_elecPFIsoDR03SumChargedHadronPt);
   fChain->SetBranchAddress("elecPFIsoDR03SumNeutralHadronPt", &elecPFIsoDR03SumNeutralHadronPt, &b_elecPFIsoDR03SumNeutralHadronPt);
   fChain->SetBranchAddress("elecPFIsoDR03SumPhotonPt", &elecPFIsoDR03SumPhotonPt, &b_elecPFIsoDR03SumPhotonPt);
   fChain->SetBranchAddress("elecGSFIp", &elecGSFIp, &b_elecGSFIp);
   fChain->SetBranchAddress("elecGSFIpError", &elecGSFIpError, &b_elecGSFIpError);
   fChain->SetBranchAddress("elecGSFBSIp", &elecGSFBSIp, &b_elecGSFBSIp);
   fChain->SetBranchAddress("elecGSFBSIpError", &elecGSFBSIpError, &b_elecGSFBSIpError);   
   fChain->SetBranchAddress("elecCTFIp", &elecCTFIp, &b_elecCTFIp);
   fChain->SetBranchAddress("elecCTFIpError", &elecCTFIpError, &b_elecCTFIpError);
   fChain->SetBranchAddress("elecMuonDistanceOfClosestApproach", &elecMuonDistanceOfClosestApproach, &b_elecMuonDistanceOfClosestApproach);
   fChain->SetBranchAddress("elecMuonDCACrossingPoint", &elecMuonDCACrossingPoint, &b_elecMuonDCACrossingPoint);
   fChain->SetBranchAddress("elecMuonTwoTrkMinDist", &elecMuonTwoTrkMinDist, &b_elecMuonTwoTrkMinDist);
   fChain->SetBranchAddress("elecMuonTwoTrkCrossingPoint", &elecMuonTwoTrkCrossingPoint, &b_elecMuonTwoTrkCrossingPoint);
   fChain->SetBranchAddress("eMuDis3D", &eMuDis3D, &b_eMuDis3D);
   fChain->SetBranchAddress("eMuDis2D", &eMuDis2D, &b_eMuDis2D);
   fChain->SetBranchAddress("elecMuonDCA3DError", &elecMuonDCA3DError, &b_elecMuonDCA3DError);
   fChain->SetBranchAddress("elecMuonDCA2DError", &elecMuonDCA2DError, &b_elecMuonDCA2DError);
   fChain->SetBranchAddress("elecMissingHits", &elecMissingHits, &b_elecMissingHits);
   fChain->SetBranchAddress("elecClass", &elecClass, &b_elecClass);
   fChain->SetBranchAddress("elecIsEcalDriven", &elecIsEcalDriven, &b_elecIsEcalDriven);
   fChain->SetBranchAddress("elecEcalDriven", &elecEcalDriven, &b_elecEcalDriven);
   fChain->SetBranchAddress("elecIsTrkDriven", &elecIsTrkDriven, &b_elecIsTrkDriven);
   fChain->SetBranchAddress("elecInEE", &elecInEE, &b_elecInEE);
   fChain->SetBranchAddress("elecInEB", &elecInEB, &b_elecInEB);
   fChain->SetBranchAddress("elecInGap", &elecInGap, &b_elecInGap);
   fChain->SetBranchAddress("heepPassedEt", &heepPassedEt, &b_heepPassedEt);
   fChain->SetBranchAddress("heepPassedPt", &heepPassedPt, &b_heepPassedPt);
   fChain->SetBranchAddress("heepPassedDetEta", &heepPassedDetEta, &b_heepPassedDetEta);
   fChain->SetBranchAddress("heepPassedCrack", &heepPassedCrack, &b_heepPassedCrack);
   fChain->SetBranchAddress("heepPassedDEtaIn", &heepPassedDEtaIn, &b_heepPassedDEtaIn);
   fChain->SetBranchAddress("heepPassedDPhiIn", &heepPassedDPhiIn, &b_heepPassedDPhiIn);
   fChain->SetBranchAddress("heepPassedHadem", &heepPassedHadem, &b_heepPassedHadem);
   fChain->SetBranchAddress("heepPassedSigmaIEtaIEta", &heepPassedSigmaIEtaIEta, &b_heepPassedSigmaIEtaIEta);
   fChain->SetBranchAddress("heepPassed2by5Over5By5", &heepPassed2by5Over5By5, &b_heepPassed2by5Over5By5);
   fChain->SetBranchAddress("heepPassedEcalHad1Iso", &heepPassedEcalHad1Iso, &b_heepPassedEcalHad1Iso);
   fChain->SetBranchAddress("heepPassedHad2Iso", &heepPassedHad2Iso, &b_heepPassedHad2Iso);
   fChain->SetBranchAddress("heepPassedTrkIso", &heepPassedTrkIso, &b_heepPassedTrkIso);
   fChain->SetBranchAddress("heepPassedEcalDriven", &heepPassedEcalDriven, &b_heepPassedEcalDriven);
   fChain->SetBranchAddress("heepPassedMissingHits", &heepPassedMissingHits, &b_heepPassedMissingHits);
   fChain->SetBranchAddress("heepPassedElecDxy", &heepPassedElecDxy, &b_heepPassedElecDxy);
   fChain->SetBranchAddress("heepPassedAllCuts", &heepPassedAllCuts, &b_heepPassedAllCuts);
   fChain->SetBranchAddress("deltaPtx", &deltaPtx, &b_deltaPtx);
   fChain->SetBranchAddress("deltaPtx", &deltaPty, &b_deltaPty);
   fChain->SetBranchAddress("deltaPtx", &deltaPtz, &b_deltaPtz);   
   fChain->SetBranchAddress("elecMVAStatus", &elecMVAStatus, &b_elecMVAStatus);
   fChain->SetBranchAddress("elecMVAOut", &elecMVAOut, &b_elecMVAOut);
   fChain->SetBranchAddress("elecMuononMass", &elecMuononMass, &b_elecMuononMass);
   fChain->SetBranchAddress("elecMuonCosDPhi", &elecMuonCosDPhi, &b_elecMuonCosDPhi);
   fChain->SetBranchAddress("elecMuonDeltaR", &elecMuonDeltaR, &b_elecMuonDeltaR);
   fChain->SetBranchAddress("elecMuonMetMass", &elecMuonMetMass, &b_elecMuonMetMass);
   fChain->SetBranchAddress("elecMuonDeltaPtMass", &elecMuonDeltaPtMass, &b_elecMuonDeltaPtMass);
   fChain->SetBranchAddress("elecMuonCollMass", &elecMuonCollMass, &b_elecMuonCollMass);
   fChain->SetBranchAddress("elecMEtMt", &elecMEtMt, &b_elecMEtMt);
   fChain->SetBranchAddress("muonMEtMt", &muonMEtMt, &b_muonMEtMt);
   fChain->SetBranchAddress("elecMuonPZetaVis", &elecMuonPZetaVis, &b_elecMuonPZetaVis);
   fChain->SetBranchAddress("elecMuonPZeta", &elecMuonPZeta, &b_elecMuonPZeta);
   fChain->SetBranchAddress("nBtagsHiEffTrkCnt", &nBtagsHiEffTrkCnt, &b_nBtagsHiEffTrkCnt);
   fChain->SetBranchAddress("nBtagsHiPurityTrkCnt", &nBtagsHiPurityTrkCnt, &b_nBtagsHiPurityTrkCnt);
   fChain->SetBranchAddress("nBTagsHiEffSimpleSecVtx", &nBTagsHiEffSimpleSecVtx, &b_nBTagsHiEffSimpleSecVtx);
   fChain->SetBranchAddress("nBTagsHiPuritySimpleSecVtx", &nBTagsHiPuritySimpleSecVtx, &b_nBTagsHiPuritySimpleSecVtx);
   fChain->SetBranchAddress("nBTagsCombSecVtxLWP", &nBTagsCombSecVtxLWP, &b_nBTagsCombSecVtxLWP);
   fChain->SetBranchAddress("jetSumEt", &jetSumEt, &b_jetSumEt);
   fChain->SetBranchAddress("jetMETSumEt", &jetMETSumEt, &b_jetMETSumEt);
   fChain->SetBranchAddress("extraTrkPtSum", &extraTrkPtSum, &b_extraTrkPtSum);
   fChain->SetBranchAddress("leptonMetCosDphi", &leptonMetCosDphi, &b_leptonMetCosDphi);
   fChain->SetBranchAddress("leptonDeltaPtCosDphi", &leptonDeltaPtCosDphi, &b_leptonDeltaPtCosDphi);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("BSx0", &BSx0, &b_BSx0);
   fChain->SetBranchAddress("BSy0", &BSy0, &b_BSy0);
   fChain->SetBranchAddress("BSz0", &BSz0, &b_BSz0);
   fChain->SetBranchAddress("beamSpot_x0", &beamSpot_x0, &b_beamSpot_x0);
   fChain->SetBranchAddress("beamSpot_y0", &beamSpot_y0, &b_beamSpot_y0);
   fChain->SetBranchAddress("beamSpot_z0", &beamSpot_z0, &b_beamSpot_z0);
   fChain->SetBranchAddress("muonVertex_x0", &muonVertex_x0, &b_muonVertex_x0);
   fChain->SetBranchAddress("muonVertex_y0", &muonVertex_y0, &b_muonVertex_y0);
   fChain->SetBranchAddress("muonVertex_z0", &muonVertex_z0, &b_muonVertex_z0);
   fChain->SetBranchAddress("muonDecayLength_x", &muonDecayLength_x, &b_muonDecayLength_x);
   fChain->SetBranchAddress("muonDecayLength_y", &muonDecayLength_y, &b_muonDecayLength_y);
   fChain->SetBranchAddress("muonDecayLength_z", &muonDecayLength_z, &b_muonDecayLength_z);
   fChain->SetBranchAddress("elecVertex_x0", &elecVertex_x0, &b_elecVertex_x0);
   fChain->SetBranchAddress("elecVertex_y0", &elecVertex_y0, &b_elecVertex_y0);
   fChain->SetBranchAddress("elecVertex_z0", &elecVertex_z0, &b_elecVertex_z0);
   fChain->SetBranchAddress("elecDecayLength_x", &elecDecayLength_x, &b_elecDecayLength_x);
   fChain->SetBranchAddress("elecDecayLength_y", &elecDecayLength_y, &b_elecDecayLength_y);
   fChain->SetBranchAddress("elecDecayLength_z", &elecDecayLength_z, &b_elecDecayLength_z);   
   fChain->SetBranchAddress("elecVx", &elecVx, &b_elecVx);
   fChain->SetBranchAddress("elecVy", &elecVy, &b_elecVy);
   fChain->SetBranchAddress("elecVz", &elecVz, &b_elecVz);
   fChain->SetBranchAddress("primaryVertex_x", &primaryVertex_x, &b_primaryVertex_x);
   fChain->SetBranchAddress("primaryVertex_y", &primaryVertex_y, &b_primaryVertex_y);
   fChain->SetBranchAddress("primaryVertex_z", &primaryVertex_z, &b_primaryVertex_z);
   fChain->SetBranchAddress("primaryVertex_xError", &primaryVertex_xError, &b_primaryVertex_xError);
   fChain->SetBranchAddress("primaryVertex_yError", &primaryVertex_yError, &b_primaryVertex_yError);
   fChain->SetBranchAddress("primaryVertex_zError", &primaryVertex_zError, &b_primaryVertex_zError);
   fChain->SetBranchAddress("genTauVtx_x", &genTauVtx_x, &b_genTauVtx_x);
   fChain->SetBranchAddress("genTauVtx_y", &genTauVtx_y, &b_genTauVtx_y);
   fChain->SetBranchAddress("genTauVtx_z", &genTauVtx_z, &b_genTauVtx_z);
   fChain->SetBranchAddress("otherGenTauVtx_x", &otherGenTauVtx_x, &b_otherGenTauVtx_x);
   fChain->SetBranchAddress("otherGenTauVtx_y", &otherGenTauVtx_y, &b_otherGenTauVtx_y);
   fChain->SetBranchAddress("otherGenTauVtx_z", &otherGenTauVtx_z, &b_otherGenTauVtx_z);  
   fChain->SetBranchAddress("genTauPx", &genTauPx, &b_genTauPx);
   fChain->SetBranchAddress("genTauPy", &genTauPy, &b_genTauPy);
   fChain->SetBranchAddress("genTauPz", &genTauPz, &b_genTauPz);
   fChain->SetBranchAddress("otherGenTauPx", &otherGenTauPx, &b_otherGenTauPx);
   fChain->SetBranchAddress("otherGenTauPy", &otherGenTauPy, &b_otherGenTauPy);
   fChain->SetBranchAddress("otherGenTauPz", &otherGenTauPz, &b_otherGenTauPz);
   fChain->SetBranchAddress("genTauCharge", &genTauCharge, &b_genTauCharge);
   fChain->SetBranchAddress("otherGenTauCharge", &otherGenTauCharge, &b_otherGenTauCharge);
 
   Notify();
}






Bool_t EMuAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EMuAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EMuAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void EMuAnalysis::GetEventWeight(){
  _eventPUWeight = 1.;
}

int EMuAnalysis::getMatchedCand(){
  int matchedCand = -1;
  if(muonMatched != 0){
    for(unsigned int mIt = 0; mIt < muonMatched->size(); ++mIt){
      if(muonMatched->at(mIt) !=  1 || elecMatched->at(mIt) != 1 || fabs(muonMotherId->at(mIt)) != 15 || fabs(eMotherId->at(mIt)) != 15)continue;
      matchedCand = mIt;
    }
  }
  return matchedCand;
}


bool EMuAnalysis::passedMatched(unsigned int theIndex){

  if (_source == "ZPrimeSSM_M-500" || _source == "ZPrimeSSM_M-750" || _source == "ZPrimeSSM_M-1000" || _source == "ZPrimeSSM_M-1250" || _source == "ZPrimeSSM_M-1500" || _source == "ZPrimeSSM_M-1750" || _source == "ZPrimeSSM_M-2000" || 
  _source == "ZPrimeSSM_M-2250" || _source == "ZPrimeSSM_M-2500") {
     if(!(elecMatched->at(theIndex) && muonMatched->at(theIndex))){
        return false;
     }
     else return true;	
  }
  else return true;
}


bool EMuAnalysis::elecExists(unsigned int theIndex){
  return elecPt->size() > 0;
}

bool EMuAnalysis::passedElecMinPt(unsigned int theIndex){
  //std::cout << "Elec pt = " << elecPt->at(theIndex) << " min elec Pt Cut = " << _minElecPtCut << "\n";
  //if (elecPt->at(theIndex) <= _minElecPtCut) std::cout << "Elec Pt failed -- pt = " << elecPt->at(theIndex) << "  cut value = " << _minElecPtCut;
  return elecPt->at(theIndex) > _minElecPtCut;
}

bool EMuAnalysis::passedElecMaxEta(unsigned int theIndex){
  //std::cout << "Elec Eta = " << fabs(elecEta->at(theIndex)) << " max ElecEta Cut = " << _maxElecEtaCut << "\n";
  return (fabs(elecSCEta->at(theIndex)) < _maxElecEtaCut) || (fabs(elecSCEta->at(theIndex)) > 1.56 &&  fabs(elecSCEta->at(theIndex)) < 2.5);
}

bool EMuAnalysis::passedElecAcc(unsigned int theIndex){
  return passedElecEt(theIndex) 
  && passedElecMaxEta(theIndex);
}

bool EMuAnalysis::elecInBarrel(unsigned int theIndex){
  return elecInEB->at(theIndex);
}
bool EMuAnalysis::elecInEndCap(unsigned int theIndex){
  return elecInEE->at(theIndex);
}
bool EMuAnalysis::elecInGaps(unsigned int theIndex){
  return elecInGap->at(theIndex);
}

bool EMuAnalysis::passedElecDeltaPhiIn(unsigned int theIndex){
  bool passed = false;
  if(elecInBarrel(theIndex)) passed = fabs(elecDeltaPhiIn->at(theIndex)) < _maxElecDeltaPhiInEBCut;
  if(elecInEndCap(theIndex)) passed = fabs(elecDeltaPhiIn->at(theIndex)) < _maxElecDeltaPhiInBBCut;
  return passed;
}
bool EMuAnalysis::passedElecDeltaEtaIn(unsigned int theIndex){
  bool passed = false;
  if(elecInBarrel(theIndex)) passed = fabs(elecDeltaEtaIn->at(theIndex)) < _maxElecDeltaEtaInEBCut;
  if(elecInEndCap(theIndex)) passed = fabs(elecDeltaEtaIn->at(theIndex)) < _maxElecDeltaEtaInBBCut;
  return passed;
}
bool EMuAnalysis::passedElecSigmaIEtaIEta(unsigned int theIndex){
  bool passed = false;
  if(elecInBarrel(theIndex)) passed = elecSigmaIEtaIEta->at(theIndex) < _maxElecSigmaIEtaIEtaEBCut;
  if(elecInEndCap(theIndex)) passed = elecSigmaIEtaIEta->at(theIndex) < _maxElecSigmaIEtaIEtaBBCut;
  return passed;
}
bool EMuAnalysis::passedElecHOverE(unsigned int theIndex){
  bool passed = false;
  if(elecInBarrel(theIndex)) passed = elecHOverEm->at(theIndex) < _maxElecHadFracEBCut;
  if(elecInEndCap(theIndex)) passed = elecHOverEm->at(theIndex) < _maxElecHadFracBBCut;
  return passed;
}
bool EMuAnalysis::passedElecEOverInvP(unsigned int theIndex){
  bool passed = false;
  if(elecInBarrel(theIndex)) passed = fabs(elecEOverInvP->at(theIndex)) < _maxOneOverEMinusOneOverPEBCut;
  if(elecInEndCap(theIndex)) passed = fabs(elecEOverInvP->at(theIndex)) < _maxOneOverEMinusOneOverPBBCut;
  return passed;
}

bool EMuAnalysis::passedElec2x5Over5x5(unsigned int theIndex){
   if (elecInBarrel(theIndex)) return (elecSCE2x5MaxOver5x5->at(theIndex) > 0.94 || elecSCE1x5Over5x5->at(theIndex) > 0.83);
   else return true;
   //return true;
}

bool EMuAnalysis::passedElecEcalHad1Iso(unsigned int theIndex){
   float ele_times_rho = 0.;
   float eta = fabs(elecSCEta->at(theIndex));
   float EcalHcal = elecEcalHcalIsoDepth1->at(theIndex); //gsfEcalIso->at(theIndex) + gsfHcalIso->at(theIndex);
   float et = gsfEnergy->at(theIndex) != 0. ? elecEt->at(theIndex) :0 ; //gsfEt->at(theIndex)/gsfEnergy->at(theIndex)*gsfCaloEnergy->at(theIndex) : 0.;
   /*
   if(eta < 1.0) ele_times_rho = 0.13;
   if(eta > 1.0 && eta < 1.479) ele_times_rho = 0.14;
   if(eta > 1.479 && eta < 2.0) ele_times_rho = 0.07;
   if(eta > 2.0 && eta < 2.2) ele_times_rho = 0.09;
   if(eta > 2.2 && eta < 2.3) ele_times_rho = 0.11;
   if(eta > 2.3 && eta < 2.4) ele_times_rho = 0.11;
   if(eta > 2.4) ele_times_rho = 0.14;
   */
   
   ele_times_rho = 0.28;
   
   
   
   if(eta < 1.479) return !(EcalHcal - ele_times_rho*rho > (2 + 0.03*et));
   else if (et < 50) return (!(EcalHcal - ele_times_rho*rho > 2.5));
   else return (!(EcalHcal - ele_times_rho*rho > (2.5 + 0.03*(et - 50))));
   
   //return true;
}

bool EMuAnalysis::passedElecRelIso(unsigned int theIndex){

  float ElecRelIso = (elecTrkIso->at(theIndex) + elecEcalIso->at(theIndex) + elecHcalIso->at(theIndex))/elecPt->at(theIndex);
  return ElecRelIso >= _minElecRelIsoCut && ElecRelIso < _maxElecRelIsoCut;
  //return true;

}

bool EMuAnalysis::passedElecTrkIso(unsigned int theIndex){
   return elecTrkIso->at(theIndex) >= _minElecTrkIsoCut && elecTrkIso->at(theIndex) < _maxElecTrkIsoCut;
}

bool EMuAnalysis::passedElecMissingHits(unsigned int theIndex){
  return elecMissingHits->at(theIndex) <= _maxElecMissingHitsCut;
}

bool EMuAnalysis::passedElecDxy(unsigned int theIndex){
  if(elecInBarrel(theIndex)) return fabs(elecGSFIp->at(theIndex)) < _maxElecDxyBarrelCut;
  else return fabs(elecGSFIp->at(theIndex)) < _maxElecDxyEndCapCut;
  //return true;
}

bool EMuAnalysis::passedElecId(unsigned int theIndex){
  return true;
}

bool EMuAnalysis::passedElecMotherId(unsigned int theIndex){
  return fabs(eMotherId->at(theIndex)) == _eMotherIdCut;
}

bool EMuAnalysis::passedMuonPt(unsigned int theIndex){
  //if(muonPt->at(theIndex) > _minMuonPtCut) std::cout << "Muon Pt Condition Passed\n";
  return muonPt->at(theIndex) > _minMuonPtCut;
}

bool EMuAnalysis::passedMuonE(unsigned int theIndex){
  return fabs(muonE->at(theIndex)) > _minMuonECut;
}

bool EMuAnalysis::passedMuonEta(unsigned int theIndex){
  return fabs(muonEta->at(theIndex)) < _maxMuonEtaCut;
}

bool EMuAnalysis::passedMuonAcc(unsigned int theIndex){
  return passedMuonPt(theIndex) && passedMuonEta(theIndex);
}

bool EMuAnalysis::passedMuonDxyVtx(unsigned int theIndex){
  return fabs(muonDxyVtx->at(theIndex)) < _maxMuonDxyVtxCut;
}

bool EMuAnalysis::passedMuonDzVtx(unsigned int theIndex){
  return fabs(muonDzVtx->at(theIndex)) < _maxMuonDzVtxCut;
}

bool EMuAnalysis::passedMuonIsGlobalMuon(unsigned int theIndex){
  return muonIsGlobalMuon->at(theIndex) >= _minMuonIsGlobalMuonCut;
}

bool EMuAnalysis::passedMuonMatched(unsigned int theIndex){
  return fabs(muonMatched->at(theIndex)) >= _minMuonMatchedCut;
}

bool EMuAnalysis::passedMuonMotherId(unsigned int theIndex){
  return fabs(muonMotherId->at(theIndex)) == _muonMotherIdCut;
}

bool EMuAnalysis::passedMuonNormChiSqrd(unsigned int theIndex){
  return fabs(muonNormChiSqrd->at(theIndex)) < _maxMuonNormChiSqrdCut;
}

bool EMuAnalysis::passedMuonValidHits(unsigned int theIndex){
  return fabs(muonValidHits->at(theIndex)) >= _minMuonValidHitsCut;
}

bool EMuAnalysis::passedMuonMatchedStations(unsigned int theIndex){
  return fabs(muonMatchedStations->at(theIndex)) >= _minMuonMatchedStationsCut;
}

bool EMuAnalysis::passedMuonPixelHits(unsigned int theIndex){
  return fabs(muonPixelHits->at(theIndex)) >= _minMuonPixelHitsCut;
}

bool EMuAnalysis::passedMuonTrkLayersWithHits(unsigned int theIndex){
  return fabs(muonTrkLayersWithHits->at(theIndex)) >= _minMuonTrkLayersWithHitsCut;
}

bool EMuAnalysis::passedMuonIsoDR03NoDB(unsigned int theIndex){
  return (muonPFIsoDR03SumChargedHadronPt->at(theIndex) + muonPFIsoDR03SumNeutralHadronPt->at(theIndex) + muonPFIsoDR03SumPhotonPt->at(theIndex))/muonPt->at(theIndex) < _maxMuonIsoDR03NoDBCut;
}

bool EMuAnalysis::passedMuonIsoDR03DB(unsigned int theIndex){
  return (muonPFIsoDR03SumChargedHadronPt->at(theIndex) + muonPFIsoDR03SumNeutralHadronPt->at(theIndex) + muonPFIsoDR03SumPhotonPt->at(theIndex) - 0.5*muonPFIsoDR03SumPUPt->at(theIndex))/muonPt->at(theIndex) < _maxMuonIsoDR03DBCut;
}

bool EMuAnalysis::passedMuonIsoDR04NoDB(unsigned int theIndex){
  return (muonPFIsoDR04SumChargedHadronPt->at(theIndex) + muonPFIsoDR04SumNeutralHadronPt->at(theIndex) + muonPFIsoDR04SumPhotonPt->at(theIndex))/muonPt->at(theIndex) < _maxMuonIsoDR04NoDBCut;
}

bool EMuAnalysis::passedMuonIsoDR04DB(unsigned int theIndex){
  
  float IsoCorrection = muonPFIsoDR04SumChargedHadronPt->at(theIndex) + std::max(0.0, muonPFIsoDR04SumNeutralHadronPt->at(theIndex) + muonPFIsoDR04SumPhotonPt->at(theIndex) - 0.5*muonPFIsoDR04SumPUPt->at(theIndex));
  float muonIso = (IsoCorrection)/(muonPt->at(theIndex));
  
  if(IsoCorrection < 0){
    std::cout << "Muon Iso Numerator = " << IsoCorrection << "\n";
  }
  
  return muonIso < _maxMuonIsoDR04DBCut && muonIso >= _minMuonIsoDR04DBCut;
}

bool EMuAnalysis::passedMuonIsHighPtMuon(unsigned int theIndex){
  return muonIsHighPtMuon->at(theIndex) >= _minMuonIsHighPtMuonCut;
}

bool EMuAnalysis::passedMuonDptOverPt(unsigned int theIndex){
  float ratio = muonTrkPtError->at(theIndex)/muonTrkPt->at(theIndex);
  if (muonTrkPt->at(theIndex) < 0) {
  	return false;
  }
  else{
   	return ratio < _maxMuonDptOverPtCut;
  }	
}

bool EMuAnalysis::passedAllMuonIDCuts(unsigned int theIndex){
  return passedMuonIsGlobalMuon(theIndex) && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
  && passedMuonValidHits(theIndex)
  && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex)
  && passedMuonTrkLayersWithHits(theIndex) && passedMuonIsoDR04DB(theIndex) && passedMuonDptOverPt(theIndex); 
}

/*
bool EMuAnalysis::passedAllMuonIDCuts(unsigned int theIndex){
  return passedMuonIsHighPtMuon(theIndex) && passedMuonIsoDR04DB(theIndex);
}
*/

bool EMuAnalysis::passedElecMatched(unsigned int theIndex){ 
  return fabs(elecMatched->at(theIndex)) >= _minElecMatchedCut;
}

bool EMuAnalysis::passedHeepPassedEt(unsigned int theIndex){
   return heepPassedEt->at(theIndex) >= _minHeepPassedEtCut;
}

bool EMuAnalysis::passedHeepPassedCrack(unsigned int theIndex){
  return fabs(heepPassedCrack->at(theIndex)) >= _minHeepPassedCrackCut;
}
bool EMuAnalysis::passedHeepPassedDEtaIn(unsigned int theIndex){
  return fabs(heepPassedDEtaIn->at(theIndex)) >= _minHeepPassedDEtaInCut;
}
bool EMuAnalysis::passedHeepPassedDPhiIn(unsigned int theIndex){
  return fabs(heepPassedDPhiIn->at(theIndex)) >= _minHeepPassedDPhiInCut;
}
bool EMuAnalysis::passedHeepPassedHadem(unsigned int theIndex){
  return fabs(heepPassedHadem->at(theIndex)) >= _minHeepPassedHademCut;
}
bool EMuAnalysis::passedHeepPassedSigmaIEtaIEta(unsigned int theIndex){
  return fabs(heepPassedSigmaIEtaIEta->at(theIndex)) >= _minHeepPassedSigmaIEtaIEtaCut;
}
bool EMuAnalysis::passedHeepPassed2by5Over5By5(unsigned int theIndex){
  return fabs(heepPassed2by5Over5By5->at(theIndex)) >= _minHeepPassed2by5Over5By5Cut;
}
bool EMuAnalysis::passedHeepPassedEcalHad1Iso(unsigned int theIndex){
  return fabs(heepPassedEcalHad1Iso->at(theIndex)) >= _minHeepPassedEcalHad1IsoCut;
}
/*
bool EMuAnalysis::passedHeepPassedHad2Iso(unsigned int theIndex){
  return fabs(heepPassedHad2Iso->at(theIndex)) > _minHeepPassedHad2IsoCut;
}
*/
bool EMuAnalysis::passedHeepPassedTrkIso(unsigned int theIndex){
  return heepPassedTrkIso->at(theIndex) >= _minHeepPassedTrkIsoCut;
}
bool EMuAnalysis::passedHeepPassedEcalDriven(unsigned int theIndex){
  return fabs(heepPassedEcalDriven->at(theIndex)) >= _minHeepPassedEcalDrivenCut;
  //return elecEcalDriven->at(theIndex) >= _minHeepPassedEcalDrivenCut;
}

bool EMuAnalysis::passedHeepPassedMissingHits(unsigned int theIndex){
   return heepPassedMissingHits->at(theIndex) == _minHeepPassedMissingHitsCut;
}

bool EMuAnalysis::passedHeepPassedElecDxy(unsigned int theIndex){
   return heepPassedElecDxy->at(theIndex) >= _minHeepPassedElecDxyCut;
}   

bool EMuAnalysis::passedElecEt(unsigned int theIndex){
  return fabs(elecEt->at(theIndex)) >= _minElecEtCut;
}


bool EMuAnalysis::passedHeepId(unsigned int theIndex){
  
  return passedHeepPassedDEtaIn(theIndex) && passedHeepPassedDPhiIn(theIndex)
  && passedHeepPassedHadem(theIndex) && passedHeepPassedSigmaIEtaIEta(theIndex)
  && passedHeepPassed2by5Over5By5(theIndex) && passedHeepPassedEcalHad1Iso(theIndex)
  && passedHeepPassedEcalDriven(theIndex) && passedHeepPassedTrkIso(theIndex)
  && passedHeepPassedMissingHits(theIndex) && passedElecDxy(theIndex);
  
  //return heepPassedAllCuts->at(theIndex) > 0;
}

bool EMuAnalysis::passedElecIpOverElecIpError(unsigned int theIndex){
  float IpOverIpErr = fabs(elecGSFBSIp->at(theIndex)/elecGSFBSIpError->at(theIndex));
  return IpOverIpErr >= _minElecIpOverElecIpErrorCut && IpOverIpErr;
}

bool EMuAnalysis::passedMuonIpOverMuonIpError(unsigned int theIndex){
  float IpOverIpErr = fabs(muonDxyBS->at(theIndex)/muonDxyBSError->at(theIndex));
  return IpOverIpErr >= _minMuonIpOverMuonIpErrorCut;
}

bool EMuAnalysis::passedCombinedIpOverIpError(unsigned int theIndex){
  float eSF = 1.0;
  float muSF = 1.0;
  if(_source == "TTJets_MADGRAPH") {
    eSF = 1.08604317;
    muSF = 1.06639004;  
  }
  float ElecIpOverElecIpErr = fabs(elecGSFBSIp->at(theIndex)/elecGSFBSIpError->at(theIndex));
  float MuonIpOverMuonIpErr = fabs(muonDxyBS->at(theIndex)/muonDxyBSError->at(theIndex));
  float ElecIp = fabs(elecGSFBSIp->at(theIndex));
  float MuonIp = fabs(muonDxyBS->at(theIndex));
  float ElecIpErr = eSF*fabs(elecGSFBSIpError->at(theIndex));
  float MuonIpErr = muSF*fabs(muonDxyBSError->at(theIndex));
  

  
  //return sqrt(ElecIpOverElecIpErr*ElecIpOverElecIpErr + MuonIpOverMuonIpErr*MuonIpOverMuonIpErr) >= _minCombinedIpOverIpErrorCut && sqrt(ElecIpOverElecIpErr*ElecIpOverElecIpErr + MuonIpOverMuonIpErr*MuonIpOverMuonIpErr) <
  //_maxCombinedIpOverIpErrorCut;
  
  return sqrt((ElecIp*ElecIp + MuonIp*MuonIp + 2*ElecIp*MuonIp)/(ElecIpErr*ElecIpErr + MuonIpErr*MuonIpErr)) >= _minCombinedIpOverIpErrorCut && 
  sqrt((ElecIp*ElecIp + MuonIp*MuonIp + 2*ElecIp*MuonIp)/(ElecIpErr*ElecIpErr + MuonIpErr*MuonIpErr)) < _maxCombinedIpOverIpErrorCut;
  
  //return ElecIpOverElecIpErr*MuonIpOverMuonIpErr >= _minCombinedIpOverIpErrorCut && ElecIpOverElecIpErr*MuonIpOverMuonIpErr < _maxCombinedIpOverIpErrorCut;
  
  
}

bool EMuAnalysis::passedElecMuonDCAOverError2D(unsigned int theIndex){
  float EMuDCAOverErr2D = elecMuonTwoTrkMinDist->at(theIndex)/elecMuonDCA2DError->at(theIndex);
  return EMuDCAOverErr2D >= _minElecMuonDCAOverError2DCut && EMuDCAOverErr2D < _maxElecMuonDCAOverError2DCut;
  
}


bool EMuAnalysis::passedElecMuonDCAOverError3D(unsigned int theIndex){
  float EMuDCAOverErr3D = elecMuonTwoTrkMinDist->at(theIndex)/elecMuonDCA3DError->at(theIndex);
  return EMuDCAOverErr3D >= _minElecMuonDCAOverError3DCut && EMuDCAOverErr3D < _maxElecMuonDCAOverError3DCut;
  
}

bool EMuAnalysis::passedLifetime(unsigned int theIndex){
  float eSF = 1.0;
  float muSF = 1.0;
  if(_source == "TTJets_MADGRAPH") {
    eSF = 1.08604317;
    muSF = 1.06639004;  
  }  
  
  float ElecIpOverElecIpErr = fabs(elecGSFIp->at(theIndex)/elecGSFIpError->at(theIndex));
  float MuonIpOverMuonIpErr = fabs(muonDxyVtx->at(theIndex)/muonDxyError->at(theIndex));
  float ElecIp = fabs(elecGSFIp->at(theIndex));
  float MuonIp = fabs(muonDxyVtx->at(theIndex));
  float ElecIpErr = eSF*fabs(elecGSFIpError->at(theIndex));
  float MuonIpErr = muSF*fabs(muonDxyError->at(theIndex));
  float CombinedIpOverErr = sqrt((ElecIp*ElecIp + MuonIp*MuonIp + 2*ElecIp*MuonIp)/(ElecIpErr*ElecIpErr + MuonIpErr*MuonIpErr));
  float EMuDCAOverErr3D = elecMuonTwoTrkMinDist->at(theIndex)/elecMuonDCA3DError->at(theIndex);
  
  return CombinedIpOverErr >= _minLifetimeCut || EMuDCAOverErr3D >= _minLifetimeCut;
  
}


bool EMuAnalysis::passedElecMuonMass(unsigned int theIndex){
  return elecMuononMass->at(theIndex) < _elecMuonMassCut;
}

bool EMuAnalysis::passedElecMuonCosDPhi(unsigned int theIndex){
  return elecMuonCosDPhi->at(theIndex) > _minElecMuonCosDPhiCut && elecMuonCosDPhi->at(theIndex) < _maxElecMuonCosDPhiCut;
}

bool EMuAnalysis::passedElecMuonDeltaR(unsigned int theIndex){
  return elecMuonDeltaR->at(theIndex) >= _elecMuonDeltaRCut;
}

bool EMuAnalysis::passedElecMuonMetMass(unsigned int theIndex){
  return elecMuonMetMass->at(theIndex) >= _minElecMuonMetMassCut && elecMuonMetMass->at(theIndex) < _maxElecMuonMetMassCut;
}

bool EMuAnalysis::passedElecMuonDeltaPtMass(unsigned int theIndex){
  return elecMuonDeltaPtMass->at(theIndex) >= _elecMuonDeltaPtMassCut;
}

bool EMuAnalysis::passedElecMetMt(unsigned int theIndex){
  return elecMEtMt->at(theIndex) < _elecMetMtCut;
}

bool EMuAnalysis::passedMuonMetMt(unsigned int theIndex){
  return (muonMEtMt->at(theIndex) > _minMuonMetMtCut && muonMEtMt->at(theIndex) < _maxMuonMetMtCut);
}

/*
bool EMuAnalysis::passedElecMuonPZetaVis(unsigned int theIndex){
  return elecMuonPZetaVis->at(theIndex) > _elecMuonPZetaVisCut;
}
*/
//Linear PZeta-PZetaVis Cut

bool EMuAnalysis::passedElecMuonPZeta(unsigned int theIndex){
  float PZeta = (elecMuonPZeta->at(theIndex) - _pZetaSlope*elecMuonPZetaVis->at(theIndex));
  return PZeta >= _minElecMuonPZetaCut && PZeta < _maxElecMuonPZetaCut;
  
  
}


//Elliptical PZeta-PZetaVis Cut
/*
bool EMuAnalysis::passedElecMuonPZeta(unsigned int theIndex){
  float x = elecMuonPZetaVis->at(theIndex);
  float x0 = _pZetaVisCentroidCut;
  float y = elecMuonPZeta->at(theIndex);
  float y0 = _pZetaCentroidCut;
  float a = _pZetaSemiMajorAxisCut;
  float b = _pZetaSemiMinorAxisCut;
  float theta = _pZetaEllipseAngleCut*3.14159/180.0;
  
  float ellipse = pow(((x - x0)*cos(theta) + (y - y0)*sin(theta)),2)/(pow(a,2)) + pow(((x - x0)*sin(theta) - (y - y0)*cos(theta)),2)/(pow(b,2));
  
  return ellipse < 1 && x > 0; 
}

*/
bool EMuAnalysis::passedNbTagsHiEffTrkCnt(unsigned int theIndex){
  return nBtagsHiEffTrkCnt->at(theIndex) <= _maxNbTagsHiEffTrkCntCut && nBtagsHiEffTrkCnt->at(theIndex) >= _minNbTagsHiEffTrkCntCut;
}

bool EMuAnalysis::passedNbTagsHiPurityTrkCnt(unsigned int theIndex){
  return nBtagsHiPurityTrkCnt->at(theIndex) <= _maxNbTagsHiPurityTrkCntCut && nBtagsHiPurityTrkCnt->at(theIndex) >= _minNbTagsHiPurityTrkCntCut;
}

bool EMuAnalysis::passedNbTagsHiEffSimpleSecVtx(unsigned int theIndex){
  return nBTagsHiEffSimpleSecVtx->at(theIndex) <= _maxNbTagsHiEffSimpleSecVtxCut && nBTagsHiEffSimpleSecVtx->at(theIndex) >= _minNbTagsHiEffSimpleSecVtxCut;
}

bool EMuAnalysis::passedNbTagsHiPuritySimpleSecVtx(unsigned int theIndex){
  return nBTagsHiPuritySimpleSecVtx->at(theIndex) <= _maxNbTagsHiPuritySimpleSecVtxCut && nBTagsHiPuritySimpleSecVtx->at(theIndex) >= _minNbTagsHiPuritySimpleSecVtxCut;
}

bool EMuAnalysis::passedNbTagsCombSecVtx(unsigned int theIndex){
  return nBTagsCombSecVtxLWP->at(theIndex) <= _maxNbTagsCombSecVtxCut && nBTagsCombSecVtxLWP->at(theIndex) >= _minNbTagsCombSecVtxCut && nJets->at(theIndex) < _nJetsCut;
}

bool EMuAnalysis::passedLeptonMetCosDPhi(unsigned int theIndex){
  return leptonMetCosDphi->at(theIndex) < _leptonMetCosDPhiCut;
}

bool EMuAnalysis::passedLeptonDeltaPtCosDPhi(unsigned int theIndex){
  return leptonDeltaPtCosDphi->at(theIndex) < _maxLeptonDeltaPtCosDPhiCut;
}

bool EMuAnalysis::passedMET(unsigned int theIndex){
  return mEt >= _metCut;
}

bool EMuAnalysis::passedThreeLepton(unsigned int theIndex){ 
 int threeLeptonVertex = eMuMuOneVertex || muEEOneVertex;
 bool passedMuonIDs = false;
 bool passedElecIDs = false;
 //std::cout << "threeLeptonMuonPt->size() = " << threeLeptonMuonPt->size() << "\n";
 //std::cout << "threeLeptonElecEt->size() = " << threeLeptonElecEt->size() << "\n";
 
 for(int i = 0; i < threeLeptonMuonPt->size(); i++){
   float threeLeptonMuonDptOverPt = threeLeptonMuonTrkPtError->at(i)/threeLeptonMuonTrkPt->at(i);
   float threeLeptonOtherMuonDptOverPt = threeLeptonOtherMuonTrkPtError->at(i)/threeLeptonOtherMuonTrkPt->at(i);
   float threeLeptonMuonIso = (threeLeptonMuonPFIsoDR04SumChargedHadronPt->at(i) + std::max(0.0, threeLeptonMuonPFIsoDR04SumCNeutralHadronPt->at(i) + threeLeptonMuonPFIsoDR04SumCPhotonPt->at(i) -
   0.5*threeLeptonMuonPFIsoDR04SumPUPt->at(i)))/threeLeptonMuonPt->at(i);
   float threeLeptonOtherMuonIso = (threeLeptonOtherMuonPFIsoDR04SumChargedHadronPt->at(i) + std::max(0.0, threeLeptonOtherMuonPFIsoDR04SumNeutralHadronPt->at(i) + threeLeptonOtherMuonPFIsoDR04SumPhotonPt->at(i) -
   0.5*threeLeptonOtherMuonPFIsoDR04SumPUPt->at(i)))/threeLeptonOtherMuonPt->at(i);   
   
   
   if(threeLeptonMuonPt->at(i) > 20. && threeLeptonOtherMuonPt->at(i) > 10. && threeLeptonIsGlobalMuon->at(i) && threeLeptonOtherMuonIsGlobalMuon->at(i) && threeLeptonMuonDxyVtx->at(i) < 0.2 && threeLeptonOtherMuonDxy->at(i) < 0.2
   && threeLeptonMuonDzVtx->at(i) < 0.5 && threeLeptonOtherMuonDzVtx->at(i) < 0.5 && threeLeptonMuonChamberHits->at(i) > 0 && threeLeptonOtherMuonChamberHits->at(i) > 0 && threeLeptonMuonPixelHits->at(i) > 0 &&
   threeLeptonOtherMuonPixelHits->at(i) > 0 && threeLeptonMuonMatchedStations->at(i) > 1 && threeLeptonOtherMuonMatchedStations->at(i) > 1 && threeLeptonMuonTrkLayersWithHits->at(i) > 5 && threeLeptonOtherMuonTrkLayersWithHits->at(i) >
   5 && threeLeptonMuonDptOverPt < 0.3 && threeLeptonOtherMuonDptOverPt < 0.3 /*&& threeLeptonMuonIso < 0.12 && threeLeptonOtherMuonIso < 0.12*/ /*&& threeLeptonMuMuDCAOverError->at(i) < 2.*/){
   /*
    std::cout << "threeLeptonMuonPt->at(i) = " <<  threeLeptonMuonPt->at(i) << "\n";
    std::cout << "threeLeptonOtherMuonPt->at(i) = " <<  threeLeptonOtherMuonPt->at(i) << "\n";
    std::cout << "threeLeptonIsGlobalMuon->at(i) = " << threeLeptonIsGlobalMuon->at(i) << "\n";
    std::cout << "threeLeptonOtherMuonIsGlobalMuon->at(i) = " << threeLeptonOtherMuonIsGlobalMuon->at(i) << "\n";
    */
   
    passedMuonIDs = true;
   }
 }
 
 for(int j = 0; j < threeLeptonElecEt->size(); j++){
   bool elecSCEta = threeLeptonElecSCEta->at(j) < 1.442 || (threeLeptonElecSCEta->at(j) > 1.56 && threeLeptonElecSCEta->at(j) < 2.5);
   bool otherElecSCEta = threeLeptonOtherElecSCEta->at(j) < 1.442 || (threeLeptonOtherElecSCEta->at(j) > 1.56 && threeLeptonOtherElecSCEta->at(j) < 2.5);
   
   if(threeLeptonElecEt->at(j) > 20. && threeLeptonOtherElecEt->at(j) > 10. && elecSCEta && otherElecSCEta && threeLeptonElecEcalDriven->at(j) && threeLeptonOtherElecEcalDriven->at(j) &&
   threeLeptonElecDEtaIn->at(j) && threeLeptonOtherElecDEtaIn->at(j) && threeLeptonElecDPhiIn->at(j) && threeLeptonOtherElecDPhiIn->at(j) && threeLeptonElecHadem->at(j) && threeLeptonOtherElecHadem->at(j) &&
   threeLeptonElecSigmaIEtaIEta->at(j) && threeLeptonOtherElecSigmaIEtaIEta->at(j) && threeLeptonElec2by5Over5By5->at(j) && threeLeptonOtherElec2by5Over5By5->at(j) && threeLeptonElecEcalHad1Iso->at(j) &&
   threeLeptonOtherElecEcalHad1Iso->at(j) /*&& threeLeptonElecTrkIso->at(j) && threeLeptonOtherElecTrkIso->at(j)*/ && threeLeptonElecMissingHits->at(j) && threeLeptonOtherElecMissingHits->at(j) &&
   threeLeptonElecDxy->at(j) && threeLeptonOtherElecDxy->at(j) /*&& elecElecDCAOverError->at(j) < 2.*/){
     
     passedElecIDs = true;
   }
     
 }
 
 
 
 //return threeLeptonVertex <= _maxThreeLeptonCut /*&& passedIDs*/;
 //return (eMuMuOneVertex <= _maxThreeLeptonCut && passedMuonIDs == true) || (muEEOneVertex <= _maxThreeLeptonCut && passedElecIDs == true);
 if((passedMuonIDs == true && eMuMuOneVertex == true) || (passedElecIDs == true && muEEOneVertex == true)){
   return 0 == _maxThreeLeptonCut;
 }
 else return true;
 
}

bool EMuAnalysis::passedNJets(unsigned int theIndex){
   
  //return (nJets->at(theIndex) < _nJetsCut && passedNbTagsCombSecVtx(theIndex)) ? true : false;
  return nJets->at(theIndex) < _nJetsCut;
}

bool EMuAnalysis::passedCharge(unsigned int theIndex){
  int combinedCharge = muonCharge->at(theIndex)*elecCharge->at(theIndex);
  return combinedCharge > _minChargeCut && combinedCharge < _maxChargeCut; 
}

bool EMuAnalysis::passedJetSumEt(unsigned int theIndex){
  return jetSumEt->at(theIndex) <= _jetSumEtCut;
}  

bool EMuAnalysis::passedElecMuonDeltaPt(unsigned int theIndex){
  float deltaPt = fabs(muonPt->at(theIndex) - elecPt->at(theIndex));
  return deltaPt >= _minElecMuonDeltaPtCut;
}

bool EMuAnalysis::passed2xPLargerLepton(unsigned int theIndex){
  float eP = fabs(elecE->at(theIndex));
  float muP = fabs(muonE->at(theIndex));
  float twoP = 2*muP;
  if(eP > muP) twoP = 2*eP;
  return twoP >= _min2xPLargerLeptonCut;
}

bool EMuAnalysis::passedElecEOverP(unsigned int theIndex){
  float IpOverIpErr = fabs(muonDxyBS->at(theIndex)/muonDxyBSError->at(theIndex));
  if(IpOverIpErr < 1.0 && elecEOverP->at(theIndex) > 2) return false;
  else return elecEOverP->at(theIndex) > _minElecEOverPCut && elecEOverP->at(theIndex) < _maxElecEOverPCut;
}

bool EMuAnalysis::passedTopology(unsigned int theIndex){
  return  passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) &&
  passedElecMuonPZeta(theIndex) && passedNbTagsCombSecVtx(theIndex) &&
  passedLeptonMetCosDPhi(theIndex) && passedLeptonDeltaPtCosDPhi(theIndex) && passedMET(theIndex) && passedCharge(theIndex) && passedJetSumEt(theIndex) && passedNJets(theIndex) && 
  passedElecIpOverElecIpError(theIndex) && passedMuonIpOverMuonIpError(theIndex) && passedCombinedIpOverIpError(theIndex) && 
  passedElecMuonDCAOverError3D(theIndex) && passedElecMuonDeltaPt(theIndex) && passed2xPLargerLepton(theIndex) && passedLifetime(theIndex);
}



//Define Control Region Cuts

bool EMuAnalysis::passedQCDCR1(unsigned int theIndex){
  float IsoCorrection = muonPFIsoDR04SumChargedHadronPt->at(theIndex) + muonPFIsoDR04SumNeutralHadronPt->at(theIndex) + muonPFIsoDR04SumPhotonPt->at(theIndex) - 0.5*muonPFIsoDR04SumPUPt->at(theIndex);
  if(IsoCorrection < 0) IsoCorrection = 0;
  float muonIso = (IsoCorrection)/(muonPt->at(theIndex));
  float ElecRelIso = (elecTrkIso->at(theIndex) + elecEcalIso->at(theIndex) + elecHcalIso->at(theIndex))/elecPt->at(theIndex);
  //return true;
  return passedElecAcc(theIndex) && passedMuonAcc(theIndex) && passedMuonIsGlobalMuon(theIndex) && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
  && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
  && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex)
  && passedMuonTrkLayersWithHits(theIndex) && muonIso > 0.1 && muonIso < 0.2 && passedHeepPassedCrack(theIndex)
  && passedHeepPassedDEtaIn(theIndex) && passedHeepPassedDPhiIn(theIndex)
  && passedHeepPassedHadem(theIndex) && passedHeepPassedSigmaIEtaIEta(theIndex)
  && passedHeepPassed2by5Over5By5(theIndex)
  && passedHeepPassedEcalDriven(theIndex) && ElecRelIso > 0. && ElecRelIso < 0.2
  && passedElecEt(theIndex) && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedJetSumEt(theIndex) && passedNJets(theIndex);
}


bool EMuAnalysis::passedQCDCR2(unsigned int theIndex){
  float IsoCorrection = muonPFIsoDR04SumChargedHadronPt->at(theIndex) + muonPFIsoDR04SumNeutralHadronPt->at(theIndex) + muonPFIsoDR04SumPhotonPt->at(theIndex) - 0.5*muonPFIsoDR04SumPUPt->at(theIndex);
  if(IsoCorrection < 0) IsoCorrection = 0;
  float muonIso = (IsoCorrection)/(muonPt->at(theIndex));  
  float ElecRelIso = (elecTrkIso->at(theIndex) + elecEcalIso->at(theIndex) + elecHcalIso->at(theIndex))/elecPt->at(theIndex);
  //return true;
  return passedElecAcc(theIndex) && passedMuonAcc(theIndex) && passedMuonIsGlobalMuon(theIndex) && passedMuonDxyVtx(theIndex) && passedMuonDzVtx(theIndex)
  && passedMuonNormChiSqrd(theIndex) && passedMuonValidHits(theIndex)
  && passedMuonMatchedStations(theIndex) && passedMuonPixelHits(theIndex)
  && passedMuonTrkLayersWithHits(theIndex) && muonIso < 0.2 && muonIso > 0. && passedHeepPassedCrack(theIndex)
  && passedHeepPassedDEtaIn(theIndex) && passedHeepPassedDPhiIn(theIndex)
  && passedHeepPassedHadem(theIndex) && passedHeepPassedSigmaIEtaIEta(theIndex)
  && passedHeepPassed2by5Over5By5(theIndex)
  && passedHeepPassedEcalDriven(theIndex) && ElecRelIso < 0.2 && ElecRelIso > 0.1
  && passedElecEt(theIndex) && passedElecMuonCosDPhi(theIndex) && passedElecMuonDeltaR(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedJetSumEt(theIndex) && passedNJets(theIndex);
}

bool EMuAnalysis::passedTTJetsCR1(unsigned int theIndex){
  return passedElecAcc(theIndex) && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex) && passedHeepId(theIndex)
&& passedElecMuonDeltaR(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedMET(theIndex) && passedCharge(theIndex);  
}

bool EMuAnalysis::passedTTJetsCR2(unsigned int theIndex){
  return passedElecAcc(theIndex) && passedHeepId(theIndex) && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex)
  && passedElecMuonDeltaR(theIndex) && passedCharge(theIndex) && mEt > 60 && muonMEtMt->at(theIndex) > 0;
}

bool EMuAnalysis::passedWJetsCR1(unsigned int theIndex){
  float ElecRelIso = (elecTrkIso->at(theIndex) + elecEcalIso->at(theIndex) + elecHcalIso->at(theIndex))/elecPt->at(theIndex);
  return passedElecAcc(theIndex) && passedHeepPassedCrack(theIndex) && passedHeepPassedDEtaIn(theIndex) && passedHeepPassedDPhiIn(theIndex) && passedHeepPassedHadem(theIndex) 
  && passedHeepPassedSigmaIEtaIEta(theIndex) && passedHeepPassed2by5Over5By5(theIndex) && passedHeepPassedEcalDriven(theIndex) && passedElecEt(theIndex) 
  && ElecRelIso > 0.1 && ElecRelIso < 0.2 && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex)
  && passedElecMuonDeltaR(theIndex) && muonMEtMt->at(theIndex) > 70 && muonMEtMt->at(theIndex) < 120 && mEt > 40 
  && passedNJets(theIndex) && passedNbTagsCombSecVtx(theIndex) && passedJetSumEt(theIndex) && passedLeptonMetCosDPhi(theIndex);
  //return true;
}

bool EMuAnalysis::passedWJetsCR2(unsigned int theIndex){
  int combinedCharge = muonCharge->at(theIndex)*elecCharge->at(theIndex);
  return passedElecAcc(theIndex) && passedHeepId(theIndex) && passedMuonAcc(theIndex) && passedAllMuonIDCuts(theIndex)
  && passedElecMuonDeltaR(theIndex) && elecMuonCosDPhi->at(theIndex) > -0.95  && elecMuonPZeta->at(theIndex) < -20 && passedNbTagsCombSecVtx(theIndex)
  && passedMET(theIndex) && passedNJets(theIndex) && combinedCharge == 1 && passedJetSumEt(theIndex) && passedLeptonMetCosDPhi(theIndex);
}

bool EMuAnalysis::passedDYToTauTauCR1(unsigned int theIndex){
	return true;
}

bool EMuAnalysis::passedDYToTauTauCR2(unsigned int theIndex){
	return true;
}


void EMuAnalysis::SetElecAccCuts(double minElecPtCut, double maxElecEtaCut){
  _minElecPtCut = minElecPtCut;
  _maxElecEtaCut = maxElecEtaCut;
}


void EMuAnalysis::SetElecMotherIdCut(double eMotherIdCut){
  _eMotherIdCut = eMotherIdCut;
}

void EMuAnalysis::SetElecIdGlobalCuts(int maxElecMissingHitsCut,
                                              double maxElecPFRelIsoCut){
  _maxElecMissingHitsCut = maxElecMissingHitsCut;
  //_elecIsoDr = elecIsoDr;
  _maxElecPFRelIsoCut = maxElecPFRelIsoCut;
}
void EMuAnalysis::SetElecIdBarrelCuts(double maxElecDeltaPhiInEBCut, double maxElecDeltaEtaInEBCut,
     					     double maxElecSigmaIEtaIEtaEBCut, double maxElecHadFracEBCut, 
					     double maxOneOverEMinusOneOverPEBCut, double maxElecDxyBarrelCut){
  _maxElecDeltaPhiInEBCut = maxElecDeltaPhiInEBCut;
  _maxElecDeltaEtaInEBCut = maxElecDeltaEtaInEBCut;
  _maxElecSigmaIEtaIEtaEBCut = maxElecSigmaIEtaIEtaEBCut;
  _maxElecHadFracEBCut = maxElecHadFracEBCut;
  _maxOneOverEMinusOneOverPEBCut = maxOneOverEMinusOneOverPEBCut;
  _maxElecDxyBarrelCut = maxElecDxyBarrelCut;
}
void EMuAnalysis::SetElecIdEndCapCuts(double maxElecDeltaPhiInBBCut, double maxElecDeltaEtaInBBCut,
                                             double maxElecSigmaIEtaIEtaBBCut, double maxElecHadFracBBCut, 
					     double maxOneOverEMinusOneOverPBBCut, double maxElecDxyEndCapCut){
  _maxElecDeltaPhiInBBCut = maxElecDeltaPhiInBBCut;
  _maxElecDeltaEtaInBBCut = maxElecDeltaEtaInBBCut;
  _maxElecSigmaIEtaIEtaBBCut = maxElecSigmaIEtaIEtaBBCut;
  _maxElecHadFracBBCut = maxElecHadFracBBCut;
  _maxOneOverEMinusOneOverPBBCut = maxOneOverEMinusOneOverPBBCut;
  _maxElecDxyEndCapCut = maxElecDxyEndCapCut;
}

void EMuAnalysis::SetMuonAccCuts(double minMuonPtCut, double maxMuonEtaCut){

  _minMuonPtCut = minMuonPtCut;
  _maxMuonEtaCut = maxMuonEtaCut;
}

void EMuAnalysis::SetMuonIdCuts(double maxMuonDxyVtxCut, double maxMuonDzVtxCut,
			      double minMuonIsGlobalMuonCut, double minMuonMatchedCut, double maxMuonNormChiSqrdCut, double minMuonValidHitsCut,
			      double minMuonMatchedStationsCut, double minMuonPixelHitsCut, double minMuonTrkLayersWithHitsCut, double minMuonIsoDR04DBCut,
			      double maxMuonIsoDR04DBCut, double minMuonIsHighPtMuonCut, double maxMuonDptOverPtCut){
			     

  _maxMuonDxyVtxCut = maxMuonDxyVtxCut;
  _maxMuonDzVtxCut = maxMuonDzVtxCut;
  _minMuonIsGlobalMuonCut = minMuonIsGlobalMuonCut;	
  _minMuonMatchedCut = minMuonMatchedCut;
  _maxMuonNormChiSqrdCut = maxMuonNormChiSqrdCut;	
  _minMuonValidHitsCut = minMuonValidHitsCut;
  _minMuonMatchedStationsCut = minMuonMatchedStationsCut;
  _minMuonPixelHitsCut = minMuonPixelHitsCut;
  _minMuonTrkLayersWithHitsCut = minMuonTrkLayersWithHitsCut;
  _minMuonIsoDR04DBCut = minMuonIsoDR04DBCut;
  _maxMuonIsoDR04DBCut = maxMuonIsoDR04DBCut;
  _minMuonIsHighPtMuonCut = minMuonIsHighPtMuonCut;
  _maxMuonDptOverPtCut = maxMuonDptOverPtCut;		     		     
  		     
}

void EMuAnalysis::SetMuonMotherIdCut(double muonMotherIdCut){
  _muonMotherIdCut = muonMotherIdCut;
}

void EMuAnalysis::SetHeepCuts(double minElecMatchedCut, double minHeepPassedEtCut, double minHeepPassedDetEtaCut, double minHeepPassedCrackCut, double minHeepPassedDEtaInCut, 
			      double minHeepPassedDPhiInCut, double minHeepPassedHademCut, double minHeepPassedSigmaIEtaIEtaCut, double minHeepPassed2by5Over5By5Cut, 
			      double minHeepPassedEcalDrivenCut, double minHeepPassedMissingHitsCut, double minHeepPassedElecDxyCut, double minElecEtCut, double minElecRelIsoCut, double maxElecRelIsoCut, double minElecTrkIsoCut, double maxElecTrkIsoCut, 
			      double minHeepPassedEcalHad1IsoCut, double minHeepPassedTrkIsoCut){

  _minElecMatchedCut = minElecMatchedCut;
  _minHeepPassedEtCut = minHeepPassedEtCut;
  _minHeepPassedDetEtaCut = minHeepPassedDetEtaCut;
  _minHeepPassedCrackCut = minHeepPassedCrackCut;
  _minHeepPassedDEtaInCut = minHeepPassedDEtaInCut;
  _minHeepPassedDPhiInCut = minHeepPassedDPhiInCut;
  _minHeepPassedHademCut = minHeepPassedHademCut;
  _minHeepPassedSigmaIEtaIEtaCut = minHeepPassedSigmaIEtaIEtaCut;
  _minHeepPassed2by5Over5By5Cut = minHeepPassed2by5Over5By5Cut;
  _minElecRelIsoCut = minElecRelIsoCut;
  _maxElecRelIsoCut = maxElecRelIsoCut;
  _minElecTrkIsoCut = minElecTrkIsoCut;
  _maxElecTrkIsoCut = maxElecTrkIsoCut;
  _minHeepPassedEcalDrivenCut = minHeepPassedEcalDrivenCut;
  _minElecEtCut = minElecEtCut; 
  _minHeepPassedEcalHad1IsoCut = minHeepPassedEcalHad1IsoCut;
  _minHeepPassedTrkIsoCut = minHeepPassedTrkIsoCut;
  _minHeepPassedMissingHitsCut = minHeepPassedMissingHitsCut;
  _minHeepPassedElecDxyCut = minHeepPassedElecDxyCut;
}

void EMuAnalysis::SetImpactParameterCuts(double minElecIpOverElecIpErrorCut, double minMuonIpOverMuonIpErrorCut, double minCombinedIpOverIpErrorCut, double maxCombinedIpOverIpErrorCut, double minElecMuonDCAOverError2DCut,
					 double maxElecMuonDCAOverError2DCut, double minElecMuonDCAOverError3DCut, double maxElecMuonDCAOverError3DCut, double minLifetimeCut){

  _minElecIpOverElecIpErrorCut = minElecIpOverElecIpErrorCut;
  _minMuonIpOverMuonIpErrorCut = minMuonIpOverMuonIpErrorCut;
  _minCombinedIpOverIpErrorCut = minCombinedIpOverIpErrorCut;
  _maxCombinedIpOverIpErrorCut = maxCombinedIpOverIpErrorCut;
  _minElecMuonDCAOverError2DCut = minElecMuonDCAOverError2DCut;
  _maxElecMuonDCAOverError2DCut = maxElecMuonDCAOverError2DCut;
  _minElecMuonDCAOverError3DCut = minElecMuonDCAOverError3DCut;
  _maxElecMuonDCAOverError3DCut = maxElecMuonDCAOverError3DCut;
  _minLifetimeCut = minLifetimeCut;

}

void EMuAnalysis::SetTopologyCuts(double elecMuonMassCut, double minElecMuonCosDPhiCut, double maxElecMuonCosDPhiCut, double elecMuonDeltaRCut, 
  		       double minElecMuonMetMassCut, double maxElecMuonMetMassCut, double elecMuonDeltaPtMassCut, double elecMetMtCut, double minMuonMetMtCut, double maxMuonMetMtCut, double elecMuonPZetaVisCut,
		       double minElecMuonPZetaCut, double maxElecMuonPZetaCut, double pZetaSlope, double minNBtagsHiEffTrkCntCut, double maxNBtagsHiEffTrkCntCut, double minNBtagsHiPurityTrkCntCut,
		       double maxNBtagsHiPurityTrkCntCut, double minNBTagsHiEffSimpleSecVtxCut, double maxNBTagsHiEffSimpleSecVtxCut, double minNBTagsHiPuritySimpleSecVtxCut, 
		       double maxNBTagsHiPuritySimpleSecVtxCut, double minNBTagsCombSecVtxCut, double maxNBTagsCombSecVtxCut,
		       double leptonMetCosDPhiCut, double maxLeptonDeltaPtCosDPhiCut, double MEtCut, double nJetsCut, double minChargeCut, double maxChargeCut, double jetSumEtCut, double minElecMuonDeltaPtCut, 
		       double min2xPLargerLeptonCut, double minElecEOverPCut, double maxElecEOverPCut, double maxThreeLeptonCut){
  _elecMuonMassCut = elecMuonMassCut;
  _minElecMuonCosDPhiCut = minElecMuonCosDPhiCut;
  _maxElecMuonCosDPhiCut = maxElecMuonCosDPhiCut;  
  _elecMuonDeltaRCut = elecMuonDeltaRCut;
  _minElecMuonMetMassCut = minElecMuonMetMassCut;
  _maxElecMuonMetMassCut = maxElecMuonMetMassCut;
  _elecMuonDeltaPtMassCut = elecMuonDeltaPtMassCut;
  _elecMetMtCut = elecMetMtCut;
  _minMuonMetMtCut = minMuonMetMtCut;
  _maxMuonMetMtCut = maxMuonMetMtCut;
  _elecMuonPZetaVisCut = elecMuonPZetaVisCut;
  _minElecMuonPZetaCut = minElecMuonPZetaCut;
  _maxElecMuonPZetaCut = maxElecMuonPZetaCut;
  _pZetaSlope = pZetaSlope;
  _minNbTagsHiEffTrkCntCut = minNBtagsHiEffTrkCntCut;
  _maxNbTagsHiEffTrkCntCut = maxNBtagsHiEffTrkCntCut;
  _minNbTagsHiPurityTrkCntCut = minNBtagsHiPurityTrkCntCut;
  _maxNbTagsHiPurityTrkCntCut = maxNBtagsHiPurityTrkCntCut;
  _minNbTagsHiEffSimpleSecVtxCut = minNBTagsHiEffSimpleSecVtxCut;
  _maxNbTagsHiEffSimpleSecVtxCut = maxNBTagsHiEffSimpleSecVtxCut;
  _minNbTagsHiPuritySimpleSecVtxCut = minNBTagsHiPuritySimpleSecVtxCut;
  _maxNbTagsHiPuritySimpleSecVtxCut = maxNBTagsHiPuritySimpleSecVtxCut;
  _minNbTagsCombSecVtxCut = minNBTagsCombSecVtxCut;
  _maxNbTagsCombSecVtxCut = maxNBTagsCombSecVtxCut;
  _leptonMetCosDPhiCut = leptonMetCosDPhiCut;
  _maxLeptonDeltaPtCosDPhiCut = maxLeptonDeltaPtCosDPhiCut;
  _nJetsCut = nJetsCut;
  _metCut = MEtCut;
  _minChargeCut = minChargeCut;
  _maxChargeCut = maxChargeCut;
  _jetSumEtCut = jetSumEtCut;
  _minElecMuonDeltaPtCut = minElecMuonDeltaPtCut;
  _min2xPLargerLeptonCut = min2xPLargerLeptonCut;
  _minElecEOverPCut = minElecEOverPCut;
  _maxElecEOverPCut = maxElecEOverPCut;
  _maxThreeLeptonCut = maxThreeLeptonCut;
  		       
	
		       
}

/*
void EMuAnalysis::SetTopologyCuts(double elecMuonMassCut, double minElecMuonCosDPhiCut, double maxElecMuonCosDPhiCut, double elecMuonDeltaRCut, 
  		       double elecMuonMetMassCut, double elecMetMtCut, double minMuonMetMtCut, double maxMuonMetMtCut, double pZetaCentroidCut,
		       double pZetaVisCentroidCut, double pZetaSemiMajorAxisCut, double pZetaSemiMinorAxisCut, double pZetaEllipseAngleCut, double minNBtagsHiEffTrkCntCut, 
		       double maxNBtagsHiEffTrkCntCut, double minNBtagsHiPurityTrkCntCut, double maxNBtagsHiPurityTrkCntCut, double minNBTagsHiEffSimpleSecVtxCut, 
		       double maxNBTagsHiEffSimpleSecVtxCut, double minNBTagsHiPuritySimpleSecVtxCut, double maxNBTagsHiPuritySimpleSecVtxCut, double minNBTagsCombSecVtxCut, 
		       double maxNBTagsCombSecVtxCut, double leptonMetCosDPhiCut, double MEtCut, double nJetsCut, double minChargeCut, double maxChargeCut, double jetSumEtCut){
		       
  _elecMuonMassCut = elecMuonMassCut;
  _minElecMuonCosDPhiCut = minElecMuonCosDPhiCut;
  _maxElecMuonCosDPhiCut = maxElecMuonCosDPhiCut;  
  _elecMuonDeltaRCut = elecMuonDeltaRCut;
  _elecMuonMetMassCut = elecMuonMetMassCut;
  _elecMetMtCut = elecMetMtCut;
  _minMuonMetMtCut = minMuonMetMtCut;
  _maxMuonMetMtCut = maxMuonMetMtCut;
  _pZetaVisCentroidCut = pZetaVisCentroidCut;
  _pZetaCentroidCut = pZetaCentroidCut;
  _pZetaSemiMajorAxisCut = pZetaSemiMajorAxisCut;
  _pZetaSemiMinorAxisCut = pZetaSemiMinorAxisCut;
  _pZetaEllipseAngleCut = pZetaEllipseAngleCut;
  _minNbTagsHiEffTrkCntCut = minNBtagsHiEffTrkCntCut;
  _maxNbTagsHiEffTrkCntCut = maxNBtagsHiEffTrkCntCut;
  _minNbTagsHiPurityTrkCntCut = minNBtagsHiPurityTrkCntCut;
  _maxNbTagsHiPurityTrkCntCut = maxNBtagsHiPurityTrkCntCut;
  _minNbTagsHiEffSimpleSecVtxCut = minNBTagsHiEffSimpleSecVtxCut;
  _maxNbTagsHiEffSimpleSecVtxCut = maxNBTagsHiEffSimpleSecVtxCut;
  _minNbTagsHiPuritySimpleSecVtxCut = minNBTagsHiPuritySimpleSecVtxCut;
  _maxNbTagsHiPuritySimpleSecVtxCut = maxNBTagsHiPuritySimpleSecVtxCut;
  _minNbTagsCombSecVtxCut = minNBTagsCombSecVtxCut;
  _maxNbTagsCombSecVtxCut = maxNBTagsCombSecVtxCut;
  _leptonMetCosDPhiCut = leptonMetCosDPhiCut;
  _nJetsCut = nJetsCut;
  _metCut = MEtCut;
  _minChargeCut = minChargeCut;
  _maxChargeCut = maxChargeCut;
  _jetSumEtCut = jetSumEtCut;
  		       
	
		       
}
*/
void EMuAnalysis::SetSignalXSection(double theXSection){
  _signalXSection = theXSection;
}

void EMuAnalysis::SetLumi(double theLumi){
  _theLumi = theLumi;
}

void EMuAnalysis::SetNEvents(int nEvents){
  _nEvents = nEvents;
}

void EMuAnalysis::SetSource(std::string theSource){
  _source = theSource;
}

void EMuAnalysis::SetOutputLogFileName(std::string theLogFileName){
  _outLogFile = theLogFileName;
}

void EMuAnalysis::SetOutputRootFileName(std::string theRootFileName){
  _outRootFileName = theRootFileName;
}

void EMuAnalysis::getEventCounters(){
  //For nonzero candidates, increment candidate by PUWeight;
  if (_matchedCounter > 0) _nMatched += _eventPUWeight;
  if (_elecExistsCounter > 0) _nElecExists += _eventPUWeight;  
  if(_elecPtCounter > 0) _nElecPt += _eventPUWeight;
  if(_elecEtaCounter > 0) _nElecEta += _eventPUWeight;
  if(_elecAccCounter > 0) _nElecAcc += _eventPUWeight;
  if (_elecDeltaPhiInCounter > 0) _nElecDeltaPhiIn += _eventPUWeight;
  if (_elecDeltaEtaInCounter > 0) _nElecDeltaEtaIn += _eventPUWeight;
  if (_elecSigmaIEtaIEtaCounter  > 0) _nElecSigmaIEtaIEta += _eventPUWeight;
  if (_elecTrkChiSqrdCounter > 0) _nElecTrkChiSqrd += _eventPUWeight;
  if (_elecHadFracCounter  > 0) _nElecHadFrac += _eventPUWeight;
  if (_elecEOverPCounter  > 0) _nElecEOverP += _eventPUWeight;
  if (_elecRelIsoCounter > 0) _nElecRelIso += _eventPUWeight;
  if (_elecTrkIsoCounter > 0) _nElecTrkIso += _eventPUWeight;
  if (_elecMissHitsCounter  > 0) _nElecMissHits += _eventPUWeight;
  //if (_elecPFRelIsoCounter > 0) _nElecPFRelIso  += _eventPUWeight;
  //if (_elecGlobalIdCounter > 0) _nElecGlobalId += _eventPUWeight;
  //if (_elecIdCounter > 0) _nElecId += _eventPUWeight;  
  if (_elecMotherIdCounter > 0) _nElecMotherId += _eventPUWeight;
  if (_elecHOverECounter > 0) _nElecHOverE += _eventPUWeight;
  if (_elec2x5Over5x5Counter > 0) _nElec2x5Over5x5 += _eventPUWeight;
  if (_elecEcalHad1IsoCounter > 0) _nElecEcalHad1Iso += _eventPUWeight;
  if (_elecDxyCounter > 0) _nElecDxy += _eventPUWeight;
  
  if (_muonPtCounter > 0) _nMuonPt += _eventPUWeight;
  if (_muonECounter > 0) _nMuonE += _eventPUWeight;
  if (_muonEtaCounter > 0) _nMuonEta += _eventPUWeight;
  if (_muonAccCounter > 0) _nMuonAcc += _eventPUWeight;  
  if (_muonDxyVtxCounter > 0) _nMuonDxyVtx += _eventPUWeight;
  if (_muonDzVtxCounter > 0) _nMuonDzVtx += _eventPUWeight;
  if (_muonIsGlobalMuonCounter > 0) _nMuonIsGlobalMuon += _eventPUWeight;
  if (_muonMatchedCounter > 0) _nMuonMatched += _eventPUWeight;
  if (_muonMotherIdCounter > 0) _nMuonMotherId += _eventPUWeight;
  if (_muonNormChiSqrdCounter > 0) _nMuonNormChiSqrd += _eventPUWeight;
  if (_muonTrkChiSqrdCounter > 0) _nMuonTrkChiSqrd += _eventPUWeight;
  if (_muonValidHitsCounter > 0) _nMuonValidHits += _eventPUWeight;
  if (_muonMatchedStationsCounter > 0) _nMuonMatchedStations += _eventPUWeight;
  if (_muonPixelHitsCounter > 0) _nMuonPixelHits += _eventPUWeight;
  if (_muonTrkLayersWithHitsCounter > 0) _nMuonTrkLayersWithHits += _eventPUWeight;
  //if (_muonIsoDR03NoDBCounter > 0) _nMuonIsoDR03NoDB += _eventPUWeight;
  //if (_muonIsoDR03DBCounter > 0) _nMuonIsoDR03DB += _eventPUWeight;
  //if (_muonIsoDR04NoDBCounter > 0) _nMuonIsoDR04NoDB += _eventPUWeight;
  if (_muonIsoDR04DBCounter > 0) _nMuonIsoDR04DB += _eventPUWeight;
  if (_muonIsHighPtMuonCounter > 0) _nMuonIsHighPtMuon += _eventPUWeight;
  if (_muonDptOverPtCounter > 0) _nMuonDptOverPt += _eventPUWeight;
  if (_muonIdCounter > 0) _nMuonId += _eventPUWeight;
  

  if (_elecMatchedCounter > 0) _nElecMatched += _eventPUWeight;
  if (_heepPassedEtCounter > 0) _nHeepPassedEt += _eventPUWeight;
  if (_heepPassedCrackCounter > 0) _nHeepPassedCrack += _eventPUWeight;
  if (_heepPassedDEtaInCounter > 0) _nHeepPassedDEtaIn += _eventPUWeight;
  if (_heepPassedDPhiInCounter > 0) _nHeepPassedDPhiIn += _eventPUWeight;
  if (_heepPassedHademCounter > 0) _nHeepPassedHadem += _eventPUWeight;
  if (_heepPassedSigmaIEtaIEtaCounter > 0) _nHeepPassedSigmaIEtaIEta += _eventPUWeight;
  if (_heepPassed2by5Over5By5Counter > 0) _nHeepPassed2by5Over5By5 += _eventPUWeight;
  if (_heepPassedEcalHad1IsoCounter > 0) _nHeepPassedEcalHad1Iso += _eventPUWeight;
  if (_heepPassedEcalDrivenCounter > 0) _nHeepPassedEcalDriven += _eventPUWeight;
  if (_heepPassedTrkIsoCounter > 0) _nHeepPassedTrkIso += _eventPUWeight;
  if (_heepPassedMissingHitsCounter > 0) _nHeepPassedMissingHits += _eventPUWeight;
  if (_heepPassedElecDxyCounter > 0) _nHeepPassedElecDxy += _eventPUWeight;
  if (_elecEtCounter > 0) _nElecEt += _eventPUWeight;    
  if (_heepIdCounter > 0) _nHeepId += _eventPUWeight;  
  
  if (_elecIpOverElecIpErrorCounter > 0) _nElecIpOverElecIpError += _eventPUWeight;
  if (_muonIpOverMuonIpErrorCounter > 0) _nMuonIpOverMuonIpError += _eventPUWeight;
  if (_combinedIpOverIpErrorCounter > 0) _nCombinedIpOverIpError += _eventPUWeight;
  if (_elecMuonDCAOverError2DCounter > 0) _nElecMuonDCAOverError2D += _eventPUWeight;
  if (_elecMuonDCAOverError3DCounter > 0) _nElecMuonDCAOverError3D += _eventPUWeight;
  if (_lifetimeCounter > 0) _nLifetime += _eventPUWeight;
  
  if (_elecMuonMassCounter > 0) _nElecMuonMass += _eventPUWeight;
  if (_elecMuonCosDPhiCounter > 0) _nElecMuonCosDPhi += _eventPUWeight;
  if (_elecMuonDeltaRCounter > 0) _nElecMuonDeltaR += _eventPUWeight;
  if (_elecMuonMetMassCounter > 0) _nElecMuonMetMass += _eventPUWeight;
  if (_elecMuonDeltaPtMassCounter > 0) _nElecMuonDeltaPtMass += _eventPUWeight; 
  if (_elecMetMtCounter > 0) _nElecMetMt += _eventPUWeight;
  if (_muonMetMtCounter > 0) _nMuonMetMt += _eventPUWeight;
  if (_elecMuonPZetaVisCounter > 0) _nElecMuonPZetaVis += _eventPUWeight;
  if (_elecMuonPZetaCounter > 0) _nElecMuonPZeta += _eventPUWeight;
  if (_nbTagsHiEffTrkCntCounter > 0) _nNbTagsHiEffTrkCnt += _eventPUWeight;
  if (_nbTagsHiPurityTrkCntCounter > 0) _nNbTagsHiPurityTrkCnt += _eventPUWeight;
  if (_nbTagsHiEffSimpleSecVtxCounter > 0) _nNbTagsHiEffSimpleSecVtx += _eventPUWeight;
  if (_nbTagsHiPuritySimpleSecVtxCounter > 0) _nNbTagsHiPuritySimpleSecVtx += _eventPUWeight;
  if (_nbTagsCombSecVtxCounter > 0) _nNbTagsCombSecVtx += _eventPUWeight;
  if (_leptonMetCosDPhiCounter > 0) _nLeptonMetCosDPhi += _eventPUWeight;
  if (_leptonDeltaPtCosDPhiCounter > 0) _nLeptonDeltaPtCosDPhi += _eventPUWeight;
  if (_metCounter > 0) _nMEt += _eventPUWeight;
  if (_threeLeptonCounter > 0) _nThreeLepton += _eventPUWeight;
  if (_nJetsCounter > 0) _nNJets += _eventPUWeight;
  if (_jetSumEtCounter > 0) _nJetSumEt += _eventPUWeight;
  if (_chargeCounter > 0) _nCharge += _eventPUWeight;
  if (_elecMuonDeltaPtCounter > 0) _nElecMuonDeltaPt += _eventPUWeight;
  if (_2xPLargerLeptonCounter > 0) _n2xPLargerLepton += _eventPUWeight;
  if (_topologyCounter > 0) _nTopology += _eventPUWeight;
  
  if (_lowPtLeptonAntiIsoCounter > 0) _nLowPtLeptonAntiIso += _eventPUWeight;
  
  if (_QCDCR1Counter > 0) _nQCDCR1 += _eventPUWeight;
  if (_QCDCR2Counter > 0) _nQCDCR2 += _eventPUWeight;
  if (_QCDOSCounter > 0) _nQCDOS += _eventPUWeight;
  if (_QCDElectronIsoCounter > 0) _nQCDElectronIso += _eventPUWeight;
  if (_QCDTopologyCounter > 0) _nQCDTopology += _eventPUWeight;
  if (_QCDMuonIsoCounter > 0) _nQCDMuonIso += _eventPUWeight;
  if (_QCDMuonAntiIsoCounter > 0) _nQCDMuonAntiIso += _eventPUWeight;
  
  
  if (_TTJetsCR1Counter > 0) _nTTJetsCR1 += _eventPUWeight;
  if (_TTJetsCR2Counter > 0) _nTTJetsCR2 += _eventPUWeight;
  if (_TTJetsJetSumEtCounter > 0) _nTTJetsJetSumEt += _eventPUWeight;
  if (_TTJetsTopologyCounter > 0) _nTTJetsTopology += _eventPUWeight;
  if (_TTJetsNoBTagsCounter > 0) _nTTJetsNoBTags += _eventPUWeight;
  if (_TTJetsAtLeastOneBTagCounter > 0) _nTTJetsAtLeastOneBTag += _eventPUWeight;
  
  if (_WJetsCR1Counter > 0) _nWJetsCR1 += _eventPUWeight;
  if (_WJetsCR2Counter > 0) _nWJetsCR2 += _eventPUWeight;
  if (_WJetsOSCounter > 0) _nWJetsOS += _eventPUWeight;
  if (_WJetsTopologyCounter > 0) _nWJetsTopology += _eventPUWeight;
  if (_WJetsElectronIsoCounter > 0) _nWJetsElectronIso += _eventPUWeight;
  if (_WJetsElectronAntiIsoCounter > 0) _nWJetsElectronAntiIso += _eventPUWeight;
  if (_WJetsMtCounter > 0) _nWJetsMt += _eventPUWeight; 
  
  if (_DYToTauTauCR1Counter > 0) _nDYToTauTauCR1 += _eventPUWeight;
  if (_DYToTauTauCR2Counter > 0) _nDYToTauTauCR2 += _eventPUWeight;  
  
  //N-1
  
  if (_nMinus1ElecPtCounter > 0) _nElecPtNMinus1 += _eventPUWeight;
  if (_nMinus1ElecEtaCounter > 0) _nElecEtaNMinus1 += _eventPUWeight; 
  if (_nMinus1ElecTrkChiSqrdCounter > 0) _nElecTrkChiSqrdNMinus1 += _eventPUWeight;   
  if (_nMinus1MuonPtCounter > 0) _nMuonPtNMinus1 += _eventPUWeight;
  if (_nMinus1MuonECounter > 0) _nMuonENMinus1 += _eventPUWeight;
  if (_nMinus1MuonEtaCounter > 0) _nMuonEtaNMinus1 += _eventPUWeight;
  if (_nMinus1MuonDxyVtxCounter > 0) _nMuonDxyVtxNMinus1 += _eventPUWeight;
  if (_nMinus1MuonDzVtxCounter > 0) _nMuonDzVtxNMinus1 += _eventPUWeight;
  if (_nMinus1MuonDzBSCounter > 0) _nMuonDzBSNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIsGlobalMuonCounter > 0) _nMuonIsGlobalMuonNMinus1 += _eventPUWeight;
  if (_nMinus1MuonMatchedCounter > 0) _nMuonMatchedNMinus1 += _eventPUWeight;
  if (_nMinus1MuonNormChiSqrdCounter > 0) _nMuonNormChiSqrdNMinus1 += _eventPUWeight;
  if (_nMinus1MuonTrkChiSqrdCounter > 0) _nMuonTrkChiSqrdNMinus1 += _eventPUWeight;
  if (_nMinus1MuonValidHitsCounter > 0) _nMuonValidHitsNMinus1 += _eventPUWeight;
  if (_nMinus1MuonMatchedStationsCounter > 0) _nMuonMatchedStationsNMinus1 += _eventPUWeight;
  if (_nMinus1MuonPixelHitsCounter > 0) _nMuonPixelHitsNMinus1 += _eventPUWeight;
  if (_nMinus1MuonTrkLayersWithHitsCounter > 0) _nMuonTrkLayersWithHitsNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIsoDR03NoDBCounter > 0) _nMuonIsoDR03NoDBNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIsoDR03DBCounter > 0) _nMuonIsoDR03DBNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIsoDR04NoDBCounter > 0) _nMuonIsoDR04NoDBNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIsoDR04DBCounter > 0) _nMuonIsoDR04DBNMinus1 += _eventPUWeight;  
  if (_nMinus1MuonIsoCounter > 0) _nMuonIsoNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIsHighPtMuonCounter > 0) _nMuonIsHighPtMuonNMinus1 += _eventPUWeight;
  if (_nMinus1MuonDptOverPtCounter > 0) _nMuonDptOverPtNMinus1 += _eventPUWeight;
  
  
  if (_nMinus1ElecMatchedCounter > 0) _nElecMatchedNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedPtCounter > 0) _nHeepPassedPtNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedDetEtaCounter > 0) _nHeepPassedDetEtaNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedCrackCounter > 0) _nHeepPassedCrackNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedDEtaInCounter > 0) _nHeepPassedDEtaInNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedDPhiInCounter > 0) _nHeepPassedDPhiInNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedHademCounter > 0) _nHeepPassedHademNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedSigmaIEtaIEtaCounter > 0) _nHeepPassedSigmaIEtaIEtaNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassed2by5Over5By5Counter > 0) _nHeepPassed2by5Over5By5NMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedEcalHad1IsoCounter > 0) _nHeepPassedEcalHad1IsoNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedTrkIsoCounter > 0) _nHeepPassedTrkIsoNMinus1 += _eventPUWeight;
  if (_nMinus1HeepPassedEcalDrivenCounter > 0) _nHeepPassedEcalDrivenNMinus1 += _eventPUWeight;
  if (_nMinus1ElecEtCounter > 0) _nElecEtNMinus1 += _eventPUWeight;  
  if (_nMinus1ElecRelIsoCounter > 0) _nElecRelIsoNMinus1 += _eventPUWeight;     
  
  if (_nMinus1ElecIpOverElecIpErrorCounter > 0) _nElecIpOverElecIpErrorNMinus1 += _eventPUWeight;
  if (_nMinus1MuonIpOverMuonIpErrorCounter > 0) _nMuonIpOverMuonIpErrorNMinus1 += _eventPUWeight;
  if (_nMinus1CombinedIpOverIpErrorCounter > 0) _nCombinedIpOverIpErrorNMinus1 += _eventPUWeight;
  if (_nMinus1ElecMuonDCAOverError2DCounter > 0) _nElecMuonDCAOverError2DNMinus1 += _eventPUWeight;
  if (_nMinus1ElecMuonDCAOverError3DCounter > 0) _nElecMuonDCAOverError3DNMinus1 += _eventPUWeight;
  
  if(_nMinus1MEtCounter > 0) _nMEtNMinus1 += _eventPUWeight;
  if(_nMinus1BTagCounter > 0) _nBTagsNMinus1 += _eventPUWeight;
  if(_nMinus1ElecMuonCosDPhiCounter > 0) _nElecMuonCosDPhiNMinus1 += _eventPUWeight;
  if(_nMinus1ElecMuonPZetaCounter > 0) _nElecMuonPZetaNMinus1 += _eventPUWeight;
  if(_nMinus1ElecMuonPZetaVisCounter > 0) _nElecMuonPZetaVisNMinus1 += _eventPUWeight;
  if(_nMinus1JetSumEtCounter > 0) _nJetSumEtNMinus1 += _eventPUWeight;
  if(_nMinus1LeptonMetCosDPhiCounter > 0) _nLeptonMetCosDPhiNMinus1 += _eventPUWeight;
  if(_nMinus1ElecMuonMetMassCounter > 0) _nElecMuonMetMassNMinus1 += _eventPUWeight;
  if(_nMinus1NJetsCounter > 0) _nNJetsNMinus1 += _eventPUWeight;
  if(_nMinus1ElecMuonMassCounter > 0) _nElecMuonMassNMinus1 += _eventPUWeight;
  if(_nMinus1MuonMetMtCounter > 0) _nMuonMetMtNMinus1 += _eventPUWeight;
  if(_nMinus1ElecMuonChargeCounter > 0) _nElecMuonChargeNMinus1 += _eventPUWeight;  
  
}

void EMuAnalysis::resetCandCounters(){
  //counters set to zero
  _elecExistsCounter = 0;
  _elecPtCounter = 0;
  _elecEtaCounter = 0;
  _elecAccCounter = 0;
  _elecDeltaPhiInCounter = 0;
  _elecDeltaEtaInCounter = 0;
  _elecSigmaIEtaIEtaCounter = 0;
  _elecTrkChiSqrdCounter = 0;
  _elecHadFracCounter =0;
  _elecEOverPCounter = 0;
  _elecRelIsoCounter = 0;
  _elecTrkIsoCounter = 0;
  _elecMissHitsCounter = 0;
  _elecPFRelIsoCounter = 0;
  _elecGlobalIdCounter = 0;
  _elecIdCounter = 0;  
  _elecMotherIdCounter = 0;
  _elecHOverECounter = 0;
  _elec2x5Over5x5Counter = 0;
  _elecEcalHad1IsoCounter = 0;
  _elecDxyCounter = 0;  

  _muonPtCounter = 0;
  _muonECounter = 0;
  _muonEtaCounter = 0;
  _muonAccCounter = 0;
  _muonDxyVtxCounter = 0;
  _muonDzVtxCounter = 0;
  _muonIsGlobalMuonCounter = 0;
  _muonMatchedCounter = 0;
  _muonMotherIdCounter = 0;
  _muonNormChiSqrdCounter = 0;
  _muonTrkChiSqrdCounter = 0;
  _muonValidHitsCounter = 0;
  _muonMatchedStationsCounter = 0;
  _muonPixelHitsCounter = 0;
  _muonTrkLayersWithHitsCounter = 0;
  _muonIsoDR03NoDBCounter = 0;
  _muonIsoDR03DBCounter = 0;
  _muonIsoDR04NoDBCounter = 0;
  _muonIsoDR04DBCounter = 0;
  _muonIsHighPtMuonCounter = 0;
  _muonDptOverPtCounter = 0;
  _muonIdCounter = 0;
  
  _elecMatchedCounter = 0;
  _heepPassedEtCounter = 0;
  _heepPassedCrackCounter = 0;
  _heepPassedDEtaInCounter = 0;
  _heepPassedDPhiInCounter = 0;
  _heepPassedHademCounter = 0;
  _heepPassedSigmaIEtaIEtaCounter = 0;
  _heepPassed2by5Over5By5Counter = 0;
  _heepPassedEcalHad1IsoCounter = 0;
  _heepPassedTrkIsoCounter = 0;
  _heepPassedEcalDrivenCounter = 0;
  _heepPassedMissingHitsCounter = 0;
  _heepPassedElecDxyCounter = 0;
  _elecEtCounter = 0;  
  _heepIdCounter = 0;
  
  _matchedCounter = 0;
  
  _elecIpOverElecIpErrorCounter = 0;
  _muonIpOverMuonIpErrorCounter = 0;
  _combinedIpOverIpErrorCounter = 0;
  _elecMuonDCAOverError2DCounter = 0;
  _elecMuonDCAOverError3DCounter = 0;
  _lifetimeCounter = 0;
  
  _elecMuonMassCounter = 0;
  _elecMuonCosDPhiCounter = 0;
  _elecMuonDeltaRCounter = 0;
  _elecMuonMetMassCounter = 0;
  _elecMuonDeltaPtMassCounter = 0;
  _elecMetMtCounter = 0;
  _muonMetMtCounter = 0;
  _elecMuonPZetaVisCounter = 0;
  _elecMuonPZetaCounter = 0;
  _nbTagsHiEffTrkCntCounter = 0;
  _nbTagsHiPurityTrkCntCounter = 0;
  _nbTagsHiEffSimpleSecVtxCounter = 0;
  _nbTagsHiPuritySimpleSecVtxCounter = 0;
  _nbTagsCombSecVtxCounter = 0;
  _leptonMetCosDPhiCounter = 0;  
  _leptonDeltaPtCosDPhiCounter = 0;
  _metCounter = 0;
  _threeLeptonCounter = 0;
  _nJetsCounter = 0;
  _jetSumEtCounter = 0;
  _chargeCounter = 0;
  _elecMuonDeltaPtCounter = 0;
  _2xPLargerLeptonCounter = 0;
  _topologyCounter = 0;
    
  _lowPtLeptonAntiIsoCounter = 0;  
    
  _QCDCR1Counter = 0;
  _QCDCR2Counter = 0;
  _QCDOSCounter = 0;
  _QCDElectronIsoCounter = 0;
  _QCDTopologyCounter = 0;
  _QCDMuonIsoCounter = 0;
  _QCDMuonAntiIsoCounter = 0;
  
  _TTJetsCR1Counter = 0;
  _TTJetsCR2Counter = 0;
  _WJetsCR1Counter = 0;
  _WJetsCR2Counter = 0;
  _DYToTauTauCR1Counter = 0;
  _DYToTauTauCR2Counter = 0;

  _TTJetsJetSumEtCounter = 0;
  _TTJetsTopologyCounter = 0;
  _TTJetsNoBTagsCounter = 0;
  _TTJetsAtLeastOneBTagCounter = 0;
  
  _WJetsOSCounter = 0;
  _WJetsTopologyCounter = 0;
  _WJetsElectronIsoCounter = 0;
  _WJetsElectronAntiIsoCounter = 0;
  _WJetsMtCounter = 0;
  
    
  //N-1 counters set to zero
  _nMinus1ElecPtCounter = 0;
  _nMinus1ElecEtaCounter = 0;
  _nMinus1ElecTrkChiSqrdCounter = 0;
  _nMinus1MuonPtCounter = 0;
  _nMinus1MuonECounter = 0;
  _nMinus1MuonEtaCounter = 0;
  _nMinus1MuonDxyVtxCounter = 0;
  _nMinus1MuonDzVtxCounter = 0;
  _nMinus1MuonDzBSCounter = 0;
  _nMinus1MuonIsGlobalMuonCounter = 0;
  _nMinus1MuonMatchedCounter = 0;
  _nMinus1MuonNormChiSqrdCounter = 0;
  _nMinus1MuonTrkChiSqrdCounter = 0;
  _nMinus1MuonValidHitsCounter = 0;
  _nMinus1MuonMatchedStationsCounter = 0;
  _nMinus1MuonPixelHitsCounter = 0;
  _nMinus1MuonTrkLayersWithHitsCounter = 0;
  _nMinus1MuonIsoDR03NoDBCounter = 0;
  _nMinus1MuonIsoDR03DBCounter = 0;
  _nMinus1MuonIsoDR04NoDBCounter = 0;
  _nMinus1MuonIsoDR04DBCounter = 0;  
  _nMinus1MuonIsHighPtMuonCounter = 0;
  _nMinus1MuonDptOverPtCounter = 0;
  _nMinus1MuonIsoCounter = 0;
  
  _nMinus1ElecMatchedCounter = 0;
  _nMinus1HeepPassedPtCounter = 0;
  _nMinus1HeepPassedDetEtaCounter = 0;
  _nMinus1HeepPassedCrackCounter = 0;
  _nMinus1HeepPassedDEtaInCounter = 0;
  _nMinus1HeepPassedDPhiInCounter = 0;
  _nMinus1HeepPassedHademCounter = 0;
  _nMinus1HeepPassedSigmaIEtaIEtaCounter = 0;
  _nMinus1HeepPassed2by5Over5By5Counter = 0;
  _nMinus1HeepPassedEcalHad1IsoCounter = 0;
  _nMinus1HeepPassedTrkIsoCounter = 0;
  _nMinus1HeepPassedEcalDrivenCounter = 0;
  _nMinus1ElecEtCounter = 0;    
  _nMinus1ElecRelIsoCounter = 0;
  
  _nMinus1ElecIpOverElecIpErrorCounter = 0;
  _nMinus1MuonIpOverMuonIpErrorCounter = 0;
  _nMinus1CombinedIpOverIpErrorCounter = 0;
  _nMinus1ElecMuonDCAOverError2DCounter = 0;
  _nMinus1ElecMuonDCAOverError3DCounter = 0;
 
  _nMinus1MEtCounter = 0;
  _nMinus1BTagCounter = 0;
  _nMinus1ElecMuonCosDPhiCounter = 0;
  _nMinus1ElecMuonPZetaCounter = 0;
  _nMinus1ElecMuonPZetaVisCounter = 0;
  _nMinus1JetSumEtCounter = 0;
  _nMinus1ElecMuonMassCounter = 0;
  _nMinus1MuonMetMtCounter = 0;
  _nMinus1LeptonMetCosDPhiCounter = 0;
  _nMinus1ElecMuonMetMassCounter = 0;
  _nMinus1NJetsCounter = 0;
  _nMinus1ElecMuonChargeCounter = 0;
  

}

void EMuAnalysis::resetHistosDefault(){
   //reset histo values to start values
  _maxElecP = -10.;
  _maxElecPt = -10.;
  _minElecEta = 10.;
  _minElecDeltaPhiInBarrel = 2.;		   
  _minElecDeltaEtaInBarrel = 1.;		   
  _minElecSigmaIEtaBarrel = 1.;
  _minElecTrkChiSqrd = 10000.;	
  _maxElecTrkNdof = -5;	   
  _minElecHadFracBarrel = 10.;		   
  _minElecEOverPBarrel = 10.;		   
  _minElecDeltaPhiInEndCap = 2.;		   
  _minElecDeltaEtaInEndCap = 1.;		   
  _minElecSigmaIEtaEndCap = 1.;		   
  _minElecHadFracEndCap = 10.;		   
  _minElecEOverPEndCap = 10.;
  //_minElecRelIso = 50.;
  _minElecTrkIso = 50.;		   		   
  _maxElecMissingHits = -10;	   
  _minPFIsoChargedHadDr03 = 1000.;	   
  _minPFIsoNeutralHadDr03 = 1000.;	   
  _minPFIsoPhotonDr03 = 1000.;		   
  _minPFRelIsoDr03 = 1000.;		   
  _minPFIsoChargedHadDr04 = 1000.;	   
  _minPFIsoNeutralHadDr04 = 1000.;	   
  _minPFIsoPhotonDr04 = 1000.;		   
  _minPFRelIsoDr04 = 1000.;		   
	
  _maxTauEnergy = -5;	
  _maxMuonPOverTauE = -5;
  _maxZPrimeEnergy = -5;
  _maxVertexDiff = -5;
  _maxVertexDiffxy = -5;
  
  _maxBeamSpotDiff = -5;
  _maxBeamSpotDiffxy = -5;
  _maxBeamSpotPrimaryzDiff = -5;
	
  _maxMuonPt = -5;
  _minMuonEta = 3;
  _maxMuonDxyVtx = -5;
  _maxMuonDzVtx = -5;
  _maxMuonNormChiSqrd = -5;
  _minMuonTrkChiSqrd = 10000.;
  _maxMuonTrkNdof = -5;
  _maxIsGlobalMuon = -2;
  _maxMuonValidHits = -5;
  _maxMuonMatchedStations = -5;
  _maxMuonPixelHits = -5;
  _maxMuonTrkLayersWithHits = -5;
  _maxMuonIsoDR04DB = -5;
  _maxMuonIsHighPtMuon = -5;
  _minMuonDptOverPt = 5;	
	
  _maxHeepPassedEcalDriven = -5;
  _maxHeepPassedCrack = -5;
  _maxHeepPassedDEtaIn = -5;
  _maxHeepPassedDPhiIn = -5;
  _maxHeepPassedHadem = -5;
  _maxHeepPassedSigmaIEtaIEta = -5;
  _maxHeepPassed2by5Over5By5 = -5;
  _maxHeepPassedEcalHad1Iso = -5;
  _maxHeepPassedTrkIso = -5;
  _maxElecEt = -5;	
  
  _maxElecSCE1x5 = -5.;
  _maxElecSCE2x5 = -5.;
  _maxElecSCE5x5 = -5.;
  _maxElecSCE1x5Over5x5 = -5.;
  _maxElecSCE2x5MaxOver5x5 = -5.;
  _maxElecSCE1x5EOverP = -5.;
  _maxElecSCE2x5EOverP = -5.;
  _maxElecSCE5x5EOverP = -5.;
  _maxMuonE = -5.;
  
  
  _maxMetMinusDeltaPt = -500;
  _maxDeltaPtx = -500;
  _maxDeltaPty = -500;
  _maxMetPx = -500;
  _maxMetPy = -500;
  
  _maxElecIpOverElecIpError = -500;
  _maxMuonIpOverMuonIpError = -500;	
  _maxCombinedIpOverIpError = -500;
  
  _maxElecMuonDistanceOfClosestApproach = -500;
  _maxElecMuonDCACrossingPoint = -500; 
  _maxElecMuonTwoTrkMinDist = -500;
  _maxElecMuonTwoTrkCrossingPoint = -500;
  _maxElecMuonDis3D = -500;
  _maxElecMuonDis2D = -500;
  _maxElecMuonDCA3DError = -500;
  _maxElecMuonDCA2DError = -500;
  
  _maxElecMuonDis3DOverError = -500;
  _maxElecMuonDis2DOverError = -500;
  
  _maxElecIp = -500;
  _maxMuonIp = -500;
  
  _maxElecIpError = -500;
  _maxMuonIpError = -500;
  
  _maxLifetime = -500;
  
  _maxElecEOverP_hiIP = -500;
  _maxElecEta_hiIP = -500;
  _maxElecEt_hiIP = -500;
  _maxElecTrkChiSqrd_hiIP = -500;
					 					           			      
  _minCharge = 5;
  _minElecMuonCosDPhi = 10;
  _maxElecMuonDeltaR = -5;
  _maxElecMuonPZeta = -10000;
  _maxElecMuonPZetaVis = -10000;
  _minLeptonMetCosDPhi = 10;	
  _minLeptonDeltaPtCosDPhi = 10;		           							      
  _maxMEt = -1;		
  _maxThreeLepton = -1;	           							      
  _minNBtagsHiEffTrkCnt = 10;
  _maxElecMuonMetMass = -10;
  _maxElecMuonDeltaPtMass = -10;
  _maxElecMuonMetMassBeginning = -10;
  _maxElecPtVsMass = -10;
  _maxMuonPtVsMass = -10;
  
  _maxBSx0 = -5;
  _maxBSy0 = -5;
  _maxBSz0 = -5;
  
  _maxMuonVertex_x0 = -5;	
  _maxMuonVertex_y0 = -5;	
  _maxMuonVertex_z0 = -5;
  _maxMuonDecayLength_x = -5;	
  _maxMuonDecayLength_y = -5;	
  _maxMuonDecayLength_z = -5;
  _maxElecVertex_x0 = -5;	
  _maxElecVertex_y0 = -5;	
  _maxElecVertex_z0 = -5;
  _maxElecDecayLength_x = -5;	
  _maxElecDecayLength_y = -5;	
  _maxElecDecayLength_z = -5; 
  _maxElecVx = -5.; 
  _maxElecVy = -5.; 
  _maxElecVz = -5.; 
  

  _maxGenTauVtx_x = -500.;
  _maxGenTauVtx_y = -500.;
  _maxGenTauVtx_z = -500.;
  _maxOtherGenTauVtx_x = -500.;
  _maxOtherGenTauVtx_y = -500.;
  _maxOtherGenTauVtx_z = -500.;  
  
  _maxElecMotherId = -5;
  _maxMuonMotherId = -5;
  
  _maxMEtVsDeltaPt = -1;
  _maxElecPtVsDeltaPt = -10;
  _maxMuonPtVsDeltaPt = -10;
  
  _maxMuonP = -5;
  
  //N-1     
  _maxElecPtNMinus1 = -10.;
  _maxElecEtNMinus1 = -10.;
  _minElecEtaNMinus1 = 10.;  
  _maxMuonPtNMinus1 = -10.;
  _minMuonIsoNMinus1 = 100;
 // _minElecRelIsoNMinus1 = 100;
  _maxMuonIsHighPtMuonNMinus1 = -5;
  _minMuonDptOverPtNMinus1 = 5;
  _maxMuonDzVtxNMinus1 = -50000;
  _maxMuonDzBSNMinus1 = -50000;
 
  _maxElecIpOverElecIpErrorNMinus1 = -500;
  _maxMuonIpOverMuonIpErrorNMinus1 = -500;
  _maxCombinedIpOverIpErrorNMinus1 = -500;
   										     					           			      
  _minChargeNMinus1 = 5;		           							      
  _maxMEtNMinus1 = -1;		
  _maxThreeLeptonNMinus1 = -10;    
  _maxThreeLeptonMassNMinus1 = -10.;       							      
  _minNBtagsHiEffTrkCntNMinus1 = 1000;
  _minNBtagsCombSecVtxNMinus1 = 20;
  _minNBtagsCombSecVtxCR2 = 20;
  _minElecMuonCosDPhiNMinus1 = 10;
  _maxElecMuonPZetaNMinus1 = -1000;
  _maxElecMuonPZetaVisNMinus1 = -1000;
  _minJetSumEtNMinus1 = 1000;
  _minLeptonMetCosDPhiNMinus1 = 1000;
  _minLeptonDeltaPtCosDPhiNMinus1 = 1000;
  _minElecMuonMetMassNMinus1 = -100000;
  _maxElecMuonDeltaPtMassNMinus1 = -10000;
  _maxNJetsNMinus1 = -10;
  _minElecMuonMassNMinus1 = 100000.;
  _minMuonMetMtNMinus1 = 1000.;
  _minElecMuonChargeNMinus1 = 10.;    	
  _maxMEtNMinus1 = -10.;		
  
  _maxElecDeltaPtxNMinus1 = -10.;	
  _maxElecDeltaPtyNMinus1 = -10.;
  _maxMuonDeltaPtNMinus1 = -10.;
  
  _maxMuonDptOverPtNMinus1 = -10.;
  _maxElecDptOverPtNMinus1 = -10.;

}


#endif // #ifdef EMuAnalysis_cxx
