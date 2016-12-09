#include "/usr/users/ajohnson/CMS/sharedBin/include/config.h"
#include "EMuAnalysis.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

using namespace std;
int main(int argc, char **argv){
 
  Config theConfig(argv[1], argv);
  
  string inRootFileName = theConfig.pString("InputRootFile");

  string source = theConfig.pString("Source");
  double signalXSection = theConfig.pDouble("SignalXSection");
  double theLumi = theConfig.pDouble("Lumi");
  int nEvents = theConfig.pInt("NEvents");

  std::cout << "inRootFileName = " << inRootFileName << "\n";
  string logName = theConfig.pString("EffLogFile");
  string rootOutName = theConfig.pString("RootOutFile");
  
  
  double minMuonPtCut = theConfig.pDouble("MinMuonPtCut");
  //double minMuonECut = theConfig.pDouble("MinMuonECut");
  double maxMuonEtaCut = theConfig.pDouble("MaxMuonEtaCut");
  double maxMuonDxyVtxCut = theConfig.pDouble("MaxMuonDxyVtxCut");
  double maxMuonDzVtxCut = theConfig.pDouble("MaxMuonDzVtxCut");
  double minMuonIsGlobalMuonCut = theConfig.pDouble("MinMuonIsGlobalMuonCut");
  double minMuonMatchedCut = theConfig.pDouble("MinMuonMatchedCut");
  double muonMotherIdCut = theConfig.pDouble("MuonMotherIdCut");
  double maxMuonNormChiSqrdCut = theConfig.pDouble("MaxMuonNormChiSqrdCut");
  double minMuonValidHitsCut = theConfig.pDouble("MinMuonValidHitsCut");
  double minMuonMatchedStationsCut = theConfig.pDouble("MinMuonMatchedStationsCut");
  double minMuonPixelHitsCut = theConfig.pDouble("MinMuonPixelHitsCut");
  double minMuonTrkLayersWithHitsCut = theConfig.pDouble("MinMuonTrkLayersWithHitsCut");
  //double maxMuonIsoDR03NoDBCut = theConfig.pDouble("MaxMuonIsoDR03NoDBCut");
  //double maxMuonIsoDR03DBCut = theConfig.pDouble("MaxMuonIsoDR03DBCut");
  //double maxMuonIsoDR04NoDBCut = theConfig.pDouble("MaxMuonIsoDR04NoDBCut");
  double minMuonIsoDR04DBCut = theConfig.pDouble("MinMuonIsoDR04DBCut");
  double maxMuonIsoDR04DBCut = theConfig.pDouble("MaxMuonIsoDR04DBCut");
  double minMuonIsHighPtMuonCut = theConfig.pDouble("MinMuonIsHighPtMuonCut");
  double maxMuonDptOverPtCut = theConfig.pDouble("MaxMuonDptOverPtCut");   
    
  double minElecPtCut = theConfig.pDouble("MinElecPtCut");
  double maxElecEtaCut = theConfig.pDouble("MaxElecEtaCut");

  double minElecMatchedCut = theConfig.pDouble("MinElecMatchedCut");
  double minHeepPassedEtCut = theConfig.pDouble("MinHeepPassedEtCut");
  double minHeepPassedDetEtaCut = theConfig.pDouble("MinHeepPassedDetEtaCut");
  double minHeepPassedCrackCut = theConfig.pDouble("MinHeepPassedCrackCut");
  double minHeepPassedDEtaInCut = theConfig.pDouble("MinHeepPassedDEtaInCut");
  double minHeepPassedDPhiInCut = theConfig.pDouble("MinHeepPassedDPhiInCut");
  double minHeepPassedHademCut = theConfig.pDouble("MinHeepPassedHademCut");
  double minHeepPassedSigmaIEtaIEtaCut = theConfig.pDouble("MinHeepPassedSigmaIEtaIEtaCut");
  double minHeepPassed2by5Over5By5Cut = theConfig.pDouble("MinHeepPassed2by5Over5By5Cut");
  double minHeepPassedEcalHad1IsoCut = theConfig.pDouble("MinHeepPassedEcalHad1IsoCut");
  double minHeepPassedTrkIsoCut = theConfig.pDouble("MinHeepPassedTrkIsoCut");
  double minHeepPassedEcalDrivenCut = theConfig.pDouble("MinHeepPassedEcalDrivenCut");
  double minHeepPassedMissingHitsCut = theConfig.pDouble("MinHeepPassedMissingHitsCut");
  double minHeepPassedElecDxyCut = theConfig.pDouble("MinHeepPassedElecDxyCut");
  double minElecEtCut = theConfig.pDouble("MinElecEtCut");     
  double eMotherIdCut = theConfig.pDouble("ElectronMotherIdCut");
  
  double minElecRelIsoCut = theConfig.pDouble("MinElecRelIsoCut");
  double maxElecRelIsoCut = theConfig.pDouble("MaxElecRelIsoCut");
  double minElecTrkIsoCut = theConfig.pDouble("MinElecTrkIsoCut");
  double maxElecTrkIsoCut = theConfig.pDouble("MaxElecTrkIsoCut");  
  int maxElecMissingHitsCut = theConfig.pInt("MaxElecMissingHits");
  //double elecIsoDr = theConfig.pDouble("ElecIsoDr");
  double maxElecPFRelIsoCut = theConfig.pDouble("MaxElecPFRelIso");
  double maxElecDeltaEtaInEBCut = theConfig.pDouble("MaxElecDeltaEtaInBarrel");
  double maxElecDeltaPhiInEBCut = theConfig.pDouble("MaxElecDeltaPhiInBarrel");
  double maxElecSigmaIEtaIEtaEBCut = theConfig.pDouble("MaxElecSigmaIEtaIEtaBarrel");
  double maxElecHadFracEBCut = theConfig.pDouble("MaxElecHadFracBarrel");
  double maxOneOverEMinusOneOverPEBCut = theConfig.pDouble("MaxEOverPBarrel");
  double maxElecDeltaEtaInBBCut = theConfig.pDouble("MaxElecDeltaEtaInEndCap");
  double maxElecDeltaPhiInBBCut = theConfig.pDouble("MaxElecDeltaPhiInEndCap");
  double maxElecSigmaIEtaIEtaBBCut = theConfig.pDouble("MaxElecSigmaIEtaIEtaEndCap");
  double maxElecHadFracBBCut = theConfig.pDouble("MaxElecHadFracEndCap");
  double maxOneOverEMinusOneOverPBBCut = theConfig.pDouble("MaxEOverPEndCap");  
  double maxElecDxyBarrelCut = theConfig.pDouble("MaxElecDxyBarrel");
  double maxElecDxyEndCapCut = theConfig.pDouble("MaxElecDxyEndCap");
  
  double minElecIpOverElecIpErrorCut = theConfig.pDouble("MinElecIpOverElecIpErrorCut");
  double minMuonIpOverMuonIpErrorCut = theConfig.pDouble("MinMuonIpOverMuonIpErrorCut");
  double minCombinedIpOverIpErrorCut = theConfig.pDouble("MinCombinedIpOverIpErrorCut");
  double maxCombinedIpOverIpErrorCut = theConfig.pDouble("MaxCombinedIpOverIpErrorCut");
  double minElecMuonDCAOverError2DCut = theConfig.pDouble("MinElecMuonDCAOverError2DCut");
  double maxElecMuonDCAOverError2DCut = theConfig.pDouble("MaxElecMuonDCAOverError2DCut");
  double minElecMuonDCAOverError3DCut = theConfig.pDouble("MinElecMuonDCAOverError3DCut");
  double maxElecMuonDCAOverError3DCut = theConfig.pDouble("MaxElecMuonDCAOverError3DCut"); 
  double elecMuonMassCut = theConfig.pDouble("ElecMuonMassCut");  
  double minElecMuonCosDPhiCut = theConfig.pDouble("MinElecMuonCosDPhiCut");
  double maxElecMuonCosDPhiCut = theConfig.pDouble("MaxElecMuonCosDPhiCut");  
  double elecMuonDeltaRCut = theConfig.pDouble("ElecMuonDeltaRCut");
  double minElecMuonMetMassCut = theConfig.pDouble("MinElecMuonMetMassCut");
  double maxElecMuonMetMassCut = theConfig.pDouble("MaxElecMuonMetMassCut");
  double elecMuonDeltaPtMassCut = theConfig.pDouble("ElecMuonDeltaPtMassCut");
  double elecMetMtCut = theConfig.pDouble("ElecMetMtCut");
  double minMuonMetMtCut = theConfig.pDouble("MinMuonMetMtCut");
  double maxMuonMetMtCut = theConfig.pDouble("MaxMuonMetMtCut");  
  
  
  double elecMuonPZetaVisCut = theConfig.pDouble("ElecMuonPZetaVisCut");
  double minElecMuonPZetaCut = theConfig.pDouble("MinElecMuonPZetaCut");
  double maxElecMuonPZetaCut = theConfig.pDouble("MaxElecMuonPZetaCut");  
  double pZetaSlope = theConfig.pDouble("PZetaSlope");
  
  double pZetaCentroidCut = theConfig.pDouble("PZetaCentroid");
  double pZetaVisCentroidCut = theConfig.pDouble("PZetaVisCentroid");
  double pZetaSemiMajorAxisCut = theConfig.pDouble("PZetaSemiMajorAxis");
  double pZetaSemiMinorAxisCut = theConfig.pDouble("PZetaSemiMinorAxis");
  double pZetaEllipseAngleCut = theConfig.pDouble("PZetaEllipseAngle");
  
  double minNbTagsHiEffTrkCntCut = theConfig.pDouble("MinNbTagsHiEffTrkCntCut");
  double maxNbTagsHiEffTrkCntCut = theConfig.pDouble("MaxNbTagsHiEffTrkCntCut"); 
  double minNbTagsHiPurityTrkCntCut = theConfig.pDouble("MinNbTagsHiPurityTrkCntCut");
  double maxNbTagsHiPurityTrkCntCut = theConfig.pDouble("MaxNbTagsHiPurityTrkCntCut");  
  double minNbTagsHiEffSimpleSecVtxCut = theConfig.pDouble("MinNbTagsHiEffSimpleSecVtxCut");
  double maxNbTagsHiEffSimpleSecVtxCut = theConfig.pDouble("MaxNbTagsHiEffSimpleSecVtxCut");  
  double minNbTagsHiPuritySimpleSecVtxCut = theConfig.pDouble("MinNbTagsHiPuritySimpleSecVtxCut");
  double maxNbTagsHiPuritySimpleSecVtxCut = theConfig.pDouble("MaxNbTagsHiPuritySimpleSecVtxCut");  
  double minNbTagsCombSecVtxCut = theConfig.pDouble("MinNbTagsCombSecVtxCut");
  double maxNbTagsCombSecVtxCut = theConfig.pDouble("MaxNbTagsCombSecVtxCut");  
  double leptonMetCosDPhiCut = theConfig.pDouble("LeptonMetCosDPhiCut");
  double maxLeptonDeltaPtCosDPhiCut = theConfig.pDouble("MaxLeptonDeltaPtCosDPhiCut");
  double metCut = theConfig.pDouble("MetCut");
  double maxThreeLeptonCut = theConfig.pDouble("MaxThreeLeptonCut");
  double nJetsCut = theConfig.pDouble("NJetsCut");
  double minChargeCut = theConfig.pDouble("MinChargeCut");
  double maxChargeCut = theConfig.pDouble("MaxChargeCut");
  double jetSumEtCut = theConfig.pDouble("JetSumEtCut");  
  
  double minElecMuonDeltaPtCut = theConfig.pDouble("MinElecMuonDeltaPtCut");
  double min2xPLargerLeptonCut = theConfig.pDouble("Min2xPLargerLeptonCut");
  
  double minElecEOverPCut = theConfig.pDouble("MinElecEOverPCut");
  double maxElecEOverPCut = theConfig.pDouble("MaxElecEOverPCut");
  
  double minLifetimeCut = theConfig.pDouble("MinLifetimeCut");
  
  
   
  /*
  
  string inRootFileName = "RootFiles/FrancescoMC_4_19_15/ZPrimeSSM_M-1500";

  string source = "ZPrimeSSM_M-1500";
  double signalXSection = 0.015;
  double theLumi = 19779;
  int nEvents = -1;

  std::cout << "inRootFileName = " << inRootFileName << "\n";
  string logName = "logFiles/%Source%_TEST.log";
  string rootOutName = "tmpRootDir/%Source%_TEST.log"; 
  
  
  double minMuonPtCut = 0.;
  double maxMuonEtaCut = 3.;
  double maxMuonDxyVtxCut = 10.;
  double maxMuonDzVtxCut = 10.;
  double minMuonIsGlobalMuonCut = 0.;
  double minMuonMatchedCut = 0.;
  double muonMotherIdCut = 24.;
  double maxMuonNormChiSqrdCut = 10.;
  double minMuonValidHitsCut = 0.;
  double minMuonMatchedStationsCut = 0.;
  double minMuonPixelHitsCut = 0.;
  double minMuonTrkLayersWithHitsCut = 0.;

  double minMuonIsoDR04DBCut = 0.;
  double maxMuonIsoDR04DBCut = 5.;
  double minMuonIsHighPtMuonCut = 0.; 
  double maxMuonDptOverPtCut = 60.;
    
  double minElecPtCut = 0.;
  double maxElecEtaCut = 60.;

  double minElecMatchedCut = 0.; 
  double minHeepPassedEtCut = 0.;

  double minHeepPassedPtCut = 0.;
  double minHeepPassedDetEtaCut = 0.;
  double minHeepPassedCrackCut = 0.;
  double minHeepPassedDEtaInCut = 0.;
  double minHeepPassedDPhiInCut = 0.;
  double minHeepPassedHademCut = 0.;
  double minHeepPassedSigmaIEtaIEtaCut = 0.;
  double minHeepPassed2by5Over5By5Cut = 0.;
  double minHeepPassedEcalDrivenCut = 0.;
  double minHeepPassedMissingHitsCut = 0.;
  double minHeepPassedElecDxyCut = 0.;
  double minHeepPassedTrkIsoCut = 0.;
  double minHeepPassedEcalHad1IsoCut = 0.;
  double minElecEtCut = 0.;
  double eMotherIdCut = 24.;
  
  double minElecRelIsoCut = 0.;
  double maxElecRelIsoCut = 7.;
  double minElecTrkIsoCut = 0.;
  double maxElecTrkIsoCut = 7.;
  int maxElecMissingHitsCut = 7.;
  double maxElecPFRelIsoCut = 7.;
  double maxElecDeltaEtaInEBCut = 7.;
  double maxElecDeltaPhiInEBCut = 7.;
  double maxElecSigmaIEtaIEtaEBCut =7.; 
  double maxElecHadFracEBCut = 8.;
  double maxOneOverEMinusOneOverPEBCut = 8.;
  double maxElecDeltaEtaInBBCut = 8.;
  double maxElecDeltaPhiInBBCut = 8.;
  double maxElecSigmaIEtaIEtaBBCut = 8.;
  double maxElecHadFracBBCut = 8.;
  double maxOneOverEMinusOneOverPBBCut = 8.;
  double maxElecDxyBarrelCut = 0.02;
  double maxElecDxyEndCapCut = 0.05;
  
  double minElecIpOverElecIpErrorCut = 0.;
  double minMuonIpOverMuonIpErrorCut = 0.;
  double minCombinedIpOverIpErrorCut = 0.;
  double maxCombinedIpOverIpErrorCut = 60.;
  double minElecMuonDCAOverError2DCut = 0.;
  double maxElecMuonDCAOverError2DCut = 1000000.;
  double minElecMuonDCAOverError3DCut = 0.;
  double maxElecMuonDCAOverError3DCut = 1000000.;  
  
  double elecMuonMassCut = 0.;
  double minElecMuonCosDPhiCut = 0.;
  double maxElecMuonCosDPhiCut = 60.;
  double elecMuonDeltaRCut = 0.;
  double minElecMuonMetMassCut = 0.;
  double maxElecMuonMetMassCut = 100000.;
  double elecMuonDeltaPtMassCut = 0.;
  double elecMetMtCut = 0.;
  double minMuonMetMtCut = 0.;
  double maxMuonMetMtCut = 6000.;
  
  
  double elecMuonPZetaVisCut = 0.;
  double minElecMuonPZetaCut = 0.;
  double maxElecMuonPZetaCut = 55.;
  double pZetaSlope = 5.;
  
  double pZetaCentroidCut = 0.;
  double pZetaVisCentroidCut = 0.;
  double pZetaSemiMajorAxisCut = 500.;
  double pZetaSemiMinorAxisCut = 500.;
  double pZetaEllipseAngleCut = 50.;
  
  double minNbTagsHiEffTrkCntCut = 0.;
  double maxNbTagsHiEffTrkCntCut =  50.;
  double minNbTagsHiPurityTrkCntCut = 0.;
  double maxNbTagsHiPurityTrkCntCut =  50.;
  double minNbTagsHiEffSimpleSecVtxCut = 0.;
  double maxNbTagsHiEffSimpleSecVtxCut = 50.;
  double minNbTagsHiPuritySimpleSecVtxCut = 0.;
  double maxNbTagsHiPuritySimpleSecVtxCut = 50.;
  double minNbTagsCombSecVtxCut = 0.;
  double maxNbTagsCombSecVtxCut = 50.; 
  double leptonMetCosDPhiCut = 5.;
  double maxLeptonDeltaPtCosDPhiCut = 5.;
  double metCut = 0.;
  double maxThreeLeptonCut = 0.;
  double nJetsCut = 0.;
  double minChargeCut = -2.;
  double maxChargeCut = 2.;
  double jetSumEtCut = 0.;
  
  double minElecMuonDeltaPtCut = 0.; 
  double min2xPLargerLeptonCut = 0.; 
  double minElecEOverPCut = 0.;
  double maxElecEOverPCut = 1000.; 
  double minLifetimeCut = 0.;
  
  
  */ 
   
  
   
   
   
   
   
    
  // instatiate the class
  TChain* theTree;
  
  //check if input file is directory or file; if dir use a chain
  struct stat sb;
  if (stat(inRootFileName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)){
    theTree = new TChain("theTree", "ChainedTrees");
    string theFiles = inRootFileName+"/*.root/analyzeHiMassTau/HMTTree";
    theTree->Add(theFiles.c_str());
  }else{
    TFile* theFile = new TFile(inRootFileName.c_str());
    theTree = (TChain*)theFile->Get("analyzeHiMassTau/HMTTree");
  }

 

  EMuAnalysis EMuAnalysis(theTree);
  
  EMuAnalysis.SetSource(source);
  EMuAnalysis.SetSignalXSection(signalXSection);
  EMuAnalysis.SetLumi(theLumi);
  EMuAnalysis.SetNEvents(nEvents);
  EMuAnalysis.SetOutputLogFileName(logName);
  EMuAnalysis.SetOutputRootFileName(rootOutName);
  
  EMuAnalysis.SetElecAccCuts(minElecPtCut, maxElecEtaCut);
  
  EMuAnalysis.SetElecIdGlobalCuts(maxElecMissingHitsCut, 
                                         maxElecPFRelIsoCut);
  EMuAnalysis.SetElecIdBarrelCuts(maxElecDeltaPhiInEBCut, maxElecDeltaEtaInEBCut, maxElecSigmaIEtaIEtaEBCut, 
                                         maxElecHadFracEBCut, maxOneOverEMinusOneOverPEBCut, maxElecDxyBarrelCut);
  EMuAnalysis.SetElecIdEndCapCuts(maxElecDeltaPhiInBBCut, maxElecDeltaEtaInBBCut, maxElecSigmaIEtaIEtaBBCut,
                                         maxElecHadFracBBCut, maxOneOverEMinusOneOverPBBCut, maxElecDxyEndCapCut);  
 
  EMuAnalysis.SetElecMotherIdCut(eMotherIdCut); 
  
  EMuAnalysis.SetMuonAccCuts(minMuonPtCut, maxMuonEtaCut);
  
  EMuAnalysis.SetMuonIdCuts(maxMuonDxyVtxCut, maxMuonDzVtxCut,
			      minMuonIsGlobalMuonCut, minMuonMatchedCut, maxMuonNormChiSqrdCut, minMuonValidHitsCut,
			      minMuonMatchedStationsCut, minMuonPixelHitsCut, minMuonTrkLayersWithHitsCut, minMuonIsoDR04DBCut,
			      maxMuonIsoDR04DBCut, minMuonIsHighPtMuonCut, maxMuonDptOverPtCut);
			  
  EMuAnalysis.SetMuonMotherIdCut(muonMotherIdCut);			      	 
			      
  EMuAnalysis.SetHeepCuts(minElecMatchedCut, minHeepPassedEtCut, minHeepPassedDetEtaCut, minHeepPassedCrackCut, minHeepPassedDEtaInCut, 
			      minHeepPassedDPhiInCut, minHeepPassedHademCut, minHeepPassedSigmaIEtaIEtaCut, minHeepPassed2by5Over5By5Cut, 
			      minHeepPassedEcalDrivenCut, minHeepPassedMissingHitsCut, minHeepPassedElecDxyCut, minElecEtCut, minElecRelIsoCut, 
			      maxElecRelIsoCut, minElecTrkIsoCut, maxElecTrkIsoCut, minHeepPassedEcalHad1IsoCut, minHeepPassedTrkIsoCut);
			      
  EMuAnalysis.SetImpactParameterCuts(minElecIpOverElecIpErrorCut, minMuonIpOverMuonIpErrorCut, minCombinedIpOverIpErrorCut, maxCombinedIpOverIpErrorCut, minElecMuonDCAOverError2DCut, maxElecMuonDCAOverError2DCut, 
  			      minElecMuonDCAOverError3DCut, maxElecMuonDCAOverError3DCut, minLifetimeCut);			      
  			      
  EMuAnalysis.SetTopologyCuts(elecMuonMassCut, minElecMuonCosDPhiCut, maxElecMuonCosDPhiCut, elecMuonDeltaRCut, minElecMuonMetMassCut,
  		        maxElecMuonMetMassCut, elecMuonDeltaPtMassCut, elecMetMtCut, minMuonMetMtCut, maxMuonMetMtCut, elecMuonPZetaVisCut,
		        minElecMuonPZetaCut, maxElecMuonPZetaCut, pZetaSlope, minNbTagsHiEffTrkCntCut, maxNbTagsHiEffTrkCntCut, 
			minNbTagsHiPurityTrkCntCut, maxNbTagsHiPurityTrkCntCut, minNbTagsHiEffSimpleSecVtxCut, maxNbTagsHiEffSimpleSecVtxCut, 
			minNbTagsHiPuritySimpleSecVtxCut, maxNbTagsHiPuritySimpleSecVtxCut, minNbTagsCombSecVtxCut, maxNbTagsCombSecVtxCut,
		        leptonMetCosDPhiCut, maxLeptonDeltaPtCosDPhiCut, metCut, nJetsCut, minChargeCut, maxChargeCut, jetSumEtCut, minElecMuonDeltaPtCut, min2xPLargerLeptonCut, minElecEOverPCut, maxElecEOverPCut, maxThreeLeptonCut);
  			
  
 /*
  EMuAnalysis.SetTopologyCuts(elecMuonMassCut, minElecMuonCosDPhiCut, maxElecMuonCosDPhiCut, elecMuonDeltaRCut, 
  		        elecMuonMetMassCut, elecMetMtCut, minMuonMetMtCut, maxMuonMetMtCut, pZetaCentroidCut,
		        pZetaVisCentroidCut, pZetaSemiMajorAxisCut, pZetaSemiMinorAxisCut, pZetaEllipseAngleCut, minNbTagsHiEffTrkCntCut, maxNbTagsHiEffTrkCntCut, 
			minNbTagsHiPurityTrkCntCut, maxNbTagsHiPurityTrkCntCut, minNbTagsHiEffSimpleSecVtxCut, maxNbTagsHiEffSimpleSecVtxCut, 
			minNbTagsHiPuritySimpleSecVtxCut, maxNbTagsHiPuritySimpleSecVtxCut, minNbTagsCombSecVtxCut, maxNbTagsCombSecVtxCut,
		        leptonMetCosDPhiCut, metCut, nJetsCut, minChargeCut, maxChargeCut, jetSumEtCut);			
					      			         

  */
  EMuAnalysis.Loop();
}
