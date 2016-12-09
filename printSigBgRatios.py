#!/usr/local/python/bin/python2.7
import os, math, sys, argparse



def main(argv):
  parser = argparse.ArgumentParser(description='Get sig/bg ratio')
  parser.add_argument('--logName', help='Variable in logfile title', required = True)
  parser.add_argument('--initValue', help='Starting value', required = True)
  parser.add_argument('--endValue', help='End value', required = True)
  parser.add_argument('--interval', help='Iteration interval', required = True)
  args = parser.parse_args()
 
  sources = ['ZPrimeSSM_M-1500', 'TTJets_MADGRAPH','WW', 'WJets', 'WZ', 'ZZ', 'QCD', 'DYJetsToLeptonLepton']
  
  ##sources = ['ZPrimeSSM_M-1500']
  
  ##sources = ['Data']
  
  workDir = os.environ['PWD']
  theInitValue = float(args.initValue)
  theEndValue = float(args.endValue)
  theInterval = float(args.interval)
  value = float(args.initValue)
    
  ratioOutFile = open('SigBgRatio.txt', 'r+')
  sigEffOutFile = open('SigEff.txt', 'r+')
  ttjetsEffOutFile = open('TTJetsEff.txt', 'r+')
  wjetsEffOutFile = open('WJetsEff.txt', 'r+')
  qcdEffOutFile = open('QCDEff.txt', 'r+')
  ztautauEffOutFile = open('DYJetsToLeptonLeptonEff.txt', 'r+')
  zeeEffOutFile = open('DYToEEEff.txt', 'r+')
  wwEffOutFile = open('WWEff.txt', 'r+')
  wzEffOutFile = open('WZEff.txt', 'r+')
  zzEffOutFile = open('ZZEff.txt', 'r+')
  dataEffOutFile = open('DataEff.txt', 'r+') 
  
  
  ratioOutFile.truncate()
  sigEffOutFile.truncate()
  ttjetsEffOutFile.truncate()
  wjetsEffOutFile.truncate()
  qcdEffOutFile.truncate()
  ztautauEffOutFile.truncate()
  zeeEffOutFile.truncate()
  wwEffOutFile.truncate()
  wzEffOutFile.truncate()
  zzEffOutFile.truncate()
  dataEffOutFile.truncate()
  
  theRange = (theEndValue - theInitValue)/theInterval
  print "The range = " + str(theRange)
  for x in range(0, int(theRange) + 1):
    
    theValue = str(value)
    if theValue.split('.')[1] == '0':
      theValue = theValue.split('.')[0]
    
    nEventsDic = {}
    theEffDic = {}
    
    #if str(value).split('.')[1] == '0':
      #if value == theInitValue:
        #value = int(theInitValue)
      #if value != theInitValue:
        #value = int(value) 
    
    for source in sources:
      #theFileName = workDir + '/ImpactParameter200GeVMassCutStandardCosDphiCut11_20_13/'+source+'_' + args.logName + '_'+theValue+'.log'
      theFileName = workDir + '/logFiles/'+source+'_' + args.logName + '_'+theValue+'.log'
      ##theFileName = workDir + '/logFiles/Data_eePairs_LifetimeStudy_noMassCut/'+source+'_' + args.logName + '_'+theValue+'.log'
      
      #print theFileNamenoPZetaNoCos
      theFile = open(theFileName, 'r') 
      for line in theFile:
	if " Lifetime" in line:
	  
	  theEff = float(line.rsplit(':')[1].split('+-')[0].strip())
	  theEffError = float(line.rsplit(':')[1].split('+-')[1].strip())
	  theEffDic[source] = (theEff, theEffError)  
   	if "@19" in line:
   	  nEvents = float(line.rsplit(':')[1].split('+-')[0].strip())
   	  nEventsError = float(line.rsplit(':')[1].split('+-')[1].strip())
	  nEventsDic[source] = (nEvents, nEventsError)	  
	  
	  print "Events for source " + str(source) + " at value " + theValue + ": " + str(nEvents)
	  
    sig = nEventsDic['ZPrimeSSM_M-1500'][0]
    ttjets = nEventsDic['TTJets_MADGRAPH'][0]
    wjets = nEventsDic['WJets'][0]
    qcd = nEventsDic['QCD'][0]
    ztautau = nEventsDic['DYJetsToLeptonLepton'][0]
    ##zee = nEventsDic['DYToEE'][0]
    ww = nEventsDic['WW'][0]
    wz = nEventsDic['WZ'][0]
    zz = nEventsDic['ZZ'][0]
    ##data = nEventsDic['Data'][0]
    
    sigErr = nEventsDic['ZPrimeSSM_M-1500'][1]
    ttjetsErr = nEventsDic['TTJets_MADGRAPH'][1]
    wjetsErr = nEventsDic['WJets'][1]
    qcdErr = nEventsDic['QCD'][1]
    ztautauErr = nEventsDic['DYJetsToLeptonLepton'][1]
    ##zeeErr = nEventsDic['DYToEE'][1]
    wwErr = nEventsDic['WW'][1]
    wzErr = nEventsDic['WZ'][1]
    zzErr = nEventsDic['ZZ'][1]
    ##dataErr = nEventsDic['Data'][1]    
    
    sigEff = theEffDic['ZPrimeSSM_M-1500'][0]
    sigEffErr = theEffDic['ZPrimeSSM_M-1500'][1]
    
    ttjetsEff = theEffDic['TTJets_MADGRAPH'][0]
    ttjetsEffErr = theEffDic['TTJets_MADGRAPH'][1]
    
    wjetsEff = theEffDic['WJets'][0]
    wjetsEffErr = theEffDic['WJets'][1]
    
    qcdEff = theEffDic['QCD'][0]
    qcdEffErr = theEffDic['QCD'][1]
    
    ztautauEff = theEffDic['DYJetsToLeptonLepton'][0]
    ztautauEffErr = theEffDic['DYJetsToLeptonLepton'][1]
    
    ##zeeEff = theEffDic['DYToEE'][0]
    ##zeeEffErr = theEffDic['DYToEE'][1]
    
    wwEff = theEffDic['WW'][0]
    wwEffErr = theEffDic['WW'][1]
    
    wzEff = theEffDic['WZ'][0]
    wzEffErr = theEffDic['WZ'][1]
    
    zzEff = theEffDic['ZZ'][0]
    zzEffErr = theEffDic['ZZ'][1]
    
    ##dataEff = theEffDic['Data'][0]
    ##dataEffErr = theEffDic['Data'][1]
    
    ##sumBg = ttjets + wjets + qcd + ztautau + ww + wz + zz
    sumBg = ttjets + ww
    sumBgErr = ttjetsErr + wjetsErr + qcdErr + ztautauErr + wwErr + wzErr + zzErr
    
    ratio = 0
    ttbarDataMCRatio = 0
    
    if math.sqrt(sig + ttjets + ww + wz + zz + wjets + qcd + ztautau) != 0:
      ##ratio = sig/math.sqrt(sig + ttjets + ww + wz + zz + wjets + qcd + ztautau)
      ratio = sig/math.sqrt(sig + ttjets + ww)
      
    #SF = (data - ztautau)/ttjets
    
    #ttjets = ttjets*SF  
    
    #if ttjets != 0:
       #ratio = dataEff/ttjetsEff  
    
    ratioErr = 0
    
    sigRatioErr = 0
    if sig != 0:
      sigRatioErr = math.pow(sigErr/sig, 2)
      
    ttjetsRatioErr = 0
    if ttjets != 0:
      ttjetsRatioErr = math.pow(ttjetsErr/ttjets, 2)      
 
    wjetsRatioErr = 0
    if wjets != 0:
      wjetsRatioErr = math.pow(wjetsErr/wjets, 2)   
      
    qcdRatioErr = 0
    if qcd != 0:
      qcdRatioErr = math.pow(qcdErr/qcd, 2)
      
    ztautauRatioErr = 0
    if ztautau != 0:
      ztautauRatioErr = math.pow(ztautauErr/ztautau, 2)                                 
    
    ##zeeRatioErr = 0
    ##if zee != 0:
      ##zeeRatioErr = math.pow(zeeErr/zee, 2)      
    
    wwRatioErr = 0
    if ww != 0:
      wwRatioErr = math.pow(wwErr/ww, 2)
   
    wzRatioErr = 0
    if wz != 0:
      wzRatioErr = math.pow(wzErr/wz, 2)   
      
    zzRatioErr = 0
    if zz != 0:
      zzRatioErr = math.pow(zzErr/zz, 2)          
            
    #ratioErr = ratio*math.sqrt(sigRatioErr + ttjetsRatioErr + wjetsRatioErr + qcdRatioErr + ztautauRatioErr + wwRatioErr + wzRatioErr + zzRatioErr)
    
    #ratioErr = ratio*math.sqrt(sigRatioErr + ttjetsRatioErr + wwRatioErr)
    
    #if sig == 0:
      #sigErr = 0
    #if ttjets == 0:
      #ttjetsErr = 0  
    #if wjets == 0:
      #wjetsErr = 0
    #if qcd == 0:
      #qcdErr = 0
    #if ztautau == 0:
      #ztautauErr = 0
    #if ww == 0:
      #wwErr = 0
    #if wz == 0:
      #wzErr = 0
    #if zz == 0:
      #zzErr = 0          
    
    #ratioErr = math.sqrt(math.pow(sigErr, 2)/sumBg + math.pow(sig*ttjetsErr,2)/(math.pow(sumBg, 3)) + math.pow(sig*wjetsErr,2)/(math.pow(sumBg, 3)) + 
    #math.pow(sig*qcdErr,2)/(math.pow(sumBg, 3)) + math.pow(sig*ztautauErr,2)/(math.pow(sumBg, 3)) + math.pow(sig*wwErr,2)/(math.pow(sumBg, 3)) + 
    #math.pow(sig*wzErr,2)/(math.pow(sumBg, 3)) + math.pow(sig*zzErr,2)/(math.pow(sumBg, 3)))
    
    ratioErr = 0
    if sumBg != 0:
      ratioErr = math.sqrt(math.pow(sigErr, 2)/sumBg + math.pow(sig*ttjetsErr,2)/(math.pow(sumBg, 3)) + math.pow(sig*wwErr,2)/(math.pow(sumBg, 3)))
    
    
    #if ttjets != 0:
       #ratioErr = math.sqrt(math.pow(dataEffErr/ttjetsEff, 2) + math.pow(dataEff*ttjetsEffErr/math.pow(ttjetsEff, 2),2))
    
    print "ratio for value " + theValue + " = " + str(ratio) + " +- " + str(ratioErr)
    print "signal efficiency for value " + theValue + " = " + str(sigEff) + " +- " + str(sigEffErr)
    
    ratioOutFile.write(theValue + " " + str(ratio) + " " + str(ratioErr) + "\n")
    sigEffOutFile.write(theValue + " " + str(sigEff) + " " + str(sigEffErr) + "\n")
    ttjetsEffOutFile.write(theValue + " " + str(ttjetsEff) + " " + str(ttjetsEffErr) + "\n")
    wjetsEffOutFile.write(str(value) + " " + str(wjetsEff) + " " + str(wjetsEffErr) + "\n")
    qcdEffOutFile.write(str(value) + " " + str(qcdEff) + " " + str(qcdEffErr) + "\n")
    ztautauEffOutFile.write(str(value) + " " + str(ztautauEff) + " " + str(ztautauEffErr) + "\n")
    ##zeeEffOutFile.write(str(value) + " " + str(zeeEff) + " " + str(zeeEffErr) + "\n")
    wwEffOutFile.write(theValue + " " + str(wwEff) + " " + str(wwEffErr) + "\n")
    wzEffOutFile.write(str(value) + " " + str(wzEff) + " " + str(wzEffErr) + "\n")
    zzEffOutFile.write(str(value) + " " + str(zzEff) + " " + str(zzEffErr) + "\n")
    ##dataEffOutFile.write(str(value) + " " + str(dataEff) + " " + str(dataEffErr) + "\n")
   
    
    value = float(value) + theInterval
    
  
   
  ratioOutFile.close()
  sigEffOutFile.close()  
  ttjetsEffOutFile.close()
  wjetsEffOutFile.close()
  qcdEffOutFile.close()
  ztautauEffOutFile.close()
  ##zeeEffOutFile.close()
  wwEffOutFile.close()
  wzEffOutFile.close()
  zzEffOutFile.close()
  ##dataEffOutFile.close()
      
if __name__ == "__main__":
  main(sys.argv[1:])
