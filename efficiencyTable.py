#!/usr/local/python/bin/python2.7
import os, math, sys, argparse



def main(argv):
  parser = argparse.ArgumentParser(description='Get sig/bg ratio')
  parser.add_argument('--variableName', help='Variable in logfile', required = True)
  parser.add_argument('--initValue', help='Starting value', required = True)
  parser.add_argument('--endValue', help='End value', required = True)
  parser.add_argument('--interval', help='Iteration interval', required = True)
  args = parser.parse_args()
 
  sources = ['ZPrime', 'WJets', 'TTJets']
  
  workDir = os.environ['PWD']
  theInitValue = float(args.initValue)
  theEndValue = float(args.endValue)
  theInterval = float(args.interval)
  value = float(args.initValue)
  
  tableOutFile = open('EfficiencyTable.txt', 'r+')

  
  tableOutFile.truncate()

  theRange = (theEndValue - theInitValue)/theInterval
  print "The range = " + str(theRange)
  for x in range(0, int(theRange) + 1):
    
    nEventsDic = {}
    theEff1Dic = {}
    theEff2Dic = {}
    
    if str(value).split('.')[1] == '0':
      if value == theInitValue:
        value = int(theInitValue)
      if value != theInitValue:
        value = int(value)
    
    for source in sources:
      theFileName = workDir + '/logFiles/'+source+'_' + args.variableName + '_'+str(value)+'.log'
      #print theFileName
      theFile = open(theFileName, 'r') 
      for line in theFile:
	if "Muon Pt" in line:
	  
	  theEff1 = float(line.rsplit(':')[1].split('+-')[0].strip())
	  theEff1Error = float(line.rsplit(':')[1].split('+-')[1].strip())
	  theEff1Dic[source] = (theEff1, theEff1Error)  
   	
	if "Elec Pt" in line:
	  theEff2 = float(line.rsplit(':')[1].split('+-')[0].strip())
	  theEff2Error = float(line.rsplit(':')[1].split('+-')[1].strip())
	  theEff2Dic[source] = (theEff2, theEff2Error)  	  
	
	if "@19600" in line:
   	  nEvents = float(line.rsplit(':')[1].split('+-')[0].strip())
   	  nEventsError = float(line.rsplit(':')[1].split('+-')[1].strip())
	  nEventsDic[source] = (nEvents, nEventsError)	  
	  
    sig = nEventsDic['ZPrime'][0]
    ttjets = nEventsDic['TTJets'][0]
    wjets = nEventsDic['WJets'][0]
    
    sigErr = nEventsDic['ZPrime'][1]
    ttjetsErr = nEventsDic['TTJets'][1]
    wjetsErr = nEventsDic['WJets'][1]
    
    sigEff1 = theEff1Dic['ZPrime'][0]
    sigEff1Err = theEff1Dic['ZPrime'][1]

    sigEff2 = theEff2Dic['ZPrime'][0]
    sigEff2Err = theEff2Dic['ZPrime'][1]
   
    ttjetsEff1 = theEff1Dic['TTJets'][0]
    ttjetsEff1Err = theEff1Dic['TTJets'][1]
   
    ttjetsEff2 = theEff2Dic['TTJets'][0]
    ttjetsEff2Err = theEff2Dic['TTJets'][1]    

    wjetsEff1 = theEff1Dic['WJets'][0]
    wjetsEff1Err = theEff1Dic['WJets'][1]
    
    wjetsEff2 = theEff2Dic['WJets'][0]
    wjetsEff2Err = theEff2Dic['WJets'][1]    
    
    
    
    print "NEvents for dual Pt cut " + str(value) + ": Signal: " + str(sig) + " +- " + str(sigErr)  + ", TTJets: " + str(ttjets) + " +- " + str(ttjetsErr) + ", W+Jets: " + str(wjets) + " +- " + str(wjetsErr) + "\n"
    
    tableOutFile.write("NEvents for dual Pt cut " + str(value) + ": Signal: " + str(sig) + " +- " + str(sigErr)  + ", TTJets: " + str(ttjets) + " +- " + str(ttjetsErr) + ", W+Jets: " + str(wjets) + " +- " + str(wjetsErr) + "\n")

    
    value = float(value) + theInterval
    
  
   
  tableOutFile.close()

      
if __name__ == "__main__":
  main(sys.argv[1:])
