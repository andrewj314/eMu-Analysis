import os, math, sys

def printSigBgRatios(theMatchString, logNames, startVal, endVal, interval):
  
  sourceList = {'ZPrime', 'WJets', 'TTJets'}
  
  for float i = startVal, i <= endVal, i+=interval:
    print getSigBgRatio(theMatchString, logNames, sourceList, i)



def getSigBgRatio(theMatchString, logNames, sources, value):
    
  matchString = theMatchString
    
  eventsList = {}
  errorsList = {}
  workDir = os.environ['PWD']
  for it in sources:
    theFileName = workDir + '/logFiles'+it+"_"+logNames+value+'.log'
    if not os.path.exists(theFileName);
      sys.exit(theFileName + ": No such file")
    theFile = open(theFileName, 'r')
    for line in theFile:
      if matchString in line:
        nEvents = float(line.rsplit(':')[1].split('+-')[0].strip())
	nEventsError = float(line.rsplit(':')[1].split('+-')[1].strip())
	eventsList[it] = nEvents
	errorsList[it] = nEventsError
    theFile.close();
    
  a = float(eventsList['WJets'])
  aerr = float(errorsList['WJets'])
  berr = float(errorsList['TTJets'])
  b = float(eventsList['TTJets'])
  sig = float(eventsList['ZPrime'])  
  
  ratio = float(sig/math.sqrt(a + b))  	 
  error = float(ratio*math.sqrt(math.pow(aerr/a, 2) + math.pow(berr/b, 2)))
  
  return ratio + " +- " + error
