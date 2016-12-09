#!/usr/local/python/bin/python2.7
import os, math, sys, argparse



def main(argv):
  parser = argparse.ArgumentParser(description='Get sig/bg ratio')
  parser.add_argument('--logName', help='Variable in logfile title', required = True)
  parser.add_argument('--initValue', help='Starting value', required = True)
  parser.add_argument('--endValue', help='End value', required = True)
  parser.add_argument('--interval', help='Iteration interval', required = True)
  args = parser.parse_args()
 
  
  sources = ['Data']
  
  workDir = os.environ['PWD']
  theInitValue = float(args.initValue)
  theEndValue = float(args.endValue)
  theInterval = float(args.interval)
  value = float(args.initValue)
    
  
  dataEffOutFile = open('DataEff.txt', 'r+') 
  
  
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
    
    theFileName = workDir + '/logFiles/Data_eePairs_LifetimeStudy_M_gt_400/Data_' + args.logName + '_'+theValue+'.log'
      
    theFile = open(theFileName, 'r') 
    for line in theFile:
      if " Lifetime" in line:
	  
	theEff = float(line.rsplit(':')[1].split('+-')[0].strip())
	theEffError = float(line.rsplit(':')[1].split('+-')[1].strip())
	theEffDic['Data'] = (theEff, theEffError)  
      if "@19" in line:
        nEvents = float(line.rsplit(':')[1].split('+-')[0].strip())
   	nEventsError = float(line.rsplit(':')[1].split('+-')[1].strip())
	nEventsDic['Data'] = (nEvents, nEventsError)	  
	  
	print "Events for data at value " + theValue + ": " + str(nEvents)
	  
    data = nEventsDic['Data'][0]
    
    
    
    dataEff = theEffDic['Data'][0]
    dataEffErr = theEffDic['Data'][1]
    
    
    
      
    
    
    
    
    
    
    dataEffOutFile.write(str(value) + " " + str(dataEff) + " " + str(dataEffErr) + "\n")
   
    
    value = float(value) + theInterval
    
  
   
  dataEffOutFile.close()
      
if __name__ == "__main__":
  main(sys.argv[1:])
