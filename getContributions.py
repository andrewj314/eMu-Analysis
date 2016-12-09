import os, math, sys

def getContributions(theMatchExp, theType, logNames, labels):
  matchExp = theMatchExp
  eventsDic = {}
  errorsDic = {}
  workDir = os.environ['PWD']
  for it in labels:
    #if it == 'data':
    #  matchExp = 'Tight'
    theFileName = workDir + '/logFiles/'+it+"_"+logNames+'.log'
    #theFileName = workDir + '/logFiles/'+it+logNames+'.log'
    if not os.path.exists(theFileName):
      sys.exit(theFileName + ": No such files")
    theFile = open(theFileName, 'r')
    for line in theFile:
      #print line
      if matchExp in line:
        print str(it) + ' Contribution = ' + line.rsplit(':')[1].split('+-')[0]
  	nEvents = float(line.rsplit(':')[1].split('+-')[0].strip())  #strip removes leading and trailing characters
  	nEventsError = math.sqrt(nEvents)
  	if(it != 'Data'):
  	  nEventsError = float(line.rsplit(':')[1].split('+-')[1].strip())
  	eventsDic[it] = nEvents
  	errorsDic[it] = nEventsError
    theFile.close()
  print nEvents
  print nEventsError
  if theType == "Values":
    return eventsDic
  if theType == "Errors":
    return errorsDic
    
    
def correctBtag(logNames,labels):
  import ROOT
  from getEff import getEfficiency, getSurvivingEvents
  ### get data histo
  theDataFile = ROOT.TFile('tmpRootDir/Data_' + logNames + '.root')
  theDataHisto = theDataFile.Get('DataHistosDirectory/_nBTagCR2Histo')
  passedBeforeCorrection = theDataHisto.GetBinContent(1)
  print 'passedBeforeCorrection = ' + str(passedBeforeCorrection)
  theDataHisto.Sumw2()
  ### get contributions
  eventContributions = getContributions("TTJets in CR2", "Values", logNames, labels)
  ### get btag control region plots for each source
  theContaminationHisto = ROOT.TH1F("theContaminationHisto", "NBtags BG", 10, 0, 10)
  for key, value in eventContributions.iteritems():
    ## open file; scale according to contribution; put them into the contamination histo
    print 'key = ' + str(key)
    if key != 'TTJets' and key != 'Data':
      theFile = ROOT.TFile('tmpRootDir/' + key + '_' + logNames + '.root')
      theHisto = theFile.Get(key + 'HistosDirectory/_nBTagCR2Histo')
      theHisto.Sumw2()
      if theHisto.Integral() > 0:
        theHisto.Scale(eventContributions[key]/theHisto.Integral())
	#theHisto.Scale(3348/theHisto.Integral())
        theContaminationHisto.Add(theHisto)
  ### subtract contamination from data (no ttbar mc)
  theDataHisto.Add(theContaminationHisto, -1.)
  
  passed = theDataHisto.GetBinContent(1)
  print 'passed = ' + str(passed)
  total = theDataHisto.Integral()
  print 'total = ' + str(total)
  theEff = getEfficiency(passed, total)[0]
  theEffError = getEfficiency(passed, total)[1] 
  oneMinusEff = 1. - theEff
  effRatio = theEff/oneMinusEff
  effRatioError = ((1./oneMinusEff + (theEff)/(oneMinusEff*oneMinusEff))*theEffError);
  
  print 'Eff\t' + str(theEff) + ' +- ' + str(theEffError)
  print 'EffRatio ' + str(effRatio) +  ' +- ' + str(effRatioError)
  
  ### recalculate bg contribution
  ### open data logfile rw; get all efficiencies from data log file
  theDataLog = 'logFiles/Data' + '_' + logNames + '.log'
  theLogFile = open(theDataLog, 'r')
  effs = {}
  effErrors = {}
  dataBTagControl = 0
  dataBTagControlError = 0
  theMatchExps = ['JetSumCR1', 'TopologyCR1', 'BTag/AntiBTag']
  for line in theLogFile:
    if 'TTJets in CR1' in line:
      dataBTagControl = float(line.rsplit(':')[1].split('+-')[0].strip())
      dataBTagControlError = math.sqrt(dataBTagControl)
    for index, value in enumerate(theMatchExps):
      if theMatchExps[index] in line:
        print theMatchExps[index]
        effs[theMatchExps[index]] = float(line.rsplit(':')[1].split('+-')[0].strip())
	effs[theMatchExps[index]] = effs[theMatchExps[index]]*0.01
	effErrors[theMatchExps[index]] = float(line.rsplit(':')[1].split('+-')[1].strip())
	effErrors[theMatchExps[index]] = effErrors[theMatchExps[index]]*0.01
  theLogFile.close()
  #eventsAfterMET = getSurvivingEvents(dataBTagControl, dataBTagControlError, effs['All MET'], effErrors['All MET'])
  eventsAfterJetSum = getSurvivingEvents(dataBTagControl, dataBTagControlError, effs['JetSumCR1'], effErrors['JetSumCR1'])
  eventsAfterTopo = getSurvivingEvents(eventsAfterJetSum[0], eventsAfterJetSum[1], effs['TopologyCR1'], effErrors['TopologyCR1'])
  eventsAfterBTag = getSurvivingEvents(eventsAfterTopo[0], eventsAfterTopo[1], effRatio, effRatioError)
  print 'eventsAfterJetSum = ' + str(eventsAfterJetSum[0]) + ' +- ' + str(eventsAfterJetSum[1]) + '\n'
  print 'eventsAfterTopo = ' + str(eventsAfterTopo[0]) + ' +- ' + str(eventsAfterTopo[1]) + '\n'
  print 'eventsAfterBTag = ' + str(eventsAfterBTag[0]) + ' +- ' + str(eventsAfterBTag[1]) + '\n'
  
  
  ### re-open file and replace new eff and final ttbar contribution
  theLogFile = open(theDataLog, 'r+')
  theOutFile = open("test.log", "w")
  for line in theLogFile:
    if 'BTagVeto' in line:
      newLine = ' BTagVeto: \t \t' + '{0:.4f}'.format(theEff) + '\t +- \t' + '{0:.4f}'.format(theEffError) +'\n'
      line = newLine
    if 'BTag/AntiBTag' in line:
      newLine = ' BTag/AntiBTag: \t' + '{0:.4f}'.format(effRatio) + '\t  +- \t ' + '{0:.4f}'.format(effRatioError) +'\n'
      line = newLine
    if 'Final TTbar' in line:
      newLine = ' Final TTbar: \t \t' + '{0:.4f}'.format(eventsAfterBTag[0]) + '\t  +- \t ' + '{0:.4f}'.format(eventsAfterBTag[1]) +'\n'
      line = newLine
    theOutFile.write(line)
  
  theLogFile.close()
  theOutFile.close()
  os.rename("test.log", theDataLog)
