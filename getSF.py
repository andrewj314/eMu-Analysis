#!/usr/bin/python
import sys, getopt, math, os
sys.path.append('/home/ajohnson/CMS/elecTauTauAnalysis/scripts/python')
from getContributions import getContributions, correctBtag

#if len(sys.argv) != 3:
#  sys.exit("Forgot the arguments: [signal] [logName]")

def main(argv):
  try:
    opt, args = getopt.getopt(sys.argv[3:], "h",["help", "hadd", "dataQCD", "dataTTbar", "dataWJets"])
  except getopt.GetoptError, err:
    print str(err)
    #usage()
    sys.exit(2)
  
  getDataQCD = False
  getDataTTbar = False
  getDataWJets = False
  runHistoAdd = False
  for o, a, in opt:
    if o in ("--hadd"):
      runHistoAdd = True
    if o in ("--dataQCD"):
      getDataQCD = True
    if o in ("--dataTTbar"):
      getDataTTbar = True
    if o in ("--dataWJets"):
      getDataWJets = True
    

  signal = sys.argv[1]
  logNames = sys.argv[2]
  sources = ['DYJetsToLeptonLepton', 'ZPrime', 'QCD', 'WJets', 'TTJets_MADGRAPH','Data', 'WW', 'WZ', 'ZZ']
  eventsDic = getContributions("Events @", "Values", logNames, sources)
  errorsDic = getContributions("Events @", "Errors", logNames, sources)
  
  #if get the contribution from data replace the values in the dictionary to calculate scale factors
  #if getDataQCD:
  #  eventsDic['qcd'] = getDataContribution('Final QCD', 'logFiles/data'+logNames+'.log')[0]
  #  errorsDic['qcd'] = getDataContribution('Final QCD', 'logFiles/data'+logNames+'.log')[1]
  if getDataTTbar:
    print 'Correcting TTBar contribution '
    correctBtag(logNames, sources)
    eventsDic['ttbar'] = getDataContribution('Final TTbar', 'logFiles/Data_'+logNames+'.log')[0]
    errorsDic['ttbar'] = getDataContribution('Final TTbar', 'logFiles/Data_'+logNames+'.log')[1]
  #if getDataWJets:
  #  eventsDic['wjets'] = getDataContribution('Final WJets', 'logFiles/data'+logNames+'.log')[0]
  #  errorsDic['wjets'] = getDataContribution('Final WJets', 'logFiles/data'+logNames+'.log')[1]
   
  bgSum = 0
  bgSigmaSqrdSum = 0
  for it in eventsDic.keys():
    if it != signal and it != 'Data':
      bgSum += eventsDic[it]
      bgSigmaSqrdSum += math.pow(errorsDic[it],2)
  controlData = eventsDic['Data'] - bgSum
  controlDataSigma = math.sqrt(math.pow(errorsDic['Data'],2) + bgSigmaSqrdSum)
  
  totalMC = eventsDic[signal] + bgSum
  totalMCSigma = math.sqrt(math.pow(errorsDic[signal],2) + bgSigmaSqrdSum)
  purity = eventsDic[signal]/totalMC
  scale = controlData/eventsDic[signal]
  scaleSigma = (1./eventsDic[signal])*math.sqrt(math.pow(controlDataSigma,2) + math.pow(scale,2)*math.pow(errorsDic[signal],2))

  print ' Total Data\t {0:.1f}'.format(eventsDic['Data'])+ ' \\pm ' + '{0:.2f}'.format(errorsDic['Data'])
  print ' Total MC\t {0:.2f}'.format(totalMC) + ' \\pm ' + '{0:.2f}'.format(totalMCSigma)
  print ' ScaleFactor\t {0:.3f}'.format(scale) + ' \\pm ' + '{0:.3f}'.format(scaleSigma)
  print ' Purity\t\t {0:.2f}'.format(purity)
  

  #runHistoAdd = raw_input('Add rootFiles? Y/N\n')
  #if runHistoAdd == 'Y':
  if runHistoAdd:
    from subprocess import call
    rootFilesDir = '/home/ajohnson/CMS/elecTauTauAnalysis/RootFiles'
    tmpRootDir = '/home/ajohnson/CMS/elecTauTauAnalysis/tmpRootDir'
    newFileName = rootFilesDir + 'dataWithMC' + logNames + '.root'
    if os.path.exists(newFileName):
      call(['hadd -f ' + newFileName + ' ' + tmpRootDir + '*'+ logNames + '.root'], shell=True)
    else:
      call(['hadd ' + newFileName + ' ' + tmpRootDir + '*'+ logNames + '.root'], shell=True)

def getDataContribution(source, dataLog):
  theFile = open(dataLog, 'r')
  thePair = []
  for line in theFile:
    if source in line:
      thePair.append(float(line.rsplit(':')[1].split('+-')[0].strip()))
      thePair.append(float(line.rsplit(':')[1].split('+-')[1].strip()))
  theFile.close()
  return thePair

if __name__ == "__main__":
  main(sys.argv[1:])
