import math, ROOT

def getEfficiency(passed, total):
  theEff = passed/total
  theEffError = math.sqrt(theEff*(1.-theEff)/total)
  if (theEff + theEffError) >= 1. or (theEff - theEffError <= 0.):
    #do one bin histo && extract bayes eff
    theNumHisto = ROOT.TH1F("theNumHisto","theNumHisto",1,0,1)
    theNumHisto.SetBinContent(1,passed)
    theNumHisto.Sumw2()
    theDenHisto = ROOT.TH1F("theDenHisto","",1,0,1)
    theDenHisto.SetBinContent(1,total)
    theDenHisto.Sumw2()
    bayesEff = ROOT.TGraphAsymmErrors()
    bayesEff.BayesDivide(theNumHisto,theDenHisto,"b")
    effErrorHigh = bayesEff.GetErrorYhigh(0)
    effErrorLow = bayesEff.GetErrorYlow(0)
    theEffError = max(effErrorHigh, effErrorLow)
  theEffAndError = [theEff, theEffError]
  return theEffAndError

def getSurvivingEvents(nEvents, nEventsError, eff, effError):
  nSurviving = nEvents*eff
  errorFirstTerm = nEvents*effError
  errorSecondTerm = eff*nEventsError
  nSurvivingError = math.sqrt((errorFirstTerm*errorFirstTerm)+(errorSecondTerm*errorSecondTerm))
  survivingEventsAndError = [nSurviving, nSurvivingError]
  #print 'NEvents\t' + str(nSurviving) + '\t+-\t' + str(nSurvivingError)  
  return survivingEventsAndError

def getCorrectedEff(eff1, eff1Err, eff2, eff2Err):
  effRatio = eff1/eff2
  effRatioSqrd = effRatio*effRatio
  
  eff1ErrSqrd = eff1Err*Eff1Err
  eff2ErrSqrd = eff2Err*Eff2Err
  
  effRatioErr = math.sqrt(eff1ErrSqrd + (effRatioSqrd*eff2ErrSqrd)/eff1Err)
  
  effRatioAndError = [effRatio,effRatioErr]
  return effRatioAndError
  
