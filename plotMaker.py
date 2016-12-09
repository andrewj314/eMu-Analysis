#!/usr/bin/python
import sys, getopt, math, os
sys.path.append('/home/ajohnson/CMS/hiMassAnalysis/elecTauTauAnalysis/scripts/python')
from getContributions import getContributions
from haddFiles import haddFiles
from setTdrStyle import setTDRStyle
from drawTCanvas import printCanvas

def main(argv):
  try:
    opt, args = getopt.getopt(sys.argv[3:], "hlr:d:",["help", "log", "rebin=", "cutLine=", "cutLine2=", "desc=", "outDir=",
    "qcdSF=", "ttbarSF=", "wjetsSF=", "wwSF=", "wzSF=", "zzSF=", "zttSF=", "zeeSF=", "zprimeSF="])
  except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit(2)

  rebinFactor = 1
  doLog = False
  drawLine = False
  lineMark = 0
  drawLine2 = False
  lineMark2 = 0
  printDesc = False
  theReq = ""
  theOutDir = "./plots"
 
  # define factors to scale plots by 
  qcdSF = 1.
  ttbarSF = 1.
  wjetsSF = 1.
  wwSF = 1.
  wzSF = 1.
  zzSF = 1.
  zttSF = 1.
  zeeSF = 1.
  zprimeSF = 1.
  #
  
  for o, a, in opt:
    if o in ("-h", "--help"):
      usage()
      exit(2)
    if o in ("-l", "--log"):
      doLog = True
    if o in ("-r", "--rebin"):
      rebinFactor = int(a)
    if o in ("--cutLine"):
      drawLine = True
      lineMark = float(a)
    if o in ("--cutLine2"):
      drawLine2 = True
      lineMark2 = float(a)
    if o in ("-d", "--desc"):
      printDesc = True
      theReq = a
    if o in ("--outDir"):
      theOutDir = a
    if o in ("--qcdSF"):
      qcdSF = float(a)
    if o in ("--ttbarSF"):
      ttbarSF = float(a)
    if o in ("--wjetsSF"):
      wjetsSF = float(a)
    if o in ("--wwSF"):
      wwSF = float(a)
    if o in ("--wzSF"):
      wzSF = float(a)
    if o in ("--zzSF"):
      zzSF = float(a)
    if o in ("--zttSF"):
      zttSF = float(a)
    if o in ("--zeeSF"):
      zeeSF = float(a)
    if o in ("--zprimeSF"):
      zprimeSF = float(a)  
      
  if len(sys.argv) < 2:
    usage()
    exit(2)
  
  import ROOT   
 
  theName = sys.argv[1]
  theHistoVar = sys.argv[2]

  sources = ['ZPrimeSSM_M-1500', 
  'TTJets_MADGRAPH', 'WJets', 'QCD', 'DYJetsToLeptonLepton', 'WW', 'WZ', 'ZZ']
  theContributions = getContributions("Events @", "Values", theName, sources)

  histoColors = {sources[0]: ROOT.kWhite, sources[1]: ROOT.kGreen, sources[2]: ROOT.kBlue, sources[3]: ROOT.kYellow, sources[4]: ROOT.kMagenta, sources[5]: ROOT.kCyan, sources[6]: ROOT.kGreen-8, sources[7]: ROOT.kBlue-8}
  legendEntry = {sources[0]: 'SSMToTauTau1500', sources[1]: 'TTJets', sources[2]: 'WJetsToLNu', sources[3]: 'QCD', sources[4]: 'DYJetsToLL_M50', sources[5]: 'WW', sources[6]: 'WZ', sources[7]: 'ZZ'}
  mcSF = {sources[0]: zprimeSF, sources[1]: ttbarSF, sources[2]: wjetsSF, sources[3]: qcdSF, sources[4]: zttSF, sources[5]: wwSF, sources[6]: wzSF, sources[7]: zzSF}

  setTDRStyle()    
  

  theCanvas = ROOT.TCanvas('theCanvas')
  theCanvas.Range(0,0,1,1)
  theCanvas.SetFillColor(0)
  theCanvas.SetBorderMode(0)
  theCanvas.SetBorderSize(2)
  theCanvas.SetFrameBorderMode(0)
  thePad = ROOT.TPad("thePad","thePad",0,0.2,1,1)
  thePad.Draw()
  thePad.cd()
  #thePad.Range(-3.742138,-40,3.754717,5677.362)
  thePad.Range(-3.742138,-40,3.754717,7000)
  
  thePad.SetFillColor(0)
  thePad.SetBorderMode(0)
  thePad.SetBorderSize(2)
  thePad.SetBottomMargin(0.005)
  thePad.SetFrameBorderMode(0)
  theRootFilesDir = './combRootFilesDir/'
  #theRootFileName = theName +".root"
  theRootFileName = haddFiles('tmpRootDir/*' + theName + '.root', theRootFilesDir + theName + '.root')
  theFile = ROOT.TFile(theRootFileName)
  #theFile = ROOT.TFile('tmpRootDir/CorrectedNBtag.root')
  
  theLegend = ROOT.TLegend(0.65, 0.525, 0.9, 0.825)
  theLegend.SetFillStyle(0)
  theLegend.SetShadowColor(0)

  ## get data histo
  #theDataHistoName = "DataHistosDirectory/" + theHistoVar
  #print theDataHistoName
  #theDataHisto = theFile.Get(theDataHistoName)
  #theDataHisto.Rebin(rebinFactor)
  #theHistoName = theDataHisto.GetTitle()
  #theDataHisto.SetLineColor(1)
  #theDataHisto.SetLineWidth(0)

  ##get signal histo - for cases when data is not available
  theSignalHistoName = "ZPrimeSSM_M-1500HistosDirectory/" + theHistoVar
  print theSignalHistoName
  theSignalHisto = theFile.Get(theSignalHistoName)
  theSignalHisto.Rebin(rebinFactor)
  theHistoName = theSignalHisto.GetTitle()
  

  ## Get Different MC sources
  sfactors = {}
  histoNames = {}
  theHistos = {}
  theHistoClones = {}
  
  theStack = ROOT.THStack("stackedHisto","")
  
  ## Make the summed Histo
  #nbins = theDataHisto.GetNbinsX()
  #lowEdge = theDataHisto.GetBinLowEdge(1)
  #hiEdge = theDataHisto.GetNbinsX() * theDataHisto.GetBinWidth(1) + lowEdge
  nbins = theSignalHisto.GetNbinsX()
  lowEdge = theSignalHisto.GetBinLowEdge(1)
  hiEdge = theSignalHisto.GetNbinsX() * theSignalHisto.GetBinWidth(1) + lowEdge
  
  ##only valid for histos with constant binning
  theHistoSum = ROOT.TH1F("theHistoSum", "", nbins, lowEdge, hiEdge) 
  theHistoSum.SetMarkerStyle(0)
  
  theHistoName = ''
  binWidth = 0
  for index, value in enumerate(sources):
    print value
    sfactors[value] = theContributions[value] * mcSF[value]
    print 'sfactor:\t {0:.1f}'.format(sfactors[value])
    theHistos[value] = theFile.Get(sources[index] + 'HistosDirectory/' + theHistoVar)
    theHistos[value].Rebin(rebinFactor)
    theHistoClones[value] = theHistos[value].Clone(sources[index] + 'HistosClone')
    theHistoClones[value].Sumw2()
    if index == 0:
      theHistoName = theHistos[value].GetTitle()
      binWidth = theHistos[value].GetBinWidth(1)
    theHistos[value].SetLineStyle(1)
    theHistos[value].SetFillColor(histoColors[value])
    if theHistos[value].Integral() > 0:
      print 'histo integral = \t {0:.1f}'.format(theHistos[value].Integral())
      #theHistos[value].Scale(sfactors[value]/theHistos[value].Integral())
      #theHistos[value].Scale(1/theHistos[value].Integral())
      
      #theHistoClones[value].Scale(sfactors[value]/theHistoClones[value].Integral())
      #theHistoClones[value].Scale(1/theHistoClones[value].Integral())
      
    theStack.Add(theHistos[value])
    theHistoSum.Add(theHistoClones[value])
    theLegend.AddEntry(theHistos[value], legendEntry[value], "f")
  
  #theMaxDataBin = theDataHisto.GetMaximum()
  theMaxDataBin = theSignalHisto.GetMaximum()
  theMaxDataBinPlusError = theMaxDataBin + math.sqrt(theMaxDataBin)
  theStack.SetMaximum(theMaxDataBinPlusError + 0.1*theMaxDataBinPlusError)
  theHistoSum.SetFillColor(ROOT.kGray + 2)
  theHistoSum.SetFillStyle(3001)
  
  #theRatioDataHisto = theDataHisto.Clone()  
  theRatioDataHisto = theSignalHisto.Clone()
  
  theRatioHisto = ROOT.TH1F("theRatio", "", nbins, lowEdge, hiEdge)
  theRatio = theRatioDataHisto
  theRatio.Divide(theHistoSum)
  theRatioHisto.Add(theRatio)
  #theRatio.SetMarkerStyle(1)
  
  #theStack.Draw()
  #theDataHisto.SetMarkerStyle(1)
  #theDataHisto.SetMarkerColor(0)
  #theDataHisto.Draw("ep same")
  theHistoSum.Draw("e2")
  theStack.Draw("same")
  theStack.GetXaxis().SetTitle(theHistoName)
  theStack.GetYaxis().SetTitle('Events/{0:.2f}'.format(binWidth))
  
  theLegend.Draw()
  
  if drawLine:
    theCutLine = ROOT.TLine(lineMark, 0.0, lineMark, theMax)
    drawTLine(theCutLine)
  if drawLine2:
    theCutLine2 = ROOT.TLine(lineMark2, 0.0, lineMark2, theMax)
    drawTLine(theCutLine2)

  cmsPre = ROOT.TLatex()
  cmsPre.SetNDC()
  cmsPre.SetTextAlign(12)
  cmsPre.DrawLatex(0.18, 0.92, "Jan 22 ReReco: 19.43 fb^{-1}")
  
  theCanvas.cd()
  
  theRatioPad = ROOT.TPad("theRatioPad","theRatioPad",0,0,1,0.2);
  theRatioPad.Draw();
  theRatioPad.cd()
  #theRatioPad.Range(-3.742138,-10.,3.754717,3.735869)
  theRatioPad.SetFillColor(0)
  theRatioPad.SetBorderMode(0)
  theRatioPad.SetBorderSize(0)
  theRatioPad.SetTopMargin(0.013)
  theRatioPad.SetBottomMargin(0.4)
  theRatioPad.SetFrameBorderMode(0)
  
  theRatioHisto.SetMarkerStyle(21)
  theRatioHisto.GetXaxis().SetTitle(theHistoVar);
  theRatioHisto.GetXaxis().SetLabelFont(63);
  theRatioHisto.GetXaxis().SetLabelSize(16);
  theRatioHisto.GetXaxis().SetTitleSize(0.1);
  theRatioHisto.GetXaxis().SetTitleFont(42);
  theRatioHisto.GetYaxis().SetLabelFont(63);
  theRatioHisto.GetYaxis().SetLabelSize(16);
  theRatioHisto.GetYaxis().SetTitleSize(0.06);
  theRatioHisto.GetYaxis().SetTitleFont(42);
  theRatioHisto.GetZaxis().SetLabelFont(42);
  theRatioHisto.GetZaxis().SetLabelSize(0.035);
  theRatioHisto.GetZaxis().SetTitleSize(0.035);
  theRatioHisto.GetZaxis().SetTitleFont(42);  
  
  theRatioHisto.Draw("ep");
  theRatioLine = ROOT.TLine(lowEdge, 1, hiEdge, 1);
  drawTLine(theRatioLine);
  theCanvas.cd();
  
  if printDesc:
    req = ROOT.TLatex()
    req.SetNDC()
    req.SetTextAlign(12)
    req.SetTextSize(0.04)
    req.SetTextColor(ROOT.kGray+3)
    req.DrawLatex(0.3, 0.85, theReq)
    
  theOutName = theHistoVar + "Stacked_" + theName
  
  printCanvas(theCanvas, theOutName + ".png", theOutDir)
  
  if doLog:
    theCanvas.SetLogy()
    theOutName = theOutName + "LogY"
    theCanvas.Update()
    printCanvas(theCanvas, theOutName + ".png", theOutDir)

def usage():
  print "Usage: test.py <name> <var>"
  print "Options: \n\t --log \n\t --rebin [n] \n\t --cutLine [cutValue] \n\t --cutLine2 [cutValue\
  \n\t --desc [Plot Description] \n\t --outDir [Output Directory (defaults to ./plots)]" 


def drawTLine(theTLine):
   theTLine.SetLineColor(2)
   theTLine.SetLineWidth(3)
   theTLine.SetLineStyle(3)
   theTLine.Draw()
  

if __name__ == "__main__":
  main(sys.argv[1:])
