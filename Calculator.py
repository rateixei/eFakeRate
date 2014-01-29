from ROOT import *
from math import *
from numpy import *
from ROOT import TGraphErrors, TAttAxis
from array import array
from ROOT import RooFit
from ROOT import RooRealVar, RooArgList, RooArgSet, RooLinkedList, RooDataHist, RooGenericPdf, RooVoigtian

def NEvtsCalc( uPars, uEPars):
	Pars 			= {}
	Pars[0]			= RooRealVar('c1', 'Exponential constant', -1, -10, 0)
	Pars[1] 		= RooRealVar('ExpY', 'Background Yield', 100, 0, 10000000)
	Pars[2]			= RooRealVar('Mean', 'Voigtian Mean' , 90.0, 20, 180.0)
	Pars[3]			= RooRealVar('Width', 'Voigtian Width' , 5.0, 0.5, 40.0)
	Pars[4]			= RooRealVar('Sigma', 'Voigtian Sigma' , 5.0, 0.5, 40.0)
	Pars[5]			= RooRealVar('VoY', 'Signal Yield', 100, 0, 10000000)
	if len(Pars) != len(uPars):
		print 'The input array has a weird number of entries...'
		return 0
	if len(uPars) != len(uEPars):
		print 'The input arrays are of different sizes...'
		return 0
	for x in xrange(0, 6):
		Pars[x].setVal(uPars[x])
		Pars[x].setError(uEPars[x])
	v			= RooRealVar('v', 'Invariant Mass (GeV)', 60, 120)
	Voigt 			= RooVoigtian('Voigt', 'Voigtian - Signal', v, Pars[2], Pars[3], Pars[4])
## Calculate integral from -2sigma to 2sigma
	VStDev_r		= Voigt.sigma(v)
	VStDev			= VStDev_r.getVal()
	v.setRange("sobRange", Pars[2].getVal() - 2.*VStDev, Pars[2].getVal() + 2.*VStDev)
	integralSig     	= Voigt.createIntegral(RooArgSet(v), RooFit.NormSet(RooArgSet(v)), RooFit.Range("sobRange"))
	FinalNumber		= integralSig.getVal()*Pars[5].getVal()
	return FinalNumber


Samples 		= ['Zep', 'Zek']
ParNames 		= ['c1', 'ExpY', 'Mean', 'Width', 'Sigma', 'VoY']
BinTypes		= ['Pt', 'MET', 'Trk']
BinN			= ['1', '2', '3', '4', '5', '6']

for s in Samples:
## Open output text file
	OutTxtName		= str(s) + '_Results.txt'
	TextFile		= open(OutTxtName, 'w')
## Open root file
	RootName		= str(s) + '_Fits.root'
	RootFile		= TFile(RootName, 'READ')
## Write Full Mass Results
	TextFile.write('FullMass \n')
	ResultName		= 'fitresult_Model_' + str(s) + '_FullMass'
	FitResult		= RootFile.Get(ResultName)
	VarList			= FitResult.floatParsFinal()
	uPars 			= {}
	uEPars 			= {}
	for x in xrange(0, 6):
		Var 		= VarList.find(str(ParNames[x]))
		uPars[x]	= Var.getVal()
		uEPars[x]	= Var.getError()
		Sentence	= str(ParNames[x]) + '\t' + str(Var.getVal()) + '\t' + str(Var.getError()) + '\n'
		TextFile.write(Sentence)
	NEvts 			= NEvtsCalc(uPars, uEPars)
	Sentence		= str('NEvts \t') + str(NEvts) + '\t' + '1 \n \n'
	TextFile.write(Sentence)
## Write Bin-wise distributions
	for b in BinTypes:
		Sentence	= str(b) + ' Histograms'
		TextFile.write(Sentence)
		for x in xrange(0,6):
			Sentence	= 'Bin: ' + str(BinN[x]) + '\n'
			TextFile.write(Sentence)
			ResultName	= str(b) + 'Bins/fitresult_Model_' + str(b) + 'Bin_' + str(BinN[x])
			FitResult 	= RootFile.Get(ResultName)
			VarList		= FitResult.floatParsFinal()
			for q in xrange(0, 6):
				Var 		= VarList.find(str(ParNames[q]))
				uPars[q]	= Var.getVal()
				uEPars[q]	= Var.getError()
				Sentence	= str(ParNames[q]) + '\t' + str(Var.getVal()) + '\t' + str(Var.getError()) + '\n'
				TextFile.write(Sentence)
			NEvts 			= NEvtsCalc(uPars, uEPars)
			Sentence		= str('NEvts \t') + str(NEvts) + '\t' + '1 \n'
			TextFile.write(Sentence)
		TextFile.write('\n')
