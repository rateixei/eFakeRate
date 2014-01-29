from ROOT import *
from math import *
from numpy import *
from ROOT import TGraphErrors, TAttAxis
from array import array
from ROOT import RooFit
from ROOT import RooRealVar, RooArgList, RooArgSet, RooLinkedList, RooDataHist, RooGenericPdf, RooVoigtian

def Drawer(data, Model, Exp, Voigt, Name):
	frame = v.frame()
	data.plotOn(frame)
	Model.plotOn(frame)
	argset = RooArgSet(Exp)
	Model.plotOn(frame, RooFit.Components(argset), RooFit.LineStyle(kDashed), RooFit.LineColor(3))
	argset = RooArgSet(Voigt)
	Model.plotOn(frame, RooFit.Components(argset), RooFit.LineColor(2))
	Model.paramOn(frame)
	c = TCanvas()
	frame.Draw()
	file_name = str('pdf/' + Name + '.pdf')
	c.SaveAs( file_name )

Samples = ['Zep', 'Zek']
BinTypes = ['Pt', 'MET', 'Trk']

## Roofit variables
v			= RooRealVar('v', 'Invariant Mass (GeV)', 60, 120)
c1 			= RooRealVar('c1', 'Exponential constant', -1, -10, 0)
Exp 			= RooGenericPdf('Exp','Exponential - Background', 'exp(v*c1)', RooArgList(v,c1))
ExpY 			= RooRealVar('ExpY', 'Background Yield', 100, 0, 10000000)

Mean 			= RooRealVar('Mean', 'Voigtian Mean' , 90.0, 20, 180.0)
Width 			= RooRealVar('Width', 'Voigtian Width' , 5.0, 0.5, 40.0)
Sigma 			= RooRealVar('Sigma', 'Voigtian Sigma' , 5.0, 0.5, 40.0)
Voigt 			= RooVoigtian('Voigt', 'Voigtian - Signal', v, Mean, Width, Sigma)
VoY 			= RooRealVar('VoY', 'Signal Yield', 100, 0, 10000000)
Model 			= RooAddPdf("Model", "Sum of Exp+Voigtian", RooArgList(Exp,Voigt), RooArgList(ExpY, VoY))

for s in Samples:
## Load root file
	FileName 	= str(s) + '_Plots.root'
	File 	 	= TFile(FileName, 'READ')
	OutName		= str(s) + '_Fits.root'
	OutFile		= TFile(OutName, 'RECREATE')
## Get Histograms
	FullMass 	= File.Get('FullMass')
	MassPt	 	= File.Get('MassPt')
	MassMET	 	= File.Get('MassMET')
	MassNTrks	= File.Get('MassNTrks')
## Fit full Mass
	FitName		= str(s) + '_FullMass'
	data 		= RooDataHist(FitName, FitName, RooArgList(v), FullMass)
	FullMassFit = Model.fitTo(data, RooFit.Save())
	FullMassFit.Write()
	FullMass.Write()
	Drawer(data, Model, Exp, Voigt, FitName)
## Fit Pt bins
	for y in BinTypes:
		HistoName 	= 'Mass' + str(y)
		if str(y) == 'Trk':
			HistoName = 'MassN' + str(y) + 's'
		Histo		= File.Get(HistoName)
		OutFile.cd()
		DirName		= str(y) + 'Bins'
		Dir		= OutFile.mkdir(DirName)
		Dir.cd()
		NBinsPt		= Histo.GetNbinsY() + 1
		for x in xrange(1, NBinsPt):
			BinFitName		= str(y) + 'Bin_' + str(x)
			BinHist			= Histo.ProjectionX(BinFitName, x, x)
			data 			= RooDataHist(BinFitName, BinFitName, RooArgList(v), BinHist)
			BinFitResult	= Model.fitTo(data, RooFit.Save())
			BinFitResult.Write()
			BinHist.Write()
			PDFName = str(s) + '_' + BinFitName
			Drawer(data, Model, Exp, Voigt, PDFName)
