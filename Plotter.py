from ROOT import *
from math import *
from numpy import *
from ROOT import TGraphErrors, TAttAxis
from array import array

#Samples = ['Zep']
Samples = ['Zep', 'Zek']
PtIntervals = [0, 30, 60, 90, 120, 200, 300]
NtrkIntervals = [0, 15, 30, 45, 60, 75, 90]

for s in Samples:
## Load root file
	print 'Opening ' + str(s) + ' file'
	FileName 	= str(s) + '.root'
	File 		= TFile(FileName, 'READ')
## Load Tree
	TreeName 	= 'T' + str(s)
	print 'Loading ' + TreeName + ' tree'
	Tree 		= File.Get(TreeName)
## Get Full Mass Histogram
	print 'Getting full mass histogram'
	HistTitle	= str(s) + ' Invariant Mass Distribution; M#_{' + str(s) + '}; N#_{Evts}'
	FullMass 	= TH1F('FullMass', HistTitle, 200, 60, 120)
	DrawVar		= str(s) + '.M() >> FullMass'
	SCutVar		= str(s) + '.M() < 120. && ' + str(s) + '.M() > 60.'
	CutVar		= TCut(SCutVar)
	Tree.Draw(DrawVar, CutVar)
## Get Mass vs Pt Histogram
	print 'Getting mass vs probe pt histogram'
	HistTitle	= str(s) + ' Invariant Mass Distribution vs Probe Pt'
	NPtBins		= len(PtIntervals) - 1
	MassPt		= TH2F('MassPt', HistTitle, 200, 60, 120, NPtBins, array('d', PtIntervals))
	DrawVar		= str(s) + '_O2.Pt():' + str(s) + '.M() >> MassPt'
	SCutVar		= str(s) + '.M() < 120. && ' + str(s) + '.M() > 60.'
	CutVar		= TCut(SCutVar)
	Tree.Draw(DrawVar, CutVar)
## Get Mass vs MET
	print 'Getting mass vs MET histogram'
	HistTitle	= str(s) + ' Invariant Mass Distribution vs Event MET'
	NPtBins		= len(PtIntervals) - 1
	MassMET		= TH2F('MassMET', HistTitle, 200, 60, 120, NPtBins, array('d', PtIntervals))
	DrawVar		= str(s) + '_MET.Mod():' + str(s) + '.M() >> MassMET'
	SCutVar		= str(s) + '.M() < 120. && ' + str(s) + '.M() > 60.'
	CutVar		= TCut(SCutVar)
	Tree.Draw(DrawVar, CutVar)
## Get Mass vs NTracks
	print 'Getting mass vs ntrks histogram'
	HistTitle	= str(s) + ' Invariant Mass Distribution vs Number of Tracks of PV'
	NNtrkInt	= len(NtrkIntervals) - 1
	MassNTrks	= TH2F('MassNTrks', HistTitle, 200, 60, 120, NNtrkInt, array('d', NtrkIntervals))
	DrawVar		= str(s) + '_ntrks:' + str(s) + '.M() >> MassNTrks'
	SCutVar		= str(s) + '.M() < 120. && ' + str(s) + '.M() > 60.'
	CutVar		= TCut(SCutVar)
	Tree.Draw(DrawVar, CutVar)
## Save Plots
	print 'Writing output file \n'
	OutFileName = str(s) + '_Plots.root'
	OutFile 	= TFile(OutFileName, 'RECREATE')
	OutFile.cd()
	FullMass.Write()
	MassPt.Write()
	MassMET.Write()
	MassNTrks.Write()
	OutFile.Write()
	OutFile.Close()
