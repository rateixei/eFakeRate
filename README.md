eFakeRate
=========

Code for calculation of fake rates of electrons faking photons

Instructions
____________

1) Run the eFakeSel code with simpleRunner.csh.
That will create skimmed ntuples with Z candidates, tag and probe, MET and ntrks information.

2) After the skimmed ntuples are created, use Plotter.py.
That will use the ntuples and create a root file with the important histograms.

3) Fitter.py
Will use the histograms created on step (2) and fit them.

4) Calculatorp.y
Will use the fits to calculate the number of events within the +2sigma and -2sigma of the Z peak.
