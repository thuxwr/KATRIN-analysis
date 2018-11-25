# Generate MC sample and fit with PyStan.
#
# Weiran, Nov.18, 2018.

import pystan, os
from ROOT import gROOT, TFile, TH1D

# Left blank for automatic cache.
#
#
#
#
#################################

# Generate MC sample, save it to ROOT file, and read it via PyROOT interface.
mass = 0.35*0.35
endpoint = 18574
gROOT.ProcessLine('.L ../../Simulation/Simulation.h')
gROOT.ProcessLine('Simulation sim')
gROOT.ProcessLine('TH1D* hsim = sim.Generate(%f,%f)'%(mass, endpoint))
gROOT.ProcessLine('hsim->SetName("simulation")')
gROOT.ProcessLine('hsim->SaveAs("_tmp.root")')
gROOT.ProcessLine('.q')

datafile = TFile("_tmp.root", "READ")
datahist = (TH1D)(datafile.Get("simulation"))

# Read input data from MC sample.
data = {}
rate = []
error = []
for nbin in range(1, datahist.GetNbinsX()+1):
    rate.append(datahist.GetBinContent(nbin))
    error.append(datahist.GetBinError(nbin))

data['nbins'] = len(rate)
data['rate'] = rate
data['error'] = error

# Configure Stan model and fitter.
sm = pystan.StanModel(file='model.stan', includes=['Pred.h'], allow_undefined=True, include_dirs=['.'], verbose=False)
fit = sm.sampling(data=data, sample_file = 'sample', iter=1000, chains=1)
#a = sm.optimizing(data=data);
print(fit.stansummary(digits_summary=4))
#print(a)

# Delete tmp files.
os.remove("_tmp.root")






