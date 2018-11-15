# Get energy loss pdf for a given scattering times.
# Elastic scattering is neglected since they can be regarded as an additional systematic shift on endpoint.
#
# Weiran, Nov.12, 2018

from ROOT import TFile, TH1D
import os

try:
    KATRINPath = os.environ['KATRIN']
except KeyError:
    raise KeyError("'KATRIN' environment variable is not defined.")

infile = TFile(KATRINPath+"/Scattering/energyloss/EnergyLoss.root", "READ")
hist = []
for i in range(3):
    hist.append(infile.Get("scat%d"%i))

def GetPdf(ScatterTimes, epsilon):
    return hist[ScatterTimes-1].GetBinContent(hist[ScatterTimes-1].FindBin(epsilon))

