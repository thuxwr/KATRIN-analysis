# Only inelastic components, in order to avoid too much calculation.

import sys

sys.path.append('../elas')
sys.path.append('../inelas')

import elas, inelas
from ROOT import RooRealVar, RooDataHist, TF1, TCanvas, RooFit, RooFFTConvPdf, RooAbsReal, RooAbsPdf, RooArgList, TFile, TH1D

#for i in range(1, 50001): #binwidth = 0.001eV
elascs = 0.29e-18
inelascs = 3.4e-18
npx = 5000000
def totalCS(epsilon):
    #return elas.EngLossPdf(epsilon[0]) * elascs + inelas.EngLossPdf(epsilon[0]) * inelascs
    return inelas.EngLossPdf(epsilon[0]) * inelascs

# Range from 0 to 500 eV to promote precision of FFT and avoid unphysical results.
fcn = TF1("pdf", totalCS, 0, 500)
fcn.SetNpx(npx)
x = RooRealVar("x","x",0,500)
x.setBins(npx)
f = RooFit.bindPdf(fcn, x)
tt = RooFFTConvPdf("tt","tt",x,f,f)
ttt = RooFFTConvPdf("ttt","ttt",x,tt,f)

tf = []
tf.append(fcn)
tf.append(tt.asTF(RooArgList(x)))
tf.append(ttt.asTF(RooArgList(x)))

pdf = []
Times = ["once", "twice", "thrice"]
wfile = TFile("../EnergyLoss.root", "RECREATE")
for i in range(len(tf)):
    tf[i].SetNpx(npx)
    histtmp = tf[i].GetHistogram()
    histtmp.Scale(1/histtmp.Integral("width"))
    # Change energy range to (0, 50)eV.
    hist = TH1D("","", int(npx/10), 0, 50)
    for nbin in range(1, int(npx/10)+1):
        hist.SetBinContent(nbin, histtmp.GetBinContent(nbin))
    hist.SetName("scat%d"%i)
    hist.GetXaxis().SetTitle("Energy loss [eV]")
    hist.GetYaxis().SetTitle("f(#varepsilon)")
    hist.SetTitle("Energy loss from inelastic scattering %s"%Times[i])
    hist.SetStats(0)
    hist.Write()

wfile.Close()
print("All energy loss pdf's generated.")
