# Energy loss function and its self convolution.
from ROOT import RooRealVar, RooDataHist, TF1, TCanvas, TH1D
import math

# Generate a function describing energy loss function.
A1 = 0.204
e1 = 12.6
w1 = 1.85
A2 = 0.0556
e2 = 14.3
w2 = 12.5
ec = 14.09
def energy_loss(x):
    if x[0]<ec:
        return A1 * math.exp(-2*((x[0]-e1)/w1)**2)
    else:
        return A2 * w2**2/(w2**2+4*(x[0]-e2)**2)

# Scale.
hist = TH1D("","",100000,0,500)
for i in range(1,100001):
    binCenter = hist.GetBinCenter(i)
    hist.SetBinContent(i, energy_loss([binCenter]))
Scale = hist.Integral("width")

def EngLossPdf(epsilon):
    return energy_loss([epsilon])/Scale
#fcn = TF1("hehe",energy_loss,0,50)
#fcn.SetNpx(500)
#c1 = TCanvas()
#fcn.Draw()
#c1.SaveAs("test.pdf")
#print(fcn.Derivative(1))


