# Generate response fcn.
#
# Weiran, Nov.12, 2018.

from ROOT import TH1D, TH2D
import sys, math

sys.path.append('../../Scattering/density')
sys.path.append('../../Scattering/energyloss')
sys.path.append('../../Scattering/prob')

import GetProb, GetEngLoss

# Detector response is insensitive to a small fluctuation in B_A, B_S, B_max, so here just fix them.
B_A = 2.68
B_S = 2.52e4
B_max = 4.2e4
m_e = 0.511e6
pi = math.pi

def gamma(energy):
    return energy/m_e + 1

def transmission(E, theta, U):
    if E*(1-math.sin(theta)**2*B_A/B_S*(gamma(E)+1)/2) > U:
        return 1
    else:
        return 0

def Response(E, U):
    if E<=U:
        return 0
    thetamax = math.asin(math.sqrt(B_S/B_max))

    response = 0
    # Deal with once, twice and three times scattering.
    # TH2D with x:theta, y:epsilon. Nbins can be adjusted according to the required precision.
    histtmp = TH2D("","",500,0,thetamax,500,0,E-U)
    for xbin in range(1,501):
        theta = histtmp.GetXaxis().GetBinCenter(xbin)
        P = []
        for nscat in range(1,4):
            P.append(GetProb.P(nscat, 0, theta))
        for ybin in range(1, 501):
            epsilon = histtmp.GetYaxis().GetBinCenter(ybin)
            EnergyLossPdf = 0
            for nscat in range(1,4):
                EnergyLossPdf += P[nscat-1] * GetEngLoss.GetPdf(nscat, epsilon)
            content = transmission(E-epsilon, theta, U) * math.sin(theta) * EnergyLossPdf
            histtmp.SetBinContent(xbin, ybin, content)
    response += histtmp.Integral("width")

    # Add non-scattered contributions.
    histtmp2 = TH1D("","",500,0,thetamax)
    for xbin in range(1,501):
        theta = histtmp2.GetBinCenter(xbin)
        histtmp2.SetBinContent(xbin, transmission(E, theta, U) * math.sin(theta) * GetProb.P(0,0,theta))
    response += histtmp2.Integral("width")
    return response

# Construct TH2D containing response matrix. x:E from (-35+E0, 10+E0), y:U from (-35+E0, 10+E0)
hist = TH2D("","",500,-35,5,500,-35,5)
for xbin in range(1,501):
    for ybin in range(1,501):
        E = hist.GetXaxis().GetBinCenter(xbin)
        U = hist.GetYaxis().GetBinCenter(ybin)
        hist.SetBinContent(xbin, ybin, Response(E, U))
        print("Finish xbin: ", xbin, " ,ybin: ", ybin)

hist.SaveAs("test.root")




