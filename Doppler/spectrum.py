# Calculate beta decay spectrum broadened by doppler effect.
# Only for test. Calculation for broadened spectrum is written in cpp.
#
# Weiran, Nov.13, 2018.

from ROOT import TH1D, RooRealVar, TF1, RooFit, RooFFTConvPdf, RooAbsReal, RooAbsPdf, RooArgList, TCanvas, kRed, kBlue, RooDataHist, RooHistPdf, RooArgSet, TFile, TGraph

import math

pi = math.pi
m_e = 0.511e6
alpha = 0.00729735257
E_0_center = 18.574e3
B_A = 2.68
B_S = 2.52e4
B_max = 4.2e4

final_state_energy = [0.053, 0.124, 0.247, 0.351, 0.442, 0.556, 0.665, 0.759, 0.850, 0.937, 1.048,
                      1.143, 1.249, 1.359, 1.451, 1.552, 1.657, 1.745, 1.834, 1.940, 2.044, 2.144,
                      2.244, 2.344, 2.510, 2.762, 3.009, 3.257, 3.507, 3.757, 4.083, 4.579, 5.132, 
                      5.647, 18.500, 19.500, 20.696, 21.658, 22.627, 23.598, 24.573, 25.550, 26.529,
                      27.510, 28.493, 29.478, 30.464, 31.455, 32.490, 33.557, 34.534, 35.492]
final_state_prob = [
    0.0069, 0.0046, 0.0233, 0.0553, 0.0457, 0.2033, 0.1649, 0.3877, 0.3808, 0.6809, 1.1214, 1.0112, 2.4406,
    3.2337, 4.0864, 6.8745, 6.6279, 5.1412, 6.5561, 5.4588, 3.7231, 2.5473, 1.6959, 1.1369, 1.6947, 1.0094,
    0.5732, 0.2806, 0.1316, 0.0623, 0.0420, 0.0080, 0.0015, 0.0000, 0.0000, 0.0000, 0.0012, 0.0113, 0.0656,
    0.2567, 0.7149, 1.4804, 2.3583, 2.9715, 3.0307, 2.5527, 1.8080, 1.1070, 0.7377, 1.0637, 1.9095, 2.2178]
# Original spectrum.
def gamma(energy):
    return energy/m_e + 1

def beta(energy):
    return math.sqrt(1.-1./gamma(energy)/gamma(energy))

def Fermi_fcn(energy, Z):
    eta = alpha * Z / beta(energy)
    return 2 * pi * eta / (1-math.exp(-2*pi*eta))

def momentum(energy):
    return m_e * beta(energy)

def epsilon(energy, V_f, E_0):
    return E_0 - energy - V_f

def w(energy):
    return (energy+m_e)/m_e

def t(energy):
    Beta = beta(energy)
    return 1/Beta * 0.5 * math.log((1+Beta)/(1-Beta)) - 1 

def radiative_correction(energy, V_f, E_0):
    W = w(energy)
    T = t(energy)
    Beta = beta(energy)
    W_0 = w(E_0-V_f)
    if W_0<=W:
        return 0
    value = math.pow(W_0-W, 2*alpha/pi * T) * (1 + 2. * alpha /pi * ( \
        T * (math.log(2.) - 1.5 + (W_0-W)/W) + 0.25 * (T+1.) * (2.*(1.+Beta*Beta) + 2.*math.log(1-Beta) + \
        math.pow(W_0-W,2)/(6.*W*W)) - 2 + 0.5 * Beta - 17./36. * Beta * Beta + 5./6. * math.pow(Beta,3) ))
    return value
        
def shape_fcn(x, ms, E_0):
    return_value = 0
    for i in range(len(final_state_energy)):
        Epsilon = epsilon(x, final_state_energy[i], E_0)
        if Epsilon<=0 or Epsilon**2<=ms:
            continue
        return_value += final_state_prob[i] * Epsilon * math.sqrt(Epsilon**2 - ms) * \
            radiative_correction(x, final_state_energy[i], E_0)
    return_value *= Fermi_fcn(x, 2) * momentum(x) * (x+m_e)
    return return_value

def decayspec(ms, E_0):
    nbins = 10000
    spec = TH1D("","",nbins,-30+E_0_center,5+E_0_center)
    for i in range(1,nbins+1):
        energy = spec.GetBinCenter(i)
        spec.SetBinContent(i, shape_fcn(energy, ms, E_0))
    spec.Scale(1/spec.Integral())
    return spec

# Broaden core, g(E_cms, E_lab).
def core(E_cms, bv=13, t=30): # bv: bulk velocity in m/s, t: temperature in K.
    k_B = 8.6173303e-5 #eV/K
    M_T2 = 6.032 * 931.4940954e6 #eV
    c = 2.99792458e8
    sigma_v = math.sqrt(k_B*t/M_T2)
    v_cms = beta(E_cms)
    cos_thetamax = math.sqrt(1-B_S/B_max)
    def fcn(x): # E_lab = E_cms + x
        E_lab = E_cms + x[0]
        v_lab = beta(E_lab)
        v_M = (v_lab-v_cms)/(1-v_lab*v_cms)
        u = bv/c
        #print(math.erf((v_M-cos_thetamax*u)/(math.sqrt(2)*sigma_v)))#-math.erf((v_M-u)/(math.sqrt(2)*sigma_v)))
        g_vM = 1/(1-cos_thetamax)/2/u * (math.erf((v_M-cos_thetamax*u)/(math.sqrt(2)*sigma_v))-math.erf((v_M-u)/(math.sqrt(2)*sigma_v)))
        g_E = g_vM/(gamma(E_cms)*m_e*v_cms) # with unit eV^-1
        return g_E
    tf = TF1("core", fcn, -17.5, 17.5)
    return tf

# Final spectrum.
# The broaden kernel is almost the same for different E_cms, so here neglect the difference.
def broadenspec(ms, E_0):
    spec = decayspec(ms, E_0)
    smear = core(E_0_center)
    newspec = TH1D("","",spec.GetNbinsX(),-17.5,17.5)
    for i in range(1, spec.GetNbinsX()+1):
        newspec.SetBinContent(i, spec.GetBinContent(i))

    #x = RooRealVar("x","x",-30+E_0_center, 5+E_0_center)
    x = RooRealVar("x","x",-17.5, 17.5)
    data = RooDataHist("","",RooArgList(x),newspec)
    specpdf = RooHistPdf("","",RooArgSet(x),data)
    #y = RooRealVar("y","y",-30, 5)
    x.setBins(10000)
    smearpdf = RooFit.bindPdf(smear, x)
    fft = RooFFTConvPdf("tt","tt", x, specpdf, smearpdf)
    #fft.setShift(0, -18574)
    #c1 = TCanvas()
    #frame = x.frame()
    #fft.plotOn(frame)
    #frame.Draw()
    tf = fft.asTF(RooArgList(x))
    tf.SetNpx(10000)
    rtf = tf.Clone()
    return rtf

tf = broadenspec(0, E_0_center)
#tf.SetNpx(10000)
#integ = tf.Integral(-30+E_0_center, 5+E_0_center)
integ = tf.Integral(-17.5, 17.5)
f = TFile("t.root", "RECREATE")
#tf.SetName("broaden")
#tf.Write()
spec = decayspec(0, E_0_center)
spec.Scale(1/spec.Integral("width") * integ)
spec.SetName("org")
specgraph = TGraph(spec)
specgraph.SetName("original")
specgraph.Write()
graph = TGraph()
for i in range(10000):
    graph.SetPoint(i, spec.GetBinCenter(i+1), tf.Eval(spec.GetBinCenter(i+1)+17.5-E_0_center-5))
graph.SetName("broaden")
graph.Write()
#print(tf.Eval(18574))
#tf.SaveAs("t.root")




