# Angular change due to elastic scattering, obtained from Fig 6.29, Grof_Stefan's thesis.

from ROOT import TGraph, TH1D
from array import array
import math, os

deg = []
pdf = []

try:
    KATRINPath = os.environ['KATRIN']
except KeyError:
    raise KeyError("'KATRIN' environment variable is not defined.")

for line in open(KATRINPath+'/Scattering/energyloss/elas/AngDist.dat','r'):
    if line[0]=='#':
        continue
    deg.append(float(line.split(',')[0]))
    pdf.append(float(line.split(',')[1]))

degarray = array('f', deg)
pdfarray = array('f', pdf)
graph = TGraph(len(deg), degarray, pdfarray)
#graph.SaveAs("test.root")

def GetAngDist(angle):
    return graph.Eval(angle)

m_e = 0.511/931.4940954 #u
m_T2 = 6.032 #u
E_0 = 18.574e3
def EnergyLoss(angle):
    return 2 * m_e / m_T2 * E_0 * (1-math.cos(angle/180. * math.pi))

# Test with bin width = 0.1 degree.
EngLossGraph = TGraph()
for i in range(200):
    angle = 0.05 + i * 0.1
    AngPdf = GetAngDist(angle)
    EngLoss = EnergyLoss(angle)
    EngLossGraph.SetPoint(i, EngLoss, AngPdf/math.sin(angle/180. * math.pi)) # Phase space changed.

# Rescale.
EngLossHist = TH1D("", "", 10000, 0, 0.2)
for i in range(1,10001):
    binCenter = EngLossHist.GetBinCenter(i)
    EngLossHist.SetBinContent(i, EngLossGraph.Eval(binCenter))
Scale = EngLossHist.Integral("width")

#EngLossGraph.SetName("hehe")
#EngLossGraph.SaveAs("test.root")
def EngLossPdf(epsilon):
    if epsilon > 0.2:
        return 0
    return EngLossGraph.Eval(epsilon)/Scale





