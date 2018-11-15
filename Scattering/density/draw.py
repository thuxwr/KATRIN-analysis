# Draw density distribution versus z. Linear interpolation under log scale.
from ROOT import TH1D, TGraph
from array import array
import math

x = []
y = []
for line in open("dens.dat","r"):
    if line[0]=='#':
        continue
    x.append(float(line.split(',')[0]))
    y.append(float(line.split(',')[1]))

z = [math.log10(y[i]) for i in range(len(y))]
xarray = array('f',x)
zarray = array('f',z)
graph = TGraph(len(x),xarray,zarray)
#graph.SaveAs("test.root")

h = TH1D("","",800,0,8)
h.SetName("density")
h.SetTitle("Relative density distribution")
for nbin in range(1, 801):
    center = h.GetBinCenter(nbin)
    h.SetBinContent(nbin, math.pow(10, graph.Eval(center)))

h.GetXaxis().SetTitle("z[m]")
h.GetYaxis().SetTitle("n/n_{0}")
h.SaveAs("density.root")



