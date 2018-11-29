from ROOT import TH1D, TGraph
from array import array
import math

x = []
y = []
error = []

for line in open("3years.txt", "r+"):
    x.append(float(line.split()[0]))
    entries = float(line.split()[1])
    time = float(line.split()[2])/86400.
    y.append(entries/time)
    error.append(math.sqrt(entries)/time)

th = TGraph()
for npoint in range(len(x)):
    th.SetPoint(npoint, x[npoint], y[npoint])
    th.Set

