from ROOT import TH2D
import math
from numpy import random

hist = TH2D("","",250, -0.25, 0.25, 250, 0, 0.25)
hist.SetName("hehe")
for i in range(250):
    truemass = float(i)/1000.
    mass = []
    lr = []
    for line in open("result/core%d.dat"%i, "r+"):
        if line[0]=='#':
            continue
        mass.append(float(line.split()[0]))
        lr.append(float(line.split()[1]))

    lr_mass = [(a,b) for a, b in zip(lr, mass)]
    lr_mass.sort()

    for j in range(math.ceil(0.95*len(mass))):
        hist.Fill(lr_mass[j][1], truemass)

hist.SaveAs("test.root")




