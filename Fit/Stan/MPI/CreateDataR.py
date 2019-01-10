import os, math

rate = []
error = []
entries = []
Time = []


try:
    KATRINpath = os.getenv("KATRIN")
except KeyError:
    raise KeyError("'KATRIN' environment variable is not defined.")

for line in open(KATRINpath+'/Data/sample.dat'):
    if line[0] == '#':
        continue
    entries.append(float(line.split()[1]))
    Time.append(float(line.split()[2]))
    rate.append(entries[-1]/Time[-1])
    error.append(math.sqrt(entries[-1])/Time[-1])

# Dump into R format
f = open("data.R", "r+")
f.write("nbins <- \n")
f.write("\t%d\n"%len(rate))
f.write("rate <- \n")
f.write("\tc(")
for i in range(len(rate)):
    f.write("%f"%(rate[i]))
    if i != len(rate)-1:
        f.write(",")

f.write(")\n")
f.write("error <- \n")
f.write("\tc(")
for i in range(len(error)):
    f.write("%f"%(error[i]))
    if i != len(error)-1:
        f.write(",")

f.write(")")

