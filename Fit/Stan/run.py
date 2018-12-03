# Generate MC sample and fit with PyStan.
#
# Weiran, Nov.18, 2018.

import pystan, os, math
import numpy

# Left blank for automatic cache.
#
#
#
#
#################################


# Read input data from MC sample.
data = {}
rate = []
error = []
entries = []
Time = []

try:
    KATRINpath = os.getenv('KATRIN')
except KeyError:
    raise KeyError("'KATRIN' environment variable is not defined.")


for line in open(KATRINpath+'/Data/sample.dat'):
    if line[0] == '#':
        continue
    entries.append(float(line.split()[1]))
    Time.append(float(line.split()[2]))
    rate.append(entries[-1]/Time[-1])
    error.append(math.sqrt(entries[-1])/Time[-1])

data['nbins'] = len(rate)
data['rate'] = rate
data['error'] = error

# Configure Stan model and fitter.
sm = pystan.StanModel(file='model.stan', includes=['Pred.h'], allow_undefined=True, include_dirs=['.'], verbose=False)
fit = sm.sampling(data=data, sample_file = 'sample', iter=4000, warmup=1000, chains=1)
print(fit.stansummary(digits_summary=4))
print("95% upper limit: ")
print(fit.summary(probs=(0.95,)))







