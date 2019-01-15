# KATRIN-analysis
A package for Bayesian analysis of KATRIN experiment.

ROOT classes used in this package are: TFile, TH1D, TH2D, TGraph, TF1, TMath, TRandom3 and TMinuit. Please make sure to include the according headers in your PyStan.

This package uses [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) for Bayesian inference, and [ROOT](https://root.cern.ch) is adopted in calculation. Please make sure to install these packages correctly in advance.

For ROOT, fftw3 and minuit2 should be turned on during installation. ROOT classes used in this package are: TFile, TH1D, TGraph, TF1, TMath, TRandom3, TMinuit and TVirtualFFT. Please link the corresponding libraries when sampling.

This package also relies on KASPER, the KATRIN analysis and simulation package. You need to install Kommon and KSC before using this package.

This package is designed for MPI multiprocessing usage on clusters, and may fail on Mac OS.

# How-to
After installing all dependencies and this package, you need to setup your compiler and environment variables. A global environment variable is ${KATRIN}, which should be the local path to this package. You should link all CmdStan, ROOT, Kasper libraries required while compiling Stan model. You may modify the Stan model in Fit/Stan/MPI, and use stanc to generate a Stan header. Then you may user mpich to compile the generated header with main.cpp together, and run Stan sampling by using mpiexec or mpirun.
