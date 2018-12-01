/*
	 Generate MC sample and fit to model using Minuit.

	 Weiran, Nov.17, 2018.
*/

#include "../../Simulation/Simulation.h"
#include "../../Configure/Configure.h"

Simulation sim;
TMinuit minuit;
KATRIN katrin;

double pmass = 0;
double pendpoint = 18574;

double* Time;
double* Voltage;
double* sample;
double* Rate;
double* error;
double* theory;
int nvoltage;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);

using namespace TMath;

void SimpleFit()
{
	/* Input time distribution and threshold. Generate or read data. One can also change the time and threshold. */
	Time = katrin.Time;
	Voltage = katrin.Voltage;
	nvoltage = katrin.Nbins;

	/* This block for simulated sample. */
	sample = sim.Generate(pmass, pendpoint, nvoltage, Voltage, Time);
	Rate = new double[nvoltage];
	error = new double[nvoltage];
	for(int i=0; i<nvoltage; i++) {
		Rate[i] = sample[i]/Time[i];
		error[i] = sqrt(sample[i])/Time[i];
	}
	/************************************/

	/* This block for real data.
	Rate = katrin.Rate;
	error = katrin.Error;
	*/

	minuit.DefineParameter(0, "mass", pmass, 0.001, 0, 0);
	minuit.DefineParameter(1, "endpoint", pendpoint, 0.01, pendpoint-5, pendpoint+5);
	minuit.DefineParameter(2, "A_sig", 1, 1e-6, 0, 0);
	minuit.DefineParameter(3, "A_bkg", 1, 1e-6, 0, 0);
	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	minuit.Migrad();

}

void fcn(int &npar, double* gin, double &f, double* par, int iflag) {
	delete theory;
	theory = sim.Asimov(par[0], par[1], nvoltage, Voltage, par[2], par[3]);
	double chi2 = 0;
	for(int i=0; i<nvoltage; i++) {
		double pred = theory[i];
		double meas = Rate[i];
		double Error = error[i];
		chi2 += pow((pred-meas)/Error, 2);
	}
	f = chi2;
}

