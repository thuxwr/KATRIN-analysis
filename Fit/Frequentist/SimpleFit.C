/*
	 Generate MC sample and fit to model using Minuit.

	 Weiran, Nov.17, 2018.
*/

#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "../../Simulation/Simulation.h"
#include "../../Configure/Configure.h"
#include <iostream>

using namespace TMath;
using namespace std;

TMinuit minuit;
KATRIN Katrin;
Simulation sim;
TStopwatch t;

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

int main(int argc, char** argv)
{
	/* Input time distribution and threshold. Generate or read data. One can also change the time and threshold. */
	Time = Katrin.Time;
	Voltage = Katrin.Voltage;
	nvoltage = Katrin.Nbins;
	sim.initialize(argc, argv);
	sim.SetSlice(50);
	sim.SetMagnetic(Katrin.B_A, Katrin.B_S, Katrin.B_max);
	sim.SetupScatParameters(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.4e-18);

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
	Rate = Katrin.Rate;
	error = Katrin.Error;
	*/

	minuit.DefineParameter(0, "mass", pmass, 0.001, 0, 0);
	minuit.DefineParameter(1, "endpoint", pendpoint, 0.01, pendpoint-5, pendpoint+5);
	minuit.DefineParameter(2, "A_sig", 1, 1e-6, 0, 0);
	minuit.DefineParameter(3, "A_bkg", 1, 1e-6, 0, 0);
	//minuit.FixParameter(0);
	//minuit.FixParameter(1);
	//minuit.FixParameter(2);
	//minuit.FixParameter(3);

	/* Nuisance parameters. */
	minuit.DefineParameter(4, "B_A", Katrin.B_A, 1e-5*Katrin.B_A, 0, 0);
	minuit.DefineParameter(5, "B_S", Katrin.B_S, 1e-5*Katrin.B_S, 0, 0);
	minuit.DefineParameter(6, "B_max", Katrin.B_max, 1e-5*Katrin.B_max, 0, 0);
	//minuit.FixParameter(4);
	//minuit.FixParameter(5);
	minuit.FixParameter(6);

	minuit.DefineParameter(7, "A1", 0.204, 0.001, 0, 0);
	minuit.DefineParameter(8, "A2", 0.0556, 0.0001, 0, 0);
	minuit.DefineParameter(9, "w1", 1.85, 0.01, 0, 0);
	minuit.DefineParameter(10, "w2", 12.5, 0.01, 0, 0);
	minuit.DefineParameter(11, "e2", 14.3, 0.01,0, 0);
	minuit.DefineParameter(12, "cs", 3.4e2, 1e1, 0, 0);
	//minuit.FixParameter(7);
	//minuit.FixParameter(8);
	//minuit.FixParameter(9);
	//minuit.FixParameter(10);
	//minuit.FixParameter(11);
	//minuit.FixParameter(12);

	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	t.Start();
	minuit.Migrad();
	t.Stop();
	t.Print();
	return 0;

}

void fcn(int &npar, double* gin, double &f, double* par, int iflag) {
	delete[] theory;
	sim.SetupScatParameters(par[7], par[8], par[9], par[10], 12.6, par[11], par[12]*1e-20);
	sim.SetMagnetic(par[4], par[5], par[6]);
	sim.SetSlice(50);
	theory = sim.Asimov(par[0], par[1], nvoltage, Voltage, par[2], par[3]);
	double chi2 = 0;
	for(int i=0; i<nvoltage; i++) {
		double pred = theory[i];
		double meas = Rate[i];
		double Error = error[i];
		chi2 += pow((pred-meas)/Error, 2);
	}

	/* Punishing. */
	chi2 += pow((par[4]-Katrin.B_A)/Katrin.B_A_sigma /2, 2);
	chi2 += pow((par[5]-Katrin.B_S)/Katrin.B_S_sigma /2, 2);
	chi2 += pow((par[6]-Katrin.B_max)/Katrin.B_max_sigma /2, 2);
	chi2 += pow((par[7]-0.204)/0.001, 2);
	chi2 += pow((par[8]-0.0556)/0.0003, 2);
	chi2 += pow((par[9]-1.85)/0.02, 2);
	chi2 += pow((par[10]-12.5)/0.1, 2);
	chi2 += pow((par[11]-14.3)/0.02, 2);
	chi2 += pow((par[12]-3.4e2)/0.07e2, 2);
	f = chi2;
}

