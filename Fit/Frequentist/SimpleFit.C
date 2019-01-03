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
	cout << "OK here " << endl;

	/* Input time distribution and threshold. Generate or read data. One can also change the time and threshold. */
	Time = Katrin.Time;
	Voltage = Katrin.Voltage;
	nvoltage = Katrin.Nbins;
	sim.initialize(argc, argv);
	sim.SetSlice(50);
	sim.SetMagnetic(Katrin.B_A, Katrin.B_S, Katrin.B_max);
	cout << "Checkpoint1" << endl;

	/* This block for simulated sample. */
	sample = sim.Generate(pmass, pendpoint, nvoltage, Voltage, Time);
	cout << "Checkpoint2" << endl;
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

	/* Nuisance parameters. */
	minuit.DefineParameter(4, "B_A", Katrin.B_A, 1e-5*Katrin.B_A, 0, 0);
	minuit.DefineParameter(5, "B_S", Katrin.B_S, 1e-5*Katrin.B_S, 0, 0);
	minuit.DefineParameter(6, "B_max", Katrin.B_max, 1e-5*Katrin.B_max, 0, 0);
	//minuit.FixParameter(4);
	//minuit.FixParameter(5);
	minuit.FixParameter(6);

	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	t.Start();
	minuit.Migrad();
	t.Stop();
	t.Print();
	return 0;

}

void fcn(int &npar, double* gin, double &f, double* par, int iflag) {
	delete theory;
	sim.SetMagnetic(par[4], par[5], par[6]);
	theory = sim.Asimov(par[0], par[1], nvoltage, Voltage, par[2], par[3]);
	double chi2 = 0;
	for(int i=0; i<nvoltage; i++) {
		double pred = theory[i];
		double meas = Rate[i];
		double Error = error[i];
		chi2 += pow((pred-meas)/Error, 2);

		/* Punishing. */
		chi2 += pow((par[4]-Katrin.B_A)/Katrin.B_A_sigma /2, 2);
		chi2 += pow((par[5]-Katrin.B_S)/Katrin.B_S_sigma /2, 2);
		chi2 += pow((par[6]-Katrin.B_max)/Katrin.B_max_sigma /2, 2);
	}
	f = chi2;
}

