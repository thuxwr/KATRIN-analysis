/*
	 Construction of Feldman Cousins confidence interval for neutrino mass.
	 A non-physical extrapolation is adopted for negative neutrino mass.

	 Weiran, Nov.30, 2018.
*/

/* I have to include all ROOT headers BEFORE including Simulation.h and Configure.h for c++ compilation. */
#include "TMinuit.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "../../../Simulation/Simulation.h"
#include "../../../Configure/Configure.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#define NExp 10
Simulation sim;
KATRIN KATRIN;
TMinuit minuit;
double* Time;
double* Voltage;
double* sample;
double* Rate;
double* error;
double* theory;
int nvoltage;

void fcn(int &npar, double* gin, double &f, double* par, int iflag);
void ResetParam(double mass);

int main(int argc, char** argv) { //Input: mass, core#.
	/* Input mass and core number for parallel calculation. */
	stringstream ss1;
	ss1 << argv[1];
	double truems;
	ss1 >> truems;

	stringstream ss2;
	ss2 << argv[2];
	int core;
	ss2 >> core;

	/* Setup output path. */
	char corenum[5];
	sprintf(corenum, "%d", core);
	string Core = corenum;
	string filepath = "./result/core" + Core + ".dat";
	ofstream fout(filepath.c_str());
	fout << "# Confidence Interval construction for true value: mass = " << truems << " eV^2." << endl;
	fout << "# mass              log likelihood ratio" << endl;

	/* KATRIN setup. */
	Time = KATRIN.Time;
	Voltage = KATRIN.Voltage;
	nvoltage = KATRIN.Nbins;

	/* Minimizer setup. */
	minuit.SetPrintLevel(-1);
	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	
	for(int iExp=0; iExp<NExp; iExp++) {
		/* Simulate a spectrum given true neutrino mass. */
		delete[] sample;
		sample = sim.Generate(truems, 18574, nvoltage, Voltage, Time);

		delete[] Rate;
		Rate = new double[nvoltage];
		delete[] error;
		error = new double[nvoltage];
		for(int i=0; i<nvoltage; i++) {
			Rate[i] = sample[i]/Time[i];
			error[i] = sqrt(sample[i])/Time[i];
		}

		/* Call minimizer to get the best fit for neutrino mass. */
		ResetParam(truems);
		minuit.Migrad();
		double val, err;
		minuit.GetParameter(0, val, err);
		if(val<0) { // Check physical boundary.
			ResetParam(0);
			minuit.FixParameter(0);
			minuit.Migrad();
		}
		double fmin, fedm, errdef;
		int npari, nparx, istat;
		minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
		double min_global = fmin;

		/* Fix neutrino mass to true value and call minimizer again. */
		ResetParam(truems);
		minuit.FixParameter(0);
		minuit.Migrad();
		minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
		double min_true = fmin;
		double lR = min_true - min_global; // log likelihood ratio.

		fout << left << setw(20) << val << setw(20) << lR << endl;
	}

	fout.close();
}










void fcn(int &npar, double* gin, double &f, double* par, int iflag) {
	delete[] theory;
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

void ResetParam(double mass) {
	minuit.mncler();
	minuit.mnrset(1);
	minuit.DefineParameter(0, "mass", mass, 0.001, 0, 0);
	minuit.DefineParameter(1, "endpoint", 18574, 0.01, 0, 0);
	minuit.DefineParameter(2, "A_sig", 1, 1e-6, 0, 0);
	minuit.DefineParameter(3, "A_bkg", 1, 1e-6, 0, 0);
}








