/*
	 Generate MC sample and fit to model using Minuit.

	 Weiran, Nov.17, 2018.
*/

#include "../../Simulation/Simulation.h"
#include <iomanip>

Simulation sim;
TMinuit minuit;

double pmass = 0.35*0.35;
double pendpoint = 18574;

int nvoltage = 41;
double voltage[41] = {18545, 18546, 18547, 18548, 18549, 18550, 18551, 18552, 18553, 18554, 18555, 18556, 18557, 18558, 18559, 18560, 18561, 18562, 18563, 18564, 18565, 18566, 18567, 18567.5, 18568, 18568.5, 18569, 18569.5, 18570, 18570.5, 18571, 18571.5, 18572, 18573, 18574, 18575, 18576, 18577, 18578, 18579, 18580};
double entries[41] = {2.43646e+07, 2.16883e+07, 1.92214e+07, 1.69654e+07, 1.48958e+07, 1.30186e+07, 1.13122e+07, 9.7757e+06, 8.38887e+06, 7.15235e+06, 6.0483e+06, 5.06647e+06, 4.20482e+06, 3.4482e+06, 2.79557e+06, 2.23058e+06, 1.75045e+06, 1.34218e+06, 1.00594e+06, 731422, 511962, 343556, 187994, 145930, 110963, 82851, 60366, 42573, 64968, 282718, 344623, 25968, 36450, 32697, 32850, 32895, 32657, 32847, 32814, 32911, 32565};
double Time[41] = {1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,902625,902625,902625,902625,902625,902625,1.96321e+06,1.21403e+07,2.06024e+07,2.00834e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06};

double* theory;
double* sample = new double[41];
double* error = new double[41];
void fcn(int &npar, double* gin, double &f, double* par, int iflag);

void Fit2()
{
	/* Generate MC sample. */
	for(int i=0; i<41; i++) {
		sample[i] = entries[i]/Time[i] * 86400.;
		error[i] = sqrt(entries[i])/Time[i] * 86400.;
	}

	minuit.DefineParameter(0, "mass", 0, 0.0001, 0, 1e4);
	minuit.FixParameter(0);
	minuit.DefineParameter(1, "endpoint", 18575, 0.01, pendpoint-5, pendpoint+5);
	//minuit.FixParameter(1);
	minuit.DefineParameter(2, "A_sig", 1, 1e-6, 0, 2);
	minuit.DefineParameter(3, "A_bkg", 1, 1e-6, 0, 2);
	minuit.SetFCN(fcn);
	minuit.SetErrorDef(1);
	minuit.Migrad();
	double endpoint, enderr;
	minuit.GetParameter(1, endpoint, enderr);
	cout << setprecision(10) << "endpoint: " << endpoint << endl;

}

void fcn(int &npar, double* gin, double &f, double* par, int iflag) {
	delete theory;
	theory = sim.Asimov(par[0], par[1], nvoltage, voltage, par[2], par[3]);
	double chi2 = 0;
	for(int i=0; i<41; i++) {
		double pred = theory[i];
		double meas = sample[i];
		double err = error[i];
		chi2 += pow((pred-meas)/err, 2);
	}
	f = chi2;
}

