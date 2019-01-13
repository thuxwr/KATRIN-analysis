/*
	 Predictions for signal rate and background rate of KATRIN experiment.

	 MPI fitter.
	 Weiran, Jan.10, 2019.
*/

#include "../../../Detect/DetectMPI.h"
#include "../../../Configure/Configure.h"

DetectMPI detect;
KATRIN Katrin;

namespace KATRIN_namespace {

int nvoltage = Katrin.Nbins;
double* voltage = Katrin.Voltage;

using namespace std;

double bkg(ostream* pstream) { return Katrin.Bkg_rate; }

template <>
inline vector<double> signal(const vector<double>& pars, ostream* pstream) {
	double mass = pars.at(0);
	double endpoint = pars.at(1);
	double B_A = pars.at(2);
	double B_S = pars.at(3);
	double B_max = pars.at(4);
	double A1 = pars.at(5);
	double A2 = pars.at(6);
	double w1 = pars.at(7);
	double w2 = pars.at(8);
	double e1 = pars.at(9);
	double e2 = pars.at(10);
	double sigma = pars.at(11);
	//double B_max = Katrin.B_max;
	vector<double> entries = {};
	detect.SetScatParams(A1, A2, w1, w2, e1, e2, sigma);
	detect.SetMagnetic(B_A, B_S, B_max);
	double* detspec = detect.DetSpec(mass, endpoint, nvoltage, voltage);
	for(int bin=0; bin<nvoltage; bin++) {
		entries.push_back(detspec[bin]);
	}
	delete[] detspec;
	return entries;
}

/* Now we force the first derivative to be zero for those nuisance parameters. */
template <>
inline vector<var> signal(const vector<var>& pars, ostream* pstream) {
	vector<double> pars_d;
	for(int i=0; i<pars.size(); i++) {
		pars_d.push_back(pars.at(i).val());
	}

	vari** operands = ChainableStack::instance().memalloc_.alloc_array<vari*>(pars.size());
	for(int i=0; i<pars.size(); i++) {
		operands[i] = pars.at(i).vi_;
	}

	int npars = 2; // Number of parameters excluding nuisance pars.
	double precision[2] = {1e-5, 1e-5}; // Precision for calculating first derivative. Only for mass and endpoint.
	vector<double> entries = signal(pars_d, pstream);
	vector<double> D1_high[npars], D1_low[npars], D2_high[npars], D2_low[npars];
	for(int i=0; i<npars; i++) {
		vector<double> new_pars[4];
		for(int j=0; j<pars.size(); j++) {
			if(i==j) {
				new_pars[0].push_back(pars_d.at(j) + precision[i]);
				new_pars[1].push_back(pars_d.at(j) - precision[i]);
				new_pars[2].push_back(pars_d.at(j) + 0.5 * precision[i]);
				new_pars[3].push_back(pars_d.at(j) - 0.5 * precision[i]);
			}
			else {
				for(int k=0; k<4; k++) new_pars[k].push_back(pars_d.at(j));
			}
		}
		D1_high[i] = signal(new_pars[0], pstream);
		D1_low[i]  = signal(new_pars[1], pstream);
		D2_high[i] = signal(new_pars[2], pstream);
		D2_low[i]  = signal(new_pars[3], pstream);
	}

	vector<var> jacob = {};
	for(int bin=0; bin<nvoltage; bin++) {
		double* grad = ChainableStack::instance().memalloc_.alloc_array<double>(pars.size());
		for(int i=0; i<pars.size(); i++) {
			if(i<npars) {
				double D_1 = (D1_high[i].at(bin)-D1_low[i].at(bin)) /2. /precision[i];
				double D_2 = (D2_high[i].at(bin)-D2_low[i].at(bin)) /precision[i];
				grad[i] = (4. * D_2 - D_1) / 3.;
			}
			else grad[i] = 0;
		}
		double value = entries.at(bin);
		jacob.push_back(var(new precomputed_gradients_vari(value, pars.size(), operands, grad)));
	}
	return jacob;
}

}
