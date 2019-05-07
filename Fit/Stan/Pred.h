/*
	 Predictions for signal rate and background rate of KATRIN experiment.

	 Weiran, Nov.18, 2018.
*/

#include "../../Detect/Detect.h"
#include "../../Configure/Configure.h"

Detect detect;
KATRIN Katrin;

int nvoltage = Katrin.Nbins;
double* voltage = Katrin.Voltage;

using namespace std;

namespace KATRIN_namespace {

double bkg(ostream* pstream) { return Katrin.Bkg_rate; }

template <>
inline vector<double> signal(const vector<double>& pars, ostream* pstream) {
	double mass = pars.at(0);
	double endpoint = pars.at(1);
	double B_A = 6.32;
	double B_S = 3.14965e4;
	double B_max = 4.2332e4;
	vector<double> entries = {};
	detect.SetScatParams(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.4e-18);
	detect.SetMagnetic(B_A, B_S, B_max);
	double* detspec = detect.DetSpec(mass, endpoint, nvoltage, voltage);
	for(int bin=0; bin<nvoltage; bin++) {
		entries.push_back(detspec[bin]);
	}
	delete[] detspec;
	return entries;
}

template <>
inline vector<stan::math::var> signal(const vector<stan::math::var>& pars, ostream* pstream) {
	vector<double> pars_d;
	for(int i=0; i<pars.size(); i++) {
		pars_d.push_back(pars.at(i).val());
	}

	vari** operands = ChainableStack::instance().memalloc_.alloc_array<vari*>(pars.size());
	for(int i=0; i<pars.size(); i++) {
		operands[i] = pars.at(i).vi_;
	}

	double precision[2] = {1e-5, 1e-5}; // precision for calculating first derivative
	vector<double> entries = signal(pars_d, pstream);
	vector<double> D1_high[pars.size()], D1_low[pars.size()], D2_high[pars.size()], D2_low[pars.size()];
	for(int i=0; i<pars.size(); i++) {
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

	vector<stan::math::var> jacob = {};
	for(int bin=0; bin<nvoltage; bin++) {
		double* grad = ChainableStack::instance().memalloc_.alloc_array<double>(pars.size());
		for(int i=0; i<pars.size(); i++) {
			double D_1 = (D1_high[i].at(bin)-D1_low[i].at(bin)) /2. /precision[i];
			double D_2 = (D2_high[i].at(bin)-D2_low[i].at(bin)) /precision[i];
			grad[i] = (4. * D_2 - D_1) / 3.;
		}
		double value = entries.at(bin);
		jacob.push_back(stan::math::var(new precomputed_gradients_vari(value, pars.size(), operands, grad)));
	}
	return jacob;
}
}
