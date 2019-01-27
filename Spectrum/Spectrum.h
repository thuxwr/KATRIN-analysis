/*
	 Theoretical calculation for beta decay spectrum. 
	 //Doppler broadening is calculated WITHOUT using RooFFTConvPdf, but by directly setting values for each bin.

	 Final state distribution should be modified to incorporate T2 purity.

	 Weiran, Nov.16, 2018.
*/

#ifndef Spectrum_h
#define Spectrum_h

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "../Configure/Configure.h"
#include "TH1D.h"

using namespace std;
using namespace Physics;

class Spectrum
{
	public:
		Spectrum() {
			/* Import final state distribution. */
			char* KATRINpath = getenv("KATRIN");
			if(KATRINpath==0) {
				cout << "Environment variable 'KATRIN' is not defined." << endl;
				exit(0);
			}

			string path = (string)KATRINpath + "/data/SSC/";

			/* T2 data. */
			{
			ifstream data((path + "FSD_DOSS_T2_rebinned.txt").c_str());
			if(!(data.is_open())) {
				cout << "Data file $KATRIN/data/SSC/FSD_DOSS_T2_rebinned.txt cannot be opened." << endl;
				exit(0);
			}

			string line;
			int nfinalstate = 0;
			while(getline(data, line)) {
				//char f = line[0];
				//if(f=='#') continue;
				istringstream iss(line);
				iss >> fsenergyT2[nfinalstate] >> fsprobT2[nfinalstate];
				nfinalstate += 1;
			}
			nfsT2 = nfinalstate;
			}

			/* DT data. */
			{
			ifstream data((path + "FSD_DOSS_DT_rebinned.txt").c_str());
			if(!(data.is_open())) {
				cout << "Data file $KATRIN/data/SSC/FSD_DOSS_DT_rebinned.txt cannot be opened." << endl;
				exit(0);
			}

			string line;
			int nfinalstate = 0;
			while(getline(data, line)) {
				//char f = line[0];
				//if(f=='#') continue;
				istringstream iss(line);
				iss >> fsenergyDT[nfinalstate] >> fsprobDT[nfinalstate];
				nfinalstate += 1;
			}
			nfsDT = nfinalstate;
			}

			/* HT data. */
			{
			ifstream data((path + "FSD_Saenz_HTmod.txt").c_str());
			if(!(data.is_open())) {
				cout << "Data file $KATRIN/data/SSC/FSD_Saenz_HTmod.txt cannot be opened." << endl;
				exit(0);
			}

			string line;
			int nfinalstate = 0;
			while(getline(data, line)) {
				//char f = line[0];
				//if(f=='#') continue;
				istringstream iss(line);
				iss >> fsenergyHT[nfinalstate] >> fsprobHT[nfinalstate];
				nfinalstate += 1;
			}
			nfsHT = nfinalstate;
			}
		}
		~Spectrum(){}

		TH1D* decayspec(double ms, double E_0) {
			TH1D* spec = new TH1D("", "", NbinsDecaySpec, LowBoundary+katrin.E_0_center, UpBoundary+katrin.E_0_center);
			for(int i=1; i<=NbinsDecaySpec; i++) {
				double energy = spec->GetBinCenter(i);
				spec->SetBinContent(i, shape_fcn(energy, ms, E_0));
			}
			spec->Scale(1/spec->Integral());
			return spec;
		}

		double GetBinWidth() {
			return (double)(UpBoundary-LowBoundary)/NbinsDecaySpec;
		}

		double GetBinCenter(int bin) {
			if(bin>NbinsDecaySpec || bin<=0) {
				cout << "Bin number out of range." << endl;
				return 0;
			}
			return LowBoundary + (bin-0.5) * GetBinWidth() + katrin.E_0_center;
		}


	private:
		double fsenergyT2[NfinalstateMAX]; //Final state energy for T2 
		double fsprobT2[NfinalstateMAX]; //Final state prob for T2
		double fsenergyDT[NfinalstateMAX];
		double fsprobDT[NfinalstateMAX];
		double fsenergyHT[NfinalstateMAX];
		double fsprobHT[NfinalstateMAX];
		int nfsT2, nfsDT, nfsHT;
		KATRIN katrin;

		/* Some useful functions. */
		double gamma(double energy) {
			return energy/m_e + 1;
		}

		double beta(double energy) {
			return sqrt(1. - 1. /gamma(energy) /gamma(energy));
		}

		double Fermi_fcn(double energy, int Z) {
			double eta = alpha * Z / beta(energy);
			return 2 * pi * eta / (1-exp(-2 * pi * eta));
		}

		double momentum(double energy) {
			return m_e * beta(energy);
		}

		double epsilon(double energy, double V_f, double E_0) {
			return E_0 - energy - V_f;
		}

		double w(double energy) {
			return (energy + m_e)/m_e;
		}

		double t(double energy) {
			double Beta = beta(energy);
			return 1./Beta * 0.5 * log((1+Beta)/(1-Beta)) - 1;
		}

		double radiative_correction(double energy, double V_f, double E_0) {
			double W = w(energy);
			double T = t(energy);
			double Beta = beta(energy);
			double W_0 = w(E_0 - V_f);
			if(W_0<=W) return 0;
			double value = pow(W_0-W, 2*alpha/pi * T) * (1 + 2. * alpha /pi * (
				T * (log(2.) - 1.5 + (W_0-W)/W) + 0.25 * (T+1.) * (2.*(1.+Beta*Beta) + 2.*log(1-Beta) + pow(W_0-W,2)/(6.*W*W))
				- 2 + 0.5 * Beta - 17./36. * Beta * Beta + 5./6. * pow(Beta,3) ));
			return value;
		}

		/* Shape. */
		double shape_fcn(double x, double ms, double E_0) {
			double shape = 0;
			/* T2 */
			for(int i=0; i<nfsT2; i++) {
				double Epsilon = epsilon(x, fsenergyT2[i], E_0);
				if(Epsilon<=0) continue; // cannot reach this final state
				if(pow(Epsilon, 2) <= ms) continue; // cannot reach this final state
				if(ms>=0)
					shape += fsprobT2[i] * Epsilon * sqrt(pow(Epsilon,2) - ms) * radiative_correction(x, fsenergyT2[i], E_0) * katrin.T2concentration;
				else {//Non-physical extrapolation.
					double mu = 0.72 * sqrt(-1 * ms);
					shape += fsprobT2[i] * (Epsilon + mu*exp(-1 - Epsilon/mu)) * sqrt(pow(Epsilon,2) - ms) * radiative_correction(x, fsenergyT2[i], E_0) * katrin.T2concentration;
				}
			}

			/* DT */
			for(int i=0; i<nfsDT; i++) {
				double Epsilon = epsilon(x, fsenergyDT[i], E_0);
				if(Epsilon<=0) continue; // cannot reach this final state
				if(pow(Epsilon, 2) <= ms) continue; // cannot reach this final state
				if(ms>=0)
					shape += fsprobDT[i] * Epsilon * sqrt(pow(Epsilon,2) - ms) * radiative_correction(x, fsenergyDT[i], E_0) * katrin.DTconcentration;
				else {//Non-physical extrapolation.
					double mu = 0.72 * sqrt(-1 * ms);
					shape += fsprobDT[i] * (Epsilon + mu*exp(-1 - Epsilon/mu)) * sqrt(pow(Epsilon,2) - ms) * radiative_correction(x, fsenergyDT[i], E_0) * katrin.DTconcentration;
				}
			}

			/* HT */
			for(int i=0; i<nfsHT; i++) {
				double Epsilon = epsilon(x, fsenergyHT[i], E_0);
				if(Epsilon<=0) continue; // cannot reach this final state
				if(pow(Epsilon, 2) <= ms) continue; // cannot reach this final state
				if(ms>=0)
					shape += fsprobHT[i] * Epsilon * sqrt(pow(Epsilon,2) - ms) * radiative_correction(x, fsenergyHT[i], E_0) * katrin.HTconcentration;
				else {//Non-physical extrapolation.
					double mu = 0.72 * sqrt(-1 * ms);
					shape += fsprobHT[i] * (Epsilon + mu*exp(-1 - Epsilon/mu)) * sqrt(pow(Epsilon,2) - ms) * radiative_correction(x, fsenergyHT[i], E_0) * katrin.HTconcentration;
				}
			}
			return shape;
		}

};

#endif
