/*
	 Simulate the detected spectrum, with y axis being decay rate.
	 Poisson distribution in each bin is approximated to be a gaussian. (Minimum entries being >10000)

	 Weiran, Nov.17, 2018.
*/

#ifndef Simulation_h
#define Simulation_h

#include <iostream>
#include "../Detect/Detect.h"
#include "../Configure/Configure.h"
#include "../Spectrum/Spectrum.h"

using namespace std;

class Simulation
{
	public:
		Simulation() {
			rndm = new TRandom3(0);
		}

		~Simulation() {
		}

		/* B_A, B_S and B_max is not allowed to change in MC generation. */
		/* Reason: Response function should be recalculated if magnetic field changes. */
		TH1D* Generate(double mass, double endpoint, double A_sig=1, double A_bkg=1, double z=0) {
			TH1D* sim = new TH1D("", "", NbinsDetSpec, DetLow+KATRIN.E_0_center, DetUp+KATRIN.E_0_center);
			TH1D* theory = detect.DetSpec(spec.decayspec(mass, endpoint), z);
			for(int bin=1; bin<=NbinsDetSpec; bin++) {
				double signal = A_sig * theory->GetBinContent(bin) * KATRIN.Time[bin-1];
				double bkg = A_bkg * KATRIN.Bkg_rate * KATRIN.Time[bin-1];
				double entries = rndm->PoissonD(signal+bkg);
				sim->SetBinContent(bin, entries/KATRIN.Time[bin-1]);
				sim->SetBinError(bin, sqrt(entries)/KATRIN.Time[bin-1]);
			}
			delete theory;
			return sim;
		}

		TH1D* Asimov(double mass, double endpoint, double A_sig=1, double A_bkg=1, double z=0) {
			TH1D* asimov = new TH1D("", "", NbinsDetSpec, DetLow+KATRIN.E_0_center, DetUp+KATRIN.E_0_center);
			TH1D* theory = detect.DetSpec(spec.decayspec(mass, endpoint), z);
			for(int bin=1; bin<=NbinsDetSpec; bin++) {
				double signal = A_sig * theory->GetBinContent(bin) * KATRIN.Time[bin-1];
				double bkg = A_bkg * KATRIN.Bkg_rate * KATRIN.Time[bin-1];
				double entries = signal+bkg;
				asimov->SetBinContent(bin, entries/KATRIN.Time[bin-1]);
				asimov->SetBinError(bin, sqrt(entries)/KATRIN.Time[bin-1]);
			}
			delete theory;
			return asimov;
		}

		
	private:
		KATRIN KATRIN;
		Detect detect;
		Spectrum spec;
		TRandom3* rndm;



};

#endif
