/*
	 A global configure file to setup all physical constants and detector parameters for KATRIN.

	 Weiran, Nov.16, 2018.
*/

#ifndef Configure_h
#define Configure_h

#define Nfinalstate 52 //Discretized final state number
#define NbinsDecaySpec 800 //Nbins for decay spectrum
#define LowBoundary -35 //Low boundary for decay spectrum respect to endpoint, in unit:eV
#define UpBoundary 5 //Up boundary for decay spectrum respect to endpoint, in unit:eV

#define NVoltageMax 70 //Nbins for detected spectrum
#define DetLow -30 //Low boundary for detected spectrum. Should be larger than boundary for decay spectrum.
#define DetUp 5 //Up boundary for detected spectrum.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
using namespace TMath;

class KATRIN
{
	public:
		KATRIN() {
			char* KATRIN_path = getenv("KATRIN");
			if(KATRIN_path==0) {
				cout << "Environment variable 'KATRIN' is not defined." << endl;
				exit(0);
			}

			ifstream data(((string)KATRIN_path + "/Data/data.dat").c_str());
			if(!(data.is_open())) {
				cout << "Data file $KATRIN/Data/data.dat cannot be opened." << endl;
				exit(0);
			}
			string line;
			Nbins = 0;
			while(getline(data, line)) {
				char f = line[0];
				if(f=='#') continue;
				if(Nbins>NVoltageMax) {
					cout << "Data size exceed maximum." << endl;
					exit(0);
				}
				istringstream iss(line);
				iss >> Voltage[Nbins] >> Count[Nbins] >> Time[Nbins];
				Nbins ++;
			}

			for(int i=0; i<Nbins; i++) {
				Rate[i] = Count[i]/Time[i];
				Error[i] = sqrt(Count[i])/Time[i];
			}
		}

		~KATRIN(){}

		/* Electromagnetic field. */
		double B_A = 2.68; // in unit:Gauss
		double B_S = 2.52e4; // in unit:Gauss
		double B_max = 4.2e4; // in unit:Gauss

		/* Source. */
		double epsilon_T = 0.95; // T2 purity
		double T_bt = 30; // temperature in WGTS, in unit:K
		double bv = 13; //weighed mean bulk velocity, in unit:m/s.

		/* Spectrum. */
		double E_0_center = 18.574e3;

		/* Measurement time distribution, in unit:day */
		/*
		double Time[NVoltageMax] = {
			5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
			5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
			5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 130, 130, 130, 130, 20, 20, 20, 20, 20, 20,
			20, 20, 20, 20, 20, 20, 20, 20, 20, 20
			};
		*/

		/* Data. */
		int Nbins;
		double Time[NVoltageMax];
		double Voltage[NVoltageMax];
		double Count[NVoltageMax];

		/* Transformed data. */
		double Rate[NVoltageMax];
		double Error[NVoltageMax];

		/* Expected signal rate and background rate, for Monte Carlo generation. */
		double Sig_rate = 27.54; // For U=18544.25V, 27.54Hz
		double Bkg_rate = 10e-3; // 10mHz

};

namespace Physics {
	double m_e = 0.511e6; // electron mass, 0.511MeV
	double alpha = 0.00729735257; // fine structure constant
	double pi = TMath::Pi(); 
	double k_B = 8.6173303e-5; // Boltzmann constant, in unit:eV/K
	double c = 2.99792458e8; // velocity of light, in unit:m/s
	double M_T2 = 6.032 * 931.4940954e6; // T2 molecular mass, in unit:eV
}

#endif
