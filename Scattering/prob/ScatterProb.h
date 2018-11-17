/*
	 Get probability for electrons to scatter n times, given z and theta.

	 Weiran, Nov.15, 2018.
*/

#ifndef prob_h
#define prob_h

#include <string>
#include <iostream>

using namespace std;

class ScatterProb
{
	public:
		ScatterProb() {
			char* KATRIN = getenv("KATRIN");
			if(KATRIN==0) {
				cout << "Environment variable 'KATRIN' is not defined." << endl;
				exit(0);
			}

			string path = (string)KATRIN + "/Scattering/density/density.root";
			file = new TFile(path.c_str(), "READ");
			dens = (TH1D*)file->Get("density");

			tmpz = -10;
			Ncalib = 5e17;
			int low = 0;
			int up = dens->FindBin(5); //Width of WGTS is 10m.
			calib = dens->Integral(low, up);
			sigma = 3.4e-18;
		}

		~ScatterProb() {
			delete dens;
			file->Close();
			delete file;
		}

		double GetProb(int s, double z, double theta) {
			double lambda = Neff(z, theta) * sigma;
			return pow(lambda, s) / TMath::Factorial(s) * exp(-1*lambda);
		}


	private:
		TFile* file;
		TH1D* dens;

		double tmpz;
		double integ;
		double calib;
		double Ncalib;
		double sigma; //Elastic: 0.29e-18; Inelastic: 3.4e-18.

		double Neff(double z, double theta) {
			if(z!=tmpz) {
				int low = dens->FindBin(abs(z));
				int up = dens->FindBin(5);
				tmpz = z;
				integ = dens->Integral(low, up);
			}
			if(z>=0) return 1/cos(theta) * integ / calib * Ncalib / 2;
			else return 1/cos(theta) * (2 - integ/calib) * Ncalib / 2;
		}
};

#endif
