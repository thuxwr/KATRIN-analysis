/*
	 Get probability for electrons to scatter n times, given z and theta.

	 Weiran, Nov.15, 2018.
*/

#ifndef prob_h
#define prob_h

#include <string>

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

			tmpz = -1;
			Ncalib = 5e17;
			calib = dens->Integral();
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
				int low = dens->FindBin(z);
				int up = dens->GetNbinsX();
				tmpz = z;
				integ = dens->Integral(low, up);
			}
			return 1/cos(theta) * integ / calib * Ncalib;
		}
};

#endif
