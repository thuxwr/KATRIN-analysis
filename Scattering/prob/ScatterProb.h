/*
	 Get probability for electrons to scatter n times, given z and theta.

	 Weiran, Nov.15, 2018.
*/

#ifndef prob_h
#define prob_h

#include <string>
#include <iostream>

#define CosThetaStep 0.0001

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
			for(int i=0; i<4; i++) delete CumulatePdf[i];
			file->Close();
			delete file;
		}

		double GetProb(int s, double z, double cos_theta) {
			double lambda = Neff(z, cos_theta) * sigma;
			return pow(lambda, s) / TMath::Factorial(s) * exp(-1*lambda);
		}

		/* Cumulated probability, integrated over theta. */
		void SetProbCumulate(double z) {
			_z = z;
			for(int i=0; i<4; i++) delete CumulatePdf[i];
			for(int s=0; s<4; s++) {
				int nbins = (int)(1./CosThetaStep);
				CumulatePdf[s] = new TH1D("", "", nbins, 0, 1);
				double cumulate = 0;
				double binwidth = 1./nbins;
				for(int bin=nbins; bin>0; bin--) {
					double bincenter = CumulatePdf[s]->GetBinCenter(bin);
					cumulate += binwidth * GetProb(s, _z, bincenter);
					CumulatePdf[s]->SetBinContent(bin, cumulate);
				}
			}
		}

		double GetProbCumulate(int s, double cosmax) { // for s<=3.
			if(s>=4) {
				cout << "Scatter times greater than three. Not supported." << endl;
				return 0;
			}
			return CumulatePdf[s]->GetBinContent(CumulatePdf[s]->FindBin(cosmax));
		}

		void SetInelasCrossSection(double cs) { 
			if(sigma==cs) return;
			sigma = cs;
			SetProbCumulate(_z);
		} // in unit: cm^2

	private:
		TFile* file;
		TH1D* dens;
		TH1D* CumulatePdf[4]; // Pdf with variable epsilon, integrated over theta.

		double tmpz;
		double integ;
		double calib;
		double Ncalib;
		double sigma; //Elastic: 0.29e-18; Inelastic: 3.4e-18.
		double _z;

		double Neff(double z, double cos_theta) {
			if(z!=tmpz) {
				int low = dens->FindBin(TMath::Abs(z));
				int up = dens->FindBin(5);
				tmpz = z;
				integ = dens->Integral(low, up);
			}
			if(z>=0) return 1/cos_theta * integ / calib * Ncalib / 2;
			else return 1/cos_theta * (2 - integ/calib) * Ncalib / 2;
		}

};

#endif
