/*
	 Read energy loss from TFile and return. Written in cpp.

	 Weiran, Nov.15, 2018.
*/

#ifndef EngLoss_h
#define EngLoss_h

#include <string>
using namespace std;

class EngLoss
{
	public:
		EngLoss() {
			char* KATRIN = getenv("KATRIN");
			if(KATRIN==0) {
				cout << "Environment variable 'KATRIN' is not defined." << endl;
				exit(0);
			}

			string path = (string)KATRIN + "/Scattering/energyloss/EnergyLoss.root";
			file = new TFile(path.c_str(), "READ");
			string name[3] = {"scat0", "scat1", "scat2"};
			for(int i=0; i<3; i++) scat[i] = (TH1D*)file->Get(name[i].c_str());
		}

		~EngLoss() {
			for(int i=0; i<3; i++) delete scat[i];
			file->Close();
			delete file;
		}

		double GetPdf(int ScatTimes, double epsilon) {
			if(ScatTimes<1) {
				cout << "Energy loss is only valid for Scatter times >= 1." << endl;
				return 0;
			}
			if(ScatTimes>3) {
				cout << "Scatter times greater than 3 is not supported." << endl;
				return 0;
			}
			int bin = scat[ScatTimes-1]->FindBin(epsilon);
			return scat[ScatTimes-1]->GetBinContent(bin);
		}
	
	private:
		TFile* file;
		TH1D* scat[3]; // We are only interested in those electrons inelastically scattered no more than 3 times.

};

#endif
