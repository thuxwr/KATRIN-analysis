/*
	 Get detector response function for given z.

	 Currently only z=0, 0.5, 1, 1.5, ..., 8 meters are calculated. 
	 z=0 is used in this toy model.

	 Weiran, Nov.15, 2018.
*/

#ifndef Response_h
#define Response_h

#include <string>

#define Nslices 1 // number of source slices
using namespace std;

class Response
{
	public:
		Response() {
			char* KATRIN = getenv("KATRIN");
			if(KATRIN==0) {
				cout << "Environment variable 'KATRIN' is not defined." << endl;
				exit(0);
			}

			string path = (string)KATRIN + "/Response/data/";
			for(int i=0; i<Nslices; i++) {
				file[i] = new TFile((path+"Core0.root").c_str(), "READ");
				response[i] = (TH2D*)file[i]->Get("response");
			}
		}

		~Response() {
			for(int i=0; i<Nslices; i++) {
				delete response[i];
				file[i]->Close();
				delete file[i];
			}
		}

		double GetResponse(double E, double U, double z=0) {
			int xbin = response[0]->GetXaxis()->FindBin(E);
			int ybin = response[0]->GetYaxis()->FindBin(U);
			return response[0]->GetBinContent(xbin, ybin);
		}


	private:
		TFile* file[Nslices];
		TH2D* response[Nslices];
		string GetSourceFile(double z); // Left blank for further usage.

};

#endif
