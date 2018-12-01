/*
	 Get detector response function for given z.

	 Currently only z=0, 0.5, 1, 1.5, ..., 8 meters are calculated. 
	 z=0 is used in this toy model.

	 Weiran, Nov.15, 2018.
*/

#ifndef Response_h
#define Response_h

#include <string>
#include <iostream>

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
				file[i] = new TFile((path+"Average.root").c_str(), "READ");
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
			if(E>=response[0]->GetXaxis()->GetXmax() || U>=response[0]->GetYaxis()->GetXmax()) return 0;
			else if(E<=response[0]->GetXaxis()->GetXmin() || U<=response[0]->GetYaxis()->GetXmin()) return 0;
			else return response[0]->Interpolate(E, U);
		}


	private:
		TFile* file[Nslices];
		TH2D* response[Nslices];
		string GetSourceFile(double z); // Left blank for further usage.

};

#endif
