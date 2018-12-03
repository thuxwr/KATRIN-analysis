/*
	 Generate sample and print to file. The purpose is to provide an interface with python.

	 Weiran, Dec.3, 2018.
*/

#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "Simulation.h"
#include "../Configure/Configure.h"
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace std;
KATRIN katrin;
Simulation sim;

int main(int argc, char** argv) {
	stringstream ss1;
	ss1 << argv[1];
	double mass;
	ss1 >> mass;
	stringstream ss2;
	ss2 << argv[2];
	double endpoint;
	ss2 >> endpoint;

	sim.SetMagnetic(katrin.B_A, katrin.B_S, katrin.B_max);
	double* sample = sim.Generate(mass, endpoint, katrin.Nbins, katrin.Voltage, katrin.Time);

	char* KATRINpath = getenv("KATRIN");
	string path = KATRINpath;
	string filepath = path + "/Data/sample.dat";
	ofstream fout(filepath.c_str());
	fout << "# Voltage/V    Count               Time/s" << endl;
	for(int i=0; i<katrin.Nbins; i++) {
		fout << left << setw(15) << katrin.Voltage[i] << setw(20) << sample[i] << setw(20) << katrin.Time[i] << endl;
	}

	fout.close();
	return 0;
}


