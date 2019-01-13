/*
	 Generate sample and print to file. The purpose is to provide an interface with python.

	 Weiran, Dec.3, 2018.

	 Add multi-slice structure.
	 Jan.11, 2019.
*/

#include "TH1D.h"
#include "TMath.h"
#include "TFile.h"
#include "Simulation.h"
#include "../Configure/Configure.h"
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace std;
KATRIN Katrin;
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

	sim.initialize(argc, argv);
	sim.SetupScatParameters(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.4e-18);
	sim.SetMagnetic(Katrin.B_A, Katrin.B_S, Katrin.B_max);
	double* sample = sim.Generate(mass, endpoint, Katrin.Nbins, Katrin.Voltage, Katrin.Time);

	char* KATRINpath = getenv("KATRIN");
	string path = KATRINpath;
	string filepath = path + "/Data/sample.dat";
	ofstream fout(filepath.c_str());
	fout << "# Voltage/V    Count               Time/s" << endl;
	for(int i=0; i<Katrin.Nbins; i++) {
		fout << left << setw(15) << Katrin.Voltage[i] << setw(20) << sample[i] << setw(20) << Katrin.Time[i] << endl;
	}

	fout.close();
	return 0;
}


