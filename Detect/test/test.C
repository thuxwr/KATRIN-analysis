#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "../Detect.h"
#include "../../Configure/Configure.h"

using namespace std;

int main(int argc, char** argv) {
	Detect det;
	KATRIN Katrin;
	det.initialize(argc, argv);
	det.SetScatParams(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.4e-18);
	//det.SetScatParams(0,0,0,0,0,0,0);
	det.SetSlice(50);
	det.SetMagnetic(Katrin.B_A, Katrin.B_S, Katrin.B_max);
	int nvoltage = 1;
	double* voltage = new double[1];
	voltage[0] = 18554;
	double* detspec = det.DetSpec(0, 18574, nvoltage, voltage);
	cout << detspec[0] << endl;
	//det.SetScatParams(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 1);
	det.SetScatParams(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.8e-18);
	det.SetSlice(50);
	det.SetMagnetic(Katrin.B_A, Katrin.B_S, Katrin.B_max);
	delete[] detspec;
	detspec = det.DetSpec(0, 18574, nvoltage, voltage);
	cout << detspec[0] << endl;
	return 0;
}


