#include "TFile.h"
#include "TGraph.h"
#include "TVirtualFFT.h"
#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "../Response.h"
#include "../../Configure/Configure.h"

using namespace std;

Response res;
KATRIN Katrin;

int main(int argc, char** argv)
{
	res.initialize(argc, argv);
	res.SetSlice(0);
	res.SetupScatParameters(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.4e-18);
	for(int i=0; i<10000; i++) {
		res.SetupResponse(Katrin.B_A, Katrin.B_S, Katrin.B_max * (1. + i * 1e-7));
		cout << "Finish: " << i << endl;
	}
	return 0;
}
	

