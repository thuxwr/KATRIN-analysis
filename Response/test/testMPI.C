#include "TFile.h"
#include "TGraph.h"
#include "TVirtualFFT.h"
#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "../ResponseMPI.h"
#include "../../Configure/Configure.h"
#include <unistd.h>
#include <iostream>

using namespace std;


int main(int argc, char** argv)
{
ResponseMPI res;
KATRIN Katrin;
	int rank;
	int time=0;
	int flag;
	int size;
	MPI_Init(&argc, &argv);
	res.initialize(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	while(true) {
		if(rank==0) 
			for(int i=1; i<size; i++) MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		else
			MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		res.SetupScatParameters(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3, 3.4e-18);
		res.SetupResponse(Katrin.B_A, Katrin.B_S, Katrin.B_max);
		if(rank==0) {
			cout << "Finish time: " << time << endl;
			time += 1;
		}
	}
	MPI_Finalize();

	//for(int i=0; i<10000; i++) {
//		res.SetupResponse(Katrin.B_A, Katrin.B_S, Katrin.B_max * (1. + i * 1e-7));
		//cout << "Finish: " << i << endl;
	//}
	return 0;
}
	

