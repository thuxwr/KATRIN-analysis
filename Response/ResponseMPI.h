/*
	 Use MPI multiprocessing to calculate detector response.

	 Weiran, Jan.8, 2019.
*/

#ifndef ResponseMPI_h
#define ResponseMPI_h

#include "Response.h"
#include "mpi.h"
#include "TGraph.h"
#include <unistd.h>
#include <iomanip>

class ResponseMPI 
{
	public:
		ResponseMPI() {
			gtotal_response = new TGraph;
		}

		~ResponseMPI() {
			delete gtotal_response;
		}

		void initialize(int argc, char** argv) { response.initialize(argc, argv); }
		bool SetupScatParameters(double A1, double A2, double w1, double w2, double e1, double e2, double InelasCS) {
			return response.SetupScatParameters(A1, A2, w1, w2, e1, e2, InelasCS);
		}
		void SetSlice(int iSlice) { response.SetSlice(iSlice); }

		void SetupResponse(double B_A, double B_S, double B_max) {
			delete gtotal_response;
			/* Get the sum of total column density. */
			int nslice = response.GetNSlices();
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			if(nslice!=size) {
				cout << "Number of cores is not in accordance with number of WGTS slices. Please use " << nslice << " cores to run again." << endl;
				exit(0);
			}
			double TotDens = 0;
			for(int i=0; i<nslice; i++) TotDens += response.GetColumnDensity(i);

			int slice;
			MPI_Comm_rank(MPI_COMM_WORLD, &slice);
			double density = response.GetColumnDensity(slice);
			double weight = density/TotDens;

			response.SetSlice(slice);
			response.SetupResponse(B_A, B_S, B_max);

			TGraph* slice_response = response.GetResponse();
			int npoints = slice_response->GetN();
			double* local_response = new double[npoints];
			double* total_response = new double[npoints];
			double* yvalue = slice_response->GetY();
			for(int i=0; i<npoints; i++) local_response[i] = yvalue[i] * weight; 

			/* Send to first core. */
			if(slice==0) {
				for(int i=0; i<npoints; i++) total_response[i] = local_response[i];
				for(int source=1; source<nslice; source++) {
					MPI_Recv(local_response, npoints, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					usleep(100000);
					for(int i=0; i<npoints; i++) total_response[i] += local_response[i];
				}
				double* xvalue = slice_response->GetX();
				gtotal_response = new TGraph(npoints, xvalue, total_response);
			}
			else {
				MPI_Send(local_response, npoints, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
				gtotal_response = new TGraph;
			}

			delete[] local_response; delete[] total_response;
		}

		double GetResponse(double E, double U) {
			if(E<U) return 0;
			else return gtotal_response->Eval(E-U);
		}

		TGraph* GetResponse() { return gtotal_response; }


	private:
		Response response;
		TGraph* gtotal_response;

};

#endif
