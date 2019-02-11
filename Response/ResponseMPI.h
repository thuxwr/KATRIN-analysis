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
#include "TMath.h"

using namespace TMath;

class ResponseMPI 
{
	public:
		ResponseMPI() {
			gtotal_response = new TGraph;
		}

		~ResponseMPI() {
			//delete gtotal_response;
		}

		void initialize(int argc, char** argv) { response.initialize(argc, argv); }
		bool SetupScatParameters(double A1, double A2, double w1, double w2, double e1, double e2, double InelasCS) {
			return response.SetupScatParameters(A1, A2, w1, w2, e1, e2, InelasCS);
		}
		void SetSlice(int iSlice) { response.SetSlice(iSlice); }

		void SetupResponse(double B_A, double B_S, double B_max) {
			/* Get the sum of total column density. */
			int nslice = response.GetNSlices();
			int size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			if(nslice!=size) {
				cout << "Number of cores is not in accordance with number of WGTS slices. Please use " << nslice << " cores to run again." << endl;
				exit(0);
			}
			double TotN = 0;
			for(int i=0; i<nslice; i++) TotN += response.GetColumnDensity(i) * response.GetTotalA(i);
			ScaleFactor = TotN;

			int slice;
			MPI_Comm_rank(MPI_COMM_WORLD, &slice);

			response.SetSlice(slice);
			response.SetupResponse(B_A, B_S, B_max);
		}

		void SetupEfficiency(double* efficiency) {
			delete gtotal_response;
			double* responsetotal = new double[response.NPoints];
			double* responselocal = new double[response.NPoints];
			for(int npoint=0; npoint<response.NPoints; npoint++) {
				responselocal[npoint] = 0;
				for(int pixel=0; pixel<NPixels; pixel++) {
					if(efficiency[pixel]<=0 || IsNaN(efficiency[pixel])) continue;
					responselocal[npoint] += response.response[pixel][npoint] * response.GetSegmentA(pixel) * efficiency[pixel];
				}
				responselocal[npoint] *= response.GetColumnDensity();
			}
			MPI_Reduce(responselocal, responsetotal, response.NPoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			/* Setup response function. */
			int slice;
			MPI_Comm_rank(MPI_COMM_WORLD, &slice);
			if(slice==0) {
				for(int npoint=0; npoint<response.NPoints; npoint++) responsetotal[npoint] /= ScaleFactor;
				gtotal_response = new TGraph(response.NPoints, response.X, responsetotal);
			}
		}

		double GetResponse(double E, double U) {
			if(E<U) return 0;
			else return gtotal_response->Eval(E-U);
		}

		TGraph* GetResponse() { return gtotal_response; }


	private:
		Response response;
		TGraph* gtotal_response;
		double ScaleFactor;

};

#endif
