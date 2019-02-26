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
			response_gather = new double[1];
			response_pixel = new double**[NVoltage];
			for(int nvoltage=0; nvoltage<NVoltage; nvoltage++) {
				response_pixel[nvoltage] = new double*[NPixels];
				for(int pixel=0; pixel<NPixels; pixel++) response_pixel[nvoltage][pixel] = new double[1];
			}
			responsetotal = new double[1];
		}

		~ResponseMPI() {
		}

		void initialize(int argc, char** argv) { response.initialize(argc, argv); }
		bool SetupScatParameters(double A1, double A2, double w1, double w2, double e1, double e2, double InelasCS) {
			return response.SetupScatParameters(A1, A2, w1, w2, e1, e2, InelasCS);
		}
		void SetSlice(int iSlice) { response.SetSlice(iSlice); }

		void SetupResponse(double B_A, double B_S, double B_max) {
			delete response_gather;
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

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);

			response.SetSlice(rank);
			response.SetupResponse(B_A, B_S, B_max);
			
			/* Send calculated response to first core. */
			double* local_response = new double[response.NPoints * NPixels * NVoltage];
			for(int nvoltage=0; nvoltage<NVoltage; nvoltage++) for(int pixel=0; pixel<NPixels; pixel++) for(int npoint=0; npoint<response.NPoints; npoint++) 
				local_response[nvoltage*NPixels*response.NPoints+pixel*response.NPoints+npoint] = response.response[nvoltage][pixel][npoint];

			response_gather = new double[size*response.NPoints*NPixels*NVoltage];
			MPI_Allgather(local_response, response.NPoints*NPixels*NVoltage, MPI_DOUBLE, response_gather, response.NPoints*NPixels*NVoltage, MPI_DOUBLE, MPI_COMM_WORLD);

			for(int nvoltage=0; nvoltage<NVoltage; nvoltage++) for(int pixel=0; pixel<NPixels; pixel++) {
				delete response_pixel[nvoltage][pixel];
				response_pixel[nvoltage][pixel] = new double[response.NPoints];
				for(int npoint=0; npoint<response.NPoints; npoint++) {
					response_pixel[nvoltage][pixel][npoint] = 0;
					for(int slice=0; slice<nslice; slice++)
						response_pixel[nvoltage][pixel][npoint] += response_gather[slice*response.NPoints*NPixels*NVoltage+nvoltage*NPixels*response.NPoints+pixel*response.NPoints+npoint] * response.GetSegmentA(pixel, slice) * response.GetColumnDensity(slice);
				}
			}
		}

		void SetupEfficiency(double* efficiency, double hv) { // Switch according to high voltage.
			delete responsetotal;
			responsetotal = new double[response.NPoints];

			int thishv = 0;
			bool IsFound = false;
			/* Allow 1eV shift. */
			for(int nhv=0; nhv<NVoltage; nhv++) {
				if(Abs(response.HV[nhv]-hv)<1) {
					thishv = nhv;
					IsFound = true;
					break;
				}
			}
			if(hv>18580) {
				for(int npoint=0; npoint<response.NPoints; npoint++) responsetotal[npoint] = 0;
				return;
			}
			if(!IsFound) {
				cout << "Response model for voltage " << hv << "eV is not calculated." << endl;
				exit(0);
			}

			for(int npoint=0; npoint<response.NPoints; npoint++) {
				responsetotal[npoint] = 0;
				for(int pixel=0; pixel<NPixels; pixel++) {
					if(efficiency[pixel]<0) continue;
					responsetotal[npoint] += response_pixel[thishv][pixel][npoint] * efficiency[pixel];
				}
				responsetotal[npoint] /= ScaleFactor;
			}
		}

		double GetResponse(double E, double U) {
			if(E<=U) return 0;
			else {
				/* Interpolate from responsetotal. */
				int low = Binary(response.NPoints, response.X, E-U);
				if(low<0) {
					cout << "Response not in range!" << endl;
					exit(0);
				}
				return (responsetotal[low+1]*(E-U-response.X[low])+responsetotal[low]*(response.X[low+1]-E+U)) / (response.X[low+1]-response.X[low]);
			}
			return -1;
		}

		//TGraph* GetResponse() { return gtotal_response; }


	private:
		Response response;
		double ScaleFactor;
		double* response_gather;
		double*** response_pixel;
		double* responsetotal;
		double weight[NPixels];

		int Binary(int total, double* array, double value) { //Binary search for the position in a sorted array.
			int low = 0;
			int up = total-1;
			if(array[up]<value || array[low]>value) return -2; // Not in the array!
			while(true) {
				int mid = (low+up)/2;
				if(array[mid]<=value) low = mid;
				else up = mid;
				if(up==1+low) break;
			}
			return low;
		}


};

#endif
