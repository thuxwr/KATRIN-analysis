/*
	 Detected spectrum is:
	 beta decay spectrum -> doppler broadening -> detector response.

	 Weiran, Nov.16, 2018.
*/

#ifndef DetectMPI_h
#define DetectMPI_h

#include <iostream>
#include "../Spectrum/Spectrum.h"
#include "../Configure/Configure.h"
#include "../Response/ResponseMPI.h"
#include "../Data/GetDataFile.h"
#include "TMath.h"
#include "mpi.h"
#include "TF1.h"

using namespace std;
using namespace Physics;
using namespace TMath;

class DetectMPI
{
	public:
		DetectMPI(double z=0) {
			detect = this;
			InitTrans(katrin.bv, katrin.T_bt);
			_mass = -1;
			_endpoint = 0;
			broadenspec = new TH1D;
			int initflag;
			MPI_Initialized(&initflag);
			if(!initflag) MPI_Init(NULL, NULL);
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &size);

			/* Determine how many subruns should be arranged for a single core. */
			njobs = CeilNint(((double)data.GetSubrunNum())/size);
			detspeclocal = new double[data.GetSubrunNum()];

		}

		~DetectMPI() {
			delete []trans;
			delete broadenspec;
		}

		void Setup(double bv, double T_bt) { // Re-setup bulk velocity and temperature.
			InitTrans(bv, T_bt);
		}

		void SetMagnetic(double B_A, double B_S, double B_max) {
			IsUpdate = IsUpdate && (_B_A==B_A) && (_B_S==B_S) && (_B_max==B_max);
			_B_A = B_A;
			_B_S = B_S;
			_B_max = B_max;
		}

		void SetScatParams(double A1, double A2, double w1, double w2, double e1, double e2, double InelasCS) {
			_A1 = A1; _A2 = A2; _w1 = w1; _w2 = w2; _e1 = e1; _e2 = e2; _InelasCS = InelasCS;
			IsUpdate = IsUpdate && response.SetupScatParameters(_A1, _A2, _w1, _w2, _e1, _e2, _InelasCS);
		}

		double* DetSpec(double mass, double endpoint) { // For discrete measurement.
			double* detspec = new double[njobs*size];
			/* Check whether this is root process(rank=0) and call MPI setup. */
			if(rank!=0) {
				cout << "This is not root process. MPI failed!" << endl;
				exit(0);
			}

			/* Send all parameters to other cores. */
			cout << "Start setting up response fcn." << endl;
			if(!IsUpdate) {
				double pars[11] = {1, _A1, _A2, _w1, _w2, _e1, _e2, _InelasCS, _B_A, _B_S, _B_max};
				MPI_Bcast(pars, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				response.SetupScatParameters(_A1, _A2, _w1, _w2, _e1, _e2, _InelasCS);
				response.SetupResponse(_B_A, _B_S, _B_max);
			}

			/* Normalization. Not suitable for FirstTritium! */
			/*
			{
				double content = 0;
				double U = 18544.25;
				for(int i=1; i<=broadenspec->GetNbinsX(); i++) {
					double E = broadenspec->GetBinCenter(i);
					content += broadenspec->GetBinContent(i) * response.GetResponse(E, U);
				}
				scale = katrin.Sig_rate/content;
			}
			*/

			/* Distributed calculation for all subruns. */
			if(true) {
				double pars[11] = {-1, mass, endpoint, 0, 0, 0, 0, 0, 0, 0, 0};
				MPI_Bcast(pars, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}

			EfficiencyDistributed(rank, mass, endpoint);
			MPI_Gather(detspeclocal, njobs, MPI_DOUBLE, detspec, njobs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			cout << "End subrun. " << endl;

			IsUpdate = true;
			return detspec;
		}

		void initialize(int argc, char** argv) { response.initialize(argc, argv); }
		void SetSlice(int iSlice) { response.SetSlice(iSlice); }

		void AssociateRun() {
			double pars[11];
			//MPI_Recv(pars, 10, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Bcast(pars, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if(pars[0]>0) { //Setup response.
				response.SetupScatParameters(pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7]);
				response.SetupResponse(pars[8], pars[9], pars[10]);
			}
			else { //Setup efficiency. pars: 1.mass, 2.endpoint.
				EfficiencyDistributed(rank, pars[1], pars[2]);
				MPI_Gather(detspeclocal, njobs, MPI_DOUBLE, detspec, njobs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
		}

	private:
		static DetectMPI* detect;
		KATRIN katrin;
		Data data;
		Spectrum spec;
		ResponseMPI response;
		double** trans; // transition matrix for i'th bin in decay spectrum and j'th bin in broadened spectrum
		int nsize;
		int size;
		int njobs;
		double core_E_cms, core_bv, core_T_bt;
		TH1D* broadenspec;
		double _mass, _endpoint;
		double _B_A, _B_S, _B_max;
		double _A1, _A2, _w1, _w2, _e1, _e2, _InelasCS;
		bool IsUpdate;
		int rank;
		double* detspeclocal;

		void EfficiencyDistributed(int nrank, double mass, double endpoint) {
			/* Setup broadened spectrum according to mass and endpoint. */
			if(mass!=_mass || endpoint!=_endpoint) {
				delete broadenspec;
				TH1D* decayspec = spec.decayspec(mass, endpoint);
				broadenspec = Broaden(decayspec);
				delete decayspec;
				_mass = mass; _endpoint = endpoint;
			}

			for(int subrun=nrank*njobs; subrun<(nrank+1)*njobs; subrun++) {
				if(subrun>=data.GetSubrunNum()) {
					detspeclocal[subrun-nrank*njobs] = 0;
					continue;
				}

				response.SetupEfficiency(efficiency[subrun]);

				double content = 0;
				double U = data.Voltage[subrun];
				for(int i=1; i<=broadenspec->GetNbinsX(); i++) {
					double E = broadenspec->GetBinCenter(i);
					content += broadenspec->GetBinContent(i) * response.GetResponse(E, U);
				}
				detspeclocal[subrun-nrank*njobs] = content;
			}
		}

		double sigma_E(double E_cms, double bv, double T_bt) {
			return Sqrt((E_cms + 2*m_e) * E_cms * k_B * T_bt / M_T2);
		}

		void InitTrans(double bv, double T_bt) { // Allow users to change default temperature and bulk velocity.
			delete []trans;
			trans = NULL;

			/* Decide size of transition matrix. */
			double binwidth = spec.GetBinWidth();
			double sigma = sigma_E(katrin.E_0_center, bv, T_bt);
			nsize = 2 * TMath::CeilNint(5 * sigma / binwidth) + 1; // Calculate convolution inside 5sigma.
			trans = new double* [NbinsDecaySpec];
			for(int i=0; i<NbinsDecaySpec; i++) {
				trans[i] = new double[nsize];
			}

			/* Calculate transition matrix. */
			for(int i=0; i<NbinsDecaySpec; i++) {
				core_E_cms = spec.GetBinCenter(i+1); core_bv = bv; core_T_bt = T_bt;
				TF1* kernel = new TF1("", core, -6*sigma, 6*sigma);
				for(int j=0; j<nsize; j++) {
					double low = (-nsize/2. + j) * spec.GetBinWidth();
					double up = low + spec.GetBinWidth();
					trans[i][j] = kernel->Integral(low, up);
				}
				delete kernel;
			}
		}

		TH1D* Broaden(TH1D* decayspec) {
			TH1D* broadenspec = (TH1D*)decayspec->Clone();
			for(int bin=1; bin<=NbinsDecaySpec; bin++) {
				if(bin <= (nsize+1)/2 || bin >= (NbinsDecaySpec - (nsize+1)/2)) {
					broadenspec->SetBinContent(bin, 0); // Spectrum near the boundary is incorrect. Neglect it.
					continue;
				}

				double content = 0;
				for(int j=0; j<nsize; j++) {
					content += decayspec->GetBinContent(bin-(nsize-1)/2+j) * trans[bin-1][j];
				}

				broadenspec->SetBinContent(bin, content);
			}
			return broadenspec;
		}

		double gamma(double energy) {
			return energy/m_e + 1;
		}

		double beta(double energy) {
			return Sqrt(1. - 1. /gamma(energy) /gamma(energy));
		}

		static double core(double* x, double* par) { // Cannot pass pars to external fcn.
			double E_cms = detect->core_E_cms;
			double bv = detect->core_bv;
			double T_bt = detect->core_T_bt;
			double E_lab = E_cms + x[0];
			double v_lab = detect->beta(E_lab);
			double v_cms = detect->beta(E_cms);
			double v_M = (v_lab - v_cms) / (1 - v_lab * v_cms);
			double u = bv / c;
			double cos_thetamax = Sqrt(1 - detect->katrin.B_S/detect->katrin.B_max);
			double sigma_v  = Sqrt(k_B*T_bt/M_T2);
			double g_vM = 1/(1-cos_thetamax)/2/u * (TMath::Erf((v_M-cos_thetamax*u)/(Sqrt(2)*sigma_v)) - TMath::Erf((v_M-u)/(Sqrt(2)*sigma_v)));
			double g_E = g_vM / (detect->gamma(E_cms) * m_e * v_cms);
			return g_E;
		}

};

DetectMPI* DetectMPI::detect; // For external invoking.

#endif
