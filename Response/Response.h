/*
	 Get detector response function for given z.

	 Currently only z=0, 0.5, 1, 1.5, ..., 8 meters are calculated. 
	 z=0 is used in this toy model.

	 Weiran, Nov.15, 2018.

	 Allow magnetic field to fluctuate, at the cost of computing time.
	 An approximation is made such that detector response only relies on E-U, but not their values respectively.
	 Weiran, Dec.1, 2018.

	 Add cyclotron radiation.
	 Weiran, Jan.2, 2019.
*/

#ifndef Response_h
#define Response_h

#include <string>
#include <iostream>
#include "../Scattering/prob/ScatterProb.h"
#include "../Scattering/energyloss/ScatEngLoss.h"
#include "../Configure/Configure.h"
#include "SSCTransmissionSynchrotron.h"
#include "SSCWGTS.h"
#include "SSCGeometry.h"
#include "XMLInitialization.h"
#include "KToolbox.h"

#define Nslices 1 // number of source slices
using namespace std;
using namespace Physics;

class Response
{
	public:
		Response() {
			/* Initialize. */
			_B_A = 0;
			_B_S = 0;
			_B_max = 0;
			tolerance = 1e-7;
			response = new TGraph;
		}

		~Response() {
		}

		void initialize(int argc, char** argv) {
			if(IsSynchrotron) {
				/* Setup synchrotron transmission from Kasper. */
				char* KASPER_path = getenv("KASPERSYS");
				if(KASPER_path==0) {
					cout << "Fatal error: Cannot find Kasper path. Maybe Kasper is not installed correctly!" << endl;
					exit(0);
				}
				if(true) { //Always use SSC setup.
					char** input = new char*[2];
					input[0] = argv[0];
					input[1] = (char*)(((string)KASPER_path + "/config/SSC/ssc.xml").data());
					Initialize(2, input);
					myWGTS = KToolbox::GetInstance().Get<SSCWGTS>("WGTS1Radial");
					myWGTS->Initialize();

					/* Cyclotron radiation only depends on z position and pitch angle. */
					trans = new SSCTransmissionSynchrotron;
					trans->Initialize();
				}
			}
		}

		/* Either use SetZ or SetSlice. */
		void SetZ(double z) {
			scat.SetProbCumulate(z);
		}

		void SetSlice(int iSlice) {
			trans->SetVoxel(&(myWGTS->GetVoxel(iSlice)));
			scat.SetProbCumulate(myWGTS->GetSlice(iSlice).GetCenterZ());
			data = &(trans->GetSynchrotronData());
		}

		double GetSyncEngLoss(double angle) {
			return data->GetEnergyLoss(angle);
		}

		void SetupResponse(double B_A, double B_S, double B_max) { // Construct a response function.
			//if(B_A==_B_A && B_S==_B_S && B_max==_B_max) return; // Avoid duplicated calculation.
			_B_A = B_A; _B_S = B_S; _B_max = B_max;
			delete response;
			response = new TGraph();
			int npoint = 0;

			/* For epsilon in filter width, more points are given. */
			double width = _B_A / _B_max * katrin.E_0_center * (gamma(katrin.E_0_center) + 1) / 2.;

			for(double x=0; x<=width; x+=0.01) { // Inside filter width, no inelastic scattering.
				double cosmax = GetCosMax(x);
				UnscatResponse = scat.GetProbCumulate(0, cosmax);
				response->SetPoint(npoint, x, UnscatResponse);
				npoint++;
			}

			for(double x=width+0.1; x<=width+8.5; x+=8.4) { // Flat.
				response->SetPoint(npoint, x, UnscatResponse);
				npoint++;
			}

			double binwidth = 0.1;
			for(double x=9.6; x<35; x+=binwidth) { // Should consider inelastic scattering.
				double ScatResponse = 0;
				for(int s=1; s<=3; s++) {
					for(double epsilon=0; epsilon<x; epsilon+=binwidth) {
						double cosmax = GetCosMax(x-epsilon-0.5*binwidth);
						ScatResponse += engloss.GetEnergyLoss(s, epsilon) * scat.GetProbCumulate(s, cosmax) * binwidth;
					}
				}
				response->SetPoint(npoint, x, UnscatResponse+ScatResponse);
				npoint++;
			}
		}

		void SetupScatParameters(double A1, double A2, double w1, double w2, double e1, double e2, double InelasCS) {
			engloss.SetupParameters(A1, A2, w1, w2, e1, e2);
			scat.SetInelasCrossSection(InelasCS);
		}

		double GetResponse(double E, double U) {
			if(E<U) return 0;
			else return response->Eval(E-U);
		}


	private:
		ScatterProb scat;
		ScatEngLoss engloss;
		KATRIN katrin;
		TGraph* response;
		SSCWGTS* myWGTS;
		SSCTransmissionSynchrotron* trans;
		SSCSynchrotronData* data;
		string GetSourceFile(double z); // Left blank for further usage.
		double _B_A, _B_S, _B_max;
		double UnscatResponse;
		double xmax;
		double tolerance;

		double gamma(double energy) {
			return energy / m_e + 1;
		}

		double transmission(double E, double theta, double U) {
			double thres = E * (1 - pow(sin(theta), 2) * _B_A / _B_S * (gamma(E)+1)/2.);
			if(thres>U) return 1;
			else return 0;
		}

		double SinSquare(double x) {
			return x / katrin.E_0_center / (gamma(katrin.E_0_center)+1) * 2 * _B_S / _B_A;
		}

		double GetCosMax(double x) {
			if(x<=0) return 1;
			double cosmaxback = sqrt(1-_B_S/_B_max); // minimum angle for reflection.
			if(IsSynchrotron) {
				double theta = acos(cosmaxback) / pi * 180;
				double syncloss = GetSyncEngLoss(theta);
				if(SinSquare(x-syncloss)>_B_S/_B_max) return cosmaxback;
			}

			double sinsquare = SinSquare(x);
			double sinsquaremax = -1;

			while(IsSynchrotron) { 
				if(sinsquare>=1) break; //Safe.
				sinsquaremax = sinsquare;
				double theta = asin(sqrt(sinsquare)) / pi * 180.;
				double syncloss = GetSyncEngLoss(theta);
				sinsquare = SinSquare(x-syncloss);
				sinsquare = (sinsquare + sinsquaremax)/2.;
				if(abs(sinsquaremax-sinsquare)<tolerance) break;
			}

			if(sinsquare>1) sinsquare = 1; // Can always be detected.
			double cosmax = sqrt(1-sinsquare);
			if(cosmax<cosmaxback) cosmax = cosmaxback;
			return cosmax;
		}


};

#endif
