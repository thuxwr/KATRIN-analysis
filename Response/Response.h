/*
	 Get detector response function for given z.

	 Currently only z=0, 0.5, 1, 1.5, ..., 8 meters are calculated. 
	 z=0 is used in this toy model.

	 Weiran, Nov.15, 2018.

	 Allow magnetic field to fluctuate, at the cost of computing time.
	 An approximation is made such that detector response only relies on E-U, but not their values respectively.
	 Weiran, Dec.1, 2018.
*/

#ifndef Response_h
#define Response_h

#include <string>
#include <iostream>
#include "../Scattering/prob/ScatterProb.h"
#include "../Scattering/energyloss/EngLoss.h"
#include "../Configure/Configure.h"

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
		}

		~Response() {
			delete response;
		}

		void SetZ(double z) {
			scat.SetProbCumulate(z);
		}

		void SetupResponse(double B_A, double B_S, double B_max) { // Construct a response function.
			if(B_A==_B_A && B_S==_B_S && B_max==_B_max) return; // Avoid duplicated calculation.
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
						ScatResponse += engloss.GetPdf(s, epsilon) * scat.GetProbCumulate(s, cosmax) * binwidth;
					}
				}
				response->SetPoint(npoint, x, UnscatResponse+ScatResponse);
				npoint++;
			}
		}

		double GetResponse(double E, double U) {
			if(E<U) return 0;
			else return response->Eval(E-U);
		}


	private:
		ScatterProb scat;
		EngLoss engloss;
		KATRIN katrin;
		TGraph* response;
		string GetSourceFile(double z); // Left blank for further usage.
		double _B_A, _B_S, _B_max;
		double UnscatResponse;
		double xmax;

		double gamma(double energy) {
			return energy / m_e + 1;
		}

		double transmission(double E, double theta, double U) {
			double thres = E * (1 - pow(sin(theta), 2) * _B_A / _B_S * (gamma(E)+1)/2.);
			if(thres>U) return 1;
			else return 0;
		}

		double GetCosMax(double x) {
			double cosmaxback = sqrt(1-_B_S/_B_max); // minimum angle for reflection.
			double sinsquare = x / katrin.E_0_center / (gamma(katrin.E_0_center)+1) * 2 * _B_S / _B_A;
			if(sinsquare>1) sinsquare = 1; // Can always be detected.
			double cosmax = sqrt(1-sinsquare);
			if(cosmax<cosmaxback) cosmax = cosmaxback;
			return cosmax;
		}


};

#endif
