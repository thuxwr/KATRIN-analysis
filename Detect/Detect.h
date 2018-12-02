/*
	 Detected spectrum is:
	 beta decay spectrum -> doppler broadening -> detector response.

	 Weiran, Nov.16, 2018.
*/

#ifndef Detect_h
#define Detect_h

#include <iostream>
#include "../Spectrum/Spectrum.h"
#include "../Configure/Configure.h"
#include "../Response/Response.h"

using namespace std;
using namespace Physics;
using namespace TMath;

class Detect
{
	public:
		Detect(double z=0) {
			detect = this;
			InitTrans(katrin.bv, katrin.T_bt);
			_mass = -1;
			_endpoint = 0;
			response.SetZ(z);
		}

		~Detect() {
			delete []trans;
			delete broadenspec;
		}

		void Setup(double bv, double T_bt) { // Re-setup bulk velocity and temperature.
			InitTrans(bv, T_bt);
		}

		void SetMagnetic(double B_A, double B_S, double B_max) {
			_B_A = B_A;
			_B_S = B_S;
			_B_max = B_max;
		}

		double* DetSpec(double mass, double endpoint, int nvoltage, double* voltage) { // For discrete measurement, with nvoltage thresholds each being voltage[i].
			if(mass!=_mass || endpoint!=_endpoint) {
				delete broadenspec;
				TH1D* decayspec = spec.decayspec(mass, endpoint);
				broadenspec = Broaden(decayspec);
				delete decayspec;
			}

			double* detspec = new double[nvoltage];
			response.SetupResponse(_B_A, _B_S, _B_max);

			/* Normalization. */
			double scale;
			{
				double content = 0;
				double U = 18544.25;
				for(int i=1; i<=broadenspec->GetNbinsX(); i++) {
					double E = broadenspec->GetBinCenter(i);
					content += broadenspec->GetBinContent(i) * response.GetResponse(E, U);
				}
				scale = katrin.Sig_rate/content;
			}

			for(int n=0; n<nvoltage; n++) {
				double content = 0;
				double U = voltage[n];
				for(int i=1; i<=broadenspec->GetNbinsX(); i++) {
					double E = broadenspec->GetBinCenter(i);
					content += broadenspec->GetBinContent(i) * response.GetResponse(E, U);
				}
				detspec[n] = content * scale;
			}

			return detspec;
		}

	private:
		static Detect* detect;
		KATRIN katrin;
		Spectrum spec;
		Response response;
		double** trans; // transition matrix for i'th bin in decay spectrum and j'th bin in broadened spectrum
		int nsize;
		double core_E_cms, core_bv, core_T_bt;
		TH1D* broadenspec;
		double _mass, _endpoint;
		double _B_A, _B_S, _B_max;

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

Detect* Detect::detect; // For external invoking.

#endif
