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

class Detect
{
	public:
		Detect() {
			detect = this;
			InitTrans(KATRIN.bv, KATRIN.T_bt);
		}

		~Detect() {
			delete []trans;
		}

		void Setup(double bv, double T_bt) { // Re-setup bulk velocity and temperature.
			InitTrans(bv, T_bt);
		}

		TH1D* DetSpec(TH1D* decayspec, double z=0) { // Get detected spectrum from decay spectrum.
			TH1D* broadenspec = Broaden(decayspec);
			TH1D* detspec = new TH1D("", "", NbinsDetSpec, DetLow+KATRIN.E_0_center, DetUp+KATRIN.E_0_center);
			for(int bin=1; bin<=NbinsDetSpec; bin++) {
				double content = 0;
				double U = detspec->GetBinCenter(bin);
				for(int i=1; i<=broadenspec->GetNbinsX(); i++) {
					double E = broadenspec->GetBinCenter(i);
					content += broadenspec->GetBinContent(i) * response.GetResponse(E, U, z);
				}
				detspec->SetBinContent(bin, content);
			}

			double scale = KATRIN.Sig_rate/detspec->GetBinContent(detspec->FindBin(18544.25));
			detspec->Scale(scale);
			delete broadenspec;
			return detspec;
		}

	private:
		static Detect* detect;
		KATRIN KATRIN;
		Spectrum spec;
		Response response;
		double** trans; // transition matrix for i'th bin in decay spectrum and j'th bin in broadened spectrum
		int nsize;
		double core_E_cms, core_bv, core_T_bt;

		double sigma_E(double E_cms, double bv, double T_bt) {
			return sqrt((E_cms + 2*m_e) * E_cms * k_B * T_bt / M_T2);
		}

		void InitTrans(double bv, double T_bt) { // Allow users to change default temperature and bulk velocity.
			delete []trans;
			trans = NULL;

			/* Decide size of transition matrix. */
			double binwidth = spec.GetBinWidth();
			double sigma = sigma_E(KATRIN.E_0_center, bv, T_bt);
			nsize = 2 * TMath::Ceil(5 * sigma / binwidth) + 1; // Calculate convolution inside 5sigma.
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
			return sqrt(1. - 1. /gamma(energy) /gamma(energy));
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
			double cos_thetamax = sqrt(1 - detect->KATRIN.B_S/detect->KATRIN.B_max);
			double sigma_v  = sqrt(k_B*T_bt/M_T2);
			double g_vM = 1/(1-cos_thetamax)/2/u * (TMath::Erf((v_M-cos_thetamax*u)/(sqrt(2)*sigma_v)) - TMath::Erf((v_M-u)/(sqrt(2)*sigma_v)));
			double g_E = g_vM / (detect->gamma(E_cms) * m_e * v_cms);
			return g_E;
		}

};

Detect* Detect::detect; // For external invoking.

#endif
