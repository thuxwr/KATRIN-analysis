/*
	 Get energy loss from inelastic scattering.

	 Weiran, Jan.3, 2019.
*/

#ifndef ScatEngLoss_h
#define ScatEngLoss_h

#include <iostream>
#include "TH1D.h"
#include "TVirtualFFT.h"

using namespace std;
using namespace TMath;

class ScatEngLoss
{
	public:
		ScatEngLoss() {
			_A1=0; _A2=0; _w1=0; _w2=0; _e1=0; _e2=0;
			ec = 14.09;
			for(int i=0; i<3; i++) pdf[i] = new TH1D;

			/* Initialize FFT. */
			npoints = 100000;
			upboundary = 50;

			fft = TVirtualFFT::FFT(1, &npoints, "R2C M K");
			fft_back = TVirtualFFT::FFT(1, &npoints, "C2R M K");
		}
		
		~ScatEngLoss() {
		}

		bool SetupParameters(double A1, double A2, double w1, double w2, double e1, double e2) {
			if(_A1==A1 && _A2==A2 && _w1==w1 && _w2==w2 && _e1==e1 && _e2==e2 ) return false;
			_A1 = A1; _A2 = A2; _w1 = w1; _w2 = w2; _e1 = e1; _e2 = e2;
			for(int i=0; i<3; i++) delete pdf[i];

			double* in = new double[npoints];
			double* x = new double[npoints];
			for(int i=0; i<npoints; i++) {
				x[i] = upboundary * 10 / npoints * (i+0.5);
				in[i] = GetBasePdf(x[i]);
			}

			double* transformed_re = new double[npoints];
			double* transformed_im = new double[npoints];

			fft->SetPoints(in);
			fft->Transform();

			fft->GetPointsComplex(transformed_re, transformed_im);
			double* scat2_re = new double[npoints]; 
			double* scat2_im = new double[npoints]; 
			double* scat3_re = new double[npoints]; 
			double* scat3_im = new double[npoints]; 
			for(int i=0; i<npoints; i++) {
				scat2_re[i] = pow(transformed_re[i], 2) - pow(transformed_im[i], 2);
				scat2_im[i] = 2. * transformed_re[i] * transformed_im[i];
				scat3_re[i] = pow(transformed_re[i], 3) - 3. * transformed_re[i] * pow(transformed_im[i], 2);
				scat3_im[i] = 3. * pow(transformed_re[i], 2) * transformed_im[i] - pow(transformed_im[i], 3);
			}

			double** pdfvalue = new double*[3];
			pdfvalue[0] = in;
			pdfvalue[1] = new double[npoints];
			pdfvalue[2] = new double[npoints];

			fft_back->SetPointsComplex(scat2_re, scat2_im);
			fft_back->Transform();
			fft_back->GetPoints(pdfvalue[1]);
			fft_back->SetPointsComplex(scat3_re, scat3_im);
			fft_back->Transform();
			fft_back->GetPoints(pdfvalue[2]);

			for(int i=0; i<3; i++) {
				double binwidth = upboundary * 10 / npoints;
				pdf[i] = new TH1D("", "", npoints, 0, upboundary*10);
				for(int bin=1; bin<=npoints; bin++) pdf[i]->SetBinContent(bin, pdfvalue[i][bin-1]);
				pdf[i]->Scale(1./pdf[i]->Integral("width"));
			}

			delete[] in; delete[] x; delete[] transformed_re; delete[] transformed_im;
			delete[] scat2_re; delete[] scat2_im; delete[] scat3_re; delete[] scat3_im;
			for(int i=1; i<3; i++) delete[] pdfvalue[i];
			delete[] pdfvalue;
			return true;
		}

		void SetNPoints(int n) { npoints = n; }
		void SetUpBoundary(double boundary) { upboundary = boundary; }

		double GetEnergyLoss(int ScatTimes, double epsilon) {
			if(ScatTimes<1) {
				cout << "Energy loss is only valid for Scatter times >= 1." << endl;
				return 0;
			}

			if(ScatTimes>3) {
				cout << "Scatter times greater than 3 is not supported yet." << endl;
				return 0;
			}

			if(epsilon>50) {
				cout << "Energy loss greater than 50eV is not calculated yet." << endl;
				return 0;
			}

			if(epsilon<0) {
				cout << "Energy loss smaller than 0eV is unphysical." << endl;
				return 0;
			}

			return pdf[ScatTimes-1]->GetBinContent(pdf[ScatTimes-1]->FindBin(epsilon));
		}

	private:
		double _A1, _A2, _w1, _w2, _e1, _e2;
		double ec;
	public:
		TH1D* pdf[3];
		int npoints;
		double upboundary;
		TVirtualFFT* fft;
		TVirtualFFT* fft_back;

		double GetBasePdf(double epsilon) {
			if(epsilon<ec) return _A1 * exp(-2. * pow((epsilon-_e1)/_w1, 2));
			else return _A2 * _w2 * _w2 / (_w2*_w2 + 4. * pow(epsilon-_e2, 2));
		}






};

#endif

