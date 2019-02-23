/*
	 Get energy loss from inelastic scattering.

	 Weiran, Jan.3, 2019.
*/

#ifndef ScatEngLoss_h
#define ScatEngLoss_h

#include <iostream>
#include <complex>
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
			for(int i=0; i<ScatTimesMax; i++) pdf[i] = new TH1D;

			/* Initialize FFT. */
			npoints = 100000;
			upboundary = 10*ScatTimesMax;

			fft = TVirtualFFT::FFT(1, &npoints, "R2C M K");
			fft_back = TVirtualFFT::FFT(1, &npoints, "C2R M K");
		}
		
		~ScatEngLoss() {
		}

		bool SetupParameters(double A1, double A2, double w1, double w2, double e1, double e2) {
			if(_A1==A1 && _A2==A2 && _w1==w1 && _w2==w2 && _e1==e1 && _e2==e2 ) return false;
			_A1 = A1; _A2 = A2; _w1 = w1; _w2 = w2; _e1 = e1; _e2 = e2;
			for(int i=0; i<ScatTimesMax; i++) delete pdf[i];

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

			for(int ScatTimes=0; ScatTimes<ScatTimesMax; ScatTimes++) {
				double* re = new double[npoints];
				double* im = new double[npoints];
				double* pdfvalue = new double[npoints];
				for(int i=0; i<npoints; i++) {
					complex<double> transformed(transformed_re[i], transformed_im[i]);
					complex<double> convoluted = pow(transformed, ScatTimes+1);
					re[i] = convoluted.real();
					im[i] = convoluted.imag();
				}
				fft_back->SetPointsComplex(re, im);
				fft_back->Transform();
				fft_back->GetPoints(pdfvalue);

				double binwidth = upboundary * 10 / npoints;
				pdf[ScatTimes] = new TH1D("", "", npoints, 0, upboundary*10);
				for(int bin=1; bin<=npoints; bin++) pdf[ScatTimes]->SetBinContent(bin, pdfvalue[bin-1]);
				pdf[ScatTimes]->Scale(1./pdf[ScatTimes]->Integral("width"));

				delete[] re; delete[] im; delete[] pdfvalue;
			}

			delete[] in; delete[] x; 
			return true;
		}

		void SetNPoints(int n) { npoints = n; }
		void SetUpBoundary(double boundary) { upboundary = boundary; }

		double GetEnergyLoss(int ScatTimes, double epsilon) {
			if(ScatTimes<1) {
				cout << "Energy loss is only valid for Scatter times >= 1." << endl;
				return 0;
			}

			if(ScatTimes>ScatTimesMax) {
				cout << "Scatter times greater than " << ScatTimesMax << " is not supported yet." << endl;
				return 0;
			}

			if(epsilon>upboundary) {
				cout << "Energy loss greater than " << upboundary << "eV is not calculated yet." << endl;
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
		TH1D* pdf[ScatTimesMax];
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

