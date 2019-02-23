/*
	 Get detector response function for given z.

	 Currently only z=0, 0.5, 1, 1.5, ..., 8 meters are calculated. 
	 z=0 is used in this toy model.

	 Weiran, Nov.15, 2018.

	 Allow magnetic field to fluctuate, at the cost of computing time.
	 An approximation is made such that detector response only relies on E-U, but not their values respectively.
	 Weiran, Dec.1, 2018.

	 Add cyclotron radiation.
	 Make sure to setup slice before any calculation!
	 Weiran, Jan.2, 2019.
*/

#ifndef Response_h
#define Response_h

#include <string>
#include <iostream>
#include "../Scattering/prob/ScatterProb.h"
#include "../Scattering/energyloss/ScatEngLoss.h"
#include "../Configure/Configure.h"
#include "../Data/GetDataFile.h"
#include "SSCTransmissionSynchrotron.h"
#include "SSCWGTS.h"
#include "SSCGeometry.h"
#include "XMLInitialization.h"
#include "KToolbox.h"
#include "TGraph.h"

#define Nslices 1 // number of source slices
#define NPointsMax 1000

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
			response = new double*[NPixels];
			for(int pixel=0; pixel<NPixels; pixel++) response[pixel] = new double[NPointsMax];
			X = new double[NPointsMax];
			IsUpdate = false;
			_slice = -1;
			//HV = new double[26] {16975, 17175, 17375, 17575, 17775, 17975, 18175, 18275, 18375, 18400, 18425, 18450, 18475, 18485, 18495, 18505, 18515, 18525, 18535, 18545, 18555, 18565, 18575, 18585, 18595, 18605};
		}

		~Response() {
		}

		void initialize(int argc, char** argv) {
			/* Setup synchrotron transmission from Kasper. */
			string KASPER_path = getenv("KASPERSYS");
			if(KASPER_path=="") {
				cout << "Fatal error: Cannot find Kasper path. Maybe Kasper is not installed correctly!" << endl;
				exit(0);
			}
			if(true) { //Always use SSC setup.
				char** input = new char*[2];
				input[0] = argv[0];
				KASPER_path += "/config/SSC/ssc.xml";
				input[1] = (char*)(KASPER_path.data());
				Initialize(2, input);
				myWGTS = KToolbox::GetInstance().Get<SSCWGTS>("WGTS1Radial");
				myWGTS->Initialize();

				/* Cyclotron radiation only depends on z position and pitch angle. */
				if(IsSynchrotron) {
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
			if(iSlice==_slice) return;
			_slice = iSlice;
			if(IsSynchrotron) {
				trans->SetVoxel(&(myWGTS->GetVoxel(iSlice)));
				data = &(trans->GetSynchrotronData());
			}
			IsUpdate = false;
		}

		double GetSyncEngLoss(double angle) {
			return data->GetEnergyLoss(angle);
		}

		void SetupResponse(double B_A, double B_S, double B_max) { // Construct a response function.
			if(CheckUpdate(B_A, B_S, B_max)) return; // Avoid duplicated calculation.
			_B_A = B_A; _B_S = B_S; _B_max = B_max;
			if(!IsUpdate) scat.SetProbCumulate(myWGTS->GetSlice(_slice).GetCenterZ());

			/* For epsilon in filter width, more points are given. */
			double width = _B_A / _B_max * katrin.E_0_center * (gamma(katrin.E_0_center) + 1) / 2.;

			/* For first tritium data, need to calculate decay spectrum for a broad range. */
			bool IsXCalculated = false;
			NPoints = 0;

			for(int npixel=0; npixel<NPixels; npixel++) {
				int npoint = 0;
				_B_A = B_A - gs.B_A[0] + gs.B_A[npixel];
				_B_S = B_S - gs.B_S[0] + gs.B_S[npixel];
				_B_max = B_max - gs.B_max[0] + gs.B_max[npixel];
				for(double x=0; x<=width; x+=0.01) { // Inside filter width, no inelastic scattering.
					if(!IsXCalculated) X[npoint] = x;
					/* Reset magnetic field for each ring. */
					double voltage = x - gs.E[0] + gs.E[npixel]; 
					double cosmax = GetCosMax(voltage);
					UnscatResponse = scat.GetProbCumulate(0, cosmax);
					response[npixel][npoint] = UnscatResponse;
					npoint++;
				}

				for(double x=width+0.1; x<=width+8.5; x+=8.4) { // Flat.
					if(!IsXCalculated) X[npoint] = x;
					response[npixel][npoint] = UnscatResponse;
					npoint++;
				}

				double binwidth = 0.2;
				for(double x=9.6; x<120; x+=binwidth) { // Should consider inelastic scattering.
					if(!IsXCalculated) X[npoint] = x;
					double ScatResponse = 0;
					for(double epsilon=0; epsilon<x; epsilon+=binwidth) {
						double voltage = x - gs.E[0] + gs.E[npixel];
						double cosmax = GetCosMax(voltage-epsilon-0.5*binwidth);
						for(int s=1; s<=ScatTimesMax; s++) {
							ScatResponse += engloss.GetEnergyLoss(s, epsilon) * scat.GetProbCumulate(s, cosmax) * binwidth;
						}
					}
					response[npixel][npoint] = UnscatResponse+ScatResponse;
					npoint++;
				}
				if(!IsXCalculated) NPoints = npoint;
				IsXCalculated = true;
			}
			_B_A = B_A;
			_B_S = B_S;
			_B_max = B_max;
			IsUpdate = true;
		}

			/*
		void SetupResponseTotal(double B_A, double B_S, double B_max) { //Multi-pixel structure, calculated with single core.
			if(CheckUpdate(B_A, B_S, B_max)) return; // Avoid duplicated calculation.
			SetSlice(0);
			SetupResponse(B_A, B_S, B_max);
			int npoints = GetResponse()->GetN();
			int nslice = GetNSlices();
			double TotDens = 0;
			for(int i=0; i<nslice; i++) TotDens += GetColumnDensity(i);

			double yvalue[npoints];
			for(int i=0; i<npoints; i++) yvalue[i] = 0;
			for(int i=0; i<nslice; i++) {
				SetSlice(i);
				SetupResponse(B_A, B_S, B_max);
				double* y = GetResponse()->GetY();
				double weight = GetColumnDensity(i)/TotDens;
				for(int j=0; j<npoints; j++) yvalue[j] += y[j] * weight;
			}
			double xvalue[npoints];
			double* x = GetResponse()->GetX();
			for(int i=0; i<npoints; i++) xvalue[i] = x[i];
			delete response;
			response = new TGraph(npoints, xvalue, yvalue);
			IsUpdate = true;
		}
		*/

		//TGraph* GetResponse() { return response[0][0]; }

		bool SetupScatParameters(double A1, double A2, double w1, double w2, double e1, double e2, double InelasCS) {
			bool IsSetup1 = engloss.SetupParameters(A1, A2, w1, w2, e1, e2);
			bool IsSetup2 = scat.SetInelasCrossSection(InelasCS);
			IsUpdate = !(IsSetup1 || IsSetup2);
			return IsUpdate;
		}

		/*
		double GetResponse(double E, double U) {
			if(E<U) return 0;
			else if(U<17075) return response[0]->Eval(E-U);
			else if(U<17275) return response[1]->Eval(E-U);
			else if(U<17475) return response[2]->Eval(E-U);
			else if(U<17675) return response[3]->Eval(E-U);
			else if(U<17875) return response[4]->Eval(E-U);
			else if(U<18075) return response[5]->Eval(E-U);
			else if(U<18225) return response[6]->Eval(E-U);
			else if(U<18325) return response[7]->Eval(E-U);
			else if(U<18387.5) return response[8]->Eval(E-U);
			else if(U<18412.5) return response[9]->Eval(E-U);
			else if(U<18437.5) return response[10]->Eval(E-U);
			else if(U<18462.5) return response[11]->Eval(E-U);
			else if(U<18480) return response[12]->Eval(E-U);
			else if(U<18490) return response[13]->Eval(E-U);
			else if(U<18500) return response[14]->Eval(E-U);
			else if(U<18510) return response[15]->Eval(E-U);
			else if(U<18520) return response[16]->Eval(E-U);
			else if(U<18530) return response[17]->Eval(E-U);
			else if(U<18540) return response[18]->Eval(E-U);
			else if(U<18550) return response[19]->Eval(E-U);
			else if(U<18560) return response[20]->Eval(E-U);
			else if(U<18570) return response[21]->Eval(E-U);
			else if(U<18580) return response[22]->Eval(E-U);
			else if(U<18590) return response[23]->Eval(E-U);
			else if(U<18600) return response[24]->Eval(E-U);
			else return response[25]->Eval(E-U);
		}
		*/

		double GetColumnDensity(int iSlice) {
			return scat.GetColumnDensity(myWGTS->GetSlice(iSlice).GetCenterZ());
		}

		double GetColumnDensity() {
			return scat.GetColumnDensity(myWGTS->GetSlice(_slice).GetCenterZ());
		}

		int GetNSlices() {
			return myWGTS->GetNSlices();
		}

		double GetTotalA(int iSlice) {
			double total = 0;
			SSCSlice& myslice = myWGTS->GetSlice(iSlice);
			for(int ring=0; ring<myslice.GetNRings(); ring++) {
				SSCRing& myring = myslice.GetRing(ring);
				total += myring.GetSegment(0).GetA() * myring.GetNSegments();
			}
			return total;
		}

		double GetRingA(int iSlice, int iRing) {
			SSCRing& myring = myWGTS->GetSlice(iSlice).GetRing(iRing);
			return myring.GetNSegments() * myring.GetSegment(0).GetA();
		}

		double GetSegmentA(int npixel, int slice=-1) {
			int iRing = 0;
			if(npixel>3) iRing = (npixel+8)/12;
			if(slice<0) slice = _slice;
			SSCRing& myring = myWGTS->GetSlice(slice).GetRing(iRing);
			return myring.GetSegment(0).GetA();
		}


		bool CheckUpdate(double B_A, double B_S, double B_max) {
			if(B_A==_B_A && B_S==_B_S && B_max==_B_max && IsUpdate) return true;
			return false;
		}

		double** response; //dimensions: pixel, E-U.
		double* X; // X axis for E-U.
		int NPoints;

	private:
		ScatterProb scat;
		ScatEngLoss engloss;
		KATRIN katrin;
		GlobalSimulation gs;
		SSCWGTS* myWGTS;
		SSCTransmissionSynchrotron* trans;
		SSCSynchrotronData* data;
		string GetSourceFile(double z); // Left blank for further usage.
		double _B_A, _B_S, _B_max;
		double UnscatResponse;
		double xmax;
		double tolerance;
		int _slice;
		bool IsUpdate;
		//double* HV;

		double gamma(double energy) {
			return energy / m_e + 1;
		}

		double transmission(double E, double theta, double U) {
			double thres = E * (1 - pow(sin(theta), 2) * _B_A / _B_S * (gamma(E)+1)/2.);
			if(thres>U) return 1;
			else return 0;
		}

		double SinSquare(double x, double hv = -1) {
			if(hv<=0) hv = katrin.E_0_center;
			return x / hv / (gamma(hv)+1) * 2 * _B_S / _B_A;
		}

		double GetCosMax(double x, double hv = -1) {
			if(x<=0) return 1;
			double cosmaxback = sqrt(1-_B_S/_B_max); // minimum angle for reflection.
			if(IsSynchrotron) {
				double theta = acos(cosmaxback) / pi * 180;
				double syncloss = GetSyncEngLoss(theta);
				if(SinSquare(x-syncloss, hv)>_B_S/_B_max) return cosmaxback;
			}

			double sinsquare = SinSquare(x, hv);
			double sinsquaremax = -1;

			while(IsSynchrotron) { 
				if(sinsquare>=1) break; //Safe.
				sinsquaremax = sinsquare;
				double theta = asin(sqrt(sinsquare)) / pi * 180.;
				double syncloss = GetSyncEngLoss(theta);
				sinsquare = SinSquare(x-syncloss, hv);
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
