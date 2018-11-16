/*
	 Detected spectrum is:
	 beta decay spectrum -> doppler broadening -> detector response.

	 Weiran, Nov.16, 2018.
*/

#ifndef Detect_h
#define Detect_h

#include <iostream>
#include "../Spectrum/Spectrum.h"

using namespace std;

class Detect
{
	public:
		Detect(){
			/* */
		}

		~Detect(){}


	private:
		Spectrum spec;
		TH1D* broadenspec;
		double* trans; // transition matrix for i'th bin in decay spectrum and j'th bin in broadened spectrum

};

#endif
