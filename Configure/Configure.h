/*
	 A global configure file to setup all physical constants and detector parameters for KATRIN.

	 Weiran, Nov.16, 2018.
*/

#ifndef Configure_h
#define Configure_h

class KATRIN
{
	public:
		/* Electromagnetic field. */
		double B_A = 2.68; // in unit:Gauss
		double B_S = 2.52e4; // in unit:Gauss
		double B_max = 4.2e4; // in unit:Gauss

		/* Source. */
		double epsilon_T = 0.95; // T2 purity

};

namespace Physics {
	double m_e = 0.511e6; // electron mass, 0.511MeV
	double alpha = 0.00729735257; // fine structure constant
	double pi = TMath::Pi(); 
}
#endif
