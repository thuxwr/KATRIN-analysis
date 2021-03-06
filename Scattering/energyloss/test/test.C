#include "../ScatEngLoss.h"
using namespace std;

void test()
{
	ScatEngLoss scat;
	scat.SetupParameters(0.204, 0.0556, 1.85, 12.5, 12.6, 14.3);
	cout << "check1" << endl;
	for(int i=0; i<1000000; i++) {
		scat.GetEnergyLoss(2, (i+0.1)*50./1000000);
	}
	cout << "check2" << endl;
	scat.pdf[2]->Draw();
}
