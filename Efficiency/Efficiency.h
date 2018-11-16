/*
	 Get detector efficiency at j's pixel.

	 Weiran, Nov.16, 2018.
*/

#ifndef Efficiency_h
#define Efficiency_h

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#define Ndet 1 // number of detector pixels

using namespace std;

class Efficiency
{
	public:
		Efficiency() {
			char* KATRIN = getenv("KATRIN");
			if(KATRIN==0) {
				cout << "Environment variable 'KATRIN' is not defined." << endl;
				exit(0);
			}

			ifstream data(((string)KATRIN + "/Efficiency/data.dat").c_str());
			if(!(data.is_open())) {
				cout << "Data file cannot be opened." << endl;
				exit(0);
			}
			string line;
			int ndet = 0;
			while(true) {
				getline(data, line);
				char f = line[0];
				if(f=='#') continue;
				istringstream iss(line);
				iss >> efficiency[ndet];
				ndet += 1;
				if(ndet >= Ndet) break;
			}
		}

		~Efficiency(){}

		double GetEfficiency(int ndet) { // Count from 1.
			if(ndet>Ndet) {
				cout << "No such detector pixel." << endl;
				return 0;
			}
			return efficiency[ndet-1];
		}

	private:
		double efficiency[Ndet];

};

#endif
