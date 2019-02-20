/*
	 Get real data via Idle. Currently support first tritium data.

	 Weiran, Feb.6, 2019.
*/

#define NSubrunMax 5000
#define NPixels 148

#ifndef GetDataFile_h
#define GetDataFile_h

#include "KIRunSummaryDocument.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "KIIdle.h"
#include "KTree.h"
#include "TMath.h"
#include "../Configure/Configure.h"

using namespace std;
using namespace katrin;
using namespace TMath;

class Data 
{
	public:
		Data() {
			char* KATRINpath = getenv("KATRIN");
			path = KATRINpath;
			path = path + "/Data";
			LiveTime = new double*[NSubrunMax];
			Efficiency = new double*[NSubrunMax];
			EventCount = new int*[NSubrunMax];

			SetDataset("FirstTritium.katrin");
			if(!GetDataFile()) {
				GetDataIdle();
				DumpData();
			}
			for(int npixel=0; npixel<NPixels; npixel++) goodlist[npixel] = true;
			Cut();
		}

		~Data(){}

		void SetDataset(string datasetname) {
			DataSet = datasetname;
		}

		void Cut() { // Data selection, and refresh dataset.
			int subruncount = 0;
			for(int subrun=0; subrun<NSubrun; subrun++) {
				/* Contain all required slow control data. */
				if(IsNaN(TritiumPurity[subrun]) || TritiumPurity[subrun]<0) continue;
				if(IsNaN(ColumnDensity[subrun]) || ColumnDensity[subrun]<0) continue;

				/* Stable gas flow. */
				if(Abs(ColumnDensity[subrun]-4.46e21)>5e18) continue;

				/* Energy in [-100, 50] eV. */
				if(Voltage[subrun]<Katrin.E_0_center-100 || Voltage[subrun]>Katrin.E_0_center+50) continue;

				/* Copy selected data. */
				TritiumPurity[subruncount] = TritiumPurity[subrun];
				ColumnDensity[subruncount] = ColumnDensity[subrun];
				Voltage[subruncount] = Voltage[subrun];
				for(int npixel=0; npixel<NPixels; npixel++) {
					Efficiency[subruncount][npixel] = Efficiency[subrun][npixel];
					if(Efficiency[subruncount][npixel]<=0) goodlist[npixel] = false;
					LiveTime[subruncount][npixel] = LiveTime[subrun][npixel];
					EventCount[subruncount][npixel] = EventCount[subrun][npixel];
				}
				subruncount ++;
			}
			NSubrun = subruncount;
		}


		int GetSubrunNum() { return NSubrun; }
		int GetSubrunCount(int subrun) {
			int count = 0;
			for(int npixel=0; npixel<NPixels; npixel++) {
				if(Efficiency[subrun][npixel]<=0) continue;
				count += EventCount[subrun][npixel];
			}
			return count;
		}

		double GetSubrunCountCorrected(int subrun) {
			double count = 0;
			for(int npixel=0; npixel<NPixels; npixel++) {
				if(Efficiency[subrun][npixel]<=0) continue;
				count += EventCount[subrun][npixel]/Efficiency[subrun][npixel];
			}
			return count;
		}

		bool IsPixelBroken(int pixel) {
			return !goodlist[pixel];
		}

	private:
		KIIdle idle;
		KATRIN Katrin;
		string DataSet;
		string path;
		int NSubrun;
		bool goodlist[NPixels];

		bool GetDataFile() {
			ifstream fin((path+"/data.dat").c_str());
			if(!fin.is_open()) return false;
			fin >> NSubrun;
			for(int run=0; run<NSubrun; run++) fin >> ColumnDensity[run];
			for(int run=0; run<NSubrun; run++) fin >> TritiumPurity[run];
			for(int run=0; run<NSubrun; run++) fin >> Voltage[run];

			for(int run=0; run<NSubrun; run++) {
				LiveTime[run] = new double[NPixels];
				for(int npixel=0; npixel<NPixels; npixel++)
					fin >> LiveTime[run][npixel];
			}
			for(int run=0; run<NSubrun; run++) {
				Efficiency[run] = new double[NPixels];
				for(int npixel=0; npixel<NPixels; npixel++)
					fin >> Efficiency[run][npixel];
			}
			for(int run=0; run<NSubrun; run++) {
				EventCount[run] = new int[NPixels];
				for(int npixel=0; npixel<NPixels; npixel++)
					fin >> EventCount[run][npixel];
			}
			fin.close();
			return true;
		}

		void GetDataIdle() {
			KTree filter = (make_tree
					("Limit", 2000)
					("FileName", "RunSummary3c")
			);
			KTabree filelist;

			auto dataset = idle.FindDataSet(DataSet);
			dataset->GetFileList(filter, filelist);
			int NRuns = filelist.NumberOfRows();

			int subrun = 0;

			/* Loop over all runs. */
			for(int run=0; run<NRuns; run++) {
				KIRunSummaryDocument doc;
				string filename = (string)filelist[run]["FileName"] + "@" + DataSet;
				doc.AddInput(filename.c_str());
				/* All subruns. */
				for(auto subRunDoc: doc.GetSubRunDocumentList()) {
					/* HV readout. */
					try {
						Voltage[subrun] = subRunDoc["Values/Transmission/K35VoltageReading"].As<double>();
						if(IsNaN(Voltage[subrun])) Voltage[subrun] = -1;
					}
					catch(KException& e) {
						Voltage[subrun] = -1;
					}

					/* Tritium purity. */
					try {
						TritiumPurity[subrun] = subRunDoc["Values/LARA/TritiumPurity"].As<double>();
						if(IsNaN(TritiumPurity[subrun])) TritiumPurity[subrun] = -1;
					}
					catch(KException& e) {
						TritiumPurity[subrun] = -1;
					}

					/* Column density. */
					try {
						ColumnDensity[subrun] = subRunDoc["Values/ColumnDensity/ColumnDensity"].As<double>();
						if(IsNaN(ColumnDensity[subrun])) ColumnDensity[subrun] = -1;
					}
					catch(KException& e) {
						ColumnDensity[subrun] = -1;
					}

					/* Event count. */
					EventCount[subrun] = new int[NPixels];
					try {
						for(int pixel=0; pixel<NPixels; pixel++) {
							EventCount[subrun][pixel] = subRunDoc["Values/Detector/EventCount"][pixel].As<int>();
						}
					}
					catch(KException& e) {
						for(int pixel=0; pixel<NPixels; pixel++) {
							EventCount[subrun][pixel] = -1;
						}
					}

					/* Livetime. */
					LiveTime[subrun] = new double[NPixels];
					try {
						for(int pixel=0; pixel<NPixels; pixel++) {
							LiveTime[subrun][pixel] = subRunDoc["Values/Detector/LiveTime"][pixel].As<double>();
						}
					}
					catch(KException& e) {
						for(int pixel=0; pixel<NPixels; pixel++) {
							LiveTime[subrun][pixel] = -1;
						}
					}

					/* Efficiency. */
					Efficiency[subrun] = new double[NPixels];
					try {
						for(int pixel=0; pixel<NPixels; pixel++) {
							Efficiency[subrun][pixel] = subRunDoc["Values/Detector/RelativeEfficiency"][pixel].As<double>();
							if(IsNaN(Efficiency[subrun][pixel])) Efficiency[subrun][pixel] = -1;
						}
					}
					catch(KException& e) {
						try {
							for(int pixel=0; pixel<NPixels; pixel++) {
								Efficiency[subrun][pixel] = subRunDoc["Values/Detector/EfficiencyCorrectedCount"][pixel].As<double>();
								if(IsNaN(Efficiency[subrun][pixel])) Efficiency[subrun][pixel] = 0;
								else Efficiency[subrun][pixel] = EventCount[subrun][pixel]/Efficiency[subrun][pixel];
							}
						}
						catch(KException& e) {
							for(int pixel=0; pixel<NPixels; pixel++) {
								Efficiency[subrun][pixel] = 1; // if no efficiency provided, then assume to be 1.
							}
						}
					}

					subrun += 1;
				}
			}

			NSubrun = subrun;
		}

		void DumpData() {
			ofstream fout((path+"/data.dat").c_str());
			fout.precision(9);

			/* First line: NSubrun. */
			fout << GetSubrunNum() << endl;

			/* Second line: column density. */
			for(int run=0; run<GetSubrunNum(); run++) 
				fout << ColumnDensity[run] << "\t";
			fout << endl;

			/* Third line: tritium purity. */
			for(int run=0; run<GetSubrunNum(); run++)
				fout << TritiumPurity[run] << "\t";
			fout << endl;

			/* Fourth line: voltage. */
			for(int run=0; run<GetSubrunNum(); run++)
				fout << Voltage[run] << "\t";
			fout << endl;

			/* 5 ~ 4+NSubrun: livetime array for each subrun. */
			for(int run=0; run<GetSubrunNum(); run++) {
				for(int npixel=0; npixel<NPixels; npixel++)
					fout << LiveTime[run][npixel] << "\t";
				fout << endl;
			}

			/* 5+NSubrun ~ 4+2NSubrun: efficiency array. */
			for(int run=0; run<GetSubrunNum(); run++) {
				for(int npixel=0; npixel<NPixels; npixel++)
					fout << Efficiency[run][npixel] << "\t";
				fout << endl;
			}

			/* 5+2NSubrun ~ 4+3NSubrun: event count. */
			for(int run=0; run<GetSubrunNum(); run++) {
				for(int npixel=0; npixel<NPixels; npixel++)
					fout << EventCount[run][npixel] << "\t";
				fout << endl;
			}

			fout.close();
		}

	public:
		double ColumnDensity[NSubrunMax];
		double TritiumPurity[NSubrunMax];
		double Voltage[NSubrunMax];
		double** LiveTime;
		double** Efficiency;
		int** EventCount;

};

class GlobalSimulation 
{
	public:
		GlobalSimulation() {
			char* KATRINpath = getenv("KATRIN");
			path = KATRINpath;
			path = path + "/Data";

			if(!GetDataFile()) {
				GetDataIdle();
				DumpData();
			}

		}

		~GlobalSimulation() {}

		double B_A[NPixels], B_S[NPixels], B_max[NPixels];
		double E[NPixels];

	private:
		KIRunSummaryDocument doc;
		string path;

		bool GetDataFile() {
			ifstream fin((path+"/field.dat").c_str());
			if(!fin.is_open()) return false;

			for(int npixel=0; npixel<NPixels; npixel++) fin >> B_A[npixel];
			for(int npixel=0; npixel<NPixels; npixel++) fin >> B_S[npixel];
			for(int npixel=0; npixel<NPixels; npixel++) fin >> B_max[npixel];
			for(int npixel=0; npixel<NPixels; npixel++) fin >> E[npixel];

			fin.close();
			return true;
		}

		void DumpData() {
			ofstream fout((path+"/field.dat").c_str());
			fout.precision(9);

			for(int npixel=0; npixel<NPixels; npixel++) 
				fout << B_A[npixel] << "\t";
			fout << endl;

			for(int npixel=0; npixel<NPixels; npixel++) 
				fout << B_S[npixel] << "\t";
			fout << endl;

			for(int npixel=0; npixel<NPixels; npixel++) 
				fout << B_max[npixel] << "\t";
			fout << endl;

			for(int npixel=0; npixel<NPixels; npixel++) 
				fout << E[npixel] << "\t";
			fout << endl;

			fout.close();
		}

		void GetDataIdle() {
			doc.AddInput("GlobalSimulation-FirstTritium-PeriodSummary_Dec2018-REF_18600V_6.0G-000002.ktf@FirstTritium.katrin");
			for(int pixel=0; pixel<NPixels; pixel++) {
				B_A[pixel] = doc["Values/AnalyzingPlane/MagneticField"][pixel].As<double>();
				E[pixel] = doc["Values/AnalyzingPlane/ElectricPotential"][pixel].As<double>();
				B_max[pixel] = doc["Values/Pinch/MagneticField"][pixel].As<double>();
				B_S[pixel] = doc["Values/PS2/MagneticField"][pixel].As<double>();
			}
		}


};

#endif
