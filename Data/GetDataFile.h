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
#include <string>
#include "KIIdle.h"
#include "KTree.h"
#include "TMath.h"

using namespace std;
using namespace katrin;
using namespace TMath;

class Data 
{
	public:
		Data() {
			LiveTime = new double*[NSubrunMax];
			Efficiency = new double*[NSubrunMax];
			EventCount = new int*[NSubrunMax];

			SetDataset("FirstTritium.katrin");
			auto dataset = idle.FindDataSet(DataSet);
			KTree filter = (make_tree
					("Limit", 2000)
					("FileName", "RunSummary3c")
			);
			dataset->GetFileList(filter, filelist);
			NRuns = filelist.NumberOfRows();
			GetData();
		}

		~Data(){}

		void SetDataset(string datasetname) {
			DataSet = datasetname;
		}

		void GetData() {
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
					}
					catch(KException& e) {
						Voltage[subrun] = -1;
					}

					/* Tritium purity. */
					try {
						TritiumPurity[subrun] = subRunDoc["Values/LARA/TritiumPurity"].As<double>();
					}
					catch(KException& e) {
						TritiumPurity[subrun] = -1;
					}

					/* Column density. */
					try {
						ColumnDensity[subrun] = subRunDoc["Values/ColumnDensity/ColumnDensity"].As<double>();
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

		int GetSubrunNum() { return NSubrun; }

	private:
		KIIdle idle;
		KTabree filelist;
		string DataSet;
		int NRuns;
		int NSubrun;

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
			doc.AddInput("GlobalSimulation-FirstTritium-PeriodSummary_Dec2018-REF_18600V_6.0G-000002.ktf@FirstTritium.katrin");
			for(int pixel=0; pixel<NPixels; pixel++) {
				B_A[pixel] = doc["Values/AnalyzingPlane/MagneticField"][pixel].As<double>();
				E[pixel] = doc["Values/AnalyzingPlane/ElectricPotential"][pixel].As<double>();
				B_max[pixel] = doc["Values/Pinch/MagneticField"][pixel].As<double>();
				B_S[pixel] = doc["Values/PS2/MagneticField"][pixel].As<double>();
			}
		}

		~GlobalSimulation() {}

		double B_A[NPixels], B_S[NPixels], B_max[NPixels];
		double E[NPixels];

	private:
		KIRunSummaryDocument doc;

};

#endif
