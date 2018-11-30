
#include "../../../Simulation/Simulation.h"

Simulation sim;

int nvoltage = 41;
double voltage[41] = {18545, 18546, 18547, 18548, 18549, 18550, 18551, 18552, 18553, 18554, 18555, 18556, 18557, 18558, 18559, 18560, 18561, 18562, 18563, 18564, 18565, 18566, 18567, 18567.5, 18568, 18568.5, 18569, 18569.5, 18570, 18570.5, 18571, 18571.5, 18572, 18573, 18574, 18575, 18576, 18577, 18578, 18579, 18580};
double entries[41] = {2.43646e+07, 2.16883e+07, 1.92214e+07, 1.69654e+07, 1.48958e+07, 1.30186e+07, 1.13122e+07, 9.7757e+06, 8.38887e+06, 7.15235e+06, 6.0483e+06, 5.06647e+06, 4.20482e+06, 3.4482e+06, 2.79557e+06, 2.23058e+06, 1.75045e+06, 1.34218e+06, 1.00594e+06, 731422, 511962, 343556, 187994, 145930, 110963, 82851, 60366, 42573, 64968, 282718, 344623, 25968, 36450, 32697, 32850, 32895, 32657, 32847, 32814, 32911, 32565};
double Time[41] = {1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,1.0436e+06,902625,902625,902625,902625,902625,902625,1.96321e+06,1.21403e+07,2.06024e+07,2.00834e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06,3.27987e+06};

double* theory;
double* sample = new double[41];
double* error = new double[41];
void Draw()
{
	for(int i=0; i<41; i++) {
		sample[i] = entries[i]/Time[i] * 86400.;
		error[i] = sqrt(entries[i])/Time[i] * 86400.;
	}

	TGraph* graph = new TGraph();
	double ex[nvoltage];
	double ey[nvoltage];
	for(int i=0; i<41; i++) {
		graph->SetPoint(i, voltage[i], sample[i]);
		ex[i] = 0;
		ey[i] = error[i];

	}
	TGraphErrors* err = new TGraphErrors(nvoltage, voltage, sample, ex, ey);

	//graph->Draw("ap");
	//err->Draw("esame");

	double x[500];
	for(int i=0; i<500; i++) x[i] = 18545 + (18580-18545)*i/499.;
	double* y = sim.Asimov(0, 18574.97447, 500, x, 0.922336, 0.999765);
	TGraph* pred = new TGraph(500, x, y);
	pred->SetLineColor(kRed);
	//pred->Draw("lsame");

	/* Draw relative difference. */
	TGraph* rela_graph = new TGraph();
	double rela_ey[nvoltage];
	double rela_y[nvoltage];
	for(int i=0; i<41; i++) {
		rela_y[i] = (sample[i]-pred->Eval(voltage[i]))/error[i];
		rela_graph->SetPoint(i, voltage[i], (sample[i]-pred->Eval(voltage[i]))/error[i]);
		rela_ey[i] = 1;
	}
	TGraphErrors* relaerr = new TGraphErrors(nvoltage, voltage, rela_y, ex, rela_ey);
	rela_graph->SetMarkerStyle(8);
	rela_graph->SetMarkerSize(0.5);
	rela_graph->GetYaxis()->SetRangeUser(-5, 5);
	rela_graph->SetTitle("(True spectrum - predicted spectrum)/error");
	rela_graph->GetXaxis()->SetTitle("Energy [eV]");
	rela_graph->GetYaxis()->SetTitle("Sigma");
	rela_graph->Draw("ap");
	relaerr->Draw("esame");
}
