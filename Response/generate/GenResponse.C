/*
	 Rewrite GenResponse.py in cpp to make it run faster.

	 Weiran, Nov.13, 2018.
*/

#include "TMath.h"
#include <string>

double B_A = 2.68;
double B_S = 2.52e4;
double B_max = 4.2e4;
double m_e = 0.511e6;
double E0 = 18.574e3;
double pi = TMath::Pi();

TH1D* density;
TH1D* scat[3];
double calib;
double Ncalib = 5e17;
double tmpz = -1;
double integ = 0;
double sigma = 3.4e-18;

double P(int s, double z, double theta);
double Neff(double z, double theta);
double gamma(double energy);
double transmission(double E, double theta, double U);
double Response(double E, double U);
double GetPdf(int ScatterTimes, double epsilon);

void GenResponse()
{
	char* KATpath = getenv("KATRIN");
	string path = KATpath;
	TFile* densfile = new TFile((path+"/Scattering/density/density.root").c_str(), "READ");
	density = (TH1D*)densfile->Get("density");
	calib = density->Integral();

	TFile* EnergyLossFile = new TFile((path+"/Scattering/energyloss/EnergyLoss.root").c_str(), "READ");
	scat[0] = (TH1D*)EnergyLossFile->Get("scat0");
	scat[1] = (TH1D*)EnergyLossFile->Get("scat1");
	scat[2] = (TH1D*)EnergyLossFile->Get("scat2");

	/* Test. */
	TH2D* hist = new TH2D("","",100,E0-35,E0+5,100,E0-35,E0+5);
	for(int xbin=1; xbin<=100; xbin++) for(int ybin=1; ybin<=100; ybin++) {
		double E = hist->GetXaxis()->GetBinCenter(xbin);
		double U = hist->GetYaxis()->GetBinCenter(ybin);
		hist->SetBinContent(xbin, ybin, Response(E, U));
		cout << "finish xbin: " << xbin << ", ybin: " << ybin << endl;
	}
	hist->SetName("haha");
	hist->SaveAs("test.root");

}

double P(int s, double z, double theta) {
	double avg = Neff(z, theta) * sigma;
	return pow(avg, s) / TMath::Factorial(s) * exp(-1.*avg);
}

double Neff(double z, double theta) {
	if(z!=tmpz) {
		int low = density->FindBin(z);
		int up = density->GetNbinsX();
		tmpz = z;
		integ = density->Integral(low, up);
	}
	return 1./cos(theta) * integ / calib * Ncalib;
}

double gamma(double energy) {
	return energy / m_e + 1;
}

double transmission(double E, double theta, double U) {
	double thres = E*(1 - pow(sin(theta), 2)*B_A/B_S*(gamma(E)+1)/2.);
	if(thres > U) return 1;
	else return 0;
}
	
double GetPdf(int ScatterTimes, double epsilon) {
	return scat[ScatterTimes-1]->GetBinContent(scat[ScatterTimes-1]->FindBin(epsilon));
}

double Response(double E, double U) {
	if(E<=U) return 0;
	double thetamax = asin(sqrt(B_S/B_max));

	double response = 0;
	TH2D* histtmp = new TH2D("","",100,0,thetamax,500,0,E-U);
	for(int xbin=1; xbin<=100; xbin++) {
		double theta = histtmp->GetXaxis()->GetBinCenter(xbin);
		double Prob[3];
		for(int nscat=1; nscat<4; nscat++) Prob[nscat-1] = P(nscat, 2, theta);
		for(int ybin=1; ybin<=500; ybin++) {
			double epsilon = histtmp->GetYaxis()->GetBinCenter(ybin);
			double EnergyLossPdf = 0;
			for(int nscat=1; nscat<4; nscat++) EnergyLossPdf += Prob[nscat-1]*GetPdf(nscat, epsilon);
			double content = transmission(E-epsilon,theta,U) * sin(theta) * EnergyLossPdf;
			histtmp->SetBinContent(xbin, ybin, content);
		}
	}
	response += histtmp->Integral("width");
	delete histtmp;

	TH1D* histtmp2 = new TH1D("","",100,0,thetamax);
	for(int xbin=1; xbin<=100; xbin++) {
		double theta = histtmp2->GetXaxis()->GetBinCenter(xbin);
		histtmp2->SetBinContent(xbin, transmission(E, theta, U) * sin(theta) * P(0,2,theta));
	}
	response += histtmp2->Integral("width");
	delete histtmp2;

	return response;
}

