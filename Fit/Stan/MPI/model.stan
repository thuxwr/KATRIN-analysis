functions {
	real[] signal(real[] pars);
	real bkg();
}

data {
	int<lower=0> nbins;
	vector[nbins] rate;
	vector[nbins] error;
}

parameters {
	real<lower=0, upper=10> mass;
	real<lower=18564, upper=18584> endpoint;
	real<lower=-9000, upper=9000> A_sig;
	real<lower=0.8, upper=1.2> A_bkg;
	real<lower=2.63, upper=2.73> B_A;
	real<lower=24950, upper=25450> B_S;
	real<lower=41500, upper=42500> B_max;
	real<lower=0.199, upper=0.208> A1;
	real<lower=0.0541, upper=0.0571> A2;
	real<lower=1.75, upper=1.95> w1;
	real<lower=12, upper=13> w2;
	real<lower=14.2, upper=14.4> e2;
	real<lower=3.05, upper=3.75> sigma;
}

model {
	real e1 = 12.6;
	real pars[12] = {mass, endpoint, B_A, B_S, B_max, A1, A2, w1, w2, e1, e2, sigma*1e-18};
	real sig[nbins] = signal(pars);
	real pred[nbins];
	//vector[nbins+2] par_std;
	for(n in 1:nbins) {
		pred[n] = sig[n] * (A_sig*0.0001 + 1) + A_bkg * bkg();
		//target += (rate[n]-pred[n]) * (pred[n]-rate[n]) / error[n] / error[n];
		//par_std[n] = (rate[n] - pred[n])/error[n];
		//print("Nbin: ", n, " Pred: ", pred[n], " rate: ", rate[n], " error: ", error[n]);
	}
	rate ~ normal(pred, error);
	B_A ~ normal(2.68, 0.01);
	B_S ~ normal(2.52e4, 50);
	B_max ~ normal(4.2e4, 100);
	A1 ~ normal(0.204, 0.001);
	A2 ~ normal(0.0556, 0.0003);
	w1 ~ normal(1.85, 0.02);
	w2 ~ normal(12.5, 0.1);
	e2 ~ normal(14.3, 0.02);
	sigma ~ normal(3.4, 0.07);
	//target += (-1)*(B_A-2.68)*(B_A-2.68)/0.0001;
	//target += (-1)*(B_S-2.52e4)*(B_S-2.52e4)/2500;
	//target += (-1)*(B_max-4.2e4)*(B_max-4.2e4)/10000;
	//par_std[nbins+1] = (B_A - 2.68)/0.01;
	//par_std[nbins+2] = (B_S - 2.52e4)/50;
	//par_std ~ normal(0,1);
}

