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
	real<lower=-10000, upper=10000> A_sig;
	real<lower=0.8, upper=1.2> A_bkg;
	real<lower=2.5, upper=2.8> B_A;
	real<lower=24600, upper=26000> B_S;
	//real<lower=4e4, upper=4.4e4> B_max;
}

model {
	real pars[4] = {mass, endpoint, B_A, B_S};
	real sig[nbins] = signal(pars);
	real pred[nbins];
	vector[nbins+2] par_std;
	for(n in 1:nbins) {
		pred[n] = sig[n] * (A_sig*0.0001 + 1) + A_bkg * bkg();
		par_std[n] = (rate[n] - pred[n])/error[n];
	}
	par_std[nbins+1] = (B_A - 2.68)/0.01;
	par_std[nbins+2] = (B_S - 2.52e4)/50;
	par_std ~ normal(0,1);
}

