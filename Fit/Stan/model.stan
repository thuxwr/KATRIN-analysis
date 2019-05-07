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
//	real<lower=0, upper=10> mass;
	real<lower=18564, upper=18584> endpoint;
	real<lower=0.9, upper=1.1> A_sig;
	real<lower=0.9, upper=1.1> A_bkg;
//	real<lower=2.5, upper=2.8> B_A;
//	real<lower=24600, upper=26000> B_S;
	//real<lower=4e4, upper=4.4e4> B_max;
}

transformed parameters {
	real dendpoint = endpoint - 18574;
}

model {
real mass=0;
//real endpoint=18574;
	real pars[2] = {mass, endpoint};
	real sig[nbins] = signal(pars);
	real pred[nbins];
	vector[nbins] par_std;
	for(n in 1:nbins) {
		pred[n] = sig[n] * (A_sig) + A_bkg * bkg();
		par_std[n] = (rate[n] - pred[n])/error[n];
	}
	par_std ~ normal(0,1);
}

