functions {
	real[] signal(real[] pars);
	real bkg();
}

data {
	int<lower=0> nbins;
	real rate[nbins];
	real error[nbins];
}

parameters {
	real<lower=0, upper=20> mass;
	real<lower=18564, upper=18584> endpoint;
	real<lower=-10000, upper=10000> A_sig;
	real<lower=0.8, upper=1.2> A_bkg;
	real<lower=2.5, upper=2.8> B_A;
	real<lower=24600, upper=26000> B_S;
	real<lower=4e4, upper=4.4e4> B_max;
}

model {
	real pars[5] = {mass, endpoint, B_A, B_S, B_max};
	real sig[nbins] = signal(pars);
	real pred[nbins];
	for(n in 1:nbins) {
		pred[n] = sig[n] * (A_sig*0.0001 + 1) + A_bkg * bkg();
	}
	rate ~ normal(pred, error);
	B_A ~ normal(2.68, 0.01);
	B_S ~ normal(2.52e4, 50);
	B_max ~ normal(4.2e4, 0.01e4);
}

