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
	real<lower=0.8, upper=1.2> A_sig;
	real<lower=0.8, upper=1.2> A_bkg;
}

model {
	real pars[2] = {mass, endpoint};
	real sig[nbins] = signal(pars);
	real pred[nbins];
	for(n in 1:nbins) {
		pred[n] = A_sig * sig[n] + A_bkg * bkg();
	}
	rate ~ normal(pred, error);
}

