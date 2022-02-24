data {
  int<lower=0> N0; //pop size with zero epg (namely group 0)
  int<lower=0> N0_4; //pop size with zero epg (group 0) & 4 slides
  int<lower=0> N1; //pop size with >0 epg (namely group 1)
  int<lower=0> N1_4; //pop size with >0 epg (group 1) & 4 slides
  int x1_4[N1_4, 4]; //egg counts from Uganda data - group 1 & 4 slides
  int x1_2[N1-N1_4, 2]; //egg counts from Uganda data - group 1 & 2 slides
  vector[N1] y1; //hb values from Uganda data - group 1
  vector[N0] y0; //hb values from Uganda data - group 0
  real alpha; //avg nr of eggs produced by each female worm
  real beta; //mean saturation level of infection
  real a; 
  real b;
  real c;
  real d;
  real e;
}
parameters {
	//all variables ending with "..1" refer to group 1, all variables ending with "..0" refer to group 0.
	vector<lower=0>[N1] w1; 
	vector<lower=0>[N0] w0; 
	vector<lower=0>[N1] ind_sui1;
	vector<lower=0>[N0] ind_sui0;
	real<upper=((5.19-a)/b)> eta; //the upper level is determined by the maximum plausible value the expected haemoglobin can take
	real<lower=0> sigma;
	real eta_w0;
	real<lower=0> phi1;
	real<lower=0> k_e;
	real<lower=0> shape; //shape and scale of the Weibull distribution for worm burdens in infected individuals
	real<lower=0> scale;
	real<lower=0, upper=1> rho;
}
transformed parameters{
	vector[N1] ec1; //effective egg counts per person - group 1
	vector[N0] ec0; //effective egg counts per person - group 0
	real log_mu0;
	real phi0; 

	//Parameters transformations
	log_mu0 = a + (b * eta);
	phi0 = c + (d * eta_w0); 

	//Functions
 	ec1 = (alpha*(w1/2)) ./ (1 + ((alpha*(w1/2)) ./ (beta * ind_sui1)));
	ec0 = (alpha*(w0/2)) ./ (1 + ((alpha*(w0/2)) ./ (beta * ind_sui0)));

}
model {
	//Declaring variables for the expected haemoglobin concentrations
	vector[N0] mu0;
	vector[N1] mu1;

	//Prior distributions for parameters shared across the entire population
	eta ~ std_normal();
	sigma ~ normal(0, 0.5);
	eta_w0 ~ std_normal();
	phi1 ~ normal(1, 3);
	k_e ~ normal(0, 2);
	rho ~ uniform(0, 1);
	shape ~ normal(0.5, 0.5); //it is likely between 0 - 1. Shape=1 reduces to an Exponential distribution
	scale ~ cauchy(2, 1); 

	//Group 1 
	target += weibull_lpdf(w1 | shape, scale);
	target += - weibull_lccdf(1 | shape, scale);
	ind_sui1 ~ gamma(e, e);
	
	mu1 = log_mu0 - log1p(exp(phi1*log(w1) - phi0));
	target += N1 * bernoulli_lpmf(0 | rho);
	y1 ~ lognormal(mu1, sigma);
	for(n in 1:N1_4){ //rows with 4 counts
		x1_4[n] ~ neg_binomial_2(ec1[n], k_e);
	}
	for(n in 1:(N1-N1_4)){ //rows with 2 counts
		x1_2[n] ~ neg_binomial_2(ec1[N1_4+n], k_e);
	}
	
	//Group 0 
	target += weibull_lpdf(w0 | shape, scale);
	target += - weibull_lccdf(1 | shape, scale);
	ind_sui0 ~ gamma(e, e);

	mu0 = log_mu0 - log1p(exp(phi1*log(w0) - phi0));
	
	for(n in 1:N0) {
		real LL0_wneg;
		real LL0_wpos;

		//Likelihood for mixture components with 0 and >=1 worms, conditional on
		//population-level parameters. N.B.: all individual-level parameters and
		//normalising constants must be included here.
		LL0_wneg = lognormal_lpdf(y0[n] | log_mu0, sigma);

		LL0_wpos = lognormal_lpdf(y0[n] | mu0[n], sigma);
		if(n <= N0_4) //rows with 4 counts
			LL0_wpos += 4 * neg_binomial_2_lpmf(0 | ec0[n], k_e); 
		if(n > N0_4) //rows with 2 counts
			LL0_wpos += 2 * neg_binomial_2_lpmf(0 | ec0[n], k_e);
	
    target += log_mix(rho, LL0_wneg, LL0_wpos);
	}
}
generated quantities {
	int x_predict[(N0+N1), 4]; //replicating egg counts
	real epg_predict[(N0+N1)]; //replicating averaged egg counts (epg)
	real y_predict[(N0+N1)]; //replicating haemoglobin concentrations
	real w_rep[(N0+N1)];
	real ec_rep[(N0+N1)];
	real indsui_rep[(N0+N1)];
	
	for(n in 1:(N0+N1)){
		w_rep[n] = weibull_rng(shape, scale);
		indsui_rep[n] = gamma_rng(e, e);
		ec_rep[n] = (alpha*(w_rep[n]/2)) ./ (1 + ((alpha*(w_rep[n]/2)) ./ (beta * indsui_rep[n])));
	}
	for(n in 1:(N0+N1)){
		if(bernoulli_rng(rho)){ //the condition would satisfy the event "zero worms"
			y_predict[n] = lognormal_rng(log_mu0, sigma);
			for(j in 1:4){
				x_predict[n, j] = 0; 
			}
		}
		else{ //>0 worms
			y_predict[n] = lognormal_rng((log_mu0 - log1p(exp(phi1*log(w_rep[n]) - phi0))), sigma);
			for(j in 1:4){
				x_predict[n, j] = neg_binomial_2_rng(ec_rep[n], k_e); 
			}
		}
	}
	for(n in 1:(N0+N1))
		epg_predict[n] = round((x_predict[n,1] + x_predict[n,2] + x_predict[n,3] + x_predict[n,4])/(0.0417*4));
}

