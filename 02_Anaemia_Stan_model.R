#################################################################
##
# Script 2: Fitting the statistical model for hookworm anaemia, to the Ugandan dataset
##
# Author: Veronica Malizia
# Created: September 2020
#
# The script contains the code to calibrate the statical anemia model using Stan software.
# The users may need to install Stan (version 2.19.1), the R package "rstan" and RTools.
# The script will call Stan code in the current repository, under the name "Stanmodel_code.stan"
#
# Required input file in Data directory:
# "Hb_eggscount_epg_adults_Uganda.csv", dataset to be used for statistical analysis
#
# Output file in Data directory:
# "Fit.posteriors.RData", the outcome of the model fit through Stan: it contains the 
# posterior distributions of the fitted parameters and the posterior predictive checks (PPC)
#################################################################

rm(list = ls())

library(tidyverse)
library(readr)
library(scales)
library(bayesplot)
library(EnvStats)

#Check if RTools are correctly installed
pkgbuild::has_build_tools(debug = TRUE)
pkgbuild::find_rtools(debug = TRUE)

library(rstan)
#library(devtools)
#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Handy function
ci <- function(X){
  ##computes the 95% confidence interval of a vector, removing NA's
  return(quantile(X, c(0.025, 0.975), na.rm=TRUE))
}

#Loading Uganda individual data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set to current source directory
data <- read.csv("Data/Hb_eggscount_epg_adults_Uganda.csv")

#Create data for running the model in Stan 
# Splitting individuals detected with positive and with zero eggs for hookworm infection
#For each of these groups, individuals are ordered to that individuals with 4 egg counts are available come first
zeros <- which(data$hookworm_epg==0)
N0 <- length(zeros)

pos <- data[-zeros, ]
pos.four.counts <- which(!is.na(pos$hook_1a)&!is.na(pos$hook_1b)&!is.na(pos$hook_2a)&!is.na(pos$hook_2b))
N1.4 <- length(pos.four.counts) # the first N1.4 rows of infected individuals have 4 counts available
pos <- bind_rows(pos[pos.four.counts,],
                 pos[-pos.four.counts,])
#Filtering egg counts only (4 KK slides or 2 based on availability)
x1_4 = pos[1:N1.4, (ncol(pos)-3):ncol(pos)]
tmp <- pos[(N1.4+1):nrow(pos), (ncol(pos)-3):ncol(pos)]
res <- c()
for(i in 1:nrow(tmp)){
  if(is.na(tmp$hook_1a[i])|is.na(tmp$hook_1b[i]))
    res <- rbind(res, c(tmp$hook_2a[i], tmp$hook_2b[i]))
  if(is.na(tmp$hook_2a[i])|is.na(tmp$hook_2b[i]))
    res <- rbind(res, c(tmp$hook_1a[i], tmp$hook_1b[i]))  
}
colnames(res) <- c("hook_1a", "hook_1b")
x1_2 = as.data.frame(res)

neg <- data[zeros, ] 
neg.four.counts <- which(!is.na(neg$hook_1a)&!is.na(neg$hook_1b)&!is.na(neg$hook_2a)&!is.na(neg$hook_2b))
N0.4 <- length(neg.four.counts) # the first N0.4 rows of uninfected individuals have 4 counts available
neg <- bind_rows(neg[neg.four.counts,],
                 neg[-neg.four.counts,])

###################

n <- nrow(data) #population size
N1 = n-N0 #positive tested population size
beta = 10000/24 #Mean saturation level of the egg production function (beta)

# Ingredients for Stan 
# The list below will be read as input in the first block "data" of the stan code 
data.stan.mixture  <- list(N0 = N0, #pop size with zero epg (namely group 0)
                           N0_4 = N0.4, #pop size with zero epg (group 0) & 4 Kato-Katz slides available
                           N1 = N1, #pop size with >0 epg (namely group 1)
                           N1_4 = N1.4, #pop size with >0 epg (group 1) & 4 Kato-Katz slides available
                           x1_4 = x1_4, #egg counts from Uganda data - group 1 with 4 KK slides
                           x1_2 = x1_2, #egg counts from Uganda data - group 1 with 2 KK slides
                           y1 = pos$hb, #hb values from Uganda data - group 1
                           y0 = neg$hb, #hb values from Uganda data - group 0
                           alpha = 200/24, #avg nr of eggs produced by each female worm
                           beta = beta, #mean saturation level for the egg production function
                           a = 4.86, #mean of log_mu0 prior distribution
                           b = 0.1,  #standard deviation of log_mu0 prior distribution
                           c = 6, #mean of phi0 prior distribution
                           d = 2, #standard deviation of phi0 prior distribution
                           e = 50 #shape of individual susceptibility to infection (ind_sus) 
)

# Compiling the code with Stan
# Make sure to remove any existing already compiled file
unlink("Stanmodel_code.rds")
sm <- stan_model("Stanmodel_code.stan")

# Launching Stan to fit the model to the Ugandan data and produce Bayesian posterior estimates of parameters 
# Launchin 4 Markov chains with 10000 iterations per each. 
# 5000 runs are used for warm-up (default)
start.time <- Sys.time()
fit.mixt <- sampling(sm, data = data.stan.mixture, 
                  chains = 4, iter = 10000, cores = 4, #refresh = 0, 
                  #init = initial, 
                  control = list(adapt_delta = 0.9,
                                 max_treedepth = 15))
end.time <- Sys.time()
end.time - start.time #It takes approximately 30 min

# Inspecting Stan results
# List of parameters we want to get posterior estimations
# Names refer to the stan code file
# These will make Table 2. of the manuscript
par <- c("log_mu0", "sigma", "phi0", "phi1",
         "shape", "scale", "k_e", "rho") 
# Print an overview of Bayesian posterior distributions and measures of Markov chain convergence
print(fit.mixt, probs = c(0.025, 0.5, 0.975), digits=3,
      pars = par)  

# Inspect the traceplots of 4 Markov chains for each parameter
traceplot(fit.mixt, pars=c("eta", "sigma",  "eta_w0", "phi1", "k_e", "shape", "scale", "rho", "lp__"), nrow=3) #, inc_warmup = T)
# Inspect the pairs plot
pairs(fit.mixt, pars = c("eta", "sigma",  "eta_w0", "phi", "k_e", "shape", "scale", "rho"),
      labels = c("eta1", "sigma",  "eta2", "phi", "k_e", "delta", "lambda")) 


## Diagnostic checks of the fit. 
#For more information look at https://mc-stan.org/rstan/reference/stan_plot_diagnostics.html
#Look at the rank plot to see how the chains differ from each other
mcmc_rank_hist(fit.mixt, pars=par)
#Look at the effective sample size
stan_ess(fit.mixt, chain=1)
#Look at the convergence parameter Rhat
stan_rhat(fit.mixt) #Gelman recommends that Rhat for each parameter be less than 1.1

diag <- rstan::get_sampler_params(fit.mixt) %>% 
       set_names(1:4) %>% 
       map_df(as_data_frame,.id = 'chain') %>% 
       group_by(chain) %>% 
       mutate(iteration = 1:length(chain)) %>% 
       mutate(warmup = iteration <= 5000)

# Check the treedepth
# Treedepth increased to 15 for better performance
diag %>%
    ggplot(aes(iteration, treedepth__, color = chain)) + 
    geom_line() + 
    geom_hline(aes(yintercept = 15), color = 'red') +
  coord_cartesian(ylim = c(0, 15))

#Check the stepsize
diag %>%
  ggplot(aes(iteration, stepsize__, color = chain)) + 
  geom_line() +
  coord_cartesian(ylim = c(0, 1))

## Saving output: 
# a data frame with the posterior distributions of the parameters listed in "par"
# the whole data frame would produce a file of 1.2 Gb
# Therefore only columns relevant for plotting and inspections are saved
# Those are listed below as "to.save"
to.save <- c("log_mu0", "sigma", "phi0", "phi1", "shape", "scale", "k_e", "rho",
             "y_predict", "epg_predict") 
bestfit <- extract(fit.mixt, pars = to.save, permuted=T) #by default warmup are not included
save(bestfit, file = "Data/Fit.posteriors.RData")
