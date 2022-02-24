#################################################################
##
# Script 4: Estimating the prevalence of anaemia attributable to hookworm as a function 
# of infection prevalence. 
##
#
#Author: Veronica Malizia
#Created: July 2021
#
#The script contains the code to load the outcomes from Stan (i.e. Bayesian posterior distributions
#of model's parameters, obtained by the model model fit to the Uganda data) and run the model to
#simulate different endimicity scenarios. 
#
#Required input file in the Data directory:
# "Fit.posteriors.RData", the outcome of the model fit through Stan: it contains the 
# posterior distributions of the fitted parameters and the posterior predictive checks (PPC)
#
# Output file in Data directory:
# "Simulations_results.RData", outcome of simulations, including prevalence measures for each endemicity scenario
#################################################################

rm(list = ls())

library(tidyverse)
library(readr)
library(scales)
library(truncnorm)
library(rstan)
library(EnvStats)
library(foreach)
library(doParallel)

#Handy functions
ci <- function(X){
  ##computes the 95% confidence interval of a vector, removing NA's
  return(quantile(X, c(0.025, 0.975), na.rm=TRUE))
}

#Load input file from Data folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set to current source directory
load("Data/Fit.posteriors.RData")

#Set parameters for simulations
N <- 1000 #Simulate 1000 individuals per each scenario
beta0 = 10000/24 #Mean saturation level for egg production

#Defining endemicity scenarios
shapes = c(0.3, 0.5, 0.7) #Varying values for the shape parameter of the distribution of worm burdens in infected individuals
mean_worms = c(1, 3*1:4) #Varying values for the average worm burden in the population
rhos = seq(0.1, 0.9, 0.2) #Varying values for the rho parameter defining the proportion of individuals having zero worms
scenarios <- expand.grid(shape = shapes, rho=rhos, mean_fworms = mean_worms)

scenarios <- scenarios %>%
  mutate(scale = mean_fworms/(rho*gamma(1+1/shape))) #scale parameter of the distribution of worm burdens in infected individuals
        
#Add the scenario resulting from the fit and approximating the Uganda dataset. 
#(For plotting purposes)
scenarios <- rbind(scenarios, c(mean(bestfit$shape),
                                mean(bestfit$rho),
                                mean(bestfit$rho)*mean(bestfit$scale)*gamma(1+1/mean(bestfit$shape)),
                                mean(bestfit$scale)))


#Defining functions for:
#computing expected egg counts given x=worm load
compute_ec <- function(x, a, b){
  ind_sui <- rgamma(1, 50, 50)
  beta <- b*ind_sui
  return((a*x) / (1 + ((a*x) /beta)))
}
#Generating worm burdens and egg counts distributions for each scenario
w <- array(dim = c(nrow(scenarios), N))
ec <- array(dim = c(nrow(scenarios), N))
for(k in 1:nrow(scenarios)){
  for(j in 1:N){
    tmp <- rweibull(1, shape=scenarios$shape[k], scale=scenarios$scale[k])*(1 - rbernoulli(1, scenarios$rho[k]))
    w[k, j] <- tmp
    ec[k, j] <- compute_ec(tmp/2, 200/24, beta0)
  }
}

#Starting simulations
s <- 5000 #sample size from the posterior

writeLines(c(""), "Sink.txt") #initiate log file

#Launching simulations. It takes approximately 
cluster <- makeCluster(min(parallel::detectCores(logical = FALSE), nrow(scenarios)))
registerDoParallel(cluster)

results <- foreach(k = 1:nrow(scenarios), #k identifies the scenario
                   .inorder = TRUE,
                   .errorhandling = "remove",
                   #.combine = bind_rows,
                   .packages = c("tidyverse")) %dopar% {
                     #for each scenario, take a random sample of size s from the posterior:
                     index <- sample(20000, s)
                     
                     hk.prev1 <- vector("numeric", s)
                     MH.prev1 <- vector("numeric", s)
                     hk.prev2 <- vector("numeric", s)
                     MH.prev2 <- vector("numeric", s)
                     hk.prev4 <- vector("numeric", s)
                     MH.prev4 <- vector("numeric", s)
                     an.prev <- vector("numeric", s)
                     b.an.prev <- vector("numeric", s)
                     
                     for(n in 1:s){
                       sink("Sink.txt", append=TRUE)
                       cat(paste(Sys.time(), ": Starting scenario", k, "post draw", n, "\n", sep = " "))
                       sink()
                       
                       i <- index[n] #i identifies the draw from the posterior
                       #extracting i-th row of the posterior distribution
                       theta <- list(k_e=bestfit$k_e[i],
                                     log_mu0=bestfit$log_mu0[i],
                                     sigma=bestfit$sigma[i],
                                     phi1=bestfit$phi1[i],
                                     phi0=bestfit$phi0[i])
                       
                       #Simulating observed egg counts & epg
                       #j identifies the individuals
                       #k identifies the scenarios
                       epg1 <- c()
                       epg2 <- c()
                       epg4 <- c()
                       reduction <- c()
                       for(j in 1:N){
                           eggs <- rnbinom(4, mu = ec[k, j], size = theta$k_e)
                           epg1 <-  c(epg1, round(mean(eggs[1])/0.0417))
                           epg2 <-  c(epg2, round(mean(eggs[1:2])/0.0417))
                           epg4 <-  c(epg4, round(mean(eggs)/0.0417))
                           if(w[k, j]==0)
                             reduction <- c(reduction, 0)
                           else
                             reduction <- c(reduction,
                                            log(1 + exp(theta$phi1*log(w[k, j]) - theta$phi0)))
                       }
                       #Saving infection prevalence
                       hk.prev1[n] <- length(which(epg1>0))/N
                       MH.prev1[n] <- length(which(epg1>2000))/N
                       hk.prev2[n] <- length(which(epg2>0))/N
                       MH.prev2[n] <- length(which(epg2>2000))/N
                       hk.prev4[n] <- length(which(epg4>0))/N
                       MH.prev4[n] <- length(which(epg4>2000))/N
                       
                       #Simulating hb concentrations in the population
                       hb <- rlnorm(N, (theta$log_mu0 - reduction), theta$sigma)
                       hb0 <- hb * exp(reduction)
                       #Saving anaemia prevalence
                       an.prev[n] <- length(which(hb<120))/N
                       b.an.prev[n] <- length(which(hb0<120))/N
                     }
                     #Collate results
                     list(hk.prev1=hk.prev1,
                          MH.prev1=MH.prev1,
                          hk.prev2=hk.prev2,
                          MH.prev2=MH.prev2,
                          hk.prev4=hk.prev4,
                          MH.prev4=MH.prev4,
                          an.prev=an.prev,
                          b.an.prev=b.an.prev)
                   }

stopCluster(cluster)

#Inspecting results
#Computing summary statistics (mean and BCI) of prevalences, for each endemicity scenario
#Hookworm infection prevalences
mean.prev <- lapply(results, function(x) sapply(x, mean))
mean.prev <- do.call(bind_rows, mean.prev)
#Prevalence of anaemia attributable to hookworm
excess.an <- lapply(results, function(x) x$an.prev - x$b.an.prev)
EX.AN <- sapply(excess.an, mean)
EX.AN.ci <- sapply(excess.an, ci)

#Preparing results for plotting
to.plot <- bind_cols(scenarios[1:4], mean.prev, EX.AN=EX.AN,
                     EX.AN.lo=EX.AN.ci[1,], EX.AN.hi=EX.AN.ci[2,])


to.plot <- to.plot %>%
  gather(scheme, hk.prev, hk.prev1:MH.prev4, factor_key = T)

#Saving simulations results to Data folder
save(to.plot, file = "Data/Simulations_results.RData")


