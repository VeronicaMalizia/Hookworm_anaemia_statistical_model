#################################################################
##
# Script 3: Plotting and inspecting the outcomes of Stan,
# after fitting the statistical model for hookworm anaemia, to the Ugandan dataset
##
# Author: Veronica Malizia
# Created: September 2020
#
# The script contains the code to access the posterior distributions obtained from Stan
# and produce plots to assess goodness of fit.
#
# Required input file in the Data directory:
# "Hb_eggscount_epg_adults_Uganda.csv", Ugandan dataset used for statistical analysis
# "Fit.posteriors.RData", the outcome of the model fit through Stan: it contains the 
# posterior distributions of the fitted parameters and the posterior predictive checks (PPC)
#
# Output file in Figures directory:
# "Fig 2" and "Fig 3" 
#################################################################


rm(list = ls())

library(tidyverse)
library(readr)
library(scales)
library(truncnorm)
library(bayesplot)
library(rstan)
library(patchwork)

#Load stan fit object
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set to current source directory
load("Data/Fit.posteriors.RData")
#load Uganda data set
data <- read.csv("Data/Hb_eggscount_epg_adults_Uganda.csv")

#Handy functions
ci <- function(X){
  ##computes the 95% confidence interval of a vector, removing NA's
  return(quantile(X, c(0.025, 0.975), na.rm=TRUE))
}
geom_mean <- function(x){
  #computes geometric mean returning zero, if the product is zero
  if(prod(x)==0)
    return(0)
  else
    return(exp(mean(log(x))))
}

#save se_mean and sd of y_predict
n <- nrow(data)
# Computing prevalences in PPC (mean and 95% prediction interval)
# Prevalence of any intensity of hookworm infection
any <- apply(bestfit$epg_predict, 1, function(x) length(which(x>0))/n)
mean(any) #54%
ci(any)
# Prevalence of moderate-to-heavy intensity of hookworm infection
mh <- apply(bestfit$epg_predict, 1, function(x) length(which(x>=2000))/n)
mean(mh) #4.5%
ci(mh)
# Prevalence of anaemia
an <- apply(bestfit$y_predict, 1, function(x) length(which(x<120))/n)
mean(an) #40%
ci(an)
# Prevalence of anaemia for the counterfactual scenario without hookworm infection
b.an <- plnorm(120, bestfit$log_mu0, bestfit$sigma)
mean(b.an)
ci(b.an)
# Prevalence of anaemia attributable to hookworm infection
mean(an - b.an)
ci(an - b.an)


# Producing Figure 2.
avg.fit.y <- apply(bestfit$y_predict, 2, mean)
avg.fit.epg <- apply(bestfit$epg_predict, 2, median)
ci.fit <- apply(bestfit$y_predict, 2, ci)
ci.fit <- t(ci.fit)
ci.fit.epg <- apply(bestfit$epg_predict, 2, ci)
ci.fit.epg <- t(ci.fit.epg)

data.to.plot <- cbind(data[, c("hb", "hookworm_epg")], 
                      avg.fit.epg, avg.fit.y, ci.fit,  
                      ci.fit.epg) 
colnames(data.to.plot) <- c("hb", "hookworm_epg", "estimated_epg", "Fitted_hb", "ci_low", "ci_hi", 
                            "ci_epg_low", "ci_epg_hi") 

f <- 
  data.to.plot %>%
  ggplot() +
  geom_point(aes(y=hb, x=hookworm_epg+3, colour="coral"), key_glyph = "rect") +
  geom_hline(yintercept = c(120), colour="darkgrey", linetype = "longdash") +
  scale_colour_discrete(name = " ") +
  scale_x_log10(name = "\n Mean egg count (epg) + 3", labels = math_format(.x)) +
  scale_y_continuous(name = "Hemoglobin concentration (g/L) \n",
                     breaks = 0 + 20*0:10,
                     limits = c(0, 200),
                     expand = c(0,0)) +
  expand_limits(x = 1) +
  theme(text = element_text(size = 12),
        panel.background = element_rect(fill = "white",
                                       colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#Plotting PPC to show the model fit to the data
  s <-20000 #the whole posterior distribution of parameters can be chosen, alternatively a random sample of it
  runs <- sample(20000, s)
  # Store PPC as a unique population of size "s"
  epg <- bestfit$epg_predict[runs,]
  hb <- bestfit$y_predict[runs,]
  # Cut the population in several categories of infection intensity
  # For each, computes geometric mean of epg, mean Hb and 95%CI of Hb
  df <- bind_cols(epg = as.vector(epg),
                  hb = as.vector(hb),
                  cat = cut(epg, breaks = c(0, 1, 50, 100, 500, 1000, 2000, 4000, 6000, 8000, 10000, 88888), right = F))
  df <- df %>%
    group_by(cat) %>%
    summarise(epg.cat = geom_mean(epg), 
              hb.cat = mean(hb),
              hb.lo = ci(hb)[1],
              hb.hi = ci(hb)[2])

  Fig2 <- f +
    geom_line(data=df, aes(x=epg.cat+3, y=hb.cat, colour="grey20"), key_glyph = "rect", size=0.8)+
    geom_line(data=df, aes(x=epg.cat+3, y=hb.lo), colour="grey", size=0.8)+
    geom_line(data=df, aes(x=epg.cat+3, y=hb.hi), colour="grey", size=0.8)+
    scale_colour_manual(name = " ",
                        values = c("coral", "grey20"),
                        labels = c("Data", "Model \n prediction"))

  #Saving Figure 2 as a tiff file  
  tiff("Figures/Fig 2.tif", width = 8, height = 6, units = 'in', res = 600, compression = "lzw")
  Fig2 
  dev.off()
  
  # Producing Figure 3: Comparison between observed data and model posterior predictive checks
  # Cumulative distribution of Hb concentrations
  hb <- ppc_ecdf_overlay(data$hb, bestfit$y_predict[runs,])+
    scale_colour_manual(name = " ",
                        values = c("blue4", "Light blue"),
                        labels = c("Data", "Predictions")) +
    scale_x_continuous(name = "\n Haemoglobin concentration (g/L)",
                       breaks = c(40*1:10)) +
    scale_y_continuous(name = "Cumulative Frequency \n",
                       breaks = c(0.2*0:5)) +
    labs(tag = "B")+
    theme_bw() +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
          legend.text=element_text(size=12))
  
# Cumulative distribution of mean egg counts (epg)
    epg <- ppc_ecdf_overlay(data$hookworm_epg+3, bestfit$epg_predict[runs,]+3, discrete=TRUE)+
    scale_x_log10(name = "\n Mean egg count (epg) + 3", 
                  #breaks = 4*c(1, 10, 100, 1000, 10000),
                  #labels = c(1, 10, 100, 1000, 10000)
                  ) +
    scale_colour_manual(name = " ",
                        values = c("blue4", "Light blue"),
                        labels = c("Data", "Predictions")) +
    scale_y_continuous(name = "Cumulative Frequency \n",
                       breaks = c(0.2*0:5)) +
    expand_limits(x = 1) +
    labs(tag = "A")+
    theme_bw()+
      theme(axis.text=element_text(size=12),axis.title=element_text(size=12), 
            legend.text=element_text(size=12))

#Saving Figure 2 as a tiff file 
  tiff("Figures/Fig 3.tif", width = 9, height = 6, units = 'in', res = 600, compression = "lzw")
  epg + hb + plot_layout(guides = "collect")
  dev.off()
  

