#################################################################
##
# Script 5: Plotting the prevalence of anaemia attributable to hookworm as a function 
# of infection prevalence. 
##
#
#Author: Veronica Malizia
#Created: July 2021
#
#The script contains the code to plot the results of simulations in script 4.
#This code will produce Fig 4 and Fig S1 of the manuscript.
#
#Required input file in the Data directory:
# "Simulations_results.RData", outcome of simulations, including prevalence measures for each endemicity scenario
#
# Output file:
# Fig 4 and S1 Fig in Figures folder
#################################################################

rm(list = ls())

library(tidyverse)
library(readr)
library(EnvStats)
library(patchwork)

#Loading simulations results from Data folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set to current source directory
load("Data/Simulations_results.RData")

#Filtering prevalence of any intensity of hookworm infection
any <- to.plot %>%
  filter(scheme == "hk.prev1" | scheme == "hk.prev2" | scheme == "hk.prev4") %>%
  mutate(int = "Any intensity prevalence") 
levels(any$scheme) <- c(rep("1 Kato-Katz slide", 2),
                        rep("2 Kato-Katz slides", 2),
                        rep("4 Kato-Katz slides", 2))

#Filtering prevalence of mod-to-heavy intensity of hookworm infection
mh <- to.plot %>%
  filter(scheme == "MH.prev1" | scheme == "MH.prev2" | scheme == "MH.prev4") %>%
  mutate(int = "M&H intensity prevalence")
levels(mh$scheme) <- c(rep("1 Kato-Katz slide", 2),
                       rep("2 Kato-Katz slides", 2),
                       rep("4 Kato-Katz slides", 2))

#Producing Fig 4 of the manuscript
p1 <- 
  any %>%
  filter(shape == 0.5) %>%
  #filter(scheme == "4 Kato-Katz slides") %>%
  ggplot(aes(x=hk.prev*100, y=EX.AN*100)) +
  geom_point(aes(color= as.factor(mean_fworms))) +
  geom_line(aes(color=as.factor(mean_fworms))) +
  geom_point(data = any[nrow(any),], aes(x=hk.prev*100, y=EX.AN*100), color = "black", shape = 17, size = 2) +
  geom_errorbar(aes(ymin=EX.AN.lo*100, ymax=EX.AN.hi*100), width=.01, alpha=0.4) +
  geom_errorbar(data = any[nrow(any),], aes(ymin=EX.AN.lo*100, ymax=EX.AN.hi*100), width=.01, alpha=0.4) +
  facet_wrap(~ scheme) +
  #tag_facets("rc", tag_levels  = c("A", "B", "C")) +
  labs(color = "Average worm load") +
  scale_x_continuous(name="\n Prevalence of any infection (%)",
                     breaks = c(20*0:5),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_y_continuous(name=" ",
                     breaks = c(5*0:20),
                     limits = c(0, 30),
                     expand = c(0, 0)) +
  expand_limits(y=0) +
  theme_bw() +
  labs(tag = "A") +
  theme(panel.spacing.x = unit(2, "lines"),
        text = element_text(size = 12),
        plot.tag.position = c(0.03, 1.03),
        plot.margin = unit(c(10,10,30,10), "pt"))

p2 <- 
  mh %>%
  filter(shape == 0.5) %>%
  #filter(scheme == "4 Kato-Katz slides") %>%
  ggplot(aes(x=hk.prev*100, y=EX.AN*100)) +
  geom_point(aes(color=as.factor(mean_fworms))) +
  geom_line(aes(color=as.factor(mean_fworms))) +
  geom_point(data = mh[nrow(mh),], aes(x=hk.prev*100, y=EX.AN*100), color = "black", shape = 17, size = 2) +
  geom_vline(xintercept=2, linetype="dashed") +
  geom_abline(slope=1, intercept = 0, linetype="dashed", color = "grey") +
  geom_errorbar(aes(ymin=EX.AN.lo*100, ymax=EX.AN.hi*100), width=.01, alpha=0.4) +
  geom_errorbar(data = mh[nrow(mh),], aes(ymin=EX.AN.lo*100, ymax=EX.AN.hi*100), width=.01, alpha=0.4) +
  facet_wrap(~ scheme) +
  labs(color = "Average worm load") +
  scale_x_continuous(name="\n Prevalence of M&HI infection (%)",
                     breaks = c(0, 2, 5*0:10),
                     limits = c(0, 55),
                     expand = c(0, 0)) +
  scale_y_continuous(name=" ",
                     breaks = c(5*0:40),
                     limits = c(0, 20),
                     expand = c(0, 0)) +
  coord_cartesian(x = c(0, 20)) +
  expand_limits(y=0) +
  theme_bw() +
  labs(tag = "B") +
  theme(panel.spacing.x = unit(2, "lines"),
        text = element_text(size = 12),
        plot.tag.position = c(0.03, 1.01),
        plot.margin = unit(c(0,0,10,10), "pt"))

ylab <- "Prevalence of anaemia due to hookworm (%)"
p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
grid::grid.draw(grid::textGrob(ylab, x = 0.012, rot = 90))

#Saving figure to a .tif file in Figures directory
tiff("Figures/Fig 4.tif", width = 7.5, height = 7, units = 'in', res = 600, compression = "lzw")
p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
grid::grid.draw(grid::textGrob(ylab, x = 0.012, rot = 90))

dev.off()

#Producing S1 Fig of the manuscript
p1 <- 
  any %>%
  filter(shape %in% c(0.3, 0.5, 0.7)) %>%
  ggplot(aes(x=hk.prev*100, y=EX.AN*100)) +
  geom_point(aes(color= as.factor(mean_fworms), shape = as.factor(shape))) +
  geom_line(aes(group=interaction(mean_fworms, shape), color = as.factor(mean_fworms)), alpha=0.5) +
  geom_errorbar(aes(ymin=EX.AN.lo*100, ymax=EX.AN.hi*100), width=.01, alpha=0.4) +
  facet_wrap(~ scheme) +
  labs(color = "Average worm load",
       shape = "Shape parameter of \n worms distribution",
       tag = "A") +
  #guides(color = F) +
  scale_x_continuous(name="\n Prevalence of any infection (%)",
                     breaks = c(20*0:5),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_y_continuous(name=" ",
                     breaks = c(5*0:20),
                     limits = c(0, 40),
                     expand = c(0, 0)) +
  expand_limits(y=0) +
  theme_bw() +
  theme(panel.spacing.x = unit(2, "lines"),
        text = element_text(size = 12),
        plot.tag.position = c(0.03, 1.03),
        plot.margin = unit(c(10,10,30,10), "pt"))

p2 <- 
  mh %>%
  filter(shape %in% c(0.3, 0.5, 0.7)) %>%
  ggplot(aes(x=hk.prev*100, y=EX.AN*100)) +
  geom_point(aes(color= as.factor(mean_fworms), shape = as.factor(shape))) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", color = "grey") +
  geom_line(aes(group=interaction(mean_fworms, shape), color = as.factor(mean_fworms)), alpha=0.5) +
  geom_vline(xintercept=2, linetype="dashed") +
  geom_errorbar(aes(ymin=EX.AN.lo*100, ymax=EX.AN.hi*100), width=.01, alpha=0.4) +
  facet_wrap(~ scheme) +
  labs(shape = "Shape parameter of \n worms distribution",
       color = "Average worm load",
       tag =  "B") +
  scale_x_continuous(name="\n Prevalence of M&HI infection (%)",
                     breaks = c(0, 2, 5*0:10),
                     limits = c(0, 55),
                     expand = c(0, 0)) +
  scale_y_continuous(name=" ",
                     breaks = c(5*0:40),
                     limits = c(0, 20),
                     expand = c(0, 0)) +
  coord_cartesian(x = c(0, 20)) +
  expand_limits(y=0) +
  theme_bw() +
  theme(panel.spacing.x = unit(2, "lines"),
        text = element_text(size = 12),
        plot.tag.position = c(0.03, 1.03),
        plot.margin = unit(c(0,0,10,10), "pt"))

ylab <- "Prevalence of anaemia due to hookworm (%)"
p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical")
grid::grid.draw(grid::textGrob(ylab, x = 0.012, rot = 90))

#Saving figure to a .tif file in Figures directory
tiff("Figures/S1 Fig.tif", width = 7.5, height = 7, units = 'in', res = 600, compression = "lzw")
p1 / p2 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.box = "vertical")
grid::grid.draw(grid::textGrob(ylab, x = 0.012, rot = 90))

dev.off()


