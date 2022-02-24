#################################################################
##
# Script 1: Process of data selection for the statistical analysis
##
# Author: Veronica Malizia
# Created: September 2020, final update: 22/02/2022
#
# The script contains the code to load the raw individual data from Uganda, to tidy and filter the dataset for the statistical analysis.
# The population size after each filter application is reported. 
# 
# Required input file in Data directory:
# "Uganda worm individual data.csv", the individual data from Uganda
#
# Output file in Data directory:
# "Hb_eggscount_epg_adults_Uganda.csv", dataset to be used for statistical analysis
#################################################################

rm(list = ls())

library(tidyverse)
library(readr)
library(data.table)
library(foreign)
library(ggplot2)
library(scales)


#Data file should be in the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set to current source directory
base.dir <- getwd()
#Load the raw dataset
count.data.file <- file.path(base.dir, "Data/Uganda worm individual data.csv")
data <- read_csv(count.data.file) #N=2037

#1. Filter individuals with hookworm data available
#Remove individuals when all the four egg counts for hookworm infection are missing
no.hk <- which(is.na(data$hook_1a)&is.na(data$hook_1b)&is.na(data$hook_2a)&is.na(data$hook_2b))
hk.data <- data[-no.hk,] #N=1925

#2. Filter individuals with haemoglobin data available
hb.data <- hk.data %>%
  filter(!is.na(hb)) #N=1840

#3. Filter individuals aged>= 20 ys old
adults <- hb.data %>%
  filter(age >=20) #N=705


#Preparing data for the model's fitting procedure, through Stan software 
## Computing mean hookworm eggs per gram (epg)
## Assigning hookworm infection intensity
## Assessing anemia 
adults <- adults %>%
  mutate(hookworm_epg4 = round((hook_1a + hook_1b + hook_2a + hook_2b)/(0.0417*4))) %>%
  mutate(hookworm_inf = if_else(hookworm_epg4==0, 0, 1))%>%
  mutate(anemia_thr = case_when(age < 5 ~ 110,
                                age >= 5 & age <= 11 ~ 115,
                                age >= 12 & age <= 14 ~ 120,
                                age >= 15 & sex == 0 ~ 120, #women
                                age >= 15 & sex == 1 ~ 130)) %>% #men
  mutate(anemia_pos = if_else(hb < anemia_thr, 1, 0)) 

## For the individuals having one or more egg count unavailable (NA), the mean hookworm eggs per gram (epg) is computed over the available counts.
adults2 <- adults %>%
  rowwise() %>%
  mutate(counts = list(c(hook_1a, hook_1b, hook_2a, hook_2b))) %>%
  mutate(detected = length(which(!is.na(counts)))) %>%
  mutate(hookworm_epg = case_when(!is.na(hookworm_epg4) ~ hookworm_epg4,
                                 is.na(hookworm_epg4) ~ round(sum(hook_1a, hook_1b, hook_2a, hook_2b, na.rm = T)/(0.0417*detected))))
  #Clean from unnecessary columns
adults <- adults2 %>%
  select(-c(hookworm_epg4, counts, detected))

#Visualise data
ggplot(adults, aes(x=hookworm_epg, y=hb))+
  geom_point()+
  geom_vline(xintercept = c(10000), linetype="dashed", color="red", size=1)+
  scale_y_continuous(name="Haemoglobin concentrations",
                     limits = c(0, 200),
                     expand = c(0, 0)) +
  scale_x_log10()

## Individuals with very high mean epg (>= 10000 epg) are assumed to be infected with Ancylostoma duodenale hookworm species

##Population prevalences (fraction) before removing Ancylostoma infections 
#Hookworm infection prevalence
pos <- length(which(adults$hookworm_epg>0)) 
pos/nrow(adults) #0.53
#Anaemia prevalence
sum(adults$anemia_pos, na.rm = T)/nrow(adults) #0.42 
#Prevalence of anaemia in individuals detected with zero eggs for hookworm
no.infected <- filter(adults, hookworm_epg == 0)
sum(no.infected$anemia_pos, na.rm = T)/nrow(adults) #0.19

#Removing individuals assumed to be infected with Ancylostoma hookworm species
final_adults <- adults[-which(adults$hookworm_epg>=10000),] #N=695

##Population prevalences (fraction) after removing Ancylostoma infections and included in the manuscript
#Hookworm infection prevalence
pos <- length(which(final_adults$hookworm_epg>0)) 
pos/nrow(final_adults) #0.524
#Prevalence of moderate-to-heavy intensity of hookworm infection
length(which(final_adults$hookworm_epg>=2000))/nrow(final_adults) #0.040
#Anaemia prevalence
sum(final_adults$anemia_pos, na.rm = T)/nrow(final_adults) #0.414 
#Prevalence of anaemia in individuals detected with zero eggs for hookworm
no.infected <- filter(final_adults, hookworm_epg == 0)
sum(no.infected$anemia_pos, na.rm = T)/nrow(final_adults) #0.194

#Save data with only relevant columns that will be used for model's fitting through Stan software
#Columns of interest: hb, hookworm_epg, hook_1a, hook_1b, hook_2a, hook_2b
write_csv(select(final_adults, c(hb, hookworm_epg, starts_with("hook_"))), 
          path = "Data/Hb_eggscount_epg_adults_Uganda.csv")
