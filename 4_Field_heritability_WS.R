

#Load packages

library(tidyverse)
library(broom)
library(lme4)



#Load data file 

Parent_offspringWS <- read.csv("Cage_data.csv", header = TRUE, sep = ',')



#Generate values needed to calculate standardised WS values for offspring


Site_meanWS_offspring <- Parent_offspringWS %>%
  group_by(Site) %>%
  summarise(Site_meanWS_offspring = mean(WS_offspring_mm, na.rm = T))

Parent_offspringWS <- merge(Parent_offspringWS, Site_meanWS_offspring, by = "Site")

Gradient_meanWS_offspring <- Parent_offspringWS %>%
  group_by(Gradient) %>%
  summarise(Gradient_meanWS_offspring = mean(WS_offspring_mm, na.rm = T))

Parent_offspringWS <- merge(Parent_offspringWS, Gradient_meanWS_offspring, by = "Gradient")

Site_SDWS_offspring <- Parent_offspringWS %>%
  group_by(Site) %>%
  summarise(Site_SDWS_offspring = sd(WS_offspring_mm, na.rm = T))

Parent_offspringWS <- merge(Parent_offspringWS, Site_SDWS_offspring, by = "Site")

Gradient_SDWS_offspring <- Parent_offspringWS %>%
  group_by(Gradient) %>%
  summarise(Gradient_SDWS_offspring = sd(WS_offspring_mm, na.rm = T))

Parent_offspringWS <- merge(Parent_offspringWS, Gradient_SDWS_offspring, by = "Gradient")

#Calculate standardised wing size for offspring

Parent_offspringWS$Z_WS_offspring_Site <- (Parent_offspringWS$WS_offspring_mm - Parent_offspringWS$Site_meanWS_offspring)/Parent_offspringWS$Site_SDWS_offspring

Parent_offspringWS$Z_WS_offspring_Gradient <- (Parent_offspringWS$WS_offspring_mm - Parent_offspringWS$Gradient_meanWS_offspring)/Parent_offspringWS$Gradient_SDWS_offspring

#Exclude sites 71 and 72 due to lack of data

Parent_offspringWS <- Parent_offspringWS %>%
  filter(Site!=71,
         Site!=72)



#Estimate field heritability of body size by running regressions of standardised offspring mean WS against standardised mean WS of lab-reared flies (before selection) at each site and overall for each gradient

by_site <- group_by(Parent_offspringWS, Site)


#Summary of models for each site

WS_fieldH2_model_summaries <- do(by_site,
                                 glance(
                                   lm(Z_WS_offspring_Site ~ Z_MeanWS, data = .)))


#Summary of estimated effects for each model term

WS_fieldH2_effects <- do(by_site,
                         tidy(
                           lm(Z_WS_offspring_Site ~ Z_MeanWS, data = .)))

#Need to add gradient and elevation to model effects and summary data frames
Site_details <- read.csv("Site_details.csv", header = TRUE, sep = ',')

WS_fieldH2_effects <- merge(WS_fieldH2_effects, Site_details, by = "Site")
WS_fieldH2_model_summaries <- merge(WS_fieldH2_model_summaries, Site_details, by = "Site")

WS_fieldH2 <- subset(WS_fieldH2_effects, term=="Z_MeanWS")


#Calculate heritability as twice the slope of the regression line ('estimate' term)

WS_fieldH2$Field_H2 <- 2*abs(WS_fieldH2$estimate)

#Cap heritability at 1 if estimate is higher than this

WS_fieldH2$Field_H2 <- ifelse(WS_fieldH2$Field_H2 > 1, 1, WS_fieldH2$Field_H2)

#Calculate p values adjusted for multiple comparisons for each gradient
WS_fieldH2 <- WS_fieldH2 %>%
  group_by(Gradient) %>%
  mutate(P_FDR = p.adjust(p.value, method = "fdr"))

#Add number of observations from model summary to WS_fieldH2 to get number of cages at each site with body size estimates for both parents and offspring

WS_fieldH2$N <- WS_fieldH2_model_summaries$nobs


#Calculate heritability overall for each gradient

#Include Site as a random effect.  Method above doesn't work for mixed models so separate by gradient then run models


Parent_offspringWS_Danbulla <- subset(Parent_offspringWS, Gradient=="Danbulla")
Parent_offspringWS_Paluma <- subset(Parent_offspringWS, Gradient=="Paluma")


WS_fieldH2_Danbulla <- lmer(Z_WS_offspring_Gradient ~ Z_MeanWS_Gradient + (1|Site), data = Parent_offspringWS_Danbulla)
WS_fieldH2_Danbulla_summary <- summary(WS_fieldH2_Danbulla)[["coefficients"]]


WS_fieldH2_Paluma <- lmer(Z_WS_offspring_Gradient ~ Z_MeanWS_Gradient + (1|Site), data = Parent_offspringWS_Paluma)
WS_fieldH2_Paluma_summary <- summary(WS_fieldH2_Paluma)[["coefficients"]]


#Test for a linear relationship between WS field heritability and elevation at each gradient, weighted by number of cages


by_gradient <- group_by(WS_fieldH2, Gradient)


#Summary of models for each site

Elevation_WS_fieldH2_model_summaries <- do(by_gradient,
                                           glance(
                                             lm(Field_H2 ~ Elevation, weights=N, data = .)))


#Summary of estimated effects for each model term

Elevation_WS_fieldH2_effects <- do(by_gradient,
                                   tidy(
                                     lm(Field_H2 ~ Elevation, weights=N, data = .)))


#Write WS_fieldH2 to file to use in estimating selection response

write.table(WS_fieldH2, "WS_fieldH2.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ',')

