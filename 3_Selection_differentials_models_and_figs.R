
#Load packages

library(tidyverse)
library(broom)
library(lme4)


#Load data file with trait means and cage fitness

Trait_fitness <- read.csv("Cage_data.csv", header = TRUE, sep = ',')




#Calculate mean emergence and survival for each site and gradient and add to data frame

Site_mean_emergence <- Trait_fitness %>%
  group_by(Site) %>%
  summarise_at(vars(Emergence), list(Site_mean_emergence = mean))


Site_mean_survival <- Trait_fitness %>%
  group_by(Site) %>%
  summarise_at(vars(Survivors), list(Site_mean_survival = mean))


Gradient_mean_emergence <- Trait_fitness %>%
  group_by(Gradient) %>%
  summarise_at(vars(Emergence), list(Gradient_mean_emergence = mean))


Gradient_mean_survival <- Trait_fitness %>%
  group_by(Gradient) %>%
  summarise_at(vars(Survivors), list(Gradient_mean_survival = mean))


Trait_fitness <- merge(Trait_fitness, Site_mean_emergence, by = "Site")
Trait_fitness <- merge(Trait_fitness, Site_mean_survival, by = "Site")
Trait_fitness <- merge(Trait_fitness, Gradient_mean_emergence, by = "Gradient")
Trait_fitness <- merge(Trait_fitness, Gradient_mean_survival, by = "Gradient")


#Calculate relative fitness (by emergence and survival) for each cage as fitness divided by site mean fitness
#Do the same for gradient

Trait_fitness$Rel_fitness_emergence <- Trait_fitness$Emergence/Trait_fitness$Site_mean_emergence
Trait_fitness$Rel_fitness_survival <- Trait_fitness$Survivors/Trait_fitness$Site_mean_survival
Trait_fitness$Rel_fitness_emergence_gradient <- Trait_fitness$Emergence/Trait_fitness$Gradient_mean_emergence
Trait_fitness$Rel_fitness_survival_gradient <- Trait_fitness$Survivors/Trait_fitness$Gradient_mean_survival



#Estimate selection differentials by running regressions of standardised trait mean against relative fitness at each site for each trait

by_site <- group_by(Trait_fitness, Site)

#Cold tolerance

#Summary of models for each site

CT_seldiff_model_summaries <- do(by_site,
                                 glance(
                                   lm(Rel_fitness_emergence ~ Z_MeanCT, data = .)))

CT_seldiff_model_summaries$Trait <- "CT"

#Summary of estimated effects for each model term
  
CT_seldiff_effects <- do(by_site,
                                 tidy(
                                   lm(Rel_fitness_emergence ~ Z_MeanCT, data = .)))

CT_seldiff_effects <- subset(CT_seldiff_effects, term == "Z_MeanCT")
  
#Rename term column 'Trait'

CT_seldiff_effects <- CT_seldiff_effects %>%
  rename(Trait = term)

#Add Nobs term from model summary

CT_seldiff_effects$Nobs <- CT_seldiff_model_summaries$nobs


#Calculate p values adjusted for multiple comparisons
CT_seldiff_effects$P_FDR <- p.adjust(CT_seldiff_effects$p.value, method = "fdr")

#Get count of lines contributing to selection estimate
#Limit to rows for which we have CT

CT <- Trait_fitness %>%
  drop_na(Z_MeanCT)

CT_LineCount <- CT %>%
  group_by(Site) %>%
  summarise(No_lines = length(unique(Line)))

#Add to CT_seldiff_effects
CT_seldiff_effects <- merge(CT_seldiff_effects, CT_LineCount, by="Site")


#Heat tolerance

#Summary of models for each site

HT_seldiff_model_summaries <- do(by_site,
                                 glance(
                                   lm(Rel_fitness_emergence ~ Z_MeanHT, data = .)))

HT_seldiff_model_summaries$Trait <- "HT"

#Summary of estimated effects for each model term

HT_seldiff_effects <- do(by_site,
                         tidy(
                           lm(Rel_fitness_emergence ~ Z_MeanHT, data = .)))

HT_seldiff_effects <- subset(HT_seldiff_effects, term == "Z_MeanHT")


#Rename term column 'Trait'

HT_seldiff_effects <- HT_seldiff_effects %>%
  rename(Trait = term)

#Add Nobs term from model summary

HT_seldiff_effects$Nobs <- HT_seldiff_model_summaries$nobs

#Calculate p values adjusted for multiple comparisons
HT_seldiff_effects$P_FDR <- p.adjust(HT_seldiff_effects$p.value, method = "fdr")

#Get count of lines contributing to selection estimate
#Limit to rows for which we have HT

HT <- Trait_fitness %>%
  drop_na(Z_MeanHT)

HT_LineCount <- HT %>%
  group_by(Site) %>%
  summarise(No_lines = length(unique(Line)))

#Add to HT_seldiff_effects
HT_seldiff_effects <- merge(HT_seldiff_effects, HT_LineCount, by="Site")



#Body size

#Summary of models for each site

WS_seldiff_model_summaries <- do(by_site,
                                 glance(
                                   lm(Rel_fitness_emergence ~ Z_MeanWS, data = .)))

WS_seldiff_model_summaries$Trait <- "WS"

#Summary of estimated effects for each model term

WS_seldiff_effects <- do(by_site,
                         tidy(
                           lm(Rel_fitness_emergence ~ Z_MeanWS, data = .)))

WS_seldiff_effects <- subset(WS_seldiff_effects, term == "Z_MeanWS")


#Rename term column 'Trait'

WS_seldiff_effects <- WS_seldiff_effects %>%
  rename(Trait = term)

#Add Nobs term from model summary

WS_seldiff_effects$Nobs <- WS_seldiff_model_summaries$nobs

#Calculate p values adjusted for multiple comparisons
WS_seldiff_effects$P_FDR <- p.adjust(WS_seldiff_effects$p.value, method = "fdr")

#Get count of lines contributing to selection estimate
#Limit to rows for which we have WS

WS <- Trait_fitness %>%
  drop_na(Z_MeanWS)

WS_LineCount <- WS %>%
  group_by(Site) %>%
  summarise(No_lines = length(unique(Line)))

#Add to WS_seldiff_effects
WS_seldiff_effects <- merge(WS_seldiff_effects, WS_LineCount, by="Site")

#Combine summaries of estimated effects for the three traits

Seldiff_model_summaries <- rbind(CT_seldiff_model_summaries, HT_seldiff_model_summaries, WS_seldiff_model_summaries)
Seldiff_effects <- rbind(CT_seldiff_effects, HT_seldiff_effects, WS_seldiff_effects)

Seldiff_effects$Trait <- recode_factor(Seldiff_effects$Trait, Z_MeanCT = "CT", Z_MeanHT = "HT", Z_MeanWS = "WS") 


#Estimate overall selection differentials for each gradient
#Include Site as a random effect.  Method above doesn't work for mixed models so separate by gradient then run models

Trait_fitness_Danbulla <- subset(Trait_fitness, Gradient=="Danbulla")
Trait_fitness_Paluma <- subset(Trait_fitness, Gradient=="Paluma")

#Cold tolerance
S_CT_Danbulla <- lmer(Rel_fitness_emergence_gradient ~ Z_MeanCT_Gradient + (1|Site), data = Trait_fitness_Danbulla)
S_CT_Danbulla_summary <- summary(S_CT_Danbulla)[["coefficients"]]

#Count cages and lines
Trait_fitness_Danbulla_CT <- Trait_fitness_Danbulla[!is.na(Trait_fitness_Danbulla$Z_MeanCT_Gradient),]
N_cages <- length(unique(Trait_fitness_Danbulla_CT$Cage))
N_lines <- length(unique(Trait_fitness_Danbulla_CT$Line))

S_CT_Paluma <- lmer(Rel_fitness_emergence_gradient ~ Z_MeanCT_Gradient + (1|Site), data = Trait_fitness_Paluma)
S_CT_Paluma_summary <- summary(S_CT_Paluma)[["coefficients"]]

#Count cages and lines
Trait_fitness_Paluma_CT <- Trait_fitness_Paluma[!is.na(Trait_fitness_Paluma$Z_MeanCT_Gradient),]
N_cages <- length(unique(Trait_fitness_Paluma_CT$Cage))
N_lines <- length(unique(Trait_fitness_Paluma_CT$Line))


#Heat tolerance
S_HT_Danbulla <- lmer(Rel_fitness_emergence_gradient ~ Z_MeanHT_Gradient + (1|Site), data = Trait_fitness_Danbulla)
S_HT_Danbulla_summary <- summary(S_HT_Danbulla)[["coefficients"]]

#Count cages and lines
Trait_fitness_Danbulla_HT <- Trait_fitness_Danbulla[!is.na(Trait_fitness_Danbulla$Z_MeanHT_Gradient),]
N_cages <- length(unique(Trait_fitness_Danbulla_HT$Cage))
N_lines <- length(unique(Trait_fitness_Danbulla_HT$Line))


S_HT_Paluma <- lmer(Rel_fitness_emergence_gradient ~ Z_MeanHT_Gradient + (1|Site), data = Trait_fitness_Paluma)
S_HT_Paluma_summary <- summary(S_HT_Paluma)[["coefficients"]]

#Count cages and lines
Trait_fitness_Paluma_HT <- Trait_fitness_Paluma[!is.na(Trait_fitness_Paluma$Z_MeanHT_Gradient),]
N_cages <- length(unique(Trait_fitness_Paluma_HT$Cage))
N_lines <- length(unique(Trait_fitness_Paluma_HT$Line))



#Body size
S_WS_Danbulla <- lmer(Rel_fitness_emergence_gradient ~ Z_MeanWS_Gradient + (1|Site), data = Trait_fitness_Danbulla)
S_WS_Danbulla_summary <- summary(S_WS_Danbulla)[["coefficients"]]

#Count cages and lines
Trait_fitness_Danbulla_WS <- Trait_fitness_Danbulla[!is.na(Trait_fitness_Danbulla$Z_MeanWS_Gradient),]
N_cages <- length(unique(Trait_fitness_Danbulla_WS$Cage))
N_lines <- length(unique(Trait_fitness_Danbulla_WS$Line))


S_WS_Paluma <- lmer(Rel_fitness_emergence_gradient ~ Z_MeanWS_Gradient + (1|Site), data = Trait_fitness_Paluma)
S_WS_Paluma_summary <- summary(S_WS_Paluma)[["coefficients"]]

#Count cages and lines
Trait_fitness_Paluma_WS <- Trait_fitness_Paluma[!is.na(Trait_fitness_Paluma$Z_MeanWS_Gradient),]
N_cages <- length(unique(Trait_fitness_Paluma_WS$Cage))
N_lines <- length(unique(Trait_fitness_Paluma_WS$Line))

#Get FDR correction for p-values for overall S estimates

p <- c(0.200, 0.126, 0.072, 0.860, 0.642, 0.457)
pFDR <- p.adjust(p, method="fdr")

#Test for a linear relationship between selection differentials and elevation for each trait at each gradient

#Need to add gradient and elevation to Seldiff_effects data frame

Site_details <- read.csv("Site_details.csv", header = TRUE, sep = ',')

Seldiff_effects <- merge(Seldiff_effects, Site_details, by = "Site")
Seldiff_model_summaries <- merge(Seldiff_model_summaries, Site_details, by = "Site")

#Group by gradient and trait to run regression models

by_gradient_trait <- Seldiff_effects %>%
  group_by(Gradient, Trait)

#Summary of models for each gradient and trait

Elevation_seldiff_model_summaries <- do(by_gradient_trait,
                                 glance(
                                   lm(estimate ~ Elevation, data = .)))

#Summary of estimated effects for each model term

Elevation_seldiff_effects <- do(by_gradient_trait,
                         tidy(
                           lm(estimate ~ Elevation, data = .)))

Elevation_seldiff_effects <- subset(Elevation_seldiff_effects, term == "Elevation")

#Calculate p values adjusted for multiple comparisons
Elevation_seldiff_effects$P_FDR <- p.adjust(Elevation_seldiff_effects$p.value, method = "fdr")

#Also test for linear relationship between the absolute value of selection differentials and elevation for each trait at each gradient

#Summary of models for each gradient and trait

Elevation_abs_seldiff_model_summaries <- do(by_gradient_trait,
                                        glance(
                                          lm(abs(estimate) ~ Elevation, data = .)))

#Summary of estimated effects for each model term

Elevation_abs_seldiff_effects <- do(by_gradient_trait,
                                tidy(
                                  lm(abs(estimate) ~ Elevation, data = .)))

Elevation_abs_seldiff_effects <- subset(Elevation_abs_seldiff_effects, term == "Elevation")

#Calculate p values adjusted for multiple comparisons
Elevation_abs_seldiff_effects$P_FDR <- p.adjust(Elevation_abs_seldiff_effects$p.value, method = "fdr")


#Make figures showing selection differentials against elevation for each trait at each gradient

# Give gradients and traits correct names
Seldiff_effects$Gradient <- recode_factor(Seldiff_effects$Gradient, Danbulla = "Mt Edith", Paluma = "Paluma") 
Seldiff_effects$Trait <- recode_factor(Seldiff_effects$Trait, CT = "Cold tolerance", HT = "Heat tolerance", WS = "Wing size") 

#Also make a dataset just for cold tolerance at Paluma so we can plot the regression line only for the significant relationship

CT_Seldiff_effects <- Seldiff_effects %>%
  filter(Gradient=="Paluma",
         Trait=="Cold tolerance")

# Plots
Elevation_seldiff_plots <- Seldiff_effects %>%
  mutate(Gradient=fct_relevel(Gradient, "Mt Edith", "Paluma")) %>%
  ggplot(aes(x=Elevation, y = estimate))+
  geom_point(aes(x=Elevation, y = estimate), size=2)+
  geom_errorbar(aes(x=Elevation, y = estimate, ymax=estimate+std.error, ymin=estimate-std.error))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid(cols=vars(Trait), rows = vars(Gradient))+
  geom_smooth(data=CT_Seldiff_effects, aes(x=Elevation, y = estimate), method = "lm", se = FALSE, color = "black")+
  xlab("Elevation")+ 
  ylab(expression(paste("Selection differential, ", italic("S"), " (\u00B1 SE)")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")

ggsave("Selection_differential_plots.png", plot = Elevation_seldiff_plots, width = 30, height = 20, units = "cm", bg = "white")

#Write Seldiff_effects to file to use in estimating selection response
#First rename 'estimate' as S

Seldiff_effects <- Seldiff_effects %>%
  rename(S = estimate)

write.table(Seldiff_effects, "Trait_selection_differentials.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ',')

