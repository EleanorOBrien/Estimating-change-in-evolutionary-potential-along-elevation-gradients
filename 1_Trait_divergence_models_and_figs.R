
#Code to analyse trait divergence between gradients and elevations and produce figures

#load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(MCMCglmm)
library(ggeffects)
library(dplyr)
library(ggplot2)


#load data

Paluma_traits <- read.csv("Paluma_Traits.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)

Danbulla_traits <- read.csv("Danbulla_Traits.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE)

Traits <- rbind(Paluma_traits, Danbulla_traits)




#Analyse variation in traits between gradients and elevations within gradients

#Ensure variable types are specified correctly
Traits$Dam <- as.factor(Traits$Mother)
Traits$Batch_CT <- as.factor(Traits$Batch_CT)
Traits$Vial_CT <- as.factor(Traits$Vial_CT)
Traits$Batch_HT <- as.factor(Traits$Batch_HT)
Traits$Vial_HT <- as.factor(Traits$Vial_HT)

#Cold tolerance

CT <- subset(Traits, Trait=="Cold_tolerance")


#Exclude zeroes

CT.1 <- subset(CT, Taken_CT=="1")

CT_model <- lmer(Offspring_CT ~ Transect*Altitude_origin + (1|Mother/Batch_CT), data = CT.1, na.action=na.exclude)

summary(CT_model)
CT.anova <- anova(CT_model)

#Calculate estimated marginal means

CT_emm <- ggemmeans(CT_model, terms=c("Transect", "Altitude_origin"))

#Significance of random effect variance

CT_model.a <- lmer(Offspring_CT ~ Transect*Altitude_origin + (1|Batch_CT), data = CT.1, na.action=na.exclude)

anova(CT_model, CT_model.a)


#Heat tolerance

HT <- subset(Traits, Trait=="Heat_tolerance")

HT_model <- lmer(TimeHT_MINS ~ Transect*Altitude_origin + (1|Mother) + (1|Batch_HT), data = HT, na.action=na.exclude)

summary(HT_model)
HT.anova <- anova(HT_model)

#Calculate estimated marginal means

HT_emm <- ggemmeans(HT_model, terms=c("Transect", "Altitude_origin"))

#Significance of random effect variance

HT_model.a <- lmer(TimeHT_MINS ~ Transect*Altitude_origin + (1|Batch_HT), data = HT, na.action=na.exclude)

anova(HT_model, HT_model.a)

HT_model.b <- lmer(TimeHT_MINS ~ Transect*Altitude_origin + (1|Mother), data = HT, na.action=na.exclude)

anova(HT_model, HT_model.b)


#Wing size

WS <- subset(Traits, Trait=="Wing_size")

#Get wing size in mm
WS$WingSize_mm <- WS$Centroid_Size*0.000861

WS_model <- lmer(WingSize_mm ~ Transect*Altitude_origin + (1|Mother) + (1|Replicate_WS), data = WS, na.action = na.exclude)

summary(WS_model)
WS.anova <- anova(WS_model)

#Calculate estimated marginal means

WS_emm <- ggemmeans(WS_model, terms=c("Transect", "Altitude_origin"))

#Significance of random effects
WS_model.a <- lmer(WingSize_mm ~ Transect*Altitude_origin + (1|Replicate_WS), data = WS, na.action = na.exclude)

anova(WS_model, WS_model.a)


#Make figures of trait variation between gradients and elevations

#Cold tolerance

#Rename transects
CT.1$Transect <- as.factor(CT.1$Transect)
CT.1$Transect <- recode_factor(CT.1$Transect, Danbulla = "Mt Edith", Paluma = "Paluma")

#Ensure estimated marginal means have same variable names
CT_emm <- rename(CT_emm, Transect = x,
                 SE = std.error,
                 Altitude_origin = group)
CT_emm <- rename(CT_emm, Offspring_CT = predicted)
CT_emm$Transect <- recode_factor(CT_emm$Transect, Danbulla = "Mt Edith", Paluma = "Paluma") 

#Make plots and overlay estimated marginal means from models
CT_figure <- CT.1 %>%
  mutate(Transect=fct_relevel(Transect, "Mt Edith", "Paluma")) %>%
  mutate(Altitude_origin=fct_relevel(Altitude_origin, "Low", "High")) %>%
  ggplot(aes(x=Altitude_origin, y = Offspring_CT, group=Altitude_origin, fill=Altitude_origin, color=Altitude_origin))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.5, binwidth=1.0)+
  facet_grid(cols=vars(Transect))+
  geom_point(data=CT_emm, aes(x=Altitude_origin, y=Offspring_CT), color="black", size=3)+
  geom_errorbar(data=CT_emm, aes(x=Altitude_origin, y=Offspring_CT, ymax=conf.high, ymin=conf.low), colour="black", width=0.05, lwd=1)+
  xlab(" ")+ ylab("No. of offspring\n   ")+
  coord_cartesian(ylim=c(0,50))+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="gray", linetype="dashed"), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(face = "bold", size=16),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18),
        strip.text.x = element_text(size = 16),
        legend.position = "none")+
  scale_color_manual(values=c("Low"="firebrick", "High"="blue4"))+
  scale_fill_discrete(name = "Source elevation")+
  ggtitle("Cold tolerance")



#Heat tolerance

#Rename transects
HT$Transect <- as.factor(HT$Transect)
HT$Transect <- recode_factor(HT$Transect, Danbulla = "Mt Edith", Paluma = "Paluma")

#Ensure estimated marginal means have same variable names
HT_emm <- rename(HT_emm, Transect = x,
                 SE = std.error,
                 Altitude_origin = group,
                 TimeHT_MINS = predicted)

HT_emm$Transect <- recode_factor(HT_emm$Transect, Danbulla = "Mt Edith", Paluma = "Paluma") 

HT_figure <- HT %>%
  mutate(Transect=fct_relevel(Transect, "Mt Edith", "Paluma")) %>%
  mutate(Altitude_origin=fct_relevel(Altitude_origin, "Low", "High")) %>%
  ggplot(aes(x=Altitude_origin, y = TimeHT_MINS, group=Altitude_origin, fill=Altitude_origin, color=Altitude_origin))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.5, binwidth=0.6)+
  facet_grid(cols=vars(Transect))+
  geom_point(data=HT_emm, aes(x=Altitude_origin, y=TimeHT_MINS), color="black", size=3)+
  geom_errorbar(data=HT_emm, aes(x=Altitude_origin, y=TimeHT_MINS, ymax=conf.high, ymin=conf.low), colour="black", width=0.05, lwd=1)+
  xlab(" ")+ ylab("Heat knockdown time\n(mins)")+
  coord_cartesian(ylim=c(0,40))+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="gray", linetype="dashed"), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(face = "bold", size=16),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18),
        strip.text.x = element_text(size = 16),
        legend.position = "none")+
  scale_color_manual(values=c("Low"="firebrick", "High"="blue4"))+
  scale_fill_discrete(name = "Source elevation")+
  ggtitle("Heat tolerance")




#Body size

#Rename transects
WS$Transect <- as.factor(WS$Transect)
WS$Transect <- recode_factor(WS$Transect, Danbulla = "Mt Edith", Paluma = "Paluma")

#Ensure estimated marginal means have same variable names
WS_emm <- rename(WS_emm, Transect = x,
                 SE = std.error,
                 Altitude_origin = group,
                 WingSize_mm = predicted)

WS_emm$Transect <- recode_factor(WS_emm$Transect, Danbulla = "Mt Edith", Paluma = "Paluma") 

WS_figure <- WS %>%
  mutate(Transect=fct_relevel(Transect, "Mt Edith", "Paluma")) %>%
  mutate(Altitude_origin=fct_relevel(Altitude_origin, "Low", "High")) %>%
  ggplot(aes(x=Altitude_origin, y = WingSize_mm, group=Altitude_origin, fill=Altitude_origin, color=Altitude_origin))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.5, binwidth=0.004)+
  facet_grid(cols=vars(Transect))+
  geom_point(data=WS_emm, aes(x=Altitude_origin, y=WingSize_mm), color= "black", size=3)+
  geom_errorbar(data=WS_emm, aes(x=Altitude_origin, y=WingSize_mm, ymax=conf.high, ymin=conf.low), colour="black", width=0.05, lwd=1)+
  xlab(" ")+ ylab("Wing centroid size\n(mm)")+
  coord_cartesian(ylim=c(1.7,2.0))+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="gray", linetype="dashed"), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(face = "bold", size=16),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18),
        strip.text.x = element_text(size = 16),
        legend.position = "none")+
  scale_color_manual(values=c("Low"="firebrick", "High"="blue4"))+
  scale_fill_discrete(name = "Source elevation")+
  ggtitle("Wing size")


#Combine figures for the 3 traits
Trait_divergence_figures <- ggarrange(CT_figure, HT_figure, WS_figure, ncol=1, nrow=3, align='v')

ggsave("Trait_divergence_plots.png", plot = Trait_divergence_figures, width = 25, height = 25, units = "cm", bg = "white")