

#Load packages

library(tidyverse)
library(lme4)
library(broom)
library(ggpubr)


#Load selection differential data 

Selection <- read.csv("Trait_selection_differentials.csv", header = TRUE, sep = ',')

#Limit to just wing size

Selection <- subset(Selection, Trait=="Wing size")

#Remove sites 71 and 72 because we don't have heritability data

Selection <- Selection %>%
  filter(Site!=71,
         Site!=72)

#Load field heritability of body size data

Field_heritability <- read.csv("WS_fieldH2.csv", header = TRUE, sep = ',')


#Make selection response data frame by combining the required columns of the selection and heritability data frames

SelectionWS <- Selection %>%
  select(Gradient, Site, Elevation, S, std.error, Nobs, No_lines) %>%
  rename(SE_S = std.error,
         Ncages_S = Nobs,
         Nlines_S = No_lines) 
  
SelectionWS$Gradient <- recode_factor(SelectionWS$Gradient, `Mt Edith` = "Danbulla",
                Paluma = "Paluma")

Field_heritabilityWS <- Field_heritability %>%
  select(Gradient, Site, Elevation, Field_H2, std.error, N) %>%
  rename(SE_H2 = std.error,
         Ncages_H2 = N)


Selection_response <- merge(SelectionWS, Field_heritabilityWS, by = c("Site", "Gradient", "Elevation"))



#Calculate selection response and standard errors

#Calculate selection response R as S*H2

Selection_response$R <- Selection_response$S * Selection_response$Field_H2


#Calculate upper and lower bounds of R using standard errors of S

Selection_response$R_upper <- (Selection_response$S + Selection_response$SE_S) * Selection_response$Field_H2
Selection_response$R_lower <- (Selection_response$S - Selection_response$SE_S) * Selection_response$Field_H2

#Also need to add columns with lab estimates of H2 for plotting

#Lab WS H2 estimates were: 0.262 at Mt Edith, 0.292 at Paluma

Selection_response_Danbulla <- Selection_response[Selection_response$Gradient=="Danbulla",]
Selection_response_Paluma <- Selection_response[Selection_response$Gradient=="Paluma",]

Selection_response_Danbulla$Lab_H2 <- 0.262
Selection_response_Paluma$Lab_H2 <- 0.292

Selection_response <- rbind(Selection_response_Danbulla, Selection_response_Paluma)



#Calculate lab selection response R_lab as S*H2

Selection_response$R_lab <- Selection_response$S * Selection_response$Lab_H2


#Calculate upper and lower bounds of R_lab using standard errors of S

Selection_response$R_lab_upper <- (Selection_response$S + Selection_response$SE_S) * Selection_response$Lab_H2
Selection_response$R_lab_lower <- (Selection_response$S - Selection_response$SE_S) * Selection_response$Lab_H2


#Test for linear relationship between magnitude of selection differential (WS) and estimated field heritability at each gradient

by_gradient <- group_by(Selection_response, Gradient)


WS_sel_heritability <- do(by_gradient,
                          glance(
                            lm(S ~ Field_H2, data = .)))
                            

WS_sel_heritability_coeffs <- do(by_gradient,
                                 tidy(
                                   lm(S ~ Field_H2, data = .)))
                                   

#Plot this
Selection_response$Gradient <- recode_factor(Selection_response$Gradient, Danbulla = "Mt Edith", Paluma = "Paluma")
Selection_response_Danbulla$Gradient <- recode_factor(Selection_response_Danbulla$Gradient, Danbulla = "Mt Edith", Paluma = "Paluma")
Selection_response_Paluma$Gradient <- recode_factor(Selection_response_Paluma$Gradient, Danbulla = "Mt Edith", Paluma = "Paluma")

WS_sel_heritability_plot <- ggplot(data = Selection_response, aes(x = Field_H2, y = S))+
  geom_errorbar(aes(x=Field_H2, ymin=S-SE_S, ymax=S+SE_S))+
  geom_errorbarh(aes(y=S, xmin=Field_H2-SE_H2, xmax=Field_H2+SE_H2))+
  geom_point(aes(x=Field_H2, y=S), size = 3, colour = "black", fill = "black")+
  facet_grid(cols=vars(Gradient))+
  coord_cartesian(xlim=c(0,1))+
  geom_smooth(data = Selection_response_Danbulla, aes(x = Field_H2, y = S), method = "lm", colour="black", se=FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.position = "none")+
  xlab(expression(paste("Field ", italic("H")^{"2"}, " (\u00B1 SE)")))+
  ylab(expression(paste("Selection differential, ", italic("S"), " (\u00B1 SE)")))

ggsave("Selection_heritability_plots.png", plot = WS_sel_heritability_plot, width = 30, height = 20, units = "cm", bg = "white")




#Test for linear relationship between field predicted selection response (WS) and elevation at each gradient

WS_field_response <- do(by_gradient,
                        glance(
                          lm(R ~ Elevation, data = .)))

WS_field_response_coeffs <- do(by_gradient,
                               tidy(
                                 lm(R ~ Elevation, data = .))) 


#Test for linear relationship between magnitude of field predicted selection response (WS) and elevation at each gradient


WS_abs_field_response <- do(by_gradient,
                        glance(
                          lm(abs(R) ~ Elevation, data = .)))

WS_abs_field_response_coeffs <- do(by_gradient,
                        tidy(
                          lm(abs(R) ~ Elevation, data = .))) 





#Plots of selection response along each gradient

#Lab & field wing size
dodge <- position_dodge(width = 2)


Plot_WS_response <-Selection_response %>%
  ggplot(aes(x=Elevation, y=R))+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_errorbar(aes(x=Elevation, ymin=R_lower, ymax=R_upper))+
  geom_point(aes(x=(Elevation), y=R), size = 4, shape = 24, fill = "gray47", colour = "black")+
  geom_errorbar(aes(x=(Elevation+5), y=R_lab, ymin=R_lab_lower, ymax=R_lab_upper))+
  geom_point(aes(x=(Elevation+5), y=R_lab), size = 4, shape = 21, fill = "gray88", colour = "black")+
  facet_grid(cols=vars(Gradient))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(0, 1050), ylim=c(-0.6,0.6))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.position = "none")+
  xlab("Elevation (m)")+
  ylab(expression(paste("Predicted response,  ", italic("R"), " (\u00B1 SE)")))


#Plot elevation against absolute value of the selection response

#Need to calculate standard errors of absolute values of R

Selection_response$abs_R_lower <- Selection_response$Field_H2 * (abs(Selection_response$S)-Selection_response$SE_S)
Selection_response$abs_R_upper <- Selection_response$Field_H2 * (abs(Selection_response$S)+Selection_response$SE_S)


Plot_abs_WS_response <-Selection_response %>%
  ggplot(aes(x=Elevation, y=abs(R)))+
  # geom_hline(aes(yintercept=0), linetype="dashed")+
  geom_errorbar(aes(x=Elevation, ymin=abs_R_lower, ymax=abs_R_upper))+
  geom_point(aes(x=Elevation, y=abs(R)), size = 4, shape = 24, fill = "gray47", colour = "black")+
  # geom_errorbar(aes(x=(Elevation+5), y=abs(R_lab), ymin=min(abs(R_lab_lower), abs(R_lab_upper)), ymax=max(abs(R_lab_lower), abs(R_lab_upper))), position = dodge)+
  # geom_point(aes(x=(Elevation+5), y=abs(R_lab)), size = 4, shape = 21, fill = "gray88", colour = "black", position = dodge)+
  facet_grid(cols=vars(Gradient))+
  geom_smooth(method = "lm", se = FALSE, color = "black", size=0.2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(0, 1050), ylim=c(0,0.6))+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_blank(),
        legend.position = "none")+
  xlab("Elevation (m)")+
  ylab(expression(paste("Absolute predicted response,  ", "|",italic("R"),"|", " (\u00B1 SE)")))

Elevation_selresponse_plots <- ggarrange(Plot_WS_response, Plot_abs_WS_response, nrow = 2, ncol = 1, align="v")

ggsave("Elevation_selresponse.png", plot = Elevation_selresponse_plots, width = 25, height = 25, units = "cm", bg = "white")