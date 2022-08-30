
This document provides descriptions of all data and code used in the analysis presented in the manuscript "Variation in the strength of selection but no trait divergence between elevational extremes in a tropical rainforest Drosophila". 

NOTE: One of the elevation gradients used in this study is referred to in the manuscript as "Mt Edith".  An alternative name that we have used previously for this gradient is "Danbulla", and this name has been used in most of the data files. In all cases, data for Danbulla relates to the Mt Edith gradient.


Description of code and data files used in each component of the study:


(1) Divergence of traits within and between elevation gradients


R Code: "1_Trait_divergence_models_and_figs.R"

This file provides the full R code used to run linear mixed models analysing trait variation among gradients and elevations and to generate figures (Figure 1).  The data files used in this analysis are also provided and are listed and described below.


Data: "Paluma_Traits.csv"; "Danbulla_Traits.csv"

These files contain data for all traits (cold tolerance, heat tolerance and wing size) and field fitness (survival and emergence) for each transect. Note that for cold tolerance, heat tolerance and wing size, there is one line per fly. Field fitness was measured at the level of the cage, therefore each line represents a transplant cage.  

Columns in data files:
ID			Individual identifier of each fly assayed for cold tolerance, heat tolerance or wing size
Mother			Individual identifier of the mother of the individual assayed
Father			Individual identifier of the father of the individual assayed	
Maternal_line		Isofemale line from which the mother was derived
Paternal_line		Isofemale line from which the father was derived
Cross_CT_HT		Identifier of the line combination making up crosses used to generate flies assayed for cold and heat tolerance
Cross_WS		Identifier of the line combination making up crosses used to generate flies assayed for wing size
Replicate		Replicate of the cross (line combination) (A - C for cold and heat tolerance; D - E for wing size) 
Transect		Transect (gradient) from which fly originated (Paluma or Danbulla)
Altitude_origin		Altitude (elevation) from which fly originated (Low or High)
Site_number		Numeric identifier of site from which fly originated (see "Site_details.csv" for elevation of each site)	
Trait			Trait assayed (Cold_tolerance, Heat_tolerance, Wing_size or Field_fitness)
Batch_CT		Batch in which cold tolerance was assayed (1 or 2)
Batch_HT		Batch in which heat tolerance was assayed (1 - 9)
Offspring_CT		Number of offspring produced following cold shock (measure of cold tolerance)
Taken_CT		Binary variable indicating whether any offspring were produced following cold shock (0 = no, 1 = yes) so zeroes can be excluded if required
TimeHT_SECS		Time taken to knockdown in response to heat shock (in seconds)	
TimeHT_MINS		Time taken to knockdown in response to heat shock (in minutes)
Centroid_Size		Wing centroid size (in pixels; see code for conversion to mm)
Log_Centroid_Size	Log of wing centroid size in pixels

# THE REMAINING COLUMNS RELATE ONLY TO LINES WHERE Trait = "Field_fitness":	
Transplant_site		Numeric identifier of site where cage was placed for assessment of field fitness (see "Site_details.csv" for elevation of each site)		
Cage			Numeric identifier of transplant cage
Emergence		Number of offspring emerging from cage
Rel_fitness_emergence	Relative fitness calculated as emergence divided by mean emergence of all cages at that transplant site
Survivors		Number of transplanted adult flies that survived to 5 days in the field (out of 10 that were originally transplanted)
Rel_fitness_survivors	Relative fitness calculated as survivors divided by mean no. survivors in all cages at that transplant site




(2) Estimating additive genetic (co)variances of traits 


R Code: "2_Trait_covariance_models.R"

This file provides the code used to estimate genetic (co)variance of traits assayed in the laboratory: cold tolerance, heat tolerance and wing size. 
We used animal models implemented in MCMCglmm (see main manuscript and Supplementary Information for full details).  

Data: 
Phenotype data: "Paluma_Traits.csv"; "Danbulla_Traits.csv" (see description in (1) above)
Pedigree data: "Paluma_pedigree.csv"; "Danbulla_pedigree.csv"

Columns in pedigree data files:
ID			Individual identifier of each fly.  Includes flies assayed as well as mothers, fathers and an identifier for isofemale line founders. Corresponds with ID in trait files
Dam			Individual identifier of the mother (dam) of each fly.  Corresponds with Mother column in trait files
Sire			Individual identifier of the father (sire) of each fly. Corresponds with Father column in trait files


Model output: "mv.Pal"; "mv.Dan" 
These are R objects containing the results of the MCMCglmm models run for each gradient to estimate genetic and maternal effect variances and covariances. 
These are provided because each model takes several hours to run so can't be replicated quickly. 
Note that results differ slightly between replicate runs of the same model due to stochastic variation in the MCMC sampling.  
We have provided the output from which the results we present were extracted.  
They can be loaded into R using the readRDS() function in base R e.g. readRDS("mv.Pal") if "mv.Pal" is in the working directory.



(3) Estimating selection on thermal tolerance and wing size along elevation gradients

R Code: "3_Selection_differentials_models_and_figs.R"

This code was used to calculate selection differentials for each trait at each transplant site using the genetic covariance of trait mean and fitness at each site.

Data: "Cage_data.csv"; "Site_details.csv"
      
This data file contains data for each cage transplanted at field sites, with one line per cage.

Columns in data files:

Cage_data.csv


Gradient		Gradient where flies were transplanted (Paluma or Danbulla). Note this is also the gradient where flies were sourced.	
Site			Numeric identifier of site where cage was transplanted (see "Site_details.csv" for elevation of each site)	
Cage			Numeric identifier of transplant cage
Line			Isofemale line from which the mother of flies in the transplant cage was derived
Emergence		Number of offspring emerging from cage	
Survivors		Number of transplanted adult flies that survived to 5 days in the field (out of 10 that were originally transplanted)	
MeanCT			Mean cold tolerance of flies from crosses with the same maternal isofemale line (reared and assayed in the lab)
MeanHT			Mean heat tolerance of flies from crosses with the same maternal isofemale line (reared and assayed in the lab)
MeanWSPixels		Mean wing size (in pixels) of flies from crosses with the same maternal isofemale line (reared and assayed in the lab)
MeanWSmm		Mean wing size (in mm) of flies from crosses with the same maternal isofemale line (reared and assayed in the lab)
SiteMeanCT		Mean cold tolerance of flies at that transplant site (obtained by averaging meanCT across all cages at the site). Used to calculate Z_MeanCT	
SiteMeanHT		Mean heat tolerance of flies at that transplant site (obtained by averaging meanHT across all cages at the site). Used to calculate Z_MeanHT
SiteMeanWS		Mean wing size (in mm) of flies at that transplant site (obtained by averaging meanWSmm across all cages at the site). Used to calculate Z_MeanWS
SiteSD_CT		Standard deviation of cold tolerance of flies at that transplant site (obtained by taking SD of meanCT across all cages at the site). Used to calculate Z_MeanCT
SiteSD_HT		Standard deviation of heat tolerance of flies at that transplant site (obtained by taking SD of meanHT across all cages at the site). Used to calculate Z_MeanHT
SiteSD_WS		Standard deviation of wing size (in mm) of flies at that transplant site (obtained by taking SD of meanWSmm across all cages at the site). Used to calculate Z_MeanWS
Z_MeanCT		Standardised cold tolerance (standardised within site so site mean = 0 and site SD = 1)	
Z_MeanHT		Standardised heat tolerance (standardised within site so site mean = 0 and site SD = 1)
Z_MeanWS		Standardised wing size (standardised within site so site mean = 0 and site SD = 1)
GradientMeanCT		Mean cold tolerance of flies at that gradient (obtained by averaging meanCT across all cages at the gradient). Used to calculate Z_MeanCT_Gradient
GradientMeanHT		Mean heat tolerance of flies at that gradient (obtained by averaging meanHT across all cages at the gradient). Used to calculate Z_MeanHT_Gradient
GradientMeanWS		Mean wing size (in mm) of flies at that gradient (obtained by averaging meanWSmm across all cages at the gradient). Used to calculate Z_MeanWS_Gradient	
GradientSD_CT		Standard deviation of cold tolerance of flies at that gradient (obtained by taking SD of meanCT across all cages at the gradient). Used to calculate Z_MeanCT_Gradient
GradientSD_HT		Standard deviation of heat tolerance of flies at that gradient (obtained by taking SD of meanHT across all cages at the gradient). Used to calculate Z_MeanHT_Gradient
GradientSD_WS		Standard deviation of wing size (in mm) of flies at that gradient (obtained by taking SD of meanWSmm across all cages at the gradient). Used to calculate Z_MeanWS_Gradient
Z_MeanCT_Gradient	Standardised cold tolerance (standardised within gradient so site mean = 0 and site SD = 1)
Z_MeanHT_Gradient	Standardised heat tolerance (standardised within gradient so site mean = 0 and site SD = 1)
Z_MeanWS_Gradient	Standardised wing size (standardised within gradient so site mean = 0 and site SD = 1)
WS_survivors		Mean wing size (in pixels) of adult flies surviving in the field cage 
WS_offspring		Mean wing size (in pixels) of offspring emerging from the field cage
WS_survivors_mm		Mean wing size (in mm) of adult flies surviving in the field cage
WS_offspring_mm		Mean wing size (in mm) of offspring emerging from the field cage

Site_details.csv

Site			Numeric identifier of sites used in transplant experiment (and where flies were sourced)
Gradient		Name of gradient (Paluma or Danbulla) where site is located
Elevation		Elevation of site (in m above sea level)


Output: "Trait_selection_differentials.csv"

File with estimates of selection differentials, along with their standard errors and p-values, for each trait at each site along each gradient (used to make Table S2 and in (5) below).



(4) Estimating field heritability of wing size along elevation gradients 


R Code: "4_field_heritability_WS.R"

This code was used to calculate site-specific heritabilities of wing size along each gradient, from regressions of wing size of female offspring emerging from cages on that of lab-reared females from the same maternal isofemale line.  

Data: "Cage_data.csv"; "Site_details.csv"

These are the same data files as were used in (3).  See descriptions above

Output: "WS_fieldH2.csv"

File with estimates of field heritability of wing size at each site along each gradient, along with standard errors and p-values (used to make Table S3 and in (5) below).



(5) Estimating response to selection on wing size along elevation gradients


R Code: "5_Selection_response.R"

This code was used to estimate the response to selection on 
Data: "Trait_selection_differentials.csv"; "WS_fieldH2.csv"

These are the output files from (3) and (4) above and were used to calculate predicted per-generation selection responses in wing size at sites along the elevation gradients
  
