rm(list=ls())

library(pavo)
library(tidyverse)
library(ggplot2)
library(stringr)

###get data
###from 300nm to 700nm for the analysis
data_smooth<-as.rspec(read.csv("spec_300_700nm.csv", sep=','))

#select only the specular data
data_smooth<-data_smooth %>% 
  select(matches("wl") | -matches("m1_")) %>% 
  select(matches("wl") | -matches("A1_")) %>% 
  select(matches("wl") | -matches("E1_")) %>% 
  select(matches("wl") | -matches("I1_")) %>%
  select(matches("wl") | -matches("P1_"))

###divide the dataset into groups
#achilles_female
data_achilles_F<-data_smooth %>% 
  select(matches("wl")|matches("_ACH_F_"))
#achilles_male
data_achilles_M<-data_smooth %>% 
  select(matches("wl")|matches("_ACH_M_"))
#helenor_female
data_helenor_F<-data_smooth %>% 
  select(matches("wl")|matches("_HEL_F_"))
#helenor_male
data_helenor_M<-data_smooth %>% 
  select(matches("wl")|matches("_HEL_M_"))
#bristowi_female
data_bristowi_F<-data_smooth %>% 
  select(matches("wl")|matches("_HB_F_"))
#bristowi_male
data_bristowi_M<-data_smooth %>% 
  select(matches("wl")|matches("_HB_M_"))
#theodorus_female
data_theodorus_F<-data_smooth %>% 
  select(matches("wl")|matches("_HT_F_"))
#theodorus_male
data_theodorus_M<-data_smooth %>% 
  select(matches("wl")|matches("_HT_M_"))

##averaging the individual duplicates
data_achilles_F_av<-aggspec(data_achilles_F, by=3, FUN=mean)
data_bristowi_F_av<-aggspec(data_bristowi_F, by=3, FUN=mean)
data_helenor_F_av<-aggspec(data_helenor_F, by=3, FUN=mean)
data_theodorus_F_av<-aggspec(data_theodorus_F, by=3, FUN=mean)
data_achilles_M_av<-aggspec(data_achilles_M, by=3, FUN=mean)
data_bristowi_M_av<-aggspec(data_bristowi_M, by=3, FUN=mean)
data_helenor_M_av<-aggspec(data_helenor_M, by=3, FUN=mean)
data_theodorus_M_av<-aggspec(data_theodorus_M, by=3, FUN=mean)

################ CREATE A MEAN SPECTRA FOR ALL POPULATIONS AT EVERY ANGLE ##################

mean_pop_spectra <- function(angle_pattern)
  #' calculates the mean reflectance spectra for each individual at 
  #' one specified angle of illumination
  #' 
  #' angle_pattern = enter angle type as written in ID, must be a character entry
  #' 
  #' ex : mean_pop_spectra('_15A_')
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, rptR, tidyr
{
  # Filter each dataset based on the given angle pattern
  data_achilles_F_0 <- data_achilles_F_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_bristowi_F_0 <- data_bristowi_F_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_helenor_F_0 <- data_helenor_F_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_theodorus_F_0 <- data_theodorus_F_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_achilles_M_0 <- data_achilles_M_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_bristowi_M_0 <- data_bristowi_M_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_helenor_M_0 <- data_helenor_M_av %>%
    select(matches("wl") | matches(angle_pattern))
  data_theodorus_M_0 <- data_theodorus_M_av %>%
    select(matches("wl") | matches(angle_pattern))
  
  #create dataframe with mean 0 spectra for each individual in each population
  a0_data_indiv<-left_join(data_achilles_F_0, data_bristowi_F_0, by="wl")
  a0_data_indiv<-left_join(a0_data_indiv, data_helenor_F_0, by="wl")
  a0_data_indiv<-left_join(a0_data_indiv, data_theodorus_F_0, by="wl")
  a0_data_indiv<-left_join(a0_data_indiv, data_achilles_M_0, by="wl")
  a0_data_indiv<-left_join(a0_data_indiv, data_bristowi_M_0, by="wl")
  a0_data_indiv<-left_join(a0_data_indiv, data_helenor_M_0, by="wl")
  a0_data_indiv<-left_join(a0_data_indiv, data_theodorus_M_0, by="wl")
  
  #rename columns
  names<-c("wl","achilles_F_1","achilles_F_2","achilles_F_3","achilles_F_4","achilles_F_5","achilles_F_6","achilles_F_7","achilles_F_8","achilles_F_9","achilles_F_10",
           "bristowi_F_1","bristowi_F_2","bristowi_F_3","bristowi_F_4","bristowi_F_5","bristowi_F_6","bristowi_F_7","bristowi_F_8","bristowi_F_9","bristowi_F_10",
           "helenor_F_1","helenor_F_2","helenor_F_3","helenor_F_4","helenor_F_5","helenor_F_6","helenor_F_7","helenor_F_8","helenor_F_9","helenor_F_10",
           "theodorus_F_1","theodorus_F_2","theodorus_F_3","theodorus_F_4","theodorus_F_5","theodorus_F_6","theodorus_F_7","theodorus_F_8","theodorus_F_9","theodorus_F_10",
           "achilles_M_1","achilles_M_2","achilles_M_3","achilles_M_4","achilles_M_5","achilles_M_6","achilles_M_7","achilles_M_8","achilles_M_9","achilles_M_10",
           "bristowi_M_1","bristowi_M_2","bristowi_M_3","bristowi_M_4","bristowi_M_5","bristowi_M_6","bristowi_M_7","bristowi_M_8","bristowi_M_9","bristowi_M_10",
           "helenor_M_1","helenor_M_2","helenor_M_3","helenor_M_4","helenor_M_5","helenor_M_6","helenor_M_7","helenor_M_8","helenor_M_9","helenor_M_10",
           "theodorus_M_1","theodorus_M_2","theodorus_M_3","theodorus_M_4","theodorus_M_5","theodorus_M_6","theodorus_M_7","theodorus_M_8","theodorus_M_9","theodorus_M_10")
  a0_data_indiv<-a0_data_indiv %>% 
    setNames(names)
  
  # Return file
  return(a0_data_indiv)
  
}

a0_data_indiv<-mean_pop_spectra("_0_")
a15A_data_indiv<-mean_pop_spectra("_15A_")
a15P_data_indiv<-mean_pop_spectra("_15P_")
a15I_data_indiv<-mean_pop_spectra("_15I_")
a15E_data_indiv<-mean_pop_spectra("_15E_")
a30A_data_indiv<-mean_pop_spectra("_30A_")
a30P_data_indiv<-mean_pop_spectra("_30P_")
a30I_data_indiv<-mean_pop_spectra("_30I_")
a30E_data_indiv<-mean_pop_spectra("_30E_")
a45A_data_indiv<-mean_pop_spectra("_45A_")
a45P_data_indiv<-mean_pop_spectra("_45P_")
a45I_data_indiv<-mean_pop_spectra("_45I_")
a45E_data_indiv<-mean_pop_spectra("_45E_")

##### MODELLING MORPHO VISION WITH R PAVO ###########

# sensmodel () models spectral sensitivities of retinas based on their peak cone sensitivity
#takes several optional arguments, but the main one is avector containing the
#peak sensitivities for the cones being modelled

morpho_model<-sensmodel(c(345,446,495,508,570,602), integrate=FALSE, beta=FALSE)
#spectral sensitivities taken from Pirih, P., Ilić, M., Meglič, A., & Belušič, G. (2022). 
#Opponent processing in the retinal mosaic of nymphalid butterflies. Philosophical Transactions of the Royal Society B, 377(1862), 20210275.
plot(morpho_model, col=c("#FF99CC","#33CCFF","forestgreen","#006633","#FFCC00","red"),
     lwd=2)


#### visual model functions #######

avian_visual_model_females <- function(mean_spectra)
  #' calculates the chromatic contrast perceived by a uv-sensitive avian visual model 
  #' for female wing reflectance only
  #' 
  #' mean_spectra = name of the file containing mean spectra information for 1 angle only
  #' 
  #' ex : avian_visual_model_females(a0_data_indiv)
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, rptR, tidyr, pavo
  {

  #apply uv-sensitive avian visual model to mean spectra at one angle of illumination
  #perceptual distance between coloured stimuli calculation
vismod_morpho_indiv<-vismodel(mean_spectra,
                              visual="avg.uv", #morpho_model=morpho vision, avg.uv=average uv-sensitive avian
                              achromatic="bt.dc", #receptor sensitivity data
                              illum="D65", #standard daylight illumination
                              relative=FALSE)


#create a grouping variable to group by females
vismod_morpho_indiv_F<-head(vismod_morpho_indiv,40)
females<-substr(rownames(vismod_morpho_indiv_F), 4, 4)

#bootstrap calculation
female_dS<-bootcoldist(vismod_morpho_indiv_F,
                       by=females,
                       n=c(1,2,2,4), #completement arbitraire,morpho model=c(1,1,1,1,1,1)
                       weber=0.1,
                       weber.achro=0.1)

female_dS<-as.data.frame(female_dS)
rownames(female_dS) <- c("helenor vs achilles","helenor vs theodorus","helenor vs bristowi",
                         "achilles vs theodorus", "achilles vs bristowi", "theodorus vs bristowi")

return(female_dS)

}

avian_visual_model_males <- function(mean_spectra)
  #' calculates the chromatic contrast perceived by a uv-sensitive avian visual model 
  #' for male wing reflectance only
  #' 
  #' mean_spectra = name of the file containing mean spectra information for 1 angle only
  #' 
  #' ex : avian_visual_model_males(a0_data_indiv)
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, rptR, tidyr, pavo
{
  
  #apply uv-sensitive avian visual model to mean spectra at one angle of illumination
  #perceptual distance between coloured stimuli calculation
  vismod_morpho_indiv<-vismodel(mean_spectra,
                                visual="avg.uv", #morpho_model=morpho vision, avg.uv=average uv-sensitive avian
                                achromatic="bt.dc", #receptor sensitivity data
                                illum="D65", #standard daylight illumination
                                relative=FALSE)
  
  
  #create a grouping variable to group by males
  vismod_morpho_indiv_M<-vismod_morpho_indiv[41:80,]
  males<-substr(rownames(vismod_morpho_indiv_M), 4, 4)
  
  #bootstrap calculation
  males_dS<-bootcoldist(vismod_morpho_indiv_M,
                        by=males,
                        n=c(1,2,2,4), #completement arbitraire, morpho model=c(1,1,1,1,1,1)
                        weber=0.1,
                        weber.achro=0.1)
  
  males_dS<-as.data.frame(males_dS)
  rownames(males_dS) <- c("helenor vs achilles","helenor vs theodorus","helenor vs bristowi",
                           "achilles vs theodorus", "achilles vs bristowi", "theodorus vs bristowi")
  
  return(males_dS)
  
}

morpho_visual_model_females <- function(mean_spectra)
  #' calculates the chromatic contrast perceived by a morpho visual model 
  #' for female wing reflectance only
  #' 
  #' mean_spectra = name of the file containing mean spectra information for 1 angle only
  #' 
  #' ex : morpho_visual_model_females(a0_data_indiv)
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, rptR, tidyr, pavo
{
  
  #apply uv-sensitive avian visual model to mean spectra at one angle of illumination
  #perceptual distance between coloured stimuli calculation
  vismod_morpho_indiv<-vismodel(mean_spectra,
                                visual=morpho_model, #morpho_model=morpho vision, avg.uv=average uv-sensitive avian
                                achromatic="bt.dc", #receptor sensitivity data
                                illum="D65", #standard daylight illumination
                                relative=FALSE)
  
  
  #create a grouping variable to group by females
  vismod_morpho_indiv_F<-head(vismod_morpho_indiv,40)
  females<-substr(rownames(vismod_morpho_indiv_F), 4, 4)
  
  #bootstrap calculation
  female_dS<-bootcoldist(vismod_morpho_indiv_F,
                         by=females,
                         n=c(1,1,1,1,1,1), #completement arbitraire,morpho model=c(1,1,1,1,1,1)
                         weber=0.1,
                         weber.achro=0.1)
  
  female_dS<-as.data.frame(female_dS)
  rownames(female_dS) <- c("helenor vs achilles","helenor vs theodorus","helenor vs bristowi",
                           "achilles vs theodorus", "achilles vs bristowi", "theodorus vs bristowi")
  
  return(female_dS)
  
}

morpho_visual_model_males <- function(mean_spectra)
  #' calculates the chromatic contrast perceived by a morpho visual model 
  #' for male wing reflectance only
  #' 
  #' mean_spectra = name of the file containing mean spectra information for 1 angle only
  #' 
  #' ex : morpho_visual_model_males(a0_data_indiv)
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, rptR, tidyr, pavo
{
  
  #apply uv-sensitive avian visual model to mean spectra at one angle of illumination
  #perceptual distance between coloured stimuli calculation
  vismod_morpho_indiv<-vismodel(mean_spectra,
                                visual=morpho_model, #morpho_model=morpho vision, avg.uv=average uv-sensitive avian
                                achromatic="bt.dc", #receptor sensitivity data
                                illum="D65", #standard daylight illumination
                                relative=FALSE)
  
  
  #create a grouping variable to group by males
  vismod_morpho_indiv_M<-vismod_morpho_indiv[41:80,]
  males<-substr(rownames(vismod_morpho_indiv_M), 4, 4)
  
  #bootstrap calculation
  males_dS<-bootcoldist(vismod_morpho_indiv_M,
                        by=males,
                        n=c(1,1,1,1,1,1), #completement arbitraire, morpho model=c(1,1,1,1,1,1)
                        weber=0.1,
                        weber.achro=0.1)
  
  males_dS<-as.data.frame(males_dS)
  rownames(males_dS) <- c("helenor vs achilles","helenor vs theodorus","helenor vs bristowi",
                          "achilles vs theodorus", "achilles vs bristowi", "theodorus vs bristowi")
  
  return(males_dS)
  
}

#a0 chromatic distances calculation
a0_morpho_female_dS<-morpho_visual_model_females(a0_data_indiv)
a0_morpho_male_dS<-morpho_visual_model_males(a0_data_indiv)
a0_avian_female_dS<-avian_visual_model_females(a0_data_indiv)
a0_avian_male_dS<-avian_visual_model_males(a0_data_indiv)

#a15A chromatic distances calculation
a15A_morpho_female_dS<-morpho_visual_model_females(a15A_data_indiv)
a15A_morpho_male_dS<-morpho_visual_model_males(a15A_data_indiv)
a15A_avian_female_dS<-avian_visual_model_females(a15A_data_indiv)
a15A_avian_male_dS<-avian_visual_model_males(a15A_data_indiv)

#a15P chromatic distances calculation
a15P_morpho_female_dS<-morpho_visual_model_females(a15P_data_indiv)
a15P_morpho_male_dS<-morpho_visual_model_males(a15P_data_indiv)
a15P_avian_female_dS<-avian_visual_model_females(a15P_data_indiv)
a15P_avian_male_dS<-avian_visual_model_males(a15P_data_indiv)

#a15I chromatic distances calculation
a15I_morpho_female_dS<-morpho_visual_model_females(a15I_data_indiv)
a15I_morpho_male_dS<-morpho_visual_model_males(a15I_data_indiv)
a15I_avian_female_dS<-avian_visual_model_females(a15I_data_indiv)
a15I_avian_male_dS<-avian_visual_model_males(a15I_data_indiv)

#a15E chromatic distances calculation
a15E_morpho_female_dS<-morpho_visual_model_females(a15E_data_indiv)
a15E_morpho_male_dS<-morpho_visual_model_males(a15E_data_indiv)
a15E_avian_female_dS<-avian_visual_model_females(a15E_data_indiv)
a15E_avian_male_dS<-avian_visual_model_males(a15E_data_indiv)

#a30A chromatic distances calculation
a30A_morpho_female_dS<-morpho_visual_model_females(a30A_data_indiv)
a30A_morpho_male_dS<-morpho_visual_model_males(a30A_data_indiv)
a30A_avian_female_dS<-avian_visual_model_females(a30A_data_indiv)
a30A_avian_male_dS<-avian_visual_model_males(a30A_data_indiv)

#a30P chromatic distances calculation
a30P_morpho_female_dS<-morpho_visual_model_females(a30P_data_indiv)
a30P_morpho_male_dS<-morpho_visual_model_males(a30P_data_indiv)
a30P_avian_female_dS<-avian_visual_model_females(a30P_data_indiv)
a30P_avian_male_dS<-avian_visual_model_males(a30P_data_indiv)

#a30I chromatic distances calculation
a30I_morpho_female_dS<-morpho_visual_model_females(a30I_data_indiv)
a30I_morpho_male_dS<-morpho_visual_model_males(a30I_data_indiv)
a30I_avian_female_dS<-avian_visual_model_females(a30I_data_indiv)
a30I_avian_male_dS<-avian_visual_model_males(a30I_data_indiv)

#a30E chromatic distances calculation
a30E_morpho_female_dS<-morpho_visual_model_females(a30E_data_indiv)
a30E_morpho_male_dS<-morpho_visual_model_males(a30E_data_indiv)
a30E_avian_female_dS<-avian_visual_model_females(a30E_data_indiv)
a30E_avian_male_dS<-avian_visual_model_males(a30E_data_indiv)

#a45A chromatic distances calculation
a45A_morpho_female_dS<-morpho_visual_model_females(a45A_data_indiv)
a45A_morpho_male_dS<-morpho_visual_model_males(a45A_data_indiv)
a45A_avian_female_dS<-avian_visual_model_females(a45A_data_indiv)
a45A_avian_male_dS<-avian_visual_model_males(a45A_data_indiv)

#a45P chromatic distances calculation
a45P_morpho_female_dS<-morpho_visual_model_females(a45P_data_indiv)
a45P_morpho_male_dS<-morpho_visual_model_males(a45P_data_indiv)
a45P_avian_female_dS<-avian_visual_model_females(a45P_data_indiv)
a45P_avian_male_dS<-avian_visual_model_males(a45P_data_indiv)

#a45I chromatic distances calculation
a45I_morpho_female_dS<-morpho_visual_model_females(a45I_data_indiv)
a45I_morpho_male_dS<-morpho_visual_model_males(a45I_data_indiv)
a45I_avian_female_dS<-avian_visual_model_females(a45I_data_indiv)
a45I_avian_male_dS<-avian_visual_model_males(a45I_data_indiv)

#a45E chromatic distances calculation
a45E_morpho_female_dS<-morpho_visual_model_females(a45E_data_indiv)
a45E_morpho_male_dS<-morpho_visual_model_males(a45E_data_indiv)
a45E_avian_female_dS<-avian_visual_model_females(a45E_data_indiv)
a45E_avian_male_dS<-avian_visual_model_males(a45E_data_indiv)


##### PLOT THE VISUAL MODEL RESULTS #######

#I have summed-up the chromatic distance measured between 
# helenor vs achilles males, helenor vs achilles females and
#bristowi vs theodorus males, bristowi vs theodorus females
#using morpho and uv-sensitive visual models
# in one file :
contrast_variation<-read.table("contrast_var_all.csv", header=TRUE, sep=";", dec=",")

#separate proximo-distal and antero-posterior data for the plot
contrast_variation_AP<-contrast_variation %>% 
  filter(grepl("P|A|0", Angle)) %>% 
  filter(!grepl("30E", Angle)) %>% 
  filter(!grepl("30I", Angle))
contrast_variation_AP$Angle<-factor(contrast_variation_AP$Angle, levels=c("45P","30P","15P","0","15A","30A","45A"))

contrast_variation_EI<-contrast_variation %>% 
  filter(grepl("E|I|0", Angle)) %>% 
  filter(!grepl("30P", Angle)) %>% 
  filter(!grepl("30A", Angle))
contrast_variation_EI$Angle<-factor(contrast_variation_EI$Angle, levels=c("45I","30I","15I","0","15E","30E","45E"))

#chromatic distances measured on the antero-posterior plane plot
ggplot(contrast_variation_AP, aes(x=Angle, y=dS, group=Type, color=Type)) +
  # Error bars with conditional transparency
  geom_errorbar(aes(ymin=dS_lwr, ymax=dS_upr, alpha = ifelse(dS < 1, 0.3, 1)), 
                width = .1, position = position_dodge(0.3)) +
  geom_line() +
  # Points with conditional transparency
  geom_point(aes(alpha = ifelse(dS < 1, 0.3, 1)), size = 4, position = position_dodge(0.3)) +
  scale_alpha_continuous(range = c(0.2, 1)) +  # Ensure proper scaling for alpha
  scale_color_manual(values = c("#FF3333", "#00CC66", "gray57")) +
  theme_classic() +
  # Horizontal line with transparency
  geom_hline(yintercept = 1, linetype = "dotted", col = "black", size = 1, alpha = 0.5) +
  ylim(0, 5.1) +
  facet_grid(Sex ~ Comparison)


#chromatic distances measured on the proximo-distal plane plot
ggplot(contrast_variation_EI, aes(x=Angle, y=dS, group=Type, color=Type)) +
  # Error bars with conditional transparency
  geom_errorbar(aes(ymin=dS_lwr, ymax=dS_upr, alpha = ifelse(dS < 1, 0.3, 1)), 
                width = .1, position = position_dodge(0.3)) +
  geom_line() +
  # Points with conditional transparency
  geom_point(aes(alpha = ifelse(dS < 1, 0.3, 1)), size = 4, position = position_dodge(0.3)) +
  scale_alpha_continuous(range = c(0.2, 1)) +  # Ensure proper scaling for alpha
  scale_color_manual(values = c("#FF3333", "#00CC66", "gray57")) +
  theme_classic() +
  # Horizontal line with transparency
  geom_hline(yintercept = 1, linetype = "dotted", col = "black", size = 1, alpha = 0.5) +
  ylim(0, 5.1) +
  facet_grid(Sex ~ Comparison)
