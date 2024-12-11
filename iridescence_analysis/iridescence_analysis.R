rm(list=ls())
###################### PACKAGES ####################################################

library(pavo)
library(tidyverse)
library(rptR)
library(lme4)
library(Matrix)
library(magrittr)
library(purrr)
library(listviewer)
library(repurrrsive)
library(tidyr)
library(tibble)
library(data.table)
library(stats)
library(reshape2)
library(cluster)
library(devtools)
library(ggfortify)
library(devtools)
library(factoextra)
library(mvnormtest)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(FactoMineR)


###################### FUNCTIONS ##################################################

extract_repet <- function(angle)
  #' calculates repeatability per individual for 1 angle 
  #' and saves it in a data frame
  #' 
  #' angle = enter angle type as written in ID, must be a character entry
  #' 
  #' ex : extract_repet('_15A_')
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, rptR, tidyr
  #' 
{ 
  angle_type=angle
  
  #1
  #create a dataframe with the selected angle
  dataX<-data_complet %>% 
    filter(grepl(angle_type, ID))
  
  #2
  #measure the repeatability per individual for the variable Chroma
  rptChroma<-rpt(Chroma ~ (1 | indiv), grname = "indiv", data = dataX, datatype = "Gaussian", nboot = 500)
  #measure the repeatability per individual for the variable Brightness
  rptBrightness<-rpt(Brightness ~ (1 | indiv), grname = "indiv", data = dataX, datatype = "Gaussian", nboot = 500)
  #measure the repeatability per individual for the variable Hue
  rptHue<-rpt(Hue ~ (1 | indiv), grname = "indiv", data = dataX, datatype = "Gaussian", nboot = 500)
  
  #3
  #extract the rptChroma list results in a Chroma_list object
  sum_Chroma = unlist(rptChroma)
  Chroma_list<-sum_Chroma[c("R.indiv","se.se","CI_emp.2.5%","CI_emp.97.5%","P.LRT_P")]
  #extract the rptChroma_ list results in a Brightness_list object
  sum_Brightness = unlist(rptBrightness)
  Brightness_list<-sum_Brightness[c("R.indiv","se.se","CI_emp.2.5%","CI_emp.97.5%","P.LRT_P")]
  #extract the rptHue list results in a Hue_list object
  sum_Hue = unlist(rptHue)
  Hue_list<-sum_Hue[c("R.indiv","se.se","CI_emp.2.5%","CI_emp.97.5%","P.LRT_P")]
  
  #4
  #set Chroma_list as dataframe and add Variable column
  Chroma_df1<-data.frame(Chroma_list)
  Chroma_df<-Chroma_df1 %>% 
    mutate(Variable='Chroma')
  #set Brightness_list as dataframe
  Brightness_df1<-data.frame(Brightness_list)
  Brightness_df<-Brightness_df1 %>% 
    mutate(Variable='Brightness')
  #set Chroma_list as dataframe
  Hue_df1<-data.frame(Hue_list)
  Hue_df<-Hue_df1 %>% 
    mutate(Variable='Hue')
  
  #5
  #merge all dataframes together
  Repet_df1<-bind_rows(Chroma_df, Brightness_df)
  Repet_df<-bind_rows(Repet_df1, Hue_df)
  #rename columns
  Repet_df_rename<-Repet_df %>% 
    rename("R"="R.indiv", 
           "SE"="se.se", 
           "Emp.2.5%"="CI_emp.2.5.", 
           "Emp.97.5%"="CI_emp.97.5.", 
           "P-value"="P.LRT_P", 
           "Variable"="Variable")
  #add angle ID column
  Repet_df_final<-Repet_df_rename %>% 
    mutate(Angle_ID=angle_type)
  
  return(Repet_df_final)
}

extract_mean <- function(df)
  #' extract mean Brightness, Hue, Chroma 
  #' and saves it in a data frame
  #' 
  #' df = dataframe with the data
  #' 
  #' ex : extract_mean(data_complet_AP)
  #' 
  #' REQUIRED PACKAGES : tidyverse, dplyr, tidyr
  #' 
{ 
  df=data_frame
  ##Brightness
  #1 Compute the global mean of the data
  global_mean_B<-data_frame %>% 
    group_by(angle_detail) %>% 
    summarize(average_B=mean(Brightness)) %>% 
    ungroup()
  #2 Compute the achilles_female mean
  data_achilles_female_B<-data_frame %>% 
    filter(sex_species=="achilles_female")
  mean_F_achilles_B<-data_achilles_female_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_achilles_female=mean(Brightness)) %>% 
    ungroup()
  #3 Compute the achilles_male mean
  data_achilles_male_B<-data_frame %>% 
    filter(sex_species=="achilles_male")
  mean_M_achilles_B<-data_achilles_male_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_achilles_male=mean(Brightness)) %>% 
    ungroup()
  #4 Compute the bristowi_female mean
  data_bristowi_female_B<-data_frame %>% 
    filter(sex_species=="bristowi_female")
  mean_F_bristowi_B<-data_bristowi_female_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_bristowi_female=mean(Brightness)) %>% 
    ungroup()
  #5 Compute the bristowi_male mean
  data_bristowi_male_B<-data_frame %>% 
    filter(sex_species=="bristowi_male")
  mean_M_bristowi_B<-data_bristowi_male_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_bristowi_male=mean(Brightness)) %>% 
    ungroup()
  #6 Compute the helenor_female mean
  data_helenor_female_B<-data_frame %>% 
    filter(sex_species=="helenor_female")
  mean_F_helenor_B<-data_helenor_female_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_helenor_female=mean(Brightness)) %>% 
    ungroup()
  #7 Compute the helenor_male mean
  data_helenor_male_B<-data_frame %>% 
    filter(sex_species=="helenor_male")
  mean_M_helenor_B<-data_helenor_male_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_helenor_male=mean(Brightness)) %>% 
    ungroup()
  #8 Compute the theodorus_female mean
  data_theodorus_female_B<-data_frame %>% 
    filter(sex_species=="theodorus_female")
  mean_F_theodorus_B<-data_theodorus_female_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_theodorus_female=mean(Brightness)) %>% 
    ungroup()
  #9 Compute the theodorus_male mean
  data_theodorus_male_B<-data_frame %>% 
    filter(sex_species=="theodorus_male")
  mean_M_theodorus_B<-data_theodorus_male_B%>% 
    group_by(angle_detail) %>% 
    summarize(B_theodorus_male=mean(Brightness)) %>% 
    ungroup()
  
  ##Hue
  #1 Compute the global mean of the data
  global_mean_H<-data_frame %>% 
    group_by(angle_detail) %>% 
    summarize(average_H=mean(Hue)) %>% 
    ungroup()
  #2 Compute the achilles_female mean
  data_achilles_female_H<-data_frame %>% 
    filter(sex_species=="achilles_female")
  mean_F_achilles_H<-data_achilles_female_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_achilles_female=mean(Hue)) %>% 
    ungroup()
  #3 Compute the achilles_male mean
  data_achilles_male_H<-data_frame %>% 
    filter(sex_species=="achilles_male")
  mean_M_achilles_H<-data_achilles_male_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_achilles_male=mean(Hue)) %>% 
    ungroup()
  #4 Compute the bristowi_female mean
  data_bristowi_female_H<-data_frame %>% 
    filter(sex_species=="bristowi_female")
  mean_F_bristowi_H<-data_bristowi_female_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_bristowi_female=mean(Hue)) %>% 
    ungroup()
  #5 Compute the bristowi_male mean
  data_bristowi_male_H<-data_frame %>% 
    filter(sex_species=="bristowi_male")
  mean_M_bristowi_H<-data_bristowi_male_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_bristowi_male=mean(Hue)) %>% 
    ungroup()
  #6 Compute the helenor_female mean
  data_helenor_female_H<-data_frame %>% 
    filter(sex_species=="helenor_female")
  mean_F_helenor_H<-data_helenor_female_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_helenor_female=mean(Hue)) %>% 
    ungroup()
  #7 Compute the helenor_male mean
  data_helenor_male_H<-data_frame %>% 
    filter(sex_species=="helenor_male")
  mean_M_helenor_H<-data_helenor_male_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_helenor_male=mean(Hue)) %>% 
    ungroup()
  #8 Compute the theodorus_female mean
  data_theodorus_female_H<-data_frame %>% 
    filter(sex_species=="theodorus_female")
  mean_F_theodorus_H<-data_theodorus_female_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_theodorus_female=mean(Hue)) %>% 
    ungroup()
  #9 Compute the theodorus_male mean
  data_theodorus_male_H<-data_frame %>% 
    filter(sex_species=="theodorus_male")
  mean_M_theodorus_H<-data_theodorus_male_H%>% 
    group_by(angle_detail) %>% 
    summarize(H_theodorus_male=mean(Hue)) %>% 
    ungroup()
  
  ##Chroma
  #1 Compute the global mean of the data
  global_mean_C<-data_frame %>% 
    group_by(angle_detail) %>% 
    summarize(average_C=mean(Chroma)) %>% 
    ungroup()
  #2 Compute the achilles_female mean
  data_achilles_female_C<-data_frame %>% 
    filter(sex_species=="achilles_female")
  mean_F_achilles_C<-data_achilles_female_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_achilles_female=mean(Chroma)) %>% 
    ungroup()
  #3 Compute the achilles_male mean
  data_achilles_male_C<-data_frame %>% 
    filter(sex_species=="achilles_male")
  mean_M_achilles_C<-data_achilles_male_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_achilles_male=mean(Chroma)) %>% 
    ungroup()
  #4 Compute the bristowi_female mean
  data_bristowi_female_C<-data_frame %>% 
    filter(sex_species=="bristowi_female")
  mean_F_bristowi_C<-data_bristowi_female_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_bristowi_female=mean(Chroma)) %>% 
    ungroup()
  #5 Compute the bristowi_male mean
  data_bristowi_male_C<-data_frame %>% 
    filter(sex_species=="bristowi_male")
  mean_M_bristowi_C<-data_bristowi_male_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_bristowi_male=mean(Chroma)) %>% 
    ungroup()
  #6 Compute the helenor_female mean
  data_helenor_female_C<-data_frame %>% 
    filter(sex_species=="helenor_female")
  mean_F_helenor_C<-data_helenor_female_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_helenor_female=mean(Chroma)) %>% 
    ungroup()
  #7 Compute the helenor_male mean
  data_helenor_male_C<-data_frame %>% 
    filter(sex_species=="helenor_male")
  mean_M_helenor_C<-data_helenor_male_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_helenor_male=mean(Chroma)) %>% 
    ungroup()
  #8 Compute the theodorus_female mean
  data_theodorus_female_C<-data_frame %>% 
    filter(sex_species=="theodorus_female")
  mean_F_theodorus_C<-data_theodorus_female_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_theodorus_female=mean(Chroma)) %>% 
    ungroup()
  #9 Compute the theodorus_male mean
  data_theodorus_male_C<-data_frame %>% 
    filter(sex_species=="theodorus_male")
  mean_M_theodorus_C<-data_theodorus_male_C%>% 
    group_by(angle_detail) %>% 
    summarize(C_theodorus_male=mean(Chroma)) %>% 
    ungroup()
  
  #merge all dataframes
  data<-global_mean_B %>% 
    left_join(mean_F_achilles_B, by="angle_detail") %>% 
    left_join(mean_M_achilles_B, by="angle_detail") %>%
    left_join(mean_F_bristowi_B, by="angle_detail") %>%
    left_join(mean_M_bristowi_B, by="angle_detail") %>%
    left_join(mean_F_helenor_B, by="angle_detail") %>%
    left_join(mean_M_helenor_B, by="angle_detail") %>%
    left_join(mean_F_theodorus_B, by="angle_detail") %>%
    left_join(mean_M_theodorus_B, by="angle_detail") %>%
    left_join(global_mean_C, by="angle_detail") %>%
    left_join(mean_F_achilles_C, by="angle_detail") %>% 
    left_join(mean_M_achilles_C, by="angle_detail") %>%
    left_join(mean_F_bristowi_C, by="angle_detail") %>%
    left_join(mean_M_bristowi_C, by="angle_detail") %>%
    left_join(mean_F_helenor_C, by="angle_detail") %>%
    left_join(mean_M_helenor_C, by="angle_detail") %>%
    left_join(mean_F_theodorus_C, by="angle_detail") %>%
    left_join(mean_M_theodorus_C, by="angle_detail") %>%
    left_join(global_mean_H, by="angle_detail") %>%
    left_join(mean_F_achilles_H, by="angle_detail") %>% 
    left_join(mean_M_achilles_H, by="angle_detail") %>%
    left_join(mean_F_bristowi_H, by="angle_detail") %>%
    left_join(mean_M_bristowi_H, by="angle_detail") %>%
    left_join(mean_F_helenor_H, by="angle_detail") %>%
    left_join(mean_M_helenor_H, by="angle_detail") %>%
    left_join(mean_F_theodorus_H, by="angle_detail") %>%
    left_join(mean_M_theodorus_H, by="angle_detail")
  
  return(data)
  
}


##############################################
############ #1 OPEN RAW DATA ################
##############################################

###get data from 300nm to 700nm for the PCA analyses
data_smooth_all<-as.rspec(read.csv("spec_300_700nm.csv", sep=','))
###from 320nm to 700nm for the analysis of Brightness, Hue, Chroma
data_smooth<-as.rspec(read.csv("spec_320_700nm.csv", sep=','))

###############################################################################################
#### #2 GET HUE, BRIGHTNESS, CHROMA AND CREATE CLEAN DATAFRAME WITH DETAILED VARIABLES ########
###############################################################################################

###extract brigthness, chroma, hue
data<-summary(data_smooth, subset=TRUE)
#first column ID
data<-tibble::rownames_to_column(data, "ID")

###create complet detailed dataframe from ID information
### variables :
#repet= repeated measures, each measure was repeated 3 times = repetitions "A", "B", or "C"
#species= species ID : "bristowi"=M.helenor bristowi, "theodorus"=M. helenor theodorus, "helenor"=M. helenor helenor, "achilles"=M. achilles achilles
#location= geographic habitat "ecuador" for allopatric populations, "french guiana" for sympatric populations
#sex= "female" or "male"
#angle= angle of illuminations "0"=0°, "15"=15°, "30"=30°, "45"=45°
#side= side of the wing illuminated "zero"=illumination at the normal of the wings, "anterieur"=on the anterior side, "interieur"=on the internal side, "posterieur"=on the posterior side, "exterieur"=on the external side
#angle_detail= angle of illumination ("0","15","30","45") + side of the wing illuminated ("A", "I", "P", "E") + discrimination between "specular" and "tilt" measurements (tilt measurements are identified with a "1" at the end of ID)
#population= species + sex ("achilles_male","achilles_female", "helenor_male", "helenor_female",....)
#inviv= each sample has a unique number from 1 to 80
data_complet<-data %>% 
  mutate (repet = case_when ( 
    grepl("_repetA", ID)~"A",
    grepl("_repetB", ID)~"B",
    grepl("_repetC", ID)~"C")) %>%
  mutate(species=case_when (
    grepl("_HB_", ID)~"bristowi",
    grepl("_HT_", ID)~"theodorus",
    grepl("_HEL_", ID)~"helenor",
    grepl("_ACH_", ID)~"achilles")) %>%
  mutate(location=case_when (
    grepl("PA23_", ID)~"ecuador",
    grepl("FG23_", ID)~"french_guiana")) %>% 
  mutate(sex=case_when (
    grepl("_F_", ID)~"female",
    grepl("_M_", ID)~"male")) %>%
  mutate(angle=case_when (
    grepl("_0_", ID)~"0",
    grepl("_15A_", ID)~"15",
    grepl("_15E_", ID)~"15",
    grepl("_15I_", ID)~"15",
    grepl("_15P_", ID)~"15",
    grepl("_30A_", ID)~"30",
    grepl("_30E_", ID)~"30",
    grepl("_30I_", ID)~"30",
    grepl("_30P_", ID)~"30",
    grepl("_30A1_", ID)~"30",
    grepl("_30E1_", ID)~"30",
    grepl("_30I1_", ID)~"30",
    grepl("_30P1_", ID)~"30",
    grepl("_30Am1_", ID)~"30",
    grepl("_30Em1_", ID)~"30",
    grepl("_30Im1_", ID)~"30",
    grepl("_30Pm1_", ID)~"30",
    grepl("_45A_", ID)~"45",
    grepl("_45E_", ID)~"45",
    grepl("_45I_", ID)~"45",
    grepl("_45P_", ID)~"45")) %>%
  mutate(angle_detail=case_when (
    grepl("_0_", ID)~"0",
    grepl("_15A_", ID)~"15A",
    grepl("_15E_", ID)~"15E",
    grepl("_15I_", ID)~"15I",
    grepl("_15P_", ID)~"15P",
    grepl("_30A_", ID)~"30A",
    grepl("_30E_", ID)~"30E",
    grepl("_30I_", ID)~"30I",
    grepl("_30P_", ID)~"30P",
    grepl("_30A1_", ID)~"30A1",
    grepl("_30E1_", ID)~"30E1",
    grepl("_30I1_", ID)~"30I1",
    grepl("_30P1_", ID)~"30P1",
    grepl("_30Am1_", ID)~"30Am1",
    grepl("_30Em1_", ID)~"30Em1",
    grepl("_30Im1_", ID)~"30Im1",
    grepl("_30Pm1_", ID)~"30Pm1",
    grepl("_45A_", ID)~"45A",
    grepl("_45E_", ID)~"45E",
    grepl("_45I_", ID)~"45I",
    grepl("_45P_", ID)~"45P")) %>%
  mutate(side=case_when (
    grepl("_0_", ID)~"zero",
    grepl("A_", ID)~"anterieur",
    grepl("A1_", ID)~"anterieur",
    grepl("Am1_", ID)~"anterieur",
    grepl("P_", ID)~"posterieur",
    grepl("P1_", ID)~"posterieur",
    grepl("Pm1_", ID)~"posterieur",
    grepl("E_", ID)~"exterieur",
    grepl("E1_", ID)~"exterieur",
    grepl("Em1_", ID)~"exterieur",
    grepl("I_", ID)~"interieur",
    grepl("I1_", ID)~"interieur",
    grepl("Im1_", ID)~"interieur")) %>%
  mutate(antero_posterieur=case_when (
    grepl("zero", side)~"yes",
    grepl("anterieur", side)~"yes",
    grepl("posterieur", side)~"yes")) %>% 
  mutate(ext_inter=case_when (
    grepl("zero", side)~"yes",
    grepl("exterieur", side)~"yes",
    grepl("interieur", side)~"yes")) %>%
  rename(Brightness = B2) %>% 
  rename(Chroma = S8) %>%
  rename(Hue = H1)

###create population variable
data_complet<-data_complet %>% unite("population", species:sex, remove=F)
###numerate individuals
y<-rep(c(1:80),each=3, times=21)
data_complet<-data_complet %>% mutate(indiv=y)
data_complet$indiv<-as.character(data_complet$indiv)

############################################################################################
################### 3) REPEATABILITY FOR HUE, BRIGHTNESS AND CHROMA ########################
############################################################################################

###use the extract_repet function created at the beginning of the file to measure repeatability for each angle
###create a list of all the angles to study
angle_list<-list('_0_', 
                 '_15A_','_15E_','_15I_','_15P_',
                 '_30A_','_30E_','_30I_','_30P_',
                 '_45A_','_45E_','_45I_','_45P_',
                 '_30A1_','_30E1_','_30I1_','_30P1_',
                 '_30Am1_','_30Em1_','_30Im1_','_30Pm1_')

###create the data frame to store the results
Repet_df<-data.frame()

###for loop to apply the extract_repet function to the angle list
for(i in angle_list){
  repet_i<-extract_repet(i)
  Repet_df<-rbind(Repet_df, repet_i)
}

######################################################################################
############## #4 PLOT HUE, BRIGHTNESS AND CHROMA ####################################
######################################################################################

### Analysis of the SPECULAR data only
### Divide the dataset into 2 planes : PROXIMODISTAL (EI) and ANTEROPOSTERIOR (AP) for SPECULAR data only
#create PROXIMODISTAL (EI) data vector
data_complet_EI<-data_complet %>% unite("sex_species", c(species,sex), remove=F) %>% 
  filter(grepl('zero|exterieur|interieur', side))
data_complet_EI$angle_detail<-as.character(data_complet_EI$angle_detail)
data_complet_EI <- data_complet_EI %>% 
  filter(angle_detail != "30E1") %>% 
  filter(angle_detail != "30Em1") %>% 
  filter(angle_detail != "30I1") %>% 
  filter(angle_detail != "30Im1")
#order on graphs
data_complet_EI$angle_detail<-factor(data_complet_EI$angle_detail, levels=c("45I","30I","15I","0","15E","30E","45E"))
#create a specific PROXIMODISTAL (EI) data vector for FRENCH GUIANA
data_complet_EI_guiana<-data_complet_EI %>% 
  filter(location=="french_guiana")
#create a specific PROXIMODISTAL (EI) data vector for ECUADOR
data_complet_EI_ecuador<-data_complet_EI %>% 
  filter(location=="ecuador")

#create ANTEROPOSTERIOR (AP) data vector
data_complet_AP<-data_complet %>% unite("sex_species", c(species,sex), remove=F) %>% 
  filter(grepl('zero|anterieur|posterieur', side))
data_complet_AP$angle_detail<-as.character(data_complet_AP$angle_detail)
data_complet_AP <- data_complet_AP %>% filter(angle_detail != "30A1") %>% 
  filter(angle_detail != "30Am1") %>% 
  filter(angle_detail != "30P1") %>% 
  filter(angle_detail != "30Pm1")

#order on graphs
data_complet_AP$angle_detail<-factor(data_complet_AP$angle_detail, levels=c("45P","30P","15P","0","15A","30A","45A"))
#create a specific ANTEROPOSTERIOR (AP) data vector for FRENCH GUIANA
data_complet_AP_guiana<-data_complet_AP %>% 
  filter(location=="french_guiana")
#create a specific ANTEROPOSTERIOR (AP) data vector for ECUADOR
data_complet_AP_ecuador<-data_complet_AP %>% 
  filter(location=="ecuador")

#calculate the mean brightness, hue and chroma value for each individual on the ANTEROPOSTERIOR plane
data_frame=data_complet_AP
mean_AP<-extract_mean(data_complet_AP)
#create specific vector for FRENCH GUIANA
mean_AP_guiana<-mean_AP %>% 
  select(contains("achilles") | contains("helenor")|contains("angle_detail"))
#create specific vector for ECUADOR
mean_AP_ecuador<-mean_AP %>% 
  select(contains("bristowi") | contains("theodorus")|contains("angle_detail"))

#calculate the mean brightness, hue and chroma value for each individual on the PROXIMODISTAL plane
data_frame=data_complet_EI
mean_EI<-extract_mean(data_complet_EI)
#create specific vector for FRENCH GUIANA
mean_EI_guiana<-mean_EI %>% 
  select(contains("achilles") | contains("helenor")|contains("angle_detail"))
#create specific vector for ECUADOR
mean_EI_ecuador<-mean_EI %>% 
  select(contains("bristowi") | contains("theodorus")|contains("angle_detail"))

#plot the boxplot with the trajectories
#change the parameters depending on the boxplots you want to generate
#example for HUE from the ANTERIOPOSTERIOR plane measured on the wings of species from FRENCH GUIANA
data_complet_AP_guiana %>% ggplot( aes(x=angle_detail, y=Hue)) +
  geom_boxplot() +
  geom_point(aes(colour = sex_species, shape=sex, size=sex, alpha=sex_species),position = position_jitterdodge())+theme_classic()+
  scale_color_manual(values=c("achilles_female"="darkolivegreen3","achilles_male"="#339933","helenor_female"="#33CCFF","helenor_male"="#3333FF"))+
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3))+
  scale_size_manual(values=c(1,2))+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_achilles_female), color="darkolivegreen3", size=3)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_achilles_female, group=1), size=1, color="darkolivegreen3")+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_achilles_male), color="#339933", size=3, shape=17)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_achilles_male, group=1), size=1, color="#339933")+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_helenor_female), color="#33CCFF", size=3)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_helenor_female, group=1), size=1, color="#33CCFF")+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_helenor_male), color="#3333FF", size=3, shape=17)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=H_helenor_male, group=1), size=1, color="#3333FF")+
  ylim(c(350,565))

#example for HUE from the ANTERIOPOSTERIOR plane measured on the wings of species from ECUADOR
data_complet_AP_ecuador %>% ggplot( aes(x=angle_detail, y=Hue)) +
  geom_boxplot() +
  geom_point(aes(colour = sex_species, shape=sex, size=sex, alpha=sex_species),position = position_jitterdodge())+theme_classic()+
  scale_color_manual(values=c("bristowi_female"="#CC99FF","bristowi_male"="#660099","theodorus_female"="#FF9966","theodorus_male"="#CC3333"))+
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3))+
  scale_size_manual(values=c(1,2))+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_bristowi_female), color="#CC99FF", size=3)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_bristowi_female, group=1), size=1, color="#CC99FF")+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_bristowi_male), color="#660099", size=3, shape=17)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_bristowi_male, group=1), size=1, color="#660099")+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_theodorus_female), color="#FF9966", size=3)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_theodorus_female, group=1), size=1, color="#FF9966")+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_theodorus_male), color="#CC3333", size=3, shape=17)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=H_theodorus_male, group=1), size=1, color="#CC3333")+
  ylim(c(350,565))

### Analysis of the TILT data only
### Divide the dataset into 2 planes : PROXIMODISTAL and ANTEROPOSTERIOR for TILT data only
#create PROXIMODISTAL (EI) data vector
data_complet_EI<-data_complet %>% unite("sex_species", c(species,sex), remove=F) %>% 
  filter(angle==30)
data_complet_EI$angle_detail<-as.character(data_complet_EI$angle_detail)
data_complet_EI <- data_complet_EI %>% 
  filter(angle_detail=="30E1"|angle_detail=="30Em1"|angle_detail=="30E"|
           angle_detail=="30I"|angle_detail=="30I1"|angle_detail=="30Im1")

#order on graphs
data_complet_EI$angle_detail<-factor(data_complet_EI$angle_detail, levels=c("30Im1","30I","30I1","30E1","30E","30Em1"))
#create a specific PROXIMODISTAL (EI) data vector for FRENCH GUIANA
data_complet_EI_guiana<-data_complet_EI %>% 
  filter(location=="french_guiana")
#create a specific PROXIMODISTAL (EI) data vector for ECUADOR
data_complet_EI_ecuador<-data_complet_EI %>% 
  filter(location=="ecuador")

#create ANTEROPOSTERIOR (AP) data vector
data_complet_AP<-data_complet %>% unite("sex_species", c(species,sex), remove=F) %>% 
  filter(angle==30)
data_complet_AP$angle_detail<-as.character(data_complet_AP$angle_detail)
data_complet_AP <- data_complet_AP %>% 
  filter(angle_detail=="30A1"|angle_detail=="30Am1"|angle_detail=="30A"|
           angle_detail=="30P"|angle_detail=="30P1"|angle_detail=="30Pm1")

#order on graphs
data_complet_AP$angle_detail<-factor(data_complet_AP$angle_detail, levels=c("30Pm1","30P","30P1","30A1","30A","30Am1"))
#create a specific ANTEROPOSTERIOR (AP) data vector for FRENCH GUIANA
data_complet_AP_guiana<-data_complet_AP %>% 
  filter(location=="french_guiana")
#create a specific ANTEROPOSTERIOR (AP) data vector for ECUADOR
data_complet_AP_ecuador<-data_complet_AP %>% 
  filter(location=="ecuador")

#calculate the mean brightness, hue and chroma value for each individual on the ANTEROPOSTERIOR plane
data_frame=data_complet_AP
mean_AP<-extract_mean(data_complet_AP)
#create specific vector for FRENCH GUIANA
mean_AP_guiana<-mean_AP %>% 
  select(contains("achilles") | contains("helenor")|contains("angle_detail"))
#create specific vector for ECUADOR
mean_AP_ecuador<-mean_AP %>% 
  select(contains("bristowi") | contains("theodorus")|contains("angle_detail"))

#calculate the mean brightness, hue and chroma value for each individual on the PROXIMODISTAL plane
data_frame=data_complet_EI
mean_EI<-extract_mean(data_complet_EI)
#create specific vector for FRENCH GUIANA
mean_EI_guiana<-mean_EI %>% 
  select(contains("achilles") | contains("helenor")|contains("angle_detail"))
#create specific vector for ECUADOR
mean_EI_ecuador<-mean_EI %>% 
  select(contains("bristowi") | contains("theodorus")|contains("angle_detail"))

#plot the boxplot with the trajectories
#change the parameters depending on the boxplots you want to generate
#example for BRIGHTNESS from the ANTERIOPOSTERIOR plane measured on the wings of species from FRENCH GUIANA
data_complet_AP_guiana %>% ggplot( aes(x=angle_detail, y=Brightness)) +
  geom_boxplot() +
  geom_point(aes(colour = sex_species, shape=sex, size=sex, alpha=sex_species),position = position_jitterdodge())+theme_classic()+
  scale_color_manual(values=c("achilles_female"="darkolivegreen3","achilles_male"="#339933","helenor_female"="#33CCFF","helenor_male"="#3333FF"))+
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3))+
  scale_size_manual(values=c(1,2))+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_achilles_female), color="darkolivegreen3", size=3)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_achilles_female, group=1), size=1, color="darkolivegreen3")+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_achilles_male), color="#339933", size=3, shape=17)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_achilles_male, group=1), size=1, color="#339933")+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_helenor_female), color="#33CCFF", size=3)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_helenor_female, group=1), size=1, color="#33CCFF")+
  geom_point(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_helenor_male), color="#3333FF", size=3, shape=17)+
  geom_line(data=mean_AP_guiana, mapping=aes(x=angle_detail, y=B_helenor_male, group=1), size=1, color="#3333FF")+
  ylim(c(10,55))
#+facet_grid(location ~ .)

#example for BRIGHTNESS from the ANTERIOPOSTERIOR plane measured on the wings of species from ECUADOR
data_complet_AP_ecuador %>% ggplot( aes(x=angle_detail, y=Brightness)) +
  geom_boxplot() +
  geom_point(aes(colour = sex_species, shape=sex, size=sex, alpha=sex_species),position = position_jitterdodge())+theme_classic()+
  scale_color_manual(values=c("bristowi_female"="#CC99FF","bristowi_male"="#660099","theodorus_female"="#FF9966","theodorus_male"="#CC3333"))+
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3))+
  scale_size_manual(values=c(1,2))+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_bristowi_female), color="#CC99FF", size=3)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_bristowi_female, group=1), size=1, color="#CC99FF")+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_bristowi_male), color="#660099", size=3, shape=17)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_bristowi_male, group=1), size=1, color="#660099")+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_theodorus_female), color="#FF9966", size=3)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_theodorus_female, group=1), size=1, color="#FF9966")+
  geom_point(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_theodorus_male), color="#CC3333", size=3, shape=17)+
  geom_line(data=mean_AP_ecuador, mapping=aes(x=angle_detail, y=B_theodorus_male, group=1), size=1, color="#CC3333")+
  ylim(c(10,55))

######################################################################################
########## #5 STATISTICAL ANALYSIS OF HUE, BRIGHTNESS AND CHROMA #####################
######################################################################################

library(tidyverse)
library(ggpubr)
library(rstatix)

#Example of PERMANOVA performed on the CHROMA data measured on the PROXIMODISTAL plane of the wings of butterflies from FRENCH GUIANA
#change the parameters depending on the variables you want to test
##check assumptions###

##normality assumptions
model<-lm(Chroma~sex*species*angle_detail, data=data_complet_EI_guiana)
#create a QQ plot of residuals
ggqqplot(residuals(model))
#compute shapiro_test
shapiro_test(residuals(model))
#homogeneity of variance assumption
data_complet_AP_guiana %>% levene_test(Hue~sex*species*angle_detail)

kruskal.test(Hue~sex*species*angle_detail, data=data_complet_AP_guiana)

#permuation ANOVA because non normal data
library(vegan)

tab<-adonis(Chroma~sex*species+angle_detail+species*angle_detail, data=data_complet_EI_ecuador, permutations=999)
summary(tab)
tab$aov.tab


################################################
########## #6 PCA OF IRIDESCENCE ###############
################################################


############ Prepare a mean wavelength per indiviual dataframe ##################

# create the mean spectra from the 3 repetitions for every individual
data_moyenne <- aggspec(data_smooth_all, by = 3, FUN = mean)

#1 
## Create a global data frame with ecological and optical variables
## of each individual

###create a data frame with each wavelengths at each angles as variables
specs_transpose<-transpose(data_moyenne)
###clean and edit new data frame
colnames(specs_transpose)<-specs_transpose[1,]
rownames(specs_transpose)<-colnames(data_moyenne)
specs_transpose<-specs_transpose[-1,]
specs_transpose<-specs_transpose %>% 
  mutate(ID=rownames(specs_transpose), .before="300")

###add to the data frame the details of the data_complet data frame
data_complet_specs<-merge(data_complet, specs_transpose, by="ID")
###and remove the Brightness, Chroma, Hue data (we are not going to use it)
data_complet_specs<-data_complet_specs %>% 
  select( -"Brightness", -"Chroma", -"Hue") %>%
  unite(indiv_var, c("indiv","repet"), remove=F) %>% 
  arrange(population)

#2
## Create a data frame with individual ID and angle details

data_complet_numerical<-data_complet_specs[,c(2,9,12:410)]
## Extend data frame : each column now stands for 1nm (between 300 and 700nm) 
## at a specific angle, for the 80 studied individuals
df.wide <- reshape(data_complet_numerical, idvar="indiv_var", timevar="angle_detail", direction="wide")
## Create a list of ecological characteristics to add to the new df.wide dataframe
## depending on the data_complet_specs topology
Population<-rep(c("achilles_french_guiana_female","achilles_french_guiana_male",
                  "bristowi_ecuador_female","bristowi_ecuador_male","helenor_french_guiana_female",
                  "helenor_french_guiana_male","theodorus_ecuador_female","theodorus_ecuador_male"), each=10,time=1)
Espece<-rep(c("achilles","bristowi","helenor","theodorus"), each=20,time=1)
Localite<-rep(c("french_guiana","ecuador"), each=20,time=2)
Sexe<-rep(c("female","male"), each=10,time=4)

## Creation of the ultimate data frame with 
## the reflectance value for each wavelength
## at each measured angle
## for the mean spectra of each individual (80)
data_acp_total<-cbind(Population, Espece, Localite, Sexe, df.wide)
Indiv<-data_acp_total %>% 
  pull(indiv_var)

#delete all residual and useless columns for the analysis
data_acp_total<-data_acp_total %>% 
  select(!c(ext_inter.0,indiv.0)) %>% 
  select(!starts_with("ext_")) %>%
  select(!starts_with("indiv."))

#####################################
##### PCA on ECUADOR data only ######
#####################################

## Select the ECUADOR data only
data_acp_ecuador<-data_acp_total %>% 
  filter(Localite == "ecuador")
## Select only the numerical data for the PCA
data_acp_ecuador_num<-data_acp_ecuador %>% 
  select_if(is.numeric)

## Compute PCA
pca_ecuador<-prcomp(data_acp_ecuador_num, scale=FALSE)
fviz_pca_var(pca_ecuador,
             axes=c(1,2))
fviz_eig(pca_ecuador)

## Plot PCA
autoplot(pca_ecuador, data = data_acp_ecuador,frame=TRUE, colour="Population",size=4, 
         scale=0)+theme_classic()+
  scale_color_manual(values=c("#CC99FF","#660099","#FF9966","#CC3333"))+
  scale_fill_manual(values=c("#CC99FF","#660099","#FF9966","#CC3333"))

##########################################
##### PCA on FRENCH GUIANA data only #####
##########################################

## Select the FRENCH GUIANA data only
data_acp_guiana<-data_acp_total %>% 
  filter(Localite == "french_guiana")
## Select only the numerical data for the PCA
data_acp_guiana_num<-data_acp_guiana %>% 
  select_if(is.numeric)

## Compute PCA
pca_guiana<-prcomp(data_acp_guiana_num, scale=FALSE)
fviz_pca_var(pca_guiana,
             axes=c(1,2))
fviz_eig(pca_guiana)
## Plot PCA
autoplot(pca_guiana, data = data_acp_guiana,frame=TRUE, colour="Population",size=4, 
         scale=0)+theme_classic()+
  scale_color_manual(values=c("darkolivegreen2","#339933","#33CCFF","#3333FF"))+
  scale_fill_manual(values=c("darkolivegreen2","#339933","#33CCFF","#3333FF"))

##############################################################
################ #7 PCA STATISTICAL ANALYSIS #################
##############################################################

####################################
### Statistics on ECUADOR data #####
####################################

## Exctract the PCA coordinates
ind <- get_pca_ind(pca_ecuador)
## Create a data frame with the coordinates per individual
coord_ecuador<-data.frame(ind$coord)
## Select the 10 first dimensions (where strongest signal lies)
coord_ecuador<-coord_ecuador %>% 
  select(Dim.1:Dim.10)
## Add the qualitative variables back
Species<-data_acp_ecuador %>% 
  pull(Espece)
Sex<-data_acp_ecuador %>% 
  pull(Sexe)
Indiv<-data_acp_ecuador %>% 
  pull(indiv_var)
coord_ecuador<-data.frame(cbind(coord_ecuador,Species,Indiv,Sex))


## Create a vector containing all the variables to test
var<-cbind(coord_ecuador$Dim.1,coord_ecuador$Dim.2,coord_ecuador$Dim.3,
           coord_ecuador$Dim.4,coord_ecuador$Dim.5,coord_ecuador$Dim.6,
           coord_ecuador$Dim.7,coord_ecuador$Dim.8,coord_ecuador$Dim.9,
           coord_ecuador$Dim.10)
## Perform a permutation MANOVA
adonis2( var~ coord_ecuador$Species*coord_ecuador$Sex, method = "euclidean",
         permutations = 999)


##########################################
### Statistics on FRENCH GUIANA data #####
##########################################

## Exctract the PCA coordinates
ind <- get_pca_ind(pca_guiana)
## Create a data frame with the coordinates per individual
coord_guiana<-data.frame(ind$coord)
## Select the 10 first dimensions (where strongest signal lies)
coord_guiana<-coord_guiana %>% 
  select(Dim.1:Dim.10)
## Add the qualitative variables back
Species<-data_acp_guiana %>% 
  pull(Espece)
Sex<-data_acp_guiana %>% 
  pull(Sexe)
Indiv<-data_acp_guiana %>% 
  pull(indiv_var)
coord_guiana<-data.frame(cbind(coord_guiana,Species,Indiv,Sex))


## Create a vector containing all the variables to test
var<-cbind(coord_guiana$Dim.1,coord_guiana$Dim.2,coord_guiana$Dim.3,
           coord_guiana$Dim.4,coord_guiana$Dim.5,coord_guiana$Dim.6,
           coord_guiana$Dim.7,coord_guiana$Dim.8,coord_guiana$Dim.9,
           coord_guiana$Dim.10)
## Perform a permutation MANOVA
adonis2( var~ coord_guiana$Species*coord_guiana$Sex, method = "euclidean",
         permutations = 999)

