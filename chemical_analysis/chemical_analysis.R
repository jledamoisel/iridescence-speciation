rm(list=ls())

library(tidyverse)
library(vegan)
library(janitor)
library(ggvegan)

################################################
############### C8-C16 ANALYSIS ################
################################################

#load the compound matrix generated with MZmine
table<-read.csv("C8_C16_compounds.csv", sep=";", dec=",")
#remove columns retention time
table<-table[,-2:-3]
#remove blanks
table<-table[,-38:-39]

############ nMDS calculation ################
table<-as.data.frame(t(table))
table<-table %>%
  row_to_names(row_number = 1)
table<-as.matrix(table)

#nmds calculation, the stress value should be less than 0.2
nmds = metaMDS(table, distance = "bray")
nmds

#points=butterflies, red crosses=compounds
plot(nmds)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data_scores = as.data.frame(scores(nmds)$sites)

species<-rep(c("achilles","achilles","helenor","helenor"),times=c(9,9,9,9))
sex<-rep(c("female","male","male","female"),times=c(9,9,9,9))
species_sex<-rep(c("achilles_female","achilles_male","helenor_male","helenor_female"),times=c(9,9,9,9))

data_scores$species <- species
data_scores$sex <- sex
data_scores$species_sex <- species_sex

#nice plot
library("ggplot2")

plot = ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, aes( shape = sex, colour = species_sex))+ 
  labs(x = "NMDS1", y = "NMDS2")  + 
  scale_colour_manual(values = c("darkolivegreen3","#339933","#33CCFF","#3333FF"))+
  scale_fill_manual(values = c("darkolivegreen3","#339933","#33CCFF","#3333FF"))+
  stat_ellipse(geom = "polygon",
               aes(fill = species_sex, colour=species_sex), 
               alpha = 0.25)+
  theme_classic()

plot

########## PERMANOVA analysis #####################
table<-as.data.frame(table)
table$species <- species
table$sex <- sex
table$species_sex <- species_sex

table_matrix<-as.matrix(table[,1:177])


adonis2( table_matrix~ species*sex, data=table, method = "bray",
         permutations = 999)

############# POST-HOC analysis ##########################
#statistical test on males only
table_males<-table %>% 
  filter(sex == "male")

table_matrix_male<-as.matrix(table_males[,1:177])
var<-c(1:177)

adonis2( table_matrix_male~ species, data=table_males, method = "bray",
         permutations = 999)

#statistical test on females only
table_females<-table %>% 
  filter(sex == "female")

table_matrix_female<-as.matrix(table_females[,1:177])
var<-c(1:177)

adonis2( table_matrix_female~ species, data=table_females, method = "bray",
         permutations = 999)


######## Indicator Value Analysis #################
library(indicspecies)

abund<-table[,1:177]
abund_male<-table[-c(1:9),1:177]
abund_male<-abund_male[-c(19:27),1:177]
abund_female<-table[-c(18:27),1:177]
abund_female<-abund_female[-c(10:17),1:177]
species_male<-rep(c("achilles","helenor"),times=c(9,9))
species_female<-rep(c("achilles","helenor"),times=c(9,9))

invf<-multipatt(abund_female, species_female, func="IndVal.g", control=how(nperm=9999))
invm<-multipatt(abund_male, species_male, func="IndVal.g", control=how(nperm=9999))
df_female<-data.frame(invf$sign)
df_male<-data.frame(invm$sign)

################################################
############### C16-C30 ANALYSIS ################
################################################

#load the compound matrix generated with MZmine
table<-read.csv("C16_C30_compounds.csv", sep=";", dec=",")
#remove columns retention time
table<-table[,-2:-3]
#remove blanks
table<-table[,-46:-50]

####### nMDS calculation #############
table<-as.data.frame(t(table))
table<-table %>%
  row_to_names(row_number = 1)
table<-as.matrix(table)

#nmds calculation, the stress value should be less than 0.2
nmds = metaMDS(table, distance = "bray")
nmds

#points=butterflies, red crosses=compounds
plot(nmds)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data_scores = as.data.frame(scores(nmds)$sites)

species<-rep(c("helenor","achilles","achilles","helenor"),times=c(10,12,11,10))
sex<-rep(c("female","female","male","male"),times=c(10,12,11,10))
species_sex<-rep(c("helenor_female","achilles_female","achilles_male","helenor_male"),times=c(10,12,11,10))

data_scores$species <- species
data_scores$sex <- sex
data_scores$species_sex <- species_sex

#nice plot
library("ggplot2")

plot = ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, aes( shape = sex, colour = species_sex))+ 
  labs(x = "NMDS1", y = "NMDS2")  + 
  scale_colour_manual(values = c("darkolivegreen3","#339933","#33CCFF","#3333FF"))+
  scale_fill_manual(values = c("darkolivegreen3","#339933","#33CCFF","#3333FF"))+
  stat_ellipse(geom = "polygon",
               aes(fill = species_sex, colour=species_sex), 
               alpha = 0.25)+
  theme_classic()

plot

########## PERMANOVA analysis #############
table<-as.data.frame(table)
table$species <- species
table$sex <- sex
table$species_sex <- species_sex

table_matrix<-as.matrix(table[,1:141])
var<-c(1:141)

adonis2( table_matrix~ species*sex, data=table, method = "bray",
         permutations = 999)

########## POST-HOC analysis #################
#statistical test on males only
table_males<-table %>% 
  filter(sex == "male")

table_matrix_male<-as.matrix(table_males[,1:141])
var<-c(1:141)

adonis2( table_matrix_male~ species, data=table_males, method = "bray",
         permutations = 999)

#statistical test on females only
table_females<-table %>% 
  filter(sex == "female")

table_matrix_female<-as.matrix(table_females[,1:141])
var<-c(1:141)

adonis2( table_matrix_female~ species, data=table_females, method = "bray",
         permutations = 999)

######## Indicator Value Analysis #################
library(indicspecies)

abund<-table[,1:141]
abund_male<-table[23:43,1:141]
abund_female<-table[1:22,1:141]
species_male<-rep(c("achilles","helenor"),times=c(11,10))
species_female<-rep(c("helenor","achilles"),times=c(10,12))

invf<-multipatt(abund_female, species_female, func="IndVal.g", control=how(nperm=9999))
invm<-multipatt(abund_male, species_male, func="IndVal.g", control=how(nperm=9999))
df_female<-data.frame(invf$sign)
df_male<-data.frame(invm$sign)
