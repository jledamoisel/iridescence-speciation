rm(list=ls())

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(reshape2)

#fonction pour mesurer des moyennes et des ecart-types par groupe
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#ouverture fichier
data<-read.table("guyane_male_choice.csv", dec=".", sep=";", header=T)
data$App_conspe<-as.numeric(data$App_conspe)
data$App_non_conspe<-as.numeric(data$App_non_conspe)

data<- data %>% 
  mutate (Motivation_approches=case_when(Zone_Achilles=='0' & Zone_Helenor=='0' ~ 'no',
                                         TRUE~'yes'))
data$Motivation_approches<-as.factor(data$Motivation_approches)
data<-data %>% 
  mutate(Sum_app=Zone_Achilles+Zone_Helenor)
data$Sum_app<-as.numeric(data$Sum_app)

data_helenor<-data  %>% 
  filter(Species == "helenor")
data_achilles<-data %>%
  filter(Species == "achilles")

sum(data_helenor$Sum_app > 0)
sum(data_achilles$Sum_app > 0)

data_proba <- data %>%
  mutate (ProbaC=(App_conspe)/(App_conspe+App_non_conspe))
data_probaNA <- data_proba %>% drop_na(ProbaC)

data_proba_helenor <- data_helenor %>%
  mutate (ProbaC=(App_conspe)/(App_conspe+App_non_conspe))
data_probaNA_helenor <- data_proba_helenor %>% drop_na(ProbaC)

data_proba_achilles <- data_achilles %>%
  mutate (ProbaC=(App_conspe)/(App_conspe+App_non_conspe))
data_probaNA_achilles <- data_proba_achilles %>% drop_na(ProbaC)

########### Proba d'approche figures ###########################
#sans moyenner par individu, data proba brutes
ggplot(data_proba,aes(x=Species,y=ProbaC, fill=Species, colour=Species))+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic()

#en moyennant les proba par individu
#helenor

dataSD_helenor<-summarySE(data_probaNA_helenor, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_helenor<-dataSD_helenor %>%
  mutate (Species="helenor")
  
#achilles

dataSD_achilles<-summarySE(data_probaNA_achilles, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_achilles<-dataSD_achilles %>%
  mutate (Species="achilles")

#fusion
dataSD<- rbind(dataSD_helenor,dataSD_achilles)

#plot
ggplot(dataSD,aes(x=Species, y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.5)+
  scale_color_manual(values=c("darkolivegreen2", "#33CCFF"))+
  scale_fill_manual(values=c("darkolivegreen2", "#33CCFF"))+theme_classic()+
  coord_flip()

########## Mean preference per species ###########################
############# achilles ###########

ggplot(dataSD_achilles,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD_achilles %>% shapiro_test(ProbaC)
ggqqplot(dataSD_achilles, x = "ProbaC") #normalite ok
stat.test<- dataSD_achilles %>% t_test(ProbaC~1, mu=0.5) #student compaison proba d'approche theodorus compar? a 0.5
mean(dataSD_achilles$ProbaC)
stat.test
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)

############# helenor ###########

ggplot(dataSD_helenor,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD_helenor %>% shapiro_test(ProbaC)
ggqqplot(dataSD_helenor, x = "ProbaC") #normalite pas ok
stat.test<- dataSD_helenor %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD_helenor$ProbaC)
dataSD_helenor %>% wilcox_test(ProbaC~1, mu=0.5)

########################## ANALYSE MOTIVATION #########################
############### Motivation approches ###########################
#pourcentage activit?
library(DHARMa)
library(lme4)
library(glmmTMB)
glmer1 <- glmer(Motivation_approches~Species+(1|Indiv)+(1|Essai), family=binomial, data=data)
summary(glmer1)

data$pred_approches <- glmer1$fitted.values

testDispersion(glmer1)
simulationOutput <- simulateResiduals(fittedModel = glmer1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

#nb d'approches
df1 <- filter(data, Sum_app>0)

sum(data_helenor$Sum_app > 0)
sum(data_achilles$Sum_app > 0)


mean_app<-summarySE(df1, measurevar="Sum_app", groupvars=c("Species"))


