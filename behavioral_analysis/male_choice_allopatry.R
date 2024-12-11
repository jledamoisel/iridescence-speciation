

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(reshape2)

#### FUNCTION ######

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
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

###### GET DATA ######

data<-read.table("ecuador_malechoice.csv", dec=".", sep=";", header=T)

#data globales

data<- data %>% 
  mutate (ID_dummy1=case_when(Paires==1~'theodorus',
                              Paires==2~'theodorus',
                              Paires==3~'bristowi_bande'))
data<- data %>% 
  mutate (ID_dummy2=case_when(Paires==1~'bristowi',
                              Paires==2~'bristowi_bande',
                              Paires==3~'bristowi'))
data<- data %>% 
  mutate (Motivation_approches=case_when(Zone_T=='0' & Zone_B=='0' ~ 'no',
                                         TRUE~'yes'))
data<-data %>% 
  mutate (Motivation_touche=case_when(Touche_T=='0' & Touche_B=='0' ~ 'no',
                                      TRUE~'yes'))
data<-data %>% 
  mutate (Sum_touche=Touche_T+Touche_B)
data<-data %>% 
  mutate(Sum_app=Zone_B+Zone_T)


data$Indiv<-as.factor(data$Indiv)
data$Paires<-as.character(data$Paires)
data$Motivation_approches<-as.factor(data$Motivation_approches)
data$Motivation_touche<-as.factor(data$Motivation_touche)
data$Sum_app<-as.numeric(data$Sum_app)
data$Sum_touche<-as.numeric(data$Sum_touche)

#data proba
data_proba <- data %>% 
  mutate (ProbaC=(App_conspe)/(App_conspe+App_non_conspe))
data_probaNA <- data_proba %>% drop_na(ProbaC)
data_proba$ProbaC[is.nan(data_proba$ProbaC)]<-0

data_proba_touche <- data %>% 
  mutate (ProbaT=(Touch_conspe)/(Touch_conspe+Touch_non_conspe))
data_proba_touche$ProbaT[is.nan(data_proba_touche$ProbaT)]<-0

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

mean_app<-summarySE(df1, measurevar="Sum_app", groupvars=c("Species"))

################ Motivation touches ##########################
library(DHARMa)
library(lme4)
library(glmmTMB)
glmer1 <- glmer(Motivation_touche~Species+(1|Indiv)+(1|Essai), family=binomial, data=data)
summary(glmer1)

testDispersion(glmer1)
simulationOutput <- simulateResiduals(fittedModel = glmer1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

#nb de touches
df1 <- filter(data, Sum_touche>0)

mean_app<-summarySE(df1, measurevar="Sum_touche", groupvars=c("Species"))


########################### ANALYSE APPROCHES ########################
######## Analyse des data dans leur totalite #########################


#creation d'une variable combinee
data$y<-cbind(data$App_conspe,data$App_non_conspe)

#glm avec effet random sur les proportions
library(DHARMa)
library(lme4)
library(glmmTMB)
glmer1 <- glmer(y~Species+(1|Indiv), family=binomial, data=data)
summary(glmer1)

testDispersion(glmer1)
simulationOutput <- simulateResiduals(fittedModel = glmer1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

data_theo<-data  %>% 
  filter(Species == "theodorus")
  
glmer1 <- glmer(y~Paires+(1|Indiv), family=binomial, data=data_theo)
summary(glmer1)

data_bristo<-data  %>% 
  filter(Species == "bristowi")

glmer2 <- glmer(y~Paires+(1|Indiv), family=binomial, data=data_bristo)
summary(glmer2)

sum(data_bristo$Sum_app > 0)
sum(data_theo$Sum_app > 0)


############# Analyse des r?sultats dans la paire 1 avec theodorus ###########
data_paire1_theo <- data_paire1 %>% 
  filter(Species == "theodorus")%>% 
  mutate (ProbaC=(Zone_T)/(Zone_T+Zone_B))
data_paire1_theo <- data_paire1_theo %>% drop_na(ProbaC)


dataSD<-summarySE(data_paire1_theo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_app_theo1<- dataSD %>%
  mutate(Species="theodorus")%>%
  mutate(Pair=1)

ggplot(dataSD_app_theo1,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5) #student compaison proba d'approche theodorus compar? a 0.5
mean(dataSD$ProbaC)
stat.test
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)

############# Analyse des r?sultats dans la paire 1 avec bristowi ###########
data_paire1_bristo <- data_paire1 %>% 
  filter(Species == "bristowi")%>% 
  mutate (ProbaC=(Zone_B)/(Zone_T+Zone_B))
data_paire1_bristo <- data_paire1_bristo %>% drop_na(ProbaC)

dataSD<-summarySE(data_paire1_bristo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_app_bristo1<- dataSD %>%
  mutate(Species="bristowi")%>%
  mutate(Pair=1)

ggplot(dataSD_app_bristo1,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite pas ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD$ProbaC)
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)

############# Analyse des r?sultats dans la paire 2 ###########
data_paire2 <- data %>% 
  filter(Paires == 2)
data_paire2 <- data_paire2 %>% 
  mutate (ProbaT=(Zone_T)/(Zone_T+Zone_B))%>% 
  mutate (ProbaC=(Zone_T)/(Zone_T+Zone_B))
data_paire2_theo <- data_paire2 %>% drop_na(ProbaT)

dataSD<-summarySE(data_paire2_theo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_app_theo2<- dataSD %>%
  mutate(Species="theodorus")%>%
  mutate(Pair=2)

ggplot(dataSD_app_theo2,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD$ProbaC)
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)


############# Analyse des r?sultats dans la paire 3 ###########
data_paire3 <- data %>% 
  filter(Paires == 3)
data_paire3 <- data_paire3 %>% 
  mutate (ProbaT=(Zone_T)/(Zone_T+Zone_B))%>% 
  mutate (ProbaC=(Zone_B)/(Zone_T+Zone_B))
data_paire3_bristo <- data_paire3 %>% drop_na(ProbaC)

dataSD<-summarySE(data_paire3_bristo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_app_bristo3<- dataSD %>%
  mutate(Species="bristowi")%>%
  mutate(Pair=3)

ggplot(dataSD_app_bristo3,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD$ProbaC)
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)


############# Graph global approche ###################################
approches_table<- rbind(dataSD_app_bristo1,dataSD_app_bristo3)
approches_table<- rbind(approches_table,dataSD_app_theo1)
approches_table<- rbind(approches_table,dataSD_app_theo2) 

ggplot(approches_table,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) + facet_wrap(~Pair)+ theme_classic() 

#approche de la paire 1 seulement
approches_table1<-rbind(dataSD_app_bristo1,dataSD_app_theo1)

ggplot(approches_table1,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.5)+
  scale_color_manual(values=c("#CC99FF", "#FF9966"))+
  scale_fill_manual(values=c("#CC99FF", "#FF9966"))+theme_classic() +
  coord_flip()

#approche de la paire 2 seulement
approches_table1<-rbind(dataSD_app_bristo1,dataSD_app_theo2)

ggplot(approches_table1,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.5)+
  scale_color_manual(values=c("#CC99FF", "#FF9966"))+
  scale_fill_manual(values=c("#CC99FF", "#FF9966"))+theme_classic() +
  coord_flip()

#approche de la paire 3 seulement
approches_table1<-rbind(dataSD_app_bristo3)

ggplot(approches_table1,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.25)+
  scale_color_manual(values=c("#CC99FF", "#FF9966"))+
  scale_fill_manual(values=c("#CC99FF", "#FF9966"))+theme_classic() +
  coord_flip()

############################ ANALYSE TOUCHES ########################

############## Analyse des resultats de la paire 1 seulement #########################

############# Analyse des r?sultats dans la paire 1 avec theodorus ###########
data_paire1_theo <- data %>% 
  filter(Species == "theodorus") %>% 
  filter(Paires == 1) %>% 
  mutate (ProbaT=(Touche_T)/(Touche_T+Touche_B)) %>% 
  mutate (ProbaC=(Touche_T)/(Touche_T+Touche_B))
data_paire1_theo <- data_paire1_theo %>% drop_na(ProbaT)

dataSD<-summarySE(data_paire1_theo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_tou_theo1<- dataSD %>%
  mutate(Species="theodorus")%>%
  mutate(Pair=1)

ggplot(dataSD_tou_theo1,aes(x=Species, y=ProbaC))+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaCSpeciesdataSD %>% shapiro_test(ProbaC))
ggqqplot(dataSD, x = "ProbaC") #normalite ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5) #student compaison proba d'approche theodorus compar? a 0
mean(dataSD$ProbaT)
stat.test

dataSD %>% wilcox_test(ProbaC~1, mu=0.5)

############# Analyse des r?sultats dans la paire 1 avec bristowi ###########
data_paire1_theo <- data %>% 
  filter(Species == "bristowi") %>% 
  filter(Paires == 1) %>% 
  mutate (ProbaT=(Touche_T)/(Touche_T+Touche_B)) %>% 
  mutate (ProbaC=(Touche_B)/(Touche_T+Touche_B))
data_paire1_theo <- data_paire1_theo %>% drop_na(ProbaT)

dataSD<-summarySE(data_paire1_theo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_tou_bristo1<- dataSD %>%
  mutate(Species="bristowi")%>%
  mutate(Pair=1)

ggplot(dataSD_tou_bristo1,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite pas ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD$ProbaC)
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)

############# Analyse des r?sultats dans la paire 2 ###########
data_paire1_theo <- data %>% 
  filter(Paires == 2) %>% 
  mutate (ProbaT=(Touche_T)/(Touche_T+Touche_B))%>% 
  mutate (ProbaC=(Touche_T)/(Touche_T+Touche_B))
data_paire1_theo <- data_paire1_theo %>% drop_na(ProbaT)

dataSD<-summarySE(data_paire1_theo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_tou_theo2<- dataSD %>%
  mutate(Species="theodorus")%>%
  mutate(Pair=2)

ggplot(dataSD_tou_theo2,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD$ProbaC)
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)


############# Analyse des r?sultats dans la paire 3 ###########
data_paire1_theo <- data %>% 
  filter(Paires == 3) %>% 
  mutate (ProbaT=(Touche_T)/(Touche_T+Touche_B))%>% 
  mutate (ProbaC=(Touche_B)/(Touche_T+Touche_B))
data_paire1_theo <- data_paire1_theo %>% drop_na(ProbaT)

dataSD<-summarySE(data_paire1_theo, measurevar="ProbaC", groupvars=c("Indiv"))
dataSD_tou_bristo3<- dataSD %>%
  mutate(Species="bristowi")%>%
  mutate(Pair=3)

ggplot(dataSD_tou_bristo3,aes(x=Species, y=ProbaC) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5)+ theme_classic() 

dataSD %>% shapiro_test(ProbaC)
ggqqplot(dataSD, x = "ProbaC") #normalite ok
stat.test<- dataSD %>% t_test(ProbaC~1, mu=0.5)
stat.test
mean(dataSD$ProbaC)
dataSD %>% wilcox_test(ProbaC~1, mu=0.5)

############# Graph global touche ###################################
approches_table<- rbind(dataSD_tou_bristo1,dataSD_tou_bristo3)
approches_table<- rbind(approches_table,dataSD_tou_theo1)
approches_table<- rbind(approches_table,dataSD_tou_theo2) 

ggplot(approches_table,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) + facet_wrap(~Pair)+ theme_classic() 

#approche de la paire 1 seulement
approches_table1<-rbind(dataSD_tou_bristo1,dataSD_tou_theo1)

ggplot(approches_table1,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.5)+
  scale_color_manual(values=c("#CC99FF", "#FF9966"))+
  scale_fill_manual(values=c("#CC99FF", "#FF9966"))+theme_classic() +
  coord_flip()

#approche de la paire 2 seulement
approches_table1<-rbind(dataSD_tou_bristo3,dataSD_tou_theo2)

ggplot(approches_table1,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.5)+
  scale_color_manual(values=c("#CC99FF", "#FF9966"))+
  scale_fill_manual(values=c("#CC99FF", "#FF9966"))+theme_classic() +
  coord_flip()

#approche de la paire 3 seulement
approches_table1<-rbind(dataSD_app_bristo3)

ggplot(approches_table1,aes(x=Species,y=ProbaC, fill=Species, colour=Species) )+
  geom_jitter(width=0.25) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5, width=0.25)+
  scale_color_manual(values=c("#CC99FF", "#FF9966"))+
  scale_fill_manual(values=c("#CC99FF", "#FF9966"))+theme_classic() +
  coord_flip()