#Mortierella Inhibition Study

#Created in MortTemp RProject
#WD="C:/Users/mrp354/OneDrive - The Pennsylvania State University/Desktop/Mortierella Projects/MortInhibit"

###############################################################################
#Load necessary packages
###############################################################################

library(dplyr)
library(tidyverse)
library(readr)
library(ScottKnottESD)
library(ScottKnott)

###############################################################################
#Load necessary files
###############################################################################

obs_data <- read.csv(file = "inputfiles/IAP.csv")

###############################################################################
#Convert to tibble, then sanity check
###############################################################################

tibbleIN <- as_tibble(obs_data)
head(tibbleIN)

###############################################################################
#Summarize Data (using dplyr), with sanity checks
###############################################################################

#Remember that some plates only had one organism so we'll see 0 values for some!
#Summarize pathogen growth by organism (use catorg column that links treatment to orgs)
pgrowthsum <- tibbleIN %>%                              
  group_by(catorg) %>% 
  summarize(min = min(pgrow),
            q1 = quantile(pgrow, 0.25),
            median = median(pgrow),
            mean = mean(pgrow),
            q3 = quantile(pgrow, 0.75),
            max = max(pgrow))
head(pgrowthsum)

#Summarize Mortierella's growth by isolate 
mgrowthsum <- tibbleIN %>%                              
  group_by(catorg) %>% 
  summarize(min = min(mgrow),
            q1 = quantile(mgrow, 0.25),
            median = median(mgrow),
            mean = mean(mgrow),
            q3 = quantile(mgrow, 0.75),
            max = max(mgrow))
head(mgrowthsum)

###############################################################################
#Reorganize data for easier calculations, particularly controls
###############################################################################

#Calculate average pathogen growth 
pavegrowth <- tibbleIN %>%
  group_by(treat, org) %>%
  mutate(Pave = mean(pgrow)) %>%
  distinct(Pave, .keep_all = TRUE) 
head(pavegrowth)

#Save pathogen control average as a value for later use
#Important: Verify this value is correct before proceeding!!
Pcon <- pavegrowth$Pave[21]

#Remove pathogen controls from pathogen growth table
FTgrowth <- pavegrowth[!(pavegrowth$treat=="C" | pavegrowth$org=="M"),]

#Get average Mortierella growth
mavegrowth <- tibbleIN %>%
  group_by(treat, org) %>%
  mutate(Mave = mean(mgrow)) %>%
  distinct(Mave, .keep_all = TRUE) 

#Separate out Mortierella control versus treatment growth 
#These need to be in their own columns for calculations!
Mcon <- mavegrowth[!(mavegrowth$treat=="F" | mavegrowth$org=="F"),]
Mtreat <- mavegrowth[!(mavegrowth$treat=="C" | mavegrowth$org=="F"),]

#Make a new tidy table with values in their own columns
m <- Mcon %>%
  select(-c(rep, round, catorg, pgrow, mgrow))
v <- cbind(m, Mtreat[8], FTgrowth[8]) 
Growthsummary <- v %>%
  group_by(org) %>%
  select(-c(treat))
Totalgrowthsummary <- Growthsummary %>% mutate_if(is.numeric, ~round(., 2))

#This renames the columns
#Fcontrol is a value, so not included in the table

colnames(Totalgrowthsummary)[2] ="MCon"
colnames(Totalgrowthsummary)[3] ="MTreat"
colnames(Totalgrowthsummary)[4] ="PTreat"
head(Totalgrowthsummary)

###############################################################################
#Calculate relative inhibition (RI)
###############################################################################

#Calculate Mortierella's growth ratio
Morat <- Totalgrowthsummary %>%
  group_by(org) %>%
  mutate(MRatio = MTreat/MCon)

#Calculate Fusarium's growth ratio
Furat <- Morat %>%
  group_by(org) %>%
  mutate(FRatio = PTreat/Pcon)

#Calculate relative inhibition 
RI <- Furat %>%
  group_by(org) %>%
  mutate(RI = FRatio/MRatio)
head(RI)
write.table(RI, file = "outputfiles/RI.csv", sep = ",", col.names = TRUE, row.names = TRUE)

###############################################################################
#Calculate percent radial growth inhibition (PRGI) of Fusarium by Mortierella
###############################################################################

#Calculate PGRI, then remove controls
PGRI_all <- tibbleIN %>%
  mutate(PGRI = ((Pcon-pgrow)/Pcon)*100)
PRGI <- PGRI_all[!(PGRI_all$treat=="C" | PGRI_all$treat=="P5"),]
head(PRGI)
write.table(PRGI, file = "outputfiles/PRGI.csv", sep = ",", col.names = TRUE, row.names = TRUE)

#Calculate average PGRI by Mortierella isolate
avePRGI <- PRGI %>%
  group_by(treat, org) %>%
  mutate(avePRGI = mean(PGRI)) %>%
  distinct(avePRGI, .keep_all = TRUE)
head(avePRGI)
write.table(avePRGI, file = "outputfiles/avePRGI.csv", sep = ",", col.names = TRUE, row.names = TRUE)

##Make a nice table
PRGItable <- avePRGI %>%
  select(-c(rep, round, catorg, pgrow, PGRI, mgrow))
head(PRGItable)

###############################################################################
#Tukey tests & Scott-Knott for growth differences by Mortierella species
###############################################################################

#Run Anova to prep for Tukey
aovPGRI <- aov(formula=PGRI ~ org, data=PRGI) 
capture.output(summary(aovPGRI), file = "outputfiles/aovPGRI.txt")
sumfit_PRGI <- summary.aov(aovPGRI)
print(aovPGRI)
print(sumfit_PRGI)

#Tukey to compare across treatments
tky_PRGI <- TukeyHSD(aovPGRI) 
capture.output(tky_PRGI, file = "outputfiles/tky_PRGI.txt")

#Identify groups based on Scott-Knott mean clustering
sk1 <- SK(PGRI ~org, 
          data = PRGI,
          which ="org", sig.level=.05,)
str(sk1)
names(sk1$clus)
print(sk1, digits = 2L, )
summary(sk1)

###############################################################################
#Plots of PRGI and Scott-Knott Clustering
###############################################################################

#This is a boxplot of PRGI, without any clustering
p<-ggplot(PRGI, aes(x=org, y=PGRI, fill=org))+ 
  labs(x= ~italic("Mortierella")~ "Isolate", y= ~italic("F. solani")~ "Percent Radial Growth Inhibition", fill="Isolate", title = ~italic("Fusarium")~ "Growth Inhibition") + 
  theme(legend.text = element_text(face = "italic")) +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",)) +
  geom_boxplot()
p

#This is the nicest Scott-Knott Boxplot I could make & I added the group letters externally
boxplot(sk1, args.legend=list(x='bottomleft'), adj=-0.5, main= ~italic("F. solani")~ "PRGI by" ~italic("Mortierella")~ "Isolate" , xlab= ~italic("Mortierella")~ "Isolate", ylab= ~italic("F. solani")~ "PRGI", margin=margin(30,0,0,0))


