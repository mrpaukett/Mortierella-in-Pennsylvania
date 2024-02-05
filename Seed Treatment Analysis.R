#Mortierella Seed Treatment Experiment

#Created in MortSeeds RProject
#WD="C:/Users/mrp354/OneDrive - The Pennsylvania State University/Desktop/Mortierella Projects/MortSeeds"

###############################################################################
#Load necessary packages
###############################################################################
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)
library(round)
library(ggpubr)   
library(gridExtra) #This is to print tables as PDFs)

###############################################################################
#Load necessary files
###############################################################################

obs_data <- read.csv(file = "inputfiles/seed_assay_finalmaster.csv")

###############################################################################
#Convert to tibble, then sanity check
###############################################################################

tibbleseed <- as_tibble(obs_data)
head(tibbleseed)

###############################################################################
#Summarize Data (using dplyr), with sanity checks
###############################################################################

#Summarize root length by treatment 
rootsum <- tibbleseed %>%                              
  group_by(treatment) %>% 
  summarize(min = min(root_length),
            q1 = quantile(root_length, 0.25),
            median = median(root_length),
            mean = mean(root_length),
            q3 = quantile(root_length, 0.75),
            max = max(root_length))
head(rootsum)

#Summarize root length by treatment 
DRsum <- tibbleseed %>%                              
  group_by(treatment) %>% 
  summarize(min = min(disease_rating),
            q1 = quantile(disease_rating, 0.25),
            median = median(disease_rating),
            mean = mean(disease_rating),
            q3 = quantile(disease_rating, 0.75),
            max = max(disease_rating))
head(DRsum)

###############################################################################
#Correlations of Root & Disease Rating
###############################################################################

#Save the two values of interest
k <- tibbleseed$disease_rating
o <- tibbleseed$root_length

#Run correlation tests with several common methods
corr <- cor.test(k, o, method=c("pearson", "kendall", "spearman"))
capture.output(corr, file = "outputfiles/corr_seeds.txt")
print(corr) 

#Visualize Spearman
s <- ggscatter(tibbleseed, x = "disease_rating", y = "root_length", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Disease Rating", ylab = "Root Length (mm)")
s

#Visualize Kendall
k <- ggscatter(tibbleseed, x = "disease_rating", y = "root_length", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "kendall",
               xlab = "Disease Rating", ylab = "Root Length (mm)")
k

#Visualize Pearson
r <- ggscatter(tibbleseed, x = "disease_rating", y = "root_length", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "kendall",
               xlab = "Disease Rating", ylab = "Root Length (mm)")
r

###############################################################################
#Count Disease Data, with sanity checks
###############################################################################

#This counts how many of each dummy variable for each treatment
countdum <- tibbleseed %>%
  group_by(treatment) %>%
  count(Dummy) 
head(countdum)

#This counts how many seeds in each treatment have each disease rating (1-4)
countDR <- tibbleseed %>%
  group_by(treatment) %>%
  count(disease_rating) 
head(countDR)

###############################################################################
#Analyze Dummy Variable Data (disease incidence), with sanity checks
###############################################################################

####
#Create a new table from counted dummy variables, with Dummy 0 & 1 in their own columns
####

#This creates the new coloumns, then binds them into a new table
x <- countdum[!(countdum$Dummy=="1" | countdum$Dummy=="1"),]
y <- countdum[!(countdum$Dummy=="0" | countdum$Dummy=="0"),]
dummytable <- cbind(x, y[3]) 

#This changes the column names to make them meaningful
colnames(dummytable)[3] ="D0"
colnames(dummytable)[4] ="D1"
head(dummytable)

####
#Calculate disease incidence (dummy ratio), based on dummy variables!
####

#Ratio calculation
DIratio <- dummytable %>%
  group_by(treatment) %>%
  mutate(DI = ((D0/D1)*100))
head(DIratio)

#Save control (p treatment) ratio as a value, for treatment change calculations
#Important: Always verify the correct value was saved before proceeding!!!
DIcont <- DIratio$DI[6]

#Calculate DI treatment difference vs control
DIdiff <- DIratio %>%  
  group_by(treatment) %>%
  mutate(DIchange = (DIcont-DI))
head(DIdiff)

###############################################################################
#Analyze Disease Rating Data (DSI), with sanity checks
###############################################################################

####
#Calculate DSI for all 3 rounds (54 total seeds, top rating was 4) 
####
dsi <- countDR  %>%
  mutate(dsiorg = (sum(disease_rating*n)/(54*4))*100)
head(dsi)

####
#Since DSI is same for each treatment, you only need one row in the table for it
#All treatments had disease ratings of 4 so keep that to create a new table
####
j <- dsi[!(dsi$disease_rating=="3" | dsi$disease_rating=="2" | dsi$disease_rating=="1" ),]

#Save control (p treatment) DSI as a value, for treatment change calculations
#Important: Always verify the correct value was saved before proceeding!!!
DSIcont <- j$dsiorg[6]

#Calculate DSI treatment difference vs control
DSIchange <- j %>%  
  group_by(treatment) %>%
  mutate(jdiff = (DSIcont-dsiorg))
head(DSIchange)

###############################################################################
#Tukey tests for differences in root lengths
###############################################################################

#Run Anova on root lengths by treatment
aovseed <- aov(formula=root_length ~ treatment, data=tibbleseed) 
capture.output(summary(aovseed), file = "outputfiles/aovseed.txt")
sumfit_seed <- summary.aov(aovseed)
print(aovseed)
print(sumfit_seed)

#Run Tukey Test to compare across treatments
tky_seed <- TukeyHSD(aovseed) 
capture.output(tky_seed, file = "outputfiles/tky_seed.txt")

#Store root length averages for use in plot
rootmeans <- aggregate(root_length ~  treatment, tibbleseed, mean)

#Make boxplot of root growth by treatment!
p<-ggplot(tibbleseed, aes(x=treatment, y=root_length, fill=treatment))+ 
  labs(x= "Treatment", y= "Root Length (mm)", fill="Treatment", title = "Root Growth") + 
  theme(legend.text = element_text(face = "italic")) +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",)) +
  geom_boxplot() + 
  geom_text(data = rootmeans, aes(label = root_length, y = root_length + 0.08))
p

###############################################################################
#Make a nice data table of DI & DSI Calculations
###############################################################################

#Get rid of unnecessary columns from DI table
cuts <- DIdiff %>% 
  subset(select = -c(Dummy,D1,D0))

#Rename columns
colnames(cuts)[1] ="Treatment" 
colnames(cuts)[2] ="DI" 
colnames(cuts)[3] ="DIchange"

#Add DSI table columns & rename
addj <- cuts %>% 
  cbind(DSIchange) %>% 
  subset(select = -c(treatment,disease_rating,n)) 
colnames(addj)[4] ="DSI" 
colnames(addj)[5] ="DSIchange"

#Round number class columns to get final table
diseasetable <- addj %>% mutate_if(is.numeric, ~round(., 1))

#Export a pdf of the nicely formatted table
pdf("outputfiles/Seed_Disease_Table.pdf")       
grid.table(diseasetable)
dev.off()

