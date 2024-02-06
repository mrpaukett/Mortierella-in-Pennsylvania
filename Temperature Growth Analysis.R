#Mortierella Temperature Study

#Created in MortTemp RProject
#WD="C:/Users/mrp354/OneDrive - The Pennsylvania State University/Desktop/Mortierella Projects/MortTemp"

###############################################################################
#Load necessary packages
###############################################################################
library(dplyr)
library(tidyverse)
library(readr)

###############################################################################
#Load necessary files
###############################################################################

obs_data <- read.csv(file = "inputfiles/temp_study_master.csv")

###############################################################################
#Convert to tibble, then sanity check
###############################################################################

tibbletemp <- as_tibble(obs_data)
head(tibbletemp)

###############################################################################
#Summarize Data (using dplyr), with sanity checks
###############################################################################

#Summarize growth by organism & temp (use catorg column that links both)
growthsum <- tibbletemp %>%                              
  group_by(catorg) %>% 
  summarize(min = min(measurement),
            q1 = quantile(measurement, 0.25),
            median = median(measurement),
            mean = mean(measurement),
            q3 = quantile(measurement, 0.75),
            max = max(measurement))
head(growthsum)

###############################################################################
#Calculate Percent Growth Reduction (PGR) by Temperature for each Organism
###############################################################################

#Calculate average growth by orgs & temp, then remove some unnecessary columns
aveTemps <- tibbletemp %>%
  group_by(temp, organism) %>%
  mutate(aveGrowth = mean(measurement)) %>%
  distinct(aveGrowth, .keep_all = TRUE) 
head(aveTemps)

#Remove control temp (25) from the table; this is the baseline for comparison
#These need to be in their own column for later calculations
aveTempTreatments <- aveTemps[!(aveTemps$treat=="C" | aveTemps$temp=="25"),]
head(aveTempTreatments)

#Put control (25) values in their own table 
#Rename column so won't cause a duplicate when bind to other table
aveTempControl <- aveTemps %>%
  subset(temp > 15) %>%
  rename(aveControlsGrowth = aveGrowth)
head(aveTempControl)

#There are two temps in the the other table for each organism
#These control values need to be duplicated so column sizes match for binding/calculations!!
replicatedavecontrolsum <- aveTempControl[rep(seq_len(nrow(aveTempControl)), each = 2), ]  
head(replicatedavecontrolsum)

#Bind the controls to the treatment table (with omitted controls)
#add the average control values by org column to the data set where the controls are omitted##
masterTemp <- cbind(aveTempTreatments, replicatedavecontrolsum[8])
head(masterTemp)

#Finally, we can calculate the growth reduction by temperature!
#Note: This is more to summarize in a small table the growth reduction. 
#This is NOT used for anova or graphs; you need a full data set for that (next section)
TempPGIsummary <- masterTemp %>%
  mutate(PGR = ((aveControlsGrowth-aveGrowth)/(aveControlsGrowth)*100)) %>%
  select(-c(treat, aveControlsGrowth, aveGrowth))
head(TempPGIsummary)

###############################################################################
#Calculate PGR by Temp for each Organism on full data set for further analysis
###############################################################################

#For the full data, set the control (25) needs to be present for analysis
#We also don't have averaged measurements per org (so 3 rows/temp/org)

#We can use the controls separated previously, BUT...
#We need 9 rows of each control value for the correct binding column size!!
replicatedavecontrols <- aveTempControl[rep(seq_len(nrow(aveTempControl)), each = 9), ]  
head(replicatedavecontrols)

#Bind the controls to the original data file (in tibble form)
masterTempfile <- cbind(tibbletemp, replicatedavecontrols[8])
head(masterTempfile)

#Calculate PGR, making it negative for the plot later
#This is the data that will be used in further analysis (Anova/Tukey) & graphs
TempPGIanalysis <- masterTempfile %>%
  mutate(PGR = (((aveControlsGrowth-measurement)/(aveControlsGrowth)*100)*-1)) %>%
  select(-c(treat, aveControlsGrowth))
head(TempPGIanalysis)

###############################################################################
#Calculate Growth Rate Per Degree Change (mm/1C increase)
###############################################################################

#Calculate the growth rate (10C change between each temp)
ratedeg <- masterTempfile %>%
  mutate(rate = ((measurement)/10)) %>%
  select(-c(treat, aveControlsGrowth))
head(ratedeg)

#Calcuate the average rate by org and temperature
averatedeg <- ratedeg  %>%
  group_by(temp, organism) %>%
  mutate(averate = mean(rate))
head(averatedeg)

#Remove replicates & some columns to get a nice small table  
addj <- averatedeg[!(averatedeg$replicate=="R2" | averatedeg$replicate=="R3"),]
m <- addj %>%
  select(-c(replicate, round, rate, measurement, catorg))
SumDegreeGrowth <- m %>% mutate_if(is.numeric, ~round(., 3))
head(SumDegreeGrowth)

#We need the different rates in thier own columns for calculations (done by row for an organism)
#So we need to separate and merge items...

#This creates the new columns, then binds them into a new table
x <- SumDegreeGrowth[!(SumDegreeGrowth$temp=="15" | SumDegreeGrowth$temp=="25"),]
y <- SumDegreeGrowth[!(SumDegreeGrowth$temp=="5" | SumDegreeGrowth$temp=="25"),]
z <- SumDegreeGrowth[!(SumDegreeGrowth$temp=="5" | SumDegreeGrowth$temp=="15"),]
v <- cbind(x, y[3], z[3]) 

#This changes the column names to make them meaningful & removes unnecessary columns
SummaryRates <- v %>%
  group_by(organism) %>%
  select(-c(temp))

colnames(SummaryRates)[2] ="five"
colnames(SummaryRates)[3] ="fifteen"
colnames(SummaryRates)[4] ="twentyfive"

#This drops duplicate rows and leaves us with the nice table (finally!)
SummaryRates <- SummaryRates[!duplicated(SummaryRates), ]
head(SummaryRates)

#Now to calculate growth change between temperatures
Ratetable <- SummaryRates %>%
  group_by(organism) %>%
  mutate(GRchange25 = twentyfive-fifteen) %>%
  mutate(GRchange5 = fifteen-five)
head(Ratetable)

###############################################################################
#Tukey tests for differences in growth 
###############################################################################

#Run Anova on PGR by treatment 
aovtemp <- aov(formula=PGR ~ catorg, data=TempPGIanalysis) 
capture.output(summary(aovtemp), file = "outputfiles/aovTemp.txt")
sumfit_temp <- summary.aov(aovtemp)
print(aovtemp)
print(sumfit_temp)

#Tukey to compare accross treatments
tky_temp <- TukeyHSD(aovtemp) 
capture.output(tky_temp, file = "outputfiles/tky_pgrtemp.txt")

#Run Anova on measurement by treatment 
aovtempmea <- aov(formula=measurement ~ catorg, data=TempPGIanalysis) 
capture.output(summary(aovtempmea), file = "outputfiles/aovTempmea.txt")
sumfit_tempmea <- summary.aov(aovtempmea)
print(aovtempmea)
print(sumfit_tempmea)

#Tukey to compare across treatments
tky_tempmea <- TukeyHSD(aovtempmea) 
capture.output(tky_tempmea, file = "outputfiles/tky_pgrtempmea.txt")

###############################################################################
#Plot Growth
###############################################################################

#Boxplot of growth reduction (standardized to 25C control)
p<-ggplot(TempPGIanalysis, aes(x=factor(catorg,level=c('25A' , '25G', '25X', '25E', '25F', '15A' , '15G', '15X', '15E', '15F', '5A' , '5G', '5X', '5E', '5F')), y=PGR, fill=organism))+ 
  labs(x= ~italic("Mortierella & Fusarium")~ "Isolate By Temperature", y= "Percent Diameter Growth Inhibition", fill="Isolate", title = "Temperature Growth Inhibition") + 
  theme(legend.text = element_text(face = "italic")) +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",)) +
  geom_boxplot()
p

#Line graph of average growth at each temperature
q <- ggplot(data=aveTemps, aes(x=temp, y=aveGrowth, group=organism)) +
  geom_line(aes(linetype=organism))+
  geom_point(aes(shape=organism))+
  labs(x= "Temperature", y= "Average Isolate Growth (mm)", title = ~italic("Mortierella & Fusarium")~ "Growth By Temperature") + 
  theme(legend.text = element_text(face = "italic")) +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust=0.5, face = "bold",))
q

