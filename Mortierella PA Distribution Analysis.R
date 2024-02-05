#Mortierella Distribution Analysis

#Created in MortDist RProject
#WD="C:/Users/mrpau/Desktop/Mortierella Projects/MortDist"

#Note: Region grouping used in analysis based on spatial coordinates
#Note: This code primarily follows the online tutorial: https://rpubs.com/an-bui/vegan-cheat-sheet

###############################################################################
#Load necessary packages
###############################################################################
library(here)      #File referencing/path-making for projects
library(tidyverse) #Data manipulation, includes dyplr & ggplot2 which need here
library(vegan)     #Community ecology, primary package for distribution analysis
library(ggvegan)   #ggplot-based plots for vegan; installed with devtools/from developer

###############################################################################
#Load necessary files
#Note: Info on paenv data is in a text file at: 'C:/Users/mrpau/Desktop/Mortierella Projects/MortDist/data'
###############################################################################

#Sites are ordered 1-30 in all files (even if no actual site column listed)
mspecies <- read.csv("inputfiles/Mortpresence.csv")  #Species abundance per site
paenv <- read.csv("inputfiles/paenviroup.csv")         #Site environmental data, includes site coordinates


###############################################################################
#Note: This code only focuses on region ('group') analysis
###############################################################################

site_type <- paenv %>% 
  select(site, group, perf, texture, county, plant)  #Make new table

###############################################################################
#Species Richness, by Region ("How specious are my communities?")
#No anova significance is fine!!
###############################################################################

sppr <- specnumber(mspecies)                                #Calculate and save number of species within each site
sppr_aov <- aov(sppr ~ group, data = site_type)             #Run anova ("Is mean species richness significantly different within region?")
summary(sppr_aov)                                           #Show results in console
capture.output(sppr_aov, file = "outputfiles/sppr_aov.txt") #Save results in project folder

###############################################################################
#Plot Species Richness, by Region 
###############################################################################

sppr_df <- sppr %>% 
  enframe()                                   #Make a data frame of species richness
sppr_df_env <- cbind(sppr_df, site_type[1:6]) #Add text environmental variables from metadata to the species data frame

plot_spprg <- ggplot(sppr_df_env, aes(x = group, y = value, fill = group)) +
  geom_boxplot()  +
  labs(x = "Region",
       y = "Number of Species per Site",
       title = "Species Richness")          #Make a boxplot of species richness
plot_spprg  #Show boxplot!

###############################################################################
#Shannon Diversity, by Region ("How diverse are my communities?")
#Shannon accounts for species abundance and evenness
###############################################################################

shannondiv <- diversity(mspecies)       #Calculate Shannon diversity
head(shannondiv)                        #Show results in console
sppdiv_aov <- aov(shannondiv ~ group, data = site_type)         #Run anova on Shannon Diversity
summary(sppdiv_aov)                                             #Show results in console
capture.output(sppdiv_aov, file = "outputfiles/sppdiv_aov.txt") #Save results in project folder

###############################################################################
#Plot Shannon Diversity, by Region 
###############################################################################

shandiv_df <- shannondiv %>% 
  enframe()                                         #Make a data frame of diversity
shandiv_df_env <- cbind(shandiv_df, site_type[1:4]) #Add text environmental variables from metadata to the species data frame

Gdiv_plot_df <- shandiv_df_env %>% 
  group_by(group) %>% 
  summarize(mean = round(mean(value), 2),
            err = sd(value)/sqrt(length(value))) %>% 
  dplyr::mutate(label = "mean") %>% 
  unite("mean_label", label, mean, sep = " = ", remove = FALSE)  #Create data frame with SE & Mean for plot labels

Gplot_shandiv <- ggplot(Gdiv_plot_df, aes(x = group, y = mean, fill = group)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err), width = 0.5) +
  geom_text(aes(x = group, y = mean + err + 0.07, label = mean_label)) +
  labs(x = "Region",
       y = "Mean Shannon Diversity",
       title = "Shannon Diversity")       #Make a boxplot of Shannon Diversity
Gplot_shandiv                             #Plot boxplot!

##############################################################################
#Run perMANOVA By Site & Region ("How different are my communities in species composition?")
#Assess differences based on dissimilarity
##############################################################################

gm_perm <- adonis(mspecies ~ group, data = paenv)         #Run perMANOVA by Region
gm_perm                                                   #Show results in console
capture.output(gm_perm, file = "outputfiles/gm_perm.txt") #Save results in project folder

###############################################################################
#PCA (Principal Components Analysis)
#Basic representation of communities (multivariate species dataset) in ordination space
#which distills the multivariate (i.e. species) dataset into two axes.
###############################################################################
mPCA <- rda(mspecies)    #Run PCA
mPCA                     #Show results in console
capture.output(mPCA, file = "outputfiles/mPCA.txt")  #Save results in project folder; use PC values from this for plot labels

PCA_biplot <- autoplot(mPCA) #Easy autoplot of PCA with ggvegan!
PCA_biplot                   #Show biplot! This can be customized with ggplot2...

PCA_fortify <- fortify(mPCA) #Makes a data frame from which elements can be extracted

PCA_fort_sites <- PCA_fortify %>% 
  filter(Score == "sites")           #Extract site coordinates (points)

PCA_fort_sites_env <- cbind(PCA_fort_sites, site_type[1:4]) #Add text environmental variables from metadata to the site coordinates data frame

PCA_fort_species <- PCA_fortify %>% 
  filter(Score == "species")         #Extract species coordinates (arrows)

GPCA_fortify_plot <- ggplot() +
  geom_point(data = PCA_fort_sites_env, aes(x = PC1, y = PC2, col = group)) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = PCA_fort_species, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = PCA_fort_species, aes(x = PC1, y = PC2, label = Label))  +
  labs(x = "PC1 (3.87%)",
       y = "PC2 (1.63%)",
       title = "Principal Components Analysis: Region") #Plot ordination by region!
##Note: PC %'s for axes were manually added so update if make changes!

GPCA_fortify_plot                                       #Show ordination plot!  

###############################################################################
#Non-metric Multidimensional Scaling (NMDS)
#Commonly used by community ecologists!
#Reduces axes for each 
#You can imagine NMDS as a reduction of axes, where all your “axes” are the species within a sample, and each sample exists relative to others on the axes. 
#Imagine an axis that describes relative abundance of ACFL - all points exist somewhere on that axis relative to the others.
#Now consider another axis describing the relative abundance of KEWA - all points still exist somewhere on that axis, but now there are two axes describing the position of your points.
#NMDS allows you to collapse all these species axes (in this case, 48) into 2 to plot in cartesian space in order to visualize the differences between samples and sites.


###############################################################################

mort_NMDS <- metaMDS(mspecies)     #Run NMDS
mort_NMDS                          #Show results in console
#This shows the stress is less than .20 which is the acceptable limit

capture.output(mort_NMDS, file = "outputfiles/mort_NMDS.txt")  #Save results in project folder

stressplot(mort_NMDS) #Plot NMDS stress 
plot(mort_NMDS)
#Looking at the stressplot for your  is an important part of evaluating how well the ordination represented the complexity in your data. The x-axis is the observed dissimilarity, and the y-axis is the ordination distance. The stressplot shows you how closely the ordination (y-axis) represents the dissimilarities calculated (x-axis). The points around the red stair steps are the communities, and the distance from the line represents the “stress”, or how they are pulled from their original position to be represented in their ordination.
#You can plot your NMDS output in base R…
#… but you can also extract elements from the output to plot in ggplot, which is much nicer.
#Things to consider about stress
#Generally, an average stress of > 0.2 is an indication that the ordination didn’t do a good job of representing community structure. However, NMDS stress increases as sample size increases, so if you have many samples (as this dataset does) your ordination will likely exhibit a lot of stress. As a quick example, I’ve subsampled 15 communities from the original set of 210. If you look at the stress from the NMDS and the stress plot, it’s within the range that’s considered “acceptable”

plot_df <- scores(mort_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") 
plot_df
k <- cbind(plot_df, site_type[1:6]) 
#rename so not two site columns
colnames(k)[4] ="Field"

#plot by region
Gplot_nmds <- ggplot(k, aes(x = NMDS1, y = NMDS2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(linewidth = 2, size = 1) +
  labs(title = "Region NMDS")
Gplot_nmds

###############################################################################
#Enfit!
#Determine the relative contribution of environmental variables to the separation of your communities along ordination axes
###############################################################################

# envfit() takes the output of metaMDS() and the species matrix you created
fit <- envfit(mort_NMDS, mspecies, perm = 999) 

fit# extract p-values for each species
fit_pvals <- fit$vectors$pvals %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::rename("pvals" = ".")

# extract coordinates for species, only keep species with p-val < 0.05***
fit_spp <- fit %>% 
  scores(., display = "vectors") %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  full_join(., fit_pvals, by = "species") %>% 
  filter(pvals < .05)

# new plot by region
Gnmds_plot_new <- ggplot(k, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() +
  geom_point(aes(color = group, shape = group), size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = group)) +
  geom_segment(data = fit_spp, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               col = "black") +
  geom_text(data = fit_spp, aes(label = species)) +
  labs(title = "Region Enfit")
Gnmds_plot_new

###############################################################################
#Constrained ordination of environmental variables
###############################################################################

#Include top variables (based on PCA) with numeric values
rdatop <- rda(mspecies~elevation+mintemp+lat+long+precip+watervol, data=paenv)
rdatop
autoplot(rdatop, arrows = TRUE)

dm_perm <- adonis(mspecies ~ elevation+temp+precip+lat+long+snow+maxtemp+mintemp+watervol+rain+temprange, data = paenv)         #Run perMANOVA by everything to see what's significant
dm_perm  

###############################################################################
#Constrained analysis!!
#I think this is what is preferred by ecologists, it connects species more to environmental factors
###############################################################################
#But this doesn't plot points! It just states rows instead.

mortCCA <- cca(mspecies ~ elevation+mintemp+lat+long+precip+watervol, data = paenv)
mortCCA
capture.output(mortCCA, file = "outputfiles/mortCCA.txt") #Save results in project folder
plot(mortCCA)

# vectors
ccavectors <- as.matrix(scores(mortCCA, display = "bp", scaling = "species")*7.627807) %>% 
  as.data.frame()

# site coordinates
site_data <- scores(mortCCA, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") 

#This is where we add site_type data
L <- cbind(site_data, site_type[1:6]) 
#rename so not two site columns
colnames(L)[4] ="Field"

# species coordinates
species_data <- scores(mortCCA, display = "species") %>% 
  as.data.frame()

# plotting by region
Gplot_cca <- ggplot(L) +
  geom_point(aes(x = CCA1, y = CCA2, color = group), shape = 19, size = 2, alpha = 0.7) +
  coord_fixed() +
  geom_segment(data = ccavectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = species_data, aes(x = CCA1, y = CCA2), shape = 17, size = 2, color = "slateblue") +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(-5, 8)) +
  geom_text(data = ccavectors, aes(x = CCA1, y = CCA2, label = rownames(ccavectors)), nudge_x = 0.3, nudge_y = 0.3) +
  labs(title = "Canonical Correspondence Analysis by Region")


Gplot_cca <- Gplot_cca + geom_text(data = species_data, aes(x = CCA1, y = CCA2, label = rownames(species_data)), nudge_x = 0.3, nudge_y = 0.3)
Gplot_cca
