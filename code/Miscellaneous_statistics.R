# The genetic landscape of hybridisation in the British flora #
# Max R. Brown, Peter M. Hollingsworth, Laura Forrest, Michelle Hart, Laura Jones, Col Ford, Tim Rich, Natasha de Vere, Alex D. Twyford #

# All code Max Brown #
# Last updated: 04.06.20 #


# Questions arising from the manuscript that need specific answers.

# Libraries needed #
library(data.table)
library(ggplot2)
library(Hmisc)
library(dplyr)

# dependent on parts 1:6 and objects therein.

##### Q1. What is the distribution of hybrid propensity? How many hybrids analysed? #####

# requires dat2_1.csv
data2.1 <- fread(file = "./data/Backups/data2.1_both200220.csv")
# add genus size
data2.1 <- data2.1[,.(Taxa, Hybrid_propensity, Genus.Size = .N), by = .(Genus)]
# adjusted the smoothing parameter, it's an exponentially decaying distribution.
ggplot(data2.1, aes(x = Hybrid_propensity))+ geom_histogram()
#+stat_density(adjust = 5)

# dist.sequ.cp4; dist.sequ.its4 needed

dist.sequ.cp4[test == TRUE] # 413
dist.sequ.its4[test == TRUE] # 279

##### Q2. How many multi-species genera do, and how many do not, have naturally occurring hybrids? #####
# for native species only #

# 123 multi species genera have naturally occurring hybrids. Largest genera (Carex, Sorbus, Juncus, Ranunculus, Potamogeton)
data2.1[, .(Genus_size = .N,
           Hybrid_propensity = sum(Hybrid_propensity)), by = "Genus"
       ][
         Genus_size > 1 & Hybrid_propensity > 0
         ][
           order(-Genus_size)]
# 121 multi species genera have no naturally occurring hybrids.
data2.1[, .(Genus_size = .N,
           Hybrid_propensity = sum(Hybrid_propensity)), by = "Genus"
       ][
         Genus_size > 1 & Hybrid_propensity == 0
         ][
           order(-Genus_size)]


##### Q3. Which large genera have no hybrids (give genera names and number of species)? #####

# see Q2. #


##### Q4. Which multi-species genera (with 5+ species) have the highest and lowest congeneric genetic divergence? #####

# PLASTID, dist.sequ.cp4 is the dataset. `test` contains whether hybrid formed or not.
# create genus size
mer <- unique(data2.1[, .(Genus1 = (Genus), Genus.Size)])
data2.2 <- mer[dist.sequ.cp4, on = "Genus1"]

# all genera with >= 5 species, hybridising or not.
data2.2[Genus.Size >=5, .(Mdist = mean(value)), by = "Genus1"][order(-Mdist)]
# only hybridising taxa
data2.2[Genus.Size >=5 & test == TRUE, .(Mdist = mean(value)), by = "Genus1"][order(-Mdist)]
# only hybridising genera
# cleaning
onlyhybs <- unique(data2.2[,.(Genus1, test)])
onlyhybs[, N := .N, by = "Genus1"]
onlyhybs[, Hybridising_genera := ifelse(test = N == 2, yes = TRUE, no = test)]
data2.3 <- unique(onlyhybs[, .(Genus1, Hybridising_genera)])[data2.2, on = "Genus1"]
# answer
data2.3[Genus.Size >=5 & Hybridising_genera == TRUE, .(Mdist = mean(value)), by = "Genus1"][order(-Mdist)]


dist.sequ.cp4[test == TRUE][,.(Mdist = mean(value)), by = .(Genus1)][order(-Mdist)]

# ITS
# merge with data created from Genetic_distance_matrices.R
data2.4 <- mer[dist.sequ.its4, on = "Genus1"]

# all genera with >= 5 species, hybridising or not.
data2.4[Genus.Size >= 5, .(Mdist = mean(value)), by = "Genus1"][order(-Mdist)]
# only hybridising taxa
data2.4[Genus.Size >=5 & test == TRUE, .(Mdist = mean(value), N=.N), by = "Genus1"][order(-Mdist)]

## THIS BIT IS IN THE MANUSCRIPT ##
# only hybridising genera
# cleaning
onlyhybs2 <- unique(data2.4[,.(Genus1, test)])
onlyhybs2[, N := .N, by = "Genus1"]
onlyhybs2[, Hybridising_genera := ifelse(test = N == 2, yes = TRUE, no = test)]
data2.5 <- unique(onlyhybs2[, .(Genus1, Hybridising_genera)])[data2.4, on = "Genus1"]
# answer(s)
its_5plus_taxa_dist_hybs <- data2.5[Genus.Size >=5 & Hybridising_genera == TRUE, .(Mdist = mean(value)), by = "Genus1"][order(-Mdist)]
data2.5[test == TRUE & Genus.Size >= 5, .(Mdist = mean(value), SE = sd(value)/sqrt(.N), N=.N), by = .(Genus1)][order(-Mdist)]

# export these results?
fwrite(its_5plus_taxa_dist_hybs, "./data/supp_its_5plus_taxa_dist_genera.csv")


##### Q4.1 Do annuals have on average lower genetic distance than perennials? #####

### ITS ###
# merge with annual perennial data
annper4 <- fread("../data/Trait_databases/annual_perennial_britf_fl_updated.csv")

test <- annper4[data2.5, on = .(Species = rn)]
annper5<- rbind(annper4, data.table(Species = unique(test[is.na(AnnPer)]$Species),
           AnnPer = c("p", "p", "p", "a", "p", "p", "p", "p", "a", "a", 
                      "a", "p", "a", "p", "p", "p", "p", "p", "p", "p", "p",
                      "a", "a", "p", "p", "p", "a", "p", "p", "p", "p",
                      "p", "a", "p")))
annpertests <- annper5[data2.5, on = .(Species = rn)]
annpertests[AnnPer == ""]$AnnPer <- "a"
annpertests[AnnPer == "b"]$AnnPer <- "p"

# annuals do have on average lower genetic distance, but the lowest distance class
# is perennial species that hybridise.
annpertests[, .(Mdist=mean(value), SEM=sd(value)/sqrt(.N)), by = .(AnnPer, test)]

### CPDNA ###

test2 <- annper5[data2.3, on = .(Species = rn)][is.na(AnnPer)]
annper6 <- rbind(annper5, data.table(Species = unique(test2[is.na(AnnPer)]$Species),
           AnnPer = c("p", "p", "p", "p", "p", "p", "p", "a", "p"))) 
annpertests2 <- annper6[data2.3, on = .(Species = rn)]
annpertests2[AnnPer == ""]$AnnPer <- "a"
annpertests2[AnnPer == "b"]$AnnPer <- "p"

annpertests2[, .(Mdist=mean(value), SEM=sd(value)/sqrt(.N)), by = .(AnnPer, test)]

##### Q5. What is the mean and SE for genetic distance for the hybridisation boxplots? Are the statistically significant differences between hybridising and non-hybridising taxa? #####

# PLASTID
# again, the dataset is dist.sequ.cp4
dist.sequ.cp4[,.(Mean_CP = mean(value), 
                 SE_CP = sd(value)/sqrt(.N)), by = "test"] %>%
  # plot
  ggplot(aes(x = test, y = Mean_CP))+
  geom_point()+
  theme_bw()+
  geom_errorbar(aes(ymin = Mean_CP -SE_CP, ymax = Mean_CP + SE_CP))

# very significantly different.
wilcox.test(value ~ test, data = dist.sequ.cp4, conf.int = TRUE)

# boxplot plastid


# ITS
# again, the dataset is dist.sequ.its4
dist.sequ.its4[,.(Mean_ITS = mean(value), 
                  SE_ITS = sd(value)/sqrt(.N)), by = "test"] %>%
  # plot
  ggplot(aes(x = test, y = Mean_ITS))+
  geom_point()+
  theme_bw()+
  geom_errorbar(aes(ymin = Mean_ITS -SE_ITS, ymax = Mean_ITS + SE_ITS))

# very significantly different.
wilcox.test(conf.int = TRUE, value ~ test, data = dist.sequ.its4)

# mean and standard error for ITS.
dist.sequ.its4[, .(Mean = mean(value),
                   SE = sd(value)/sqrt(.N)), by = "test"]

##### Q6. Which are good examples of closely related taxa that do not hybridise? (What is their mean genetic distance and number of species per genus) #####

# PLASTID
# ALL GENERA
data2.2[Genus.Size >=2, .(Mdist = mean(value)), by = "Genus1"][order(Mdist)][1:10]
# ONLY NON-HYBRIDISING TAXA
data2.2[Genus.Size >=2 & test == FALSE, .(Mdist = mean(value)), by = "Genus1"][order(Mdist)][1:10]
# ONLY NON-HYBRIDISING GENERA
data2.3[Genus.Size >=2 & Hybridising_genera == FALSE, .(Mdist = mean(value)), by = "Genus1"][order(Mdist)][1:10]

# ITS
# ALL GENERA (used in manuscript.)
data2.4[Genus.Size >=2, .(Mdist = mean(value)), by = "Genus1"][order(Mdist)][1:10]
# ONLY NON-HYBRIDISING TAXA
data2.4[Genus.Size >=2 & test == FALSE, .(Mdist = mean(value)), by = "Genus1"][order(Mdist)][1:10]
# ONLY NON-HYBRIDISING GENERA
data2.5[Genus.Size >=2 & Hybridising_genera == FALSE, .(Mdist = mean(value)), by = "Genus1"][order(Mdist)][1:10]

