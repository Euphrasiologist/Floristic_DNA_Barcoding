# Floristic DNA Barcoding reveals the landscape of hybridisation in the British flora #
# Max Brown, Jarrod Hadfield, Peter M. Hollingsworth, Alex D. Twyford #

# All code Max Brown: 20.02.20 #

# Libraries needed throughout # 
library(ggplot2)
library(data.table)
devtools::install_github("Euphrasiologist/VCVglmm")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)
library(MCMCglmm)
library(lme4)
library(plyr)
library(dplyr)
library(stringr)
library(reshape2)
library(phangorn)
library(phytools)
library(coda)
library(ape)
devtools::install_github("https://github.com/YuLab-SMU/ggtree")
library(ggtree)
library(ggstance)
library(brranching)
library(truncreg)
library(Taxonstand)
library(circlize)
library(phylogram)
library(treeio)
library(seqinr)
library(Hmisc)
library(numbers)
library(stringr)

# import hybriddata
hybriddata <- fread("../data/Hybrid_flora_of_the_British_Isles/hybriddata_Stace4_removed.csv")

# the columns of interest are parent A and parent B. Save these.
parenta<-hybriddata$`Parent A`
parentb<-hybriddata$`Parent B`
# concatenate the columns
catparents<-c(as.character(parenta), as.character(parentb))
# table creates the frequency of each parent, sort orders species alphabetically
hybridpropensity<-data.table(sort(table(catparents)))
# set the names
setnames(x = hybridpropensity, old = colnames(hybridpropensity), new = c("Taxa", "Hybrid_propensity"))
# create the genus factor
hybridpropensity$Genus <- sub(" .*", "", hybridpropensity$Taxa)

# finaldata is the hybrid propensity of each species, **now purely for plotting purposes**
finaldata <- hybridpropensity

# remove hybrids
finaldata <- finaldata[!grepl(pattern = " x ", x = finaldata$Taxa),]

# pre merging with tree!
finaldata[Hybrid_propensity >0,] 
# species level hybridisation propensity that is scaled (sensu Whitney et al. 2010)
finaldata2 <- finaldata[, .(N = .N), by = "Genus"][finaldata, on = "Genus"][, Scaled.HP := 100*Hybrid_propensity/((N^2-N)/2)]
finaldata$Scaled.HP <- ifelse(test = is.nan(finaldata2$Scaled.HP) | is.infinite(finaldata2$Scaled.HP), yes = 0, no = finaldata2$Scaled.HP)

# write this 
fwrite(x = finaldata, file = "../data/Hybrid_flora_of_the_British_Isles/finaldata_updated210220.csv")
