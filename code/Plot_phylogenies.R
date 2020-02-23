# Floristic DNA Barcoding reveals the landscape of hybridisation in the British flora #
# Max Brown, Jarrod Hadfield, Peter M. Hollingsworth, Alex D. Twyford #

# All code Max Brown: 20.02.20 #

# Libraries needed # 
library(ggplot2)
library(data.table)
library(ape)
library(phangorn)
library(phytools)
library(ggtree)
library(circlize)
library(phylogram)
library(taxonlookup)
library(treeio)
library(dplyr)
library(Taxonstand)
library(RColorBrewer)

##### Part 1: Import and plot phylogeny ##### 

##### Part 1.1: Cleaning phylogeny #####

# phylogeny has been generated using IQTREE and three concatenated genes, two plastid
# rbcL and matK, and one nuclear, ITS. The ITS was aligned only between species in the
# same genus and the plastid genes were used to resolve deeper nodes in the phylogeny.
# The tree was also constrained in its topology, using the latest (APG4) angiosperm 
# phylogeny to make sure that deep relationships were constrained to be correct.

##### Part 1.2: Plotting phylogeny #####

# read in tree file
tree.3 <- read.tree(file = "../data/Tree_files/tree.3.stace4")

# use this to determine which clade to extract, should we want to
ggtree(tree.3)+geom_tiplab(cex = 0.9)+geom_text2(aes(subset=!isTip, label=node), hjust=-.5)

# For example Rosaceae at node 1518. The Sorbi, although polytomous within section,
# resolves the sections Soraria, Aria and Tormaria
ggtree(extract.clade(tree.3, node = 1531))+geom_tiplab()

# Apiaceae (+ Araliaceae...)
ggtree(extract.clade(tree.3, node = 1945))+geom_tiplab()
# Orobanchaceae, observe complexity of Euprhasia
ggtree(extract.clade(tree.3, node = 2036 + 17))+geom_tiplab()
# Cyperaceae, Carex is most species rich genus
ggtree(extract.clade(tree.3, node = 2559))+geom_tiplab()
# Poaceae
ggtree(extract.clade(tree.3, node = 2418 + 18))+geom_tiplab()
# Asteraceae. The split to Gnaphalium and Filago really fit the tree here!
ggtree(extract.clade(tree.3, node = 1792 + 16))+geom_tiplab()
# make tiplabs ugly again (must do)
tree.3$tip.label <- gsub(tree.3$tip.label, pattern = " ", replacement = "_", fixed = TRUE)

##### Part 2: Make sure that names match data and phylogeny #####

# create the data that will plot the hybrid propensity.
finaldata <- fread("../data/Hybrid_flora_of_the_British_Isles/finaldata_updated210220.csv")[,-"V1"]

data <- data.table(
  factor = factor(rep(1, 1406)),
  x = tree.3$tip.label,
  Taxa = gsub(pattern = "_", replacement = " ", x = tree.3$tip.label, fixed = TRUE)
)
data2 <- finaldata[data, on = c("Taxa")]
# if hybrid propensity is NA, change it to zero
data2[is.na(data2$Hybrid_propensity)]$Hybrid_propensity <- 0
# resurrect the genus column
data2$Genus <- gsub(pattern = " .*", replacement = "", x = data2$Taxa)

# create the data that will plot the posterior modes of the phylogenetic model
# need the tree that was used in the final model to be loaded here.
load("../data/Backups/tree.vcv5.RData")

data_modes <- data.frame(
  factor = factor(rep(1, 1110)),
  x = tree.vcv5$tip.label,
  Taxa = gsub(pattern = "_", replacement = " ", x = tree.vcv5$tip.label, fixed = TRUE)
)
# make it a data table
setDT(data_modes); data_modes
# resurrect the genus column
data_modes$Genus <- gsub(pattern = " .*", replacement = "", x = data_modes$Taxa)

##### Part 3: Circular Plot of Phylogeny #####

# set colours 
palette(RColorBrewer::brewer.pal(8, name = "Greens"))

# make the dendrogram
dend_modes <- phylogram::as.dendrogram(phytools::force.ultrametric(tree.vcv5))
# extract the height of the dendrogram
dend_modes_height <- attr(dend_modes, "height")
# make sure row order is preserved (before or after merge?)
# model data.
data2_modes <- data_modes[data.table(x = factor(labels(dend_modes), levels = unique(labels(dend_modes)))), on = "x"]
data2_modes$pos <- 1:nrow(data2_modes)

# add higher taxonomic levels
higher_orders <- setDT(lookup_table(species_list = as.character(data2_modes$Taxa)))
setnames(x = higher_orders, old = c("genus", "family", "order", "group"), new = c("Genus", "Family", "Order", "Group"))
# merge with phylogeny data
data2_modes <- data2_modes[higher_orders, on = "Genus"]

# 17.1.20 just had a thought that the circular phylogeny could have the
#posterior modes of the probability of each taxon hybridising
load("../data/Backups/mcmc.vcv5.RData") # or another run of mcmc.vcv5. 

tree_modes <- VCVglmm::Solapply(mcmc.vcv5)[Group == "Sp1"][order(-Grouped_Value)]
tree_modesp <- tree_modes[!grepl("Node", Variable)][, Taxa := substring(Variable, 5)]


data2.1_modes <- tree_modesp[, .(Taxa, Mode = Grouped_Value)][data2_modes, on = .(Taxa)]

data2.1_both <- data2[, .(Taxa, Hybrid_propensity)][data2.1_modes, on = .(Taxa)]
data2.1_both[, Hybrid_propensity := ifelse(is.na(Hybrid_propensity), 0, Hybrid_propensity)]

# remove duplicates
data2.1_both <- data2.1_both[!duplicated(pos)]

treevar.vcv5 <- VCVglmm::Bbar(phylo = tree.vcv5, type = "Species")*mcmc.vcv5$VCV[,"Sp1+Sp2"]
# make the estimates of the mode on the probability scale, averaged over random effects
marg_post_mode <- pnorm(
  summary(mcmc.vcv5)$solutions[1, 1] + # intercept
    summary(mcmc.vcv5)$solutions[3, 1] * mean(vcv5$Hectads_shared) + # at mean hectad sharing
    summary(mcmc.vcv5)$solutions[6, 1] * mean(vcv5$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5)$solutions[4, 1] + #annual - perennial hybridisation
    summary(mcmc.vcv5)$solutions[2, 1] * mean(vcv5$dist2) +  # effect of branch length
    data2.1_both$Mode,
  sd = sqrt(mean(
    mcmc.vcv5$VCV[, "units"] + 2 * mcmc.vcv5$VCV[, "T1+T2"] + treevar.vcv5 # accounting for phylogeny
  ))
) 
data2.1_both[, marg_post_mode := marg_post_mode]


# General parameters
circos.clear()
# original parameters
#circos.par("track.height" = 0.4, "gap.degree" = 10, "track.margin" = c(0,0))
# adding in scaled hp
circos.par("track.height" = 0.2, "gap.degree" = 8, "track.margin"=c(0,0))

## attempt at plotting posterior modes.
# Initialize chart
circos.initialize(factors = data2.1_both$factor, x = data2.1_both$pos)
circos.trackPlotRegion(factors = data2.1_both$factor, y=data2.1_both$Mode, panel.fun = function(x, y) {
  # plot the 10 most species rich orders
  dat <- data2.1_modes[,.(Min = min(pos),
                    Max = max(pos),
                    Order = unique(Order),
                    N =.N), by = "Order"][, -4][order(-N)][1:20,]
  dat <- as.data.frame(lapply(dat, unlist))
  for(i in 1:10){
    circos.lines(x = c(dat[i,2], dat[i, 3]), y = c(1.9, 1.9), col = 5, lwd = 4) # line
    circos.text(x = mean(c(dat[i,2], dat[i, 3])), y = 2.8, 
                labels = paste(as.character(dat[i,1])), 
                facing = "outside", 
                niceFacing = TRUE, cex = 3) # label 
  }
  # plot the 5 genera with the highest mode sum
  dat2 <- data2.1_modes[,.(Min = min(pos),
                     Max = max(pos),
                     Genus = unique(Genus),
                     N =.N,
                     Mode = sum(Mode)), by = "Genus"][, -4][order(-Mode)][1:20,]
  
  dat2 <- as.data.frame(lapply(dat2, unlist))
  dat2 <- dat2[1:5,]
  dat2$ID <- c(1:5)
  for(i in 1:5){
    circos.lines(x = c(dat2[i,2], dat2[i, 3]), y = c(1.5, 1.5), col = 1, lwd = 1, lty = 2, area = TRUE) # line
    circos.text(x = mean(c(dat2[i,2], dat2[i, 3])), y = -1.5, 
                labels = paste(as.character(dat2[i,6])), 
                facing = "outside", 
                niceFacing = TRUE, cex = 1.5) 
  }
  
  # plot the y axis
  circos.yaxis(side = "right", labels.niceFacing = TRUE, labels.cex = 1.5)
  #circos.axis(labels = data2.1$Taxa, labels.cex=0.5, labels.font=1, lwd=0.8, h="top", direction="outside", labels.niceFacing = T, labels.facing = "clockwise", major.at = 1:1408)
})
# add hybridisation lines across the phylogeny
# zero line
circos.trackLines(factors = data2.1_both$factor, x = data2.1_both$pos, y = rep(0, length(data2.1_both$pos)), col = scales::alpha("tomato", 0.4), lwd=3)
# order messed up otherwise
circos.trackLines(factors = data2.1_both[order(pos)]$factor, x = data2.1_both[order(pos)]$pos, y = data2.1_both[order(pos)]$Mode, col = "lightskyblue3", lwd=4.5)

circos.trackPlotRegion(factors = data2.1_both$factor, y=data2.1_both$Hybrid_propensity,ylim = c(0,15), panel.fun = function(x,y){
  circos.yaxis(side = "right", labels.niceFacing = TRUE, labels.cex = 1.5)
})
circos.trackLines(factors = data2.1_modes$factor, x = data2.1_modes$pos, y = data2.1_both$Hybrid_propensity, col = 5, lwd=4.5)

# add ultrametric dendrogram. slow!
circos.track(ylim = c(0, dend_modes_height), bg.border = NA, 
             track.height = 0.4, panel.fun = function(x, y) {
               circos.dendrogram(dend_modes)
             })

# back up dat2.1_both
write.csv(x = data2.1_both, file = "../data/Backups/data2.1_both200220.csv")
