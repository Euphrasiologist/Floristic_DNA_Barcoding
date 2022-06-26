# The genetic landscape of hybridisation in the British flora #
# Max R. Brown, Peter M. Hollingsworth, Laura Forrest, Michelle Hart, Laura Jones, Col Ford, Tim Rich, Natasha de Vere, Alex D. Twyford #

# All code Max Brown #
# Last updated: 04.06.20 #

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

# this first chunk of code takes you from the output of the IQTREE software to a phylogeny
# that has had taxa removed and taxonomy updated. Go to tree.3 (~ line 200) to skip the 
# tree cleaning.

# load the tree
tree <- read.tree(file = "../data/Tree_files/ITS_correct.treefile") # 1421 tips
# root the tree with the most basal taxon(?)
tree.2 <- root(tree, outgroup = "Taxus_baccata", resolve.root = TRUE)
# remove Frankenia laevis (branch length too long), Euphrasia officinalis (an aggregate taxon), 
# Berberis vulgaris (not in constraint tree, grouped with Euphrasia...),
# Cerastium arcticum (syn of Cerastium nigrescens), Triglochin maritima (duplicate of Triglochin maritimum, spelled wrong), 
# Alchemilla cairnwellensis (taxon doesnt exist?), Rubus arcticus (extinct since the 1990's at least), Salix fragilis
# is now known as Salix x fragilis (alba x euxina). Callitriche hamulata is subsp of C. brutia. Lythrum hyssopifolium
# is duplicated

# remove apomictic taxa for phylogenetic analysis: Taraxacum agg, Hieracium agg and Rubus fruticosus
tree.2 <- ape::drop.tip(phy = tree.2, tip = c("Frankenia_laevis", "Berberis_vulgaris", "Euphrasia_officinalis", 
                                         "Cerastium_arcticum", "Triglochin_maritima", "Alchemilla_cairnwellensis",
                                         "Rubus_arcticus", "Taraxacum_agg", "Hieracium_agg", "Rubus_fruticosus",
                                         "Salix_fragilis", "Callitriche_hamulata", "Lythrum_hyssopifolium")) # 1408 tips

# make tiplabs pretty
tree.2$tip.label <- gsub(tree.2$tip.label, pattern = "_", replacement = " ", fixed = TRUE)
# save this version of the tree
tree.2.1 <- tree.2
# this might mess things up... but it is ideal to standardise the taxa names in the phylogeny and external sources.
dashes <- lapply(X = strsplit(tree.2$tip.label, " "), FUN = function(x) ifelse(length(x) > 2, paste(x[1], paste(x[2], x[3], sep = "-"), sep = " "), paste(x[1], x[2], sep = " ")))
dashes <- sapply(dashes, "[[", 1)
# get the updated names from Plant List.
updated_species <- Taxonstand::TPL(dashes)
# isolate the new names
updated_species2 <- setDT(updated_species)[, .(paste(New.Genus, New.Species, sep = " "))][[1]]
# sort out Dactylorhiza
updated_species2[grepl("Dactylorhiza majalis", updated_species2)] <- c("Dactylorhiza trausteinerioides", "Dactylorhiza kerryensis")
# sort out Euphrasia, the officinalis generated here must be E. anglica
updated_species2[grepl("Euphrasia campbellae", updated_species2)] <- "Euphrasia campbelliae"
updated_species2[grepl("Euphrasia officinalis", updated_species2)] <- "Euphrasia anglica"
# note that Carex lepidocarpa, demissa and oederi are treated as one taxon by the tree (viridula)
# but Stace treats them as three species.

# save this object.
#save(updated_species2, file = "./hybridpropensity/data/finaltreenames.RData")

  ### Keep this here for now, but these names cannot all be fixed. UPDATE THIS BIT
# taco4 is from Part_Four_Models.R. It tells us which names are wrong... or that are in our data here, 
# but not in the distribution database
correct_names <- setdiff(setDT(updated_species)[, .(paste(New.Genus, New.Species, sep = " "))][[1]], taco4$Species)
  ###

# update the species names on the tree
#tree.2$tip.label <- setDT(updated_species)[, .(paste(New.Genus, New.Species, sep = " "))][[1]]
tree.2$tip.label <- updated_species2


# now clean up tree to include Stace 4 names
# make sure to do the ones changed in hybriddata
# I know this is crazy, but wanted to manually curate the names
tree.2$tip.label[grepl("Senecio erucifolius", tree.2$tip.label)] <- "Jacobaea erucifolia"
tree.2$tip.label[grepl("Sedum telephium", tree.2$tip.label)] <- "Hylotelephium telephium"
tree.2$tip.label[grepl("Elymus farctus", tree.2$tip.label)] <- "Elymus junceiformis"
tree.2$tip.label[grepl("Rosa obtusifolia", tree.2$tip.label)]  <- "Rosa tomentella"
tree.2$tip.label[grepl("Erodium lebelii", tree.2$tip.label)] <- "Erodium aethiopicum"
# now onto the other names
tree.2$tip.label[grepl("Papaver hybridum", tree.2$tip.label)]  <- "Roemeria hispida"
tree.2$tip.label[grepl("Papaver argemone", tree.2$tip.label)]  <- "Roemeria argemone"
tree.2$tip.label[grepl("Saxifraga nivalis", tree.2$tip.label)]  <- "Micranthes nivalis"
tree.2$tip.label[grepl("Saxifraga stellaris", tree.2$tip.label)]  <- "Micranthes stellaris"
tree.2$tip.label[grepl("Sedum rosea", tree.2$tip.label)]  <- "Rhodiola rosea"
tree.2$tip.label[grepl("Sedum forsterianum", tree.2$tip.label)]  <- "Petrosedum forsterianum"
tree.2$tip.label[grepl("Vicia sylvatica", tree.2$tip.label)]  <- "Ervilia sylvatica"
tree.2$tip.label[grepl("Vicia hirsuta", tree.2$tip.label)]  <- "Ervilia hirsuta"
tree.2$tip.label[grepl("Vicia parviflora", tree.2$tip.label)]  <- "Ervum gracile"
tree.2$tip.label[grepl("Vicia tetrasperma", tree.2$tip.label)]  <- "Ervum tetraspermum"
tree.2$tip.label[grepl("Rosa obtusifolia", tree.2$tip.label)]  <- "Rosa tomentella"
tree.2$tip.label[grepl("Rosa xanthina", tree.2$tip.label)]  <- "Rosa spinosissima"
tree.2$tip.label[grepl("Viola persicifolia", tree.2$tip.label)] <- "Viola stagnina"
tree.2$tip.label[grepl("Radiola linoides", tree.2$tip.label)] <- "Linum radiola"
tree.2$tip.label[grepl("Draba muralis", tree.2$tip.label)] <- "Drabella muralis"
tree.2$tip.label[grepl("Persicaria bistorta", tree.2$tip.label)] <- "Bistorta officinalis"
tree.2$tip.label[grepl("Persicaria vivipara", tree.2$tip.label)] <- "Bistorta vivipara"
tree.2$tip.label[grepl("Polygonum arenastrum", tree.2$tip.label)] <- "Polygonum depressum"
tree.2$tip.label[grepl("Minuartia hybrida", tree.2$tip.label)] <- "Sabulina tenuifolia"
tree.2$tip.label[grepl("Minuartia rubella", tree.2$tip.label)] <- "Sabulina rubella"
tree.2$tip.label[grepl("Minuartia verna", tree.2$tip.label)] <- "Sabulina verna"
tree.2$tip.label[grepl("Minuartia stricta", tree.2$tip.label)] <- "Sabulina stricta"
tree.2$tip.label[grepl("Minuartia sedoides", tree.2$tip.label)] <- "Cherleria sedoides"
tree.2$tip.label[grepl("Myosoton", tree.2$tip.label)] 
tree.2$tip.label[grepl("Chenopodium polyspermum", tree.2$tip.label)]  <- "Lipandra polysperma"
tree.2$tip.label[grepl("Chenopodium hybridum", tree.2$tip.label)]  <- "Chenopodiastrum hybridum"
tree.2$tip.label[grepl("Chenopodium murale", tree.2$tip.label)]  <- "Chenopodiastrum murale"
tree.2$tip.label[grepl("Chenopodium chenopodioides", tree.2$tip.label)] <- "Oxybasis chenopodioides"
tree.2$tip.label[grepl("Chenopodium glaucum", tree.2$tip.label)]  <- "Oxybasis glauca"
tree.2$tip.label[grepl("Chenopodium rubrum", tree.2$tip.label)]  <- "Oxybasis rubra"
tree.2$tip.label[grepl("Chenopodium urbicum", tree.2$tip.label)]  <- "Oxybasis urbica"
tree.2$tip.label[grepl("Chenopodium bonus-henricus", tree.2$tip.label)] <- "Blitum bonus-henricus"
tree.2$tip.label[grepl("Salicornia pusilla", tree.2$tip.label)] <- "Salicornia disarticulata"
tree.2$tip.label[grepl("Anagallis tenella", tree.2$tip.label)] <- "Lysimachia tenella"
tree.2$tip.label[grepl("Anagallis arvensis", tree.2$tip.label)] <- "Lysimachia arvensis"
tree.2$tip.label[grepl("Anagallis minima", tree.2$tip.label)] <- "Lysimachia minima"
tree.2$tip.label[grepl("Centaurium scilloides", tree.2$tip.label)] <- "Centaurium portense"
tree.2$tip.label[grepl("Gentianella anglica", tree.2$tip.label)] <- "Gentianella amarella" # there are subspecies now. This will be a duplicate, so remove later.
tree.2$tip.label[grepl("Lithospermum purpureocoeruleum", tree.2$tip.label)] <- "Aegonychon purpureocaeruleum" # misspelt anyway oops
tree.2$tip.label[grepl("Anchusa arvensis", tree.2$tip.label)] <- "Lycopsis arvensis"
tree.2$tip.label[grepl("Thymus praecox", tree.2$tip.label)] <- "Thymus drucei"
tree.2$tip.label[grepl("Orobanche purpurea", tree.2$tip.label)] <- "Phelipanche purpurea"
tree.2$tip.label[grepl("Filago minima", tree.2$tip.label)] <- "Logfia minima"
tree.2$tip.label[grepl("Filago gallica", tree.2$tip.label)] <- "Logfia gallica"
tree.2$tip.label[grepl("Gnaphalium norvegicum", tree.2$tip.label)] <- "Omalotheca norvegica"
tree.2$tip.label[grepl("Gnaphalium sylvaticum", tree.2$tip.label)] <- "Omalotheca sylvatica"
tree.2$tip.label[grepl("Gnaphalium supinum", tree.2$tip.label)] <- "Omalotheca supina"
tree.2$tip.label[grepl("Erigeron acer", tree.2$tip.label)] <- "Erigeron acris"
tree.2$tip.label[grepl("Apium nodiflorum", tree.2$tip.label)] <- "Helosciadium nodiflorum"
tree.2$tip.label[grepl("Apium repens", tree.2$tip.label)] <- "Helosciadium repens"
tree.2$tip.label[grepl("Apium inundatum", tree.2$tip.label)] <- "Helosciadium inundatum"
tree.2$tip.label[grepl("Petroselinum segetum", tree.2$tip.label)] <- "Sison segetum"
tree.2$tip.label[grepl("Carum verticillatum", tree.2$tip.label)] <- "Trocdaris verticillata"
tree.2$tip.label[grepl("Peucedanum palustre", tree.2$tip.label)] <- "Thysselinum palustre"
tree.2$tip.label[grepl("Ruppia cirrhosa", tree.2$tip.label)] <- "Ruppia spiralis"
tree.2$tip.label[grepl("Kobresia simpliciuscula", tree.2$tip.label)] <- "Carex simpliciuscula"
tree.2$tip.label[grepl("Festuca altissima", tree.2$tip.label)] <- "Drymochloa sylvatica"
tree.2$tip.label[grepl("Festuca arundinacea", tree.2$tip.label)] <- "Schedonorus arundinaceus"
tree.2$tip.label[grepl("Festuca gigantea", tree.2$tip.label)] <- "Schedonorus giganteus"
tree.2$tip.label[grepl("Festuca pratensis", tree.2$tip.label)] <- "Schedonorus pratensis"
tree.2$tip.label[grepl("Pseudosclerochloa rupestris", tree.2$tip.label)] <- "Puccinellia rupestris"
tree.2$tip.label[grepl("Bromus sterilis", tree.2$tip.label)] <- "Anisantha sterilis"
tree.2$tip.label[grepl("Bromus benekenii", tree.2$tip.label)] <- "Bromopsis benekenii"
tree.2$tip.label[grepl("Bromus ramosus", tree.2$tip.label)] <- "Bromopsis ramosa"
tree.2$tip.label[grepl("Bromus erectus", tree.2$tip.label)] <- "Bromopsis erecta"
tree.2$tip.label[grepl("Bromus commutatus", tree.2$tip.label)] <- "Bromus racemosus"
tree.2$tip.label[grepl("Deschampsia flexuosa", tree.2$tip.label)] <- "Avenella flexuosa"
tree.2$tip.label[grepl("Deschampsia setacea", tree.2$tip.label)] <- "Aristavena setacea"
tree.2$tip.label[tree.2$tip.label == "Triglochin concinna"] <- "Triglochin maritima"
tree.2$tip.label[tree.2$tip.label == "Dactylorhiza viridis"] <- "Coeloglossum viride"
tree.2$tip.label[tree.2$tip.label == "Carex cuprina"] <- "Carex otrubae"
tree.2$tip.label[tree.2$tip.label == "Juncus hybridus"] <- "Juncus ranarius"

# the updated names have generated duplicate species names. Christ.
tree.2$tip.label[duplicated(tree.2$tip.label)]

# try numbers and drop tips corresponding to the duplicates
# interestingly, not all are sister taxa
# the above taxonomic changes introduced two duplciates, Gentianella amarella and Bromus racemosus
numtabs <- data.table(Num = 1:length(tree.2$tip.label), tree.2$tip.label)[V2 %in% tree.2$tip.label[duplicated(tree.2$tip.label)]]
# this shows the original names
data.table(Num = 1:length(tree.2.1$tip.label), tree.2.1$tip.label)[Num %in% numtabs$Num]
# double check duplicated names
tree.2$tip.label[tree.2$tip.label == "Ranunculus peltatus"][2] <- "Ranunculus baudotii"
tree.2$tip.label[tree.2$tip.label == "Lamium purpureum"][2] <- "Lamium hybridum"
tree.2$tip.label[tree.2$tip.label == "Festuca rubra"][1] <- "Festuca arenaria"
tree.2$tip.label[tree.2$tip.label == "Brachypodium pinnatum"][2] <- "Brachypodium rupestre"
tree.2$tip.label[tree.2$tip.label == "Asparagus officinalis"][2] <- "Asparagus prostratus"
tree.2$tip.label[tree.2$tip.label == "Ulmus minor"] <- c("Ulmus minor", "Ulmus plotii", "Ulmus procera")
tree.2$tip.label[tree.2$tip.label == "Trifolium repens"][1] <- "Trifolium occidentale"
tree.2$tip.label[tree.2$tip.label == "Ononis spinosa"][1] <- "Ononis repens"
tree.2$tip.label[tree.2$tip.label == "Erophila verna"] <- c("Erophila glabrescens", "Erophila verna", "Erophila majuscula")
# I have issue with Euphorbia peplis only I think? I need a goddamn taxonomist...
data.table(Num = 1:length(tree.2.1$tip.label), tree.2.1$tip.label)[grepl("Euphorbia", V2)] # #1251
tree.2$tip.label[1251] <- "Euphorbia peplis"

numtabs2 <- data.table(Num = 1:length(tree.2$tip.label), tree.2$tip.label)[V2 %in% tree.2$tip.label[duplicated(tree.2$tip.label)]]
numtabs2 
# tips to delete, numerically.
# create vector manually...
vec.tree <- c(354, 737)

# create a data frame of old and new names
name_changes <- data.table(old = updated_species2,
                           new = tree.2$tip.label)

#write.csv(name_changes, "./hybridpropensity/data/name_changes.csv")

tree.3 <- ape::drop.tip(phy = tree.2, tip = vec.tree)

##### Part 1.2: Plotting phylogeny #####

# read it in again
#tree.3 <- ape::read.tree(file = "./data/Tree_files/tree.3.stace4")
tree.3 <- ape::read.tree(file = "./Barcoding_Phylogeny_Dating/DNA_Barcoding.dated.treefile")

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
finaldata <- fread("./data/Hybrid_flora_of_the_British_Isles/finaldata_updated210220.csv")[,-"V1"]

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
  x = tree.vcv5d$tip.label,
  Taxa = gsub(pattern = "_", replacement = " ", x = tree.vcv5d$tip.label, fixed = TRUE)
)
# make it a data table
setDT(data_modes); data_modes
# resurrect the genus column
data_modes$Genus <- gsub(pattern = " .*", replacement = "", x = data_modes$Taxa)

##### Part 3: Circular Plot of Phylogeny #####

# set colours 
palette(RColorBrewer::brewer.pal(8, name = "Greens"))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# make the dendrogram
dend_modes <- phylogram::as.dendrogram(phytools::force.ultrametric(tree.vcv5d))
# extract the height of the dendrogram
dend_modes_height <- attr(dend_modes, "height")
# make sure row order is preserved (before or after merge?)
# model data.
data2_modes <- data_modes[data.table(x = factor(labels(dend_modes), levels = unique(labels(dend_modes)))), on = "x"]
data2_modes$pos <- 1:nrow(data2_modes)

# add higher taxonomic levels
higher_orders <- setDT(taxonlookup::lookup_table(species_list = as.character(data2_modes$Taxa)))
setnames(x = higher_orders, old = c("genus", "family", "order", "group"), new = c("Genus", "Family", "Order", "Group"))
# merge with phylogeny data
data2_modes <- data2_modes[higher_orders, on = "Genus"]

# 17.1.20 just had a thought that the circular phylogeny could have the
#posterior modes of the probability of each taxon hybridising
load("../data/Backups/mcmc.vcv5.RData") # or another run of mcmc.vcv5. May not upload as it's 100+Mbs

tree_modes <- VCVglmm::Solapply(mcmc.vcv5d, FUN = mean)[Group == "Sp1"][order(-Grouped_Value)]
tree_modesp <- tree_modes[!grepl("Node", Variable)][, Taxa := substring(Variable, 5)]


data2.1_modes <- tree_modesp[, .(Taxa, Mean = Grouped_Value)][data2_modes, on = .(Taxa)]

data2.1_both <- data2[, .(Taxa, Hybrid_propensity)][data2.1_modes, on = .(Taxa)]
#data2.1_both[, Hybrid_propensity := ifelse(is.na(Hybrid_propensity), 0, Hybrid_propensity)]

# create scaled hybrid propensity.
# get genus size from vcv5d

data2.1_both <- data2.1_both[unique(vcv5d[, .(Genus = Genus1, Genus.Size)]), on = .(Genus)]
# remove NA
data2.1_both <- data2.1_both[!is.na(Taxa)]

# remove duplicates
data2.1_both <- data2.1_both[!duplicated(pos)]

# order by position
data2.1_both <- data2.1_both[order(pos)]

# hybridisation propensity
data2.1_both$Hybrid_propensity <- ifelse(is.na(data2.1_both$Hybrid_propensity),
       0, data2.1_both$Hybrid_propensity)
data2.1_both$Hybrid_propensity <- data2.1_both$Hybrid_propensity/data2.1_both$Genus.Size
data2.1_both$Hybrid_propensity <- ifelse(is.na(data2.1_both$Hybrid_propensity),
                                         0, data2.1_both$Hybrid_propensity)

treevar.vcv5d <- VCVglmm::Bbar(phylo = tree.vcv5d, type = "Species")*mcmc.vcv5d$VCV[,"Sp1+Sp2"]
# make the estimates of the mode on the probability scale, averaged over random effects
marg_post_mode <- pnorm(
  summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
    summary(mcmc.vcv5d)$solutions[3, 1] * mean(vcv5d$Hectads_shared) + # at mean hectad sharing
    summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[4, 1] + #annual - perennial hybridisation
    summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) +  # effect of branch length
    data2.1_both$Mode,
  sd = sqrt(mean(
    mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d # accounting for phylogeny
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
## ylim = c(-3.9, 3.2) if CI's wanted.
circos.trackPlotRegion(factors = data2.1_both$factor, y=data2.1_both$Mean, panel.fun = function(x, y) {
  # plot the 10 most species rich orders
  dat <- data2.1_modes[,.(Min = min(pos),
                    Max = max(pos),
                    Order = unique(Order),
                    N =.N), by = "Order"][, -4][order(-N)][1:20,]
  dat <- as.data.frame(lapply(dat, unlist))
  for(i in 1:10){
    circos.lines(x = c(dat[i,2], dat[i, 3]), y = c(2.9, 2.9), col = 5, lwd = 4) # line
    circos.text(x = mean(c(dat[i,2], dat[i, 3])), y = 3.8, 
                labels = paste(as.character(dat[i,1])), 
                facing = "outside", 
                niceFacing = TRUE, cex = 3) # label 
  }
  # plot the 5 genera with the highest mode sum
  dat2 <- data2.1_modes[,.(Min = min(pos),
                     Max = max(pos),
                     Genus = unique(Genus),
                     N =.N,
                     Mean = sum(Mean)), by = "Genus"][, -4][order(-Mean)][1:20,]
  
  dat2 <- as.data.frame(lapply(dat2, unlist))
  dat2 <- dat2[1:5,]
  dat2$ID <- c(1:5)
  for(i in 1:5){
    circos.lines(x = c(dat2[i,2], dat2[i, 3]), y = c(2, 2), col = 1, lwd = 1, lty = 2, area = TRUE) # line
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
circos.trackLines(factors = data2.1_both[order(pos)]$factor, x = data2.1_both[order(pos)]$pos, y = data2.1_both[order(pos)]$Mean, col = cbPalette[3], lwd=4.5)
# hybrid propensity
circos.trackPlotRegion(factors = data2.1_both$factor, y=data2.1_both$Hybrid_propensity,ylim = c(0,6), panel.fun = function(x,y){
  circos.yaxis(side = "right", labels.niceFacing = TRUE, labels.cex = 1.5)
})
circos.trackLines(factors = data2.1_modes$factor, x = data2.1_modes$pos, y = data2.1_both$Hybrid_propensity, col = 6, lwd=4.5)

# add ultrametric dendrogram. slow!
circos.track(ylim = c(0, dend_modes_height), bg.border = NA, 
             track.height = 0.4, panel.fun = function(x, y) {
               circos.dendrogram(dend_modes)
             })

# back up dat2.1_both
write.csv(x = data2.1_both, file = "../data/Backups/data2.1_both200220.csv")
