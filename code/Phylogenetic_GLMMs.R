# The genetic landscape of hybridisation in the British flora #
# Max R. Brown, Peter M. Hollingsworth, Laura Forrest, Michelle Hart, Laura Jones, Col Ford, Tim Rich, Natasha de Vere, Alex D. Twyford #

# All code Max Brown #
# Last updated: 04.06.20 #

# Libraries needed #
#install_github("Euphrasiologist/VCVglmm")
library(ape)
library(phytools)
library(MCMCglmm)
library(data.table)
library(Taxonstand)
library(pbapply)
library(ggplot2)
library(devtools)
library(VCVglmm)
library(gridExtra)
library(lattice)

# functions specified #

# useful little function that splits a character vector and then sorts it
# alphabetically after the split.
# the pbapply can be replaced with apply.
sort_names <- function(x, split = " ", collapse = ""){
  sort1 <- pbapply::pblapply(strsplit(x = x, split = split), sort)
  sort2 <- pbapply::pblapply(sort1, function(x) paste(x, collapse = collapse))
  return(unlist(sort2))
}

# thanks https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##### Part 1: Phylogenetic generalised linear mixed effect models. #####

# The simplest possible model here, with no explanatory variables.
# This is because we are simply interested in variance components
tree.3 <- read.tree(file = "../data/Tree_files/tree.3.stace4")
tree.3$tip.label <- gsub(pattern = "_", replacement = " ", tree.3$tip.label)
tree.3 # 1406 tips!


# otherwise the inverseA function fails.
tree.3$node.label <- NULL

# multi membership model
# y <- vector of 0's and 1's as to whether a species pair gave rise to a successful hybrid or not
# Sp1 <- vector of species 1
# Sp2 <- vector of species 2
# T1 <- copy of Sp1
# T2 <- copy of Sp2
# dist <- variance covariance matrix of phylogeny
# cooc <- co-occurrence of parental taxa (would be good if I could get the data!!)

# Sp1, Sp2 and the VCV
# just to generate the pairwise comparisons, vcv values discarded later.
v <- vcv(tree.3)

# create the data table object
vcv1 <- data.table(as.data.frame(as.matrix(v)), keep.rownames = TRUE)
vcv2 <- melt.data.table(vcv1)
setnames(x = vcv2, old = c("rn", "variable", "value"), new = c("Sp1", "Sp2", "dist"))

# create column which we can sort on later (expensive call)
vcv2[, On := sort_names(paste(Sp1, Sp2, sep = " "))]
# remove dist, updated dist2 below
vcv2 <- vcv2[,-"dist"]
# genus columns
vcv2[, c("Genus1", "Genus2") := list(gsub(" .*", "", Sp1),
                                     gsub(" .*", "", Sp2))]
# these are all of the hybrids.
hybriddata <- fread("../data/Hybrid_flora_of_the_British_Isles/hybriddata_Stace4_removed.csv")
merg <- hybriddata[, On := sort_names(paste(`Parent A`, `Parent B`, sep = " "))][, .(On = On,
                                                                                 y = 1)]
# merge the data sets
vcv3 <- merg[vcv2, on = "On"]

# let's add in the annual perennial data.
# need to download from the bsbi website and standardise names.

# Are all hybrids represented? No... what's going on???
# it's because the phylogeny is not entirely sequenced/correct.
# Nothing that can be done about this right now.
setdiff(merg$On[grepl(pattern = "Euphrasia", merg$On)], 
        vcv3[Genus1 == "Euphrasia" & y == 1 & !duplicated(On)]$On)

# add in data on annual and perennial species. From https://database.bsbi.org/checklists.php
annper <- fread("../data/Trait_databases/annual_perennial_britf_fl.csv")
name_changes <- fread("../data/Hybrid_flora_of_the_British_Isles/name_changes.csv")
name_changes <- name_changes[-1, .(Species = V2, New.Species = V3)]
# standardise the names
annpernames <- Taxonstand::TPL(annper$prefTaxonName, silent = FALSE)
# set as a data table
setDT(annpernames) 
# bind the updated names to the annual perennial columns
annper2 <- cbind(annpernames[,.(Species = paste(New.Genus, New.Species, sep = " "))], annper[,.(AnnPer = dataValue)])

# NA values will be introduced. Unfortunately add manually.
annper2.1 <- unique(annper2[name_changes, on = .(Species)][, .(Species = New.Species, AnnPer)])
# to see NA's
annper2.1[is.na(AnnPer)]
annual.perennial <- c("p", "p", "a", "a", "a", "p", "a", "a", "p", "a", "p", "p", "p", "p", "p", "p", "p", "p", rep("p", 9), "p", "a")
annper2.1[is.na(AnnPer), AnnPer := annual.perennial]


# merge with the names in the phylogeny
annper3 <- unique(annper2.1)[unique(vcv3[,.(Species = Sp1)]), on = "Species"]
# why are there 1409 in the new annper, not 1406...
# remove the table from below
annper3[Species %in% annper3[duplicated(Species)]$Species][c(1,3,5)]
setkey(annper3, Species, AnnPer)
# and now subset.
annper4 <- annper3[!.(annper3[Species %in% annper3[duplicated(Species)]$Species][c(1,3,5)])]

fwrite(annper4, "../data/Trait_databases/annual_perennial_britf_fl_updated.csv")

# merge with vcv3
vcv3.1 <- vcv3[annper4, on = .(Sp1 = Species)]
vcv3.2 <- vcv3.1[annper4, on = .(Sp2 = Species)]
# make biennial perennial
vcv3.2[, AnnPer := ifelse(AnnPer == "b", yes = "p", no = AnnPer)]
vcv3.2[, i.AnnPer := ifelse(i.AnnPer == "b", yes = "p", no = i.AnnPer)]
# make annual perennial parental combination
vcv3.2[, Annual_Perennial := paste(AnnPer, i.AnnPer, sep = "-")]
# make p-a equal to a-p
vcv3.2[, Annual_Perennial := ifelse(Annual_Perennial == "p-a", yes = "a-p", no = Annual_Perennial)]
# remove AnnPer and i.AnnPer
vcv3.2 <- vcv3.2[,-c("AnnPer", "i.AnnPer")]

# filter for congeneric
vcv4 <- vcv3.2[Genus1 == Genus2]
# make zero values for y
vcv4[, y := ifelse(test = is.na(y), yes = 0, no = y)]
# duplicate Sp1, Sp2 for species specific variance estimates
vcv4[, c("T1", "T2") := list(Sp1, Sp2)]
# hold ye horses with removing On. Might need it yet.
# vcv4[, On := NULL]

# Jarrod Hadfield's help here! (http://jarrod.bio.ed.ac.uk/)
vcv4 <- subset(vcv4,Sp1 != Sp2)
# make the relevant columns factors
as_factor <- names(vcv4)[-c(1,2)]
for(names in as_factor){
  set(x = vcv4, j = names, value = as.factor(vcv4[[names]]))
}

# tree manipulation
tree.4 <- tree.3
# remove tips that are not in both Sp1 and Sp2 (can Jarrod explain this please?)
tree.5 <- ape::drop.tip(tree.4, tree.4$tip.label[!tree.4$tip.label %in% union(vcv4$Sp1, vcv4$Sp2)])
tree.5$node.label <- NULL
# Distance!
D<-ape::cophenetic.phylo(x = tree.5)
# cool little mapply function here
vcv4$dist2<-mapply(as.character(vcv4$Sp1), as.character(vcv4$Sp2), FUN=function(x,y){D[x,y]})

# pull out the inverse matrix... 
Ainv1<-inverseA(tree.5, scale=FALSE)$Ainv
# make factor levels the same
vcv4$T1<-factor(vcv4$T1, levels=levels(vcv4$T2))
vcv4$Sp1<-factor(vcv4$Sp1, levels=levels(vcv4$Sp2))

load(file = "../data/Backups/vcv4.RData")

## NOT USED IN THE MANUSCRIPT (but here for completeness.) ##

prior.h1 <-list(R=list(V=diag(1), fix=1), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*100),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*100)))

# gonna be a long run
mcmc.h1 <- MCMCglmm(y ~ dist2 + Annual_Perennial, 
                    family = "threshold",
                    random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
                    ginverse = list(Sp1 = Ainv1,
                                    Sp2 = Ainv1),
                    data= vcv4,
                    prior=prior.h1,
                    pr=TRUE, 
                    nitt = 13000*100,
                    thin = 10*100,
                    burnin = 3000*100)

# 6110 pairwise zero combinations
# 558 pairwise one combinations

##### Part 2: Add in taxon co-occurrence #####

# this model and the ploidy one below are the two models described in the manuscript

taco <- fread("../data/Trait_databases/results20180326172536.csv")
taco$prefTaxonName <- gsub(x = taco$prefTaxonName, replacement = "", pattern = " [[:upper:]].*|\\(.*|=.*")
# remove hybrids
taco2 <- taco[!grepl(pattern = " x | x$", x = taco$prefTaxonName),]
# this destroys half of the hybrid name but tha's okay because we dont care about hybrids
taco2$key <- gsub(x = taco2$prefTaxonName, pattern = " $", replacement = "")
# probably would be a good time to do this here. 
#updat_sp <- TPL(splist = taco2$key)
load(file = "../data/Backups/updat_sp.RData")

taco3 <- taco2[, .(Species = setDT(updat_sp)[, .(Species = paste(New.Genus, New.Species))]$Species,
                   Hectads = dataValue)]
# hectads are now a list
taco3[, Hectads := strsplit(Hectads, split = ",")]
# merge duplicates
taco4 <- unique(taco3[, .(unlist(Hectads)), , by = "Species"])[, .(paste(V1, collapse = ",")), by = "Species"]
# turn hectads back to list!
taco4 <- taco4[, Hectads := strsplit(V1, split = ",")][,-"V1"]

# update names in the hectad dataset

name_changes[New.Species == "Gentianella amarella", which = TRUE]
name_changes[New.Species == "Bromus racemosus", which = TRUE]

name_changes2 <- name_changes[-c(354, 740)]
# ta-da
taco5 <- taco4[name_changes2, on = .(Species)]

# make all combinations
# match hectads of combinations
# setdiff on the match

# New.Species <- Species
taco5 <- taco5[, .(Species = New.Species, Hectads)]

# make all of the combinations of 2 taxa in the taco5 dataset (lots)
combs <- data.table(Species = t(combn(taco5$Species, 2)))
# join to one column of the combinations, then the other 
fjoin <- taco5[combs, on = .(Species = Species.V1)]
sjoin <- taco5[fjoin, on = .(Species = Species.V2)]
# modify sjoin to calculate the number of hectads shared between each species. 
# This takes about three minutes to run on this dataset.
all_combs <- sjoin[, .(Species1 = Species, Species2 = i.Species,
                   Hectads_shared = pbapply::pblapply(1:dim(sjoin)[1], 
                                                      function(x) length(intersect(sjoin$Hectads[[x]], sjoin$i.Hectads[[x]]))))]
# make the 'On' column
### 17.10.19 something is causing NA values to appear in Hectads_shared... ###
all_combs[, On := sort_names(paste(Species1, Species2, sep = " "))]
# join to vcv4
vcv5 <- all_combs[,.(Hectads_shared, On)][vcv4, on = "On"]
# turn list to vector.
vcv5[, Hectads_shared := unlist(as.numeric(as.character(vcv5$Hectads_shared)))]
# vcv4 now updated to vcv5 with taxon-hectad co-occurrence. vcv5 for main model!

### MAIN MODEL ###

# make the levels of each of T1 and Sp1 the same as T2 and Sp2.
vcv5$T1<-factor(vcv5$T1, levels=levels(vcv5$T2))
vcv5$Sp1<-factor(vcv5$Sp1, levels=levels(vcv5$Sp2))
vcv5<-subset(vcv5, !is.na(Hectads_shared))

## add in genus size ##
# rowbind each of the species and genera, then find unique
gs <- unique(rbind(vcv5[,.(Species = Sp1, Genus1 = Genus1)], vcv5[,.(Species = Sp2, Genus1 = Genus2)]))
# add the genus size
gs[, Genus.Size := .N, by = .(Genus1)]

vcv5 <- vcv5[unique(gs[, .(Genus1, Genus.Size)]), on = .(Genus1)]

load(file = "../data/Backups/vcv5.RData")

### test run on ploidy dataset... ###
tree.vcv5 <- ape::drop.tip(tree.4, tree.4$tip.label[!tree.4$tip.label %in% union(vcv5$Sp1, vcv5$Sp2)])
tree.vcv5 # 1110 species
# pull out the inverse matrix
Ainv.vcv5<-inverseA(tree.vcv5, scale=FALSE)$Ainv

# basic stats
# of all pairwise combinations, how many produce hybrids? 7.7%
vcv5[duplicated(On)][y == 1, .(N=.N)]/vcv5[duplicated(On)][,.(N=.N)]
# how many multi species genera are there?
mult.sp<-unique(rbind(vcv5[,.(Genus = Genus1, Genus.Size, y)], vcv5[,.(Genus = Genus2, Genus.Size, y)]))[order(Genus.Size)]
# looks like 96 in this dataset.
mult.sp[, .(y = sum(y)), by = .(Genus)][y > 0]
# what are the top hybridising genera in this analysis?
tops <- unique(vcv5[y == 1 & duplicated(On),.(.N, Genus.Size), by = .(Genus1)][order(-N)])
# order by genus size
ggplot(tops, aes(x = Genus.Size, y = N))+geom_point() +ggrepel::geom_label_repel(aes(label=Genus1[c(T,NA, NA, NA)]))
tops[order(-Genus.Size)][1:20]
# and how many hybrid combinations are there?
unique(vcv5[,.(y, Sp1, Sp2, On)])[!duplicated(On)][y == 1]


# of the parental species, how many are perennial? 68%
table(unique(vcv5[y == 1, .(Sp1, AnnPer = gsub("(.){2}$", "", Annual_Perennial))])$AnnPer)
# numbers of plants forming hybrids or not broken down by parental species life history
vcv5[duplicated(On), .(n=.N), by = .(y, Annual_Perennial)][order(y)]
vcv5[Annual_Perennial == "a-p" & duplicated(On) & y == 1]
vcv5[Annual_Perennial == "a-a" & duplicated(On) & y == 1]

# are perennial only genera larger than annual only genera? Only 2 of the top twenty largest genera with unique life histories are annual.
# in this study... or of the species we looked at.
unique(vcv5[,.(Genus = Genus1, Genus.Size, Annual_Perennial)])[order(-Genus.Size)][, N :=.N, by = .(Genus)][N == 1][1:20]


prior.vcv5 <-list(R=list(V=diag(1), fix=1), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# this takes 3+ hours...
# XXX genera
# XXX species 
mcmc.vcv5 <- MCMCglmm(y ~ dist2 + Hectads_shared + Annual_Perennial + Genus.Size, 
                       family = "threshold",
                       random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
                       ginverse = list(Sp1 = Ainv.vcv5,
                                       Sp2 = Ainv.vcv5),
                       data= vcv5,
                       prior=prior.vcv5,
                       pr=TRUE, 
                       nitt = 13000*100,
                       thin = 10*100,
                       burnin = 3000*100)
summary(mcmc.vcv5)

# relative importance of coefficients
sddist2 <- sd(vcv5$dist2)
sdHectads_shared <- sd(vcv5$Hectads_shared)
sdGenus <- sd(vcv5$Genus.Size)

VCVglmm::posterior.sd(mcmc.vcv5)

# seems like genetic distance mode important [not sure about the annual-perennial categorical covariates]
write.csv(x = data.table(covar = rownames(summary(mcmc.vcv5)$solutions)[c(2,3,6)],
           post.mean = summary(mcmc.vcv5)$solutions[c(2,3,6)],
           sd = c(sddist2, sdHectads_shared, sdGenus),
           scaled.post.mean = summary(mcmc.vcv5)$solutions[c(2,3,6),1]*c(sddist2, sdHectads_shared, sdGenus)), 
          file = "../data/Model_outputs/mcmc.vcv.5.standardised.covars.csv")

write.csv(x = specify_decimal(summary(mcmc.vcv5)$solutions, 4),
          file = "../data/Model_outputs/mcmc.vsv5.summary.csv")

write.csv(x = cbind(c("Chi2", "df", "p"), specify_decimal(as.data.table(VCVglmm::Wald.test.auto(mcmc.vcv5)), 4)), 
          file = "../data/Model_outputs/mcmc.vsv5.wald.tests.csv")

# effect of annual and perennial
# phylogenetic components here
treevar.vcv5 <- VCVglmm::Bbar(phylo = tree.vcv5, type = "Species")*mcmc.vcv5$VCV[,"Sp1+Sp2"]

# keeping in mind, with Annual Perennial added...
a_a <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
                       summary(mcmc.vcv5)$solutions[2,1]*mean(vcv5$dist2) + # at mean distance
                       summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$dist2) + # at mean genus size
                       summary(mcmc.vcv5)$solutions[3,1]*mean(vcv5$Hectads_shared), # at mean hectads
                     sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
p_p <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
                        summary(mcmc.vcv5)$solutions[2,1]*mean(vcv5$dist2) + # at mean distance
                        summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$dist2) + # at mean genus size
                        summary(mcmc.vcv5)$solutions[5,1] + # perennial perennial effect
                        summary(mcmc.vcv5)$solutions[3,1]*mean(vcv5$Hectads_shared), # at mean hectads
                      sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
# 1.32
p_p/a_a

# hectad sharing for perennial plants.
lower <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
               summary(mcmc.vcv5)$solutions[2,1]*mean(vcv5$dist2) + # at mean distance
                 summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$dist2) + # at mean genus size
               summary(mcmc.vcv5)$solutions[5,1] + # perennial perennial effect
               summary(mcmc.vcv5)$solutions[3,1]*0, # at lowest hectad sharing
             sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
higher <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
               summary(mcmc.vcv5)$solutions[2,1]*mean(vcv5$dist2) + # at mean distance
                 summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$dist2) + # at mean genus size
               summary(mcmc.vcv5)$solutions[5,1] + # perennial perennial effect
               summary(mcmc.vcv5)$solutions[3,1]*3755, # at mean hectads
             sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
# 2.51 times
higher/lower

# phylogenetic variance explained/signal
p1 <- posterior.mode(treevar.vcv5/rowSums(cbind(treevar.vcv5, 
                                     2*mcmc.vcv5$VCV[,"T1+T2"],
                                     mcmc.vcv5$VCV[,"units"])))
hpd1 <- HPDinterval(treevar.vcv5/rowSums(cbind(treevar.vcv5, 
                                  2*mcmc.vcv5$VCV[,"T1+T2"],
                                  mcmc.vcv5$VCV[,"units"])))
# species variance
p2 <- posterior.mode(2*mcmc.vcv5$VCV[,"T1+T2"]/rowSums(cbind(treevar.vcv5,
                                                     2*mcmc.vcv5$VCV[,"T1+T2"],
                                                     mcmc.vcv5$VCV[,"units"])))
hpd2 <- HPDinterval(2*mcmc.vcv5$VCV[,"T1+T2"]/rowSums(cbind(treevar.vcv5,
                                                  2*mcmc.vcv5$VCV[,"T1+T2"],
                                                  mcmc.vcv5$VCV[,"units"])))

# heatmap of the joint distribution of hectad sharing and branch length.

x_test <- seq(0,0.3,0.01)
y_test <- seq(0, 3000, 100)

vec <- matrix(nrow = length(x_test), ncol = length(y_test))
for(i in 1:length(x_test)){
  for(j in 1:length(y_test)){
    vec[i,j] <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
                        summary(mcmc.vcv5)$solutions[3,1]*y_test[j]+ # at mean hectad sharing
                        summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$Genus.Size) + # at mean genus size
                        summary(mcmc.vcv5)$solutions[4,1] + #annual - perennial hybridisation
                        summary(mcmc.vcv5)$solutions[2,1]*x_test[i], # effect of branch length
                      sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5)))
    
    
  }
}
colnames(vec)<- y_test
rownames(vec)<- x_test

col_fun <- colorRampPalette(colors = c(cbPalette[1], cbPalette[2]))

pdf(file = "../figures/Supplementary/Correlation_genetic_geographic_mod_output_vcv5.pdf")
lattice::levelplot(t(vec), col.regions=col_fun(1000), 
                   ylab = "Branch length", xlab = "Hectad sharing",
                   scales = list(tck = c(1,0), x = list(rot=45)))

  
# the effect of branch length
mcmc.vcv5.plot <- function(){
  x <- seq(0,0.3,0.0001)
  y <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
               summary(mcmc.vcv5)$solutions[3,1]*mean(vcv5$Hectads_shared) + # at mean hectad sharing
               summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$Genus.Size) + # at mean genus size
               summary(mcmc.vcv5)$solutions[4,1] + #annual - perennial hybridisation
               summary(mcmc.vcv5)$solutions[2,1]*x, # effect of branch length
             sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
  yup <- pnorm(summary(mcmc.vcv5)$solutions[1,3] + # intercept
                 summary(mcmc.vcv5)$solutions[3,3]*mean(vcv5$Hectads_shared) + # at mean hectad sharing
                 summary(mcmc.vcv5)$solutions[6,3]*mean(vcv5$Genus.Size) + # at mean genus size
                 summary(mcmc.vcv5)$solutions[4,3] + #annual - perennial hybridisation
                 summary(mcmc.vcv5)$solutions[2,3]*x, # effect of branch length
               sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
  ylo <- pnorm(summary(mcmc.vcv5)$solutions[1,2] + # intercept
                 summary(mcmc.vcv5)$solutions[3,2]*mean(vcv5$Hectads_shared) + # at mean hectad sharing
                 summary(mcmc.vcv5)$solutions[6,2]*mean(vcv5$Genus.Size) + # at mean genus size
                 summary(mcmc.vcv5)$solutions[4,2] + #annual - perennial hybridisation
                 summary(mcmc.vcv5)$solutions[2,2]*x, # effect of branch length
               sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
  ggplot(data.table(x=x, y=y, yup=yup), aes(x = x)) + 
    geom_line(size = 2, aes(y = y)) +
    geom_line(size = 1, lty = 2, aes(y = yup)) +
    geom_line(size = 1, lty = 2, aes(y = ylo))+
    geom_vline(xintercept = vcv5[, .(mean(dist2)), by = .(Genus1)][, .(mean(V1))][1]$V1, lty =2, col = cbPalette[7], size = 2) +
    xlab(label = "Branch length distance between two parental species") +
    ylab(label = "Probability of forming a hybrid")+
    theme_bw()+
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.text.y =  element_text(size = 30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    xlim(0,0.3)
} 
mcmc.vcv5.plot()

ggsave(filename = "../figures/Supplementary/Distance_hectads_branch_model1.pdf", plot = mcmc.vcv5.plot(), width = 12, height = 9, device = "pdf")

# the effect of hectad sharing

x <- seq(0,3000,1)
y <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
             summary(mcmc.vcv5)$solutions[2,1]*mean(vcv5$dist2) + # mean branch length
             summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$Genus.Size) + # at mean genus size
             summary(mcmc.vcv5)$solutions[3,1]*x, # hectads shared
           sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny

c("Min" = min(y), "Max" = max(y), "Mean", mean(y))

mcmc.vcv5.plot2 <- function(){
  x <- seq(0,3000,1)
  y <- pnorm(summary(mcmc.vcv5)$solutions[1,1] + # intercept
               summary(mcmc.vcv5)$solutions[2,1]*mean(vcv5$dist2) + # mean branch length
               summary(mcmc.vcv5)$solutions[6,1]*mean(vcv5$Genus.Size) + # at mean genus size
               summary(mcmc.vcv5)$solutions[3,1]*x, # hectads shared
             sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
  yup <- pnorm(summary(mcmc.vcv5)$solutions[1,3] + # intercept
                 summary(mcmc.vcv5)$solutions[2,3]*mean(vcv5$dist2) + # mean branch length
                 summary(mcmc.vcv5)$solutions[6,3]*mean(vcv5$Genus.Size) + # at mean genus size
                 summary(mcmc.vcv5)$solutions[3,3]*x, # hectads shared
               sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
  ylo <- pnorm(summary(mcmc.vcv5)$solutions[1,2] + # intercept
                 summary(mcmc.vcv5)$solutions[2,2]*mean(vcv5$dist2) + # mean branch length
                 summary(mcmc.vcv5)$solutions[6,2]*mean(vcv5$Genus.Size) + # at mean genus size
                 summary(mcmc.vcv5)$solutions[3,2]*x, # hectads shared
               sd = sqrt(mean(mcmc.vcv5$VCV[,"units"] + 2*mcmc.vcv5$VCV[,"T1+T2"] + treevar.vcv5))) # accounting for phylogeny
  ggplot(data.table(x=x, y=y), aes(x = x, y = y)) + 
    geom_line(size = 2) +
    geom_line(size = 1, lty = 2, aes(y = yup)) +
    geom_line(size = 1, lty = 2, aes(y = ylo))+
    geom_vline(xintercept = mean(vcv5$Hectads_shared), lty =2, col = "red", size = 2) +
    xlab(label = expression(paste("Pairwise overlap in distribution", " km"^2))) +
    ylab(label = "Probability of forming a hybrid")+
    theme_bw()+
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.text.y =  element_text(size = 30))
}
mcmc.vcv5.plot2()
## not in supplementary materials ##
ggsave(filename = "../figures/Supplementary/Pairwise_overlap_model1.pdf", plot = mcmc.vcv5.plot2(), width = 12, height = 9, device = "pdf")

# extract the posterior modes, which are the greatest in value?
tree_modes <- VCVglmm::Solapply(mcmc.vcv5, mean)[Group == "Sp1"][order(-Grouped_Value)]
# can we subset the important nodes? 1053, 68 (top), 450, 225
tree_modes2 <- tree_modes[grepl("Node", Variable)][, Variable := substring(Variable, 9)]
tree_modes3 <- tree_modes2[Variable %in% c(1059, 75, 207, 457)]

pdf(file = "../figures/Supplementary/two_phyogenies.pdf", width = 16, height = 16)
# two clades. Many nodes will be subtrees of other trees.
par(mfrow = c(1,2), mai = c(0.05, 0.05, 0.05, 0.05))
ape::plot.phylo(treeio::tree_subset(node = 1059, tree = tree.vcv5), font = 4, edge.width = 3)
ape::plot.phylo(treeio::tree_subset(node = 207, tree = tree.vcv5), font = 4, edge.width = 3)
dev.off()

### PLOIDY BELOW ###

# now add into the mix the ploidy level data. 
# *be careful: note to self* after doing the intial add
# manually check against the hybrid flora to check the hybridising taxa ploidy levels are correct
# this can be done by filtering on y == 1 and unique(On)!
# added a few more species from the Kew C-value data base 21.11.19
chromosome_data6 <- fread("../data/Trait_databases/chromosome_data_updated181119.csv")

# update species names
# name changes 2 from above
# so 865 species in the phylogeny for which ploidy is known with certainty
chromosome_data6[name_changes2, on = .(Species)][!is.na(Ploidy)] %>%
  ggplot(aes(x = Ploidy)) + geom_histogram()
# merge with new names
chromosome_data6 <- chromosome_data6[name_changes2, on = .(Species)][!is.na(Ploidy)][, .(Species = New.Species, Ploidy, Chromosome_number, Base, N)]
# make all combinations again
combs2 <- data.table(Species = t(combn(chromosome_data6$Species, 2)))
# join to the chromosome data
fjoin2 <- chromosome_data6[combs2, on = .(Species = Species.V1)]
sjoin2 <- chromosome_data6[fjoin2, on = .(Species = Species.V2)]
# another function to determine if a hybrid/potential hybrid is cross-ploid or not
# only species are included where only one ploidy is known in the UK
# apart from those species that may have multiple ploidy levels, but are recorded in the hybrid flora
# as having been formed from a specific cross... does that make sense?
all_combs2 <- sjoin2[, .(Species1 = Species, Species2 = i.Species,
                       Cross_Ploid = pbapply::pblapply(1:dim(sjoin2)[1], 
                                                       function(x) ifelse(test = sjoin2$Ploidy[[x]] == sjoin2$i.Ploidy[[x]],
                                                                          yes = "Homoploid", 
                                                                          no = "Heteroploid")))]
# make the 'On' column
all_combs2[, On := sort_names(paste(Species1, Species2, sep = " "))]
# join to main data
vcv6 <- all_combs2[,.(Cross_Ploid, On)][vcv5, on = "On"]
# turn list to vector.
vcv6[, Cross_Ploid := unlist(as.factor(as.character(vcv6$Cross_Ploid)))]
vcv6 <- vcv6[Cross_Ploid != "NULL"]
# always check structure
str(vcv6)
# change these four columns to factor columns (perhaps this is the problem, as.factor, factor or as.character?)
#for_factor <- c("Sp1", "Sp2", "T1", "T2")
#for(col in for_factor){
#  set(vcv6, j=col, value=as.character(vcv6[[col]]))
#}
# make the levels of each of T1 and Sp1 the same as T2 and Sp2.
vcv6$T1<-factor(vcv6$T1, levels=levels(vcv6$T2))
vcv6$Sp1<-factor(vcv6$Sp1, levels=levels(vcv6$Sp2))
vcv6<-subset(vcv6, !is.na(Hectads_shared))
# relevel so homoploid is baseline?
vcv6$Cross_Ploid <- relevel(x = vcv6$Cross_Ploid, ref = "Homoploid")

### test run on ploidy dataset... ###
# from 568 species in the old run, to 684 species.
tree.ploid <- ape::drop.tip(tree.4, tree.4$tip.label[!tree.4$tip.label %in% union(vcv6$Sp1, vcv6$Sp2)])
tree.ploid
# pull out the inverse matrix
Ainv.ploid<-inverseA(tree.ploid, scale=FALSE)$Ainv

prior.h1 <-list(R=list(V=diag(1), fix=1), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# gonna be a long run (it was.)
# 165 genera, up from 141!
# 684 species, up from 568 (nearly half British Flora!)
mcmc.ploid <- MCMCglmm(y ~ dist2 + Cross_Ploid + Hectads_shared + Annual_Perennial + Genus.Size, 
                    family = "threshold",
                    random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
                    ginverse = list(Sp1 = Ainv.ploid,
                                    Sp2 = Ainv.ploid),
                    data= vcv6,
                    prior=prior.h1,
                    pr=TRUE, 
                    nitt = 13000*100,
                    thin = 10*100,
                    burnin = 3000*100)
summary(mcmc.ploid)

write.csv(x = specify_decimal(summary(mcmc.ploid)$solutions, 4),
          file = "../data/Model_outputs/mcmc.ploid.summary.csv")

write.csv(x = cbind(c("Chi2", "df", "p"), specify_decimal(as.data.table(VCVglmm::Wald.test.auto(mcmc.ploid)), 4)), 
          file = "../data/Model_outputs/mcmc.ploid.wald.tests.csv")

VCVglmm::MCMCfixplot(mod = mcmc.ploid)

# What is the difference between the same ploidy and cross ploidy levels?
# as tree wasnt scaled, we need to multiply the phylogenetic variance by the tree variance.
treevar.ploid <- VCVglmm::Bbar(phylo = tree.ploid, type = "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]

# keeping in mind, with Annual Perennial added...
same_ploidy <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
                             summary(mcmc.ploid)$solutions[2,1]*mean(vcv6$dist2) + # at mean distance
                             summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size)+ # mean genus size
                             summary(mcmc.ploid)$solutions[4,1]*mean(vcv6$Hectads_shared), # at mean hectads
                           sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + treevar.ploid))) # accounting for phylogeny
cross_ploidy <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
                       summary(mcmc.ploid)$solutions[2,1]*mean(vcv6$dist2) + # at mean distance
                        summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size)+ # mean genus size
                       summary(mcmc.ploid)$solutions[3,1] + # cross ploid effect
                       summary(mcmc.ploid)$solutions[4,1]*mean(vcv6$Hectads_shared), # at mean hectads
                     sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + treevar.ploid))) # accounting for phylogeny
# 1.308
same_ploidy/cross_ploidy

# what about the effect of annual perennial, at mean genetic distance
# at mean hectad sharing at same ploidy level
annual_effect <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept (annual-annual parents)
                       summary(mcmc.ploid)$solutions[2,1]*mean(vcv6$dist2) + # at mean distance
                       summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size)+ # mean genus size
                       summary(mcmc.ploid)$solutions[4,1]*mean(vcv6$Hectads_shared), # at mean hectads
                     sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + treevar.ploid))) # accounting for phylogeny
perennial_effect <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
                        summary(mcmc.ploid)$solutions[6,1]*mean(vcv6$dist2) + # perennial-perennial parents
                        summary(mcmc.ploid)$solutions[2,1]*mean(vcv6$dist2) + # at mean distance
                        summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size)+ # mean genus size
                        summary(mcmc.ploid)$solutions[4,1]*mean(vcv6$Hectads_shared), # at mean hectads
                      sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + treevar.ploid))) # accounting for phylogeny
# perennial-perennial effect is 1.02 times higher probability of creating successful hybrid
# than annual-annual effect. This is not significant.
perennial_effect/annual_effect


# visualise the fixed effects from raw data.
mcmc.ploidvis <- vcv6[, .(Prob = mean(y),
         SE.Prob = sd(y)/sqrt(.N)), by = .(Cross_Ploid)]
ggplot(mcmc.ploidvis, aes(y = Prob, x = Cross_Ploid))+
  geom_point()+
  geom_errorbar(aes(ymin = Prob - SE.Prob, ymax = Prob + SE.Prob), width = 0.3) +
  xlab(label = "Parental Ploidy difference") +
  ylab(label = "Probability of forming a hybrid")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y =  element_text(size = 20))
  
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# predicted model outputs
ploid_plot <- function(){
  x <- seq(0,3000,1)
  y <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
               summary(mcmc.ploid)$solutions[2,1]*mean(vcv6$dist2) + # mean of the genetic distance
               summary(mcmc.ploid)$solutions[5,1] + # annual-perennial
               summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size) + # mean genus size
               summary(mcmc.ploid)$solutions[4,1]*x, # for a range of hectads
             sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 
                              2*mcmc.ploid$VCV[,"T1+T2"] + 
                              Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # homoploid
  
  yup <-  pnorm(summary(mcmc.ploid)$solutions[1,3] + # intercept
                  summary(mcmc.ploid)$solutions[2,3]*mean(vcv6$dist2) + # mean of the genetic distance
                  summary(mcmc.ploid)$solutions[5,3] + # annual-perennial
                  summary(mcmc.ploid)$solutions[7,3]*mean(vcv6$Genus.Size) + # mean genus size
                  summary(mcmc.ploid)$solutions[4,3]*x, # for a range of hectads
                sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 
                                 2*mcmc.ploid$VCV[,"T1+T2"] + 
                                 Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # homoploid
  ylo <- pnorm(summary(mcmc.ploid)$solutions[1,2] + # intercept
                 summary(mcmc.ploid)$solutions[2,2]*mean(vcv6$dist2) + # mean of the genetic distance
                 summary(mcmc.ploid)$solutions[7,2]*mean(vcv6$Genus.Size)+ # mean genus size
                 summary(mcmc.ploid)$solutions[5,2] + # annual-perennial
                 summary(mcmc.ploid)$solutions[4,2]*x, # for a range of hectads
               sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 
                                2*mcmc.ploid$VCV[,"T1+T2"] + 
                                Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # homoploid
  
  y2<- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
               summary(mcmc.ploid)$solutions[2,1]*mean(vcv6$dist2) + # mean of the genetic distance
               summary(mcmc.ploid)$solutions[3,1] + # effect of heteroploidy
               summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size) + # mean genus size
               summary(mcmc.ploid)$solutions[5,1] + # annual perennial
               summary(mcmc.ploid)$solutions[4,1]*x, # for a range of hectads
             sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # heteroploid
  y2up<- pnorm(summary(mcmc.ploid)$solutions[1,3] + # intercept
               summary(mcmc.ploid)$solutions[2,3]*mean(vcv6$dist2) + # mean of the genetic distance
                 summary(mcmc.ploid)$solutions[7,3]*mean(vcv6$Genus.Size) + # mean genus size
               summary(mcmc.ploid)$solutions[3,3] + # effect of heteroploidy
                 summary(mcmc.ploid)$solutions[5,3] + # annual-perennial
               summary(mcmc.ploid)$solutions[4,3]*x, # for a range of hectads
             sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # heteroploid
  y2lo<- pnorm(summary(mcmc.ploid)$solutions[1,2] + # intercept
               summary(mcmc.ploid)$solutions[2,2]*mean(vcv6$dist2) + # mean of the genetic distance
                 summary(mcmc.ploid)$solutions[7,2]*mean(vcv6$Genus.Size) + # mean genus size
               summary(mcmc.ploid)$solutions[3,2] + # effect of heteroploidy
                 summary(mcmc.ploid)$solutions[5,2] + # annual-perennial
               summary(mcmc.ploid)$solutions[4,2]*x, # for a range of hectads
             sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # heteroploid
  
  data <- melt(data.table(x=x, Homoploid=y, Heteroploid=y2), measure.vars = c("Homoploid", "Heteroploid"), 
               variable.name = "Parental Ploidy")
  data$lCI <- c(ylo, y2lo)
  data$uCI <- c(yup, y2up)
  
  ggplot(data, aes(x = x, y = value, colour = `Parental Ploidy`)) + 
    geom_line(size = 2) +
    geom_line(aes(y = lCI), lty = 2, size =1.2) +
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = `Parental Ploidy`), alpha = 0.15) +
    geom_line(aes(y = uCI), lty = 2, size =1.2) +
    geom_vline(xintercept = mean(vcv6$Hectads_shared), lty =2, col = cbPalette[7], size = 2) +
    xlab(label = expression(paste("Pairwise overlap in distribution", " km"^2))) +
    ylab(label = "Probability of forming a hybrid")+
    theme_bw()+
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.text.y =  element_text(size = 30),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_manual(values = c(cbPalette[1], cbPalette[2]))+
    scale_color_manual(values = c(cbPalette[1], cbPalette[2]))
}
ploid_plot()
ggsave(filename = "../figures/Supplementary/Cross_ploid_hectads.pdf", plot = ploid_plot(), width = 12, height = 9, device = "pdf")

ploid_plot_2 <- function(){
  x <- seq(0, 0.3,0.01)
  y <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
               summary(mcmc.ploid)$solutions[2,1]*x + # genetic distance
               summary(mcmc.ploid)$solutions[5,1] + # annual-perennial
               summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size) + # mean genus size
               summary(mcmc.ploid)$solutions[4,1]*mean(vcv6$Hectads_shared), # mean of hectads
             sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 
                              2*mcmc.ploid$VCV[,"T1+T2"] + 
                              Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # homoploid
  
  yup <-  pnorm(summary(mcmc.ploid)$solutions[1,3] + # intercept
                  summary(mcmc.ploid)$solutions[2,3]*x + #  genetic distance
                  summary(mcmc.ploid)$solutions[5,3] + # annual-perennial
                  summary(mcmc.ploid)$solutions[7,3]*mean(vcv6$Genus.Size) + # mean genus size
                  summary(mcmc.ploid)$solutions[4,3]*mean(vcv6$Hectads_shared), # for a range of hectads
                sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 
                                 2*mcmc.ploid$VCV[,"T1+T2"] + 
                                 Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # homoploid
  ylo <- pnorm(summary(mcmc.ploid)$solutions[1,2] + # intercept
                 summary(mcmc.ploid)$solutions[2,2]*x + #  genetic distance
                 summary(mcmc.ploid)$solutions[7,2]*mean(vcv6$Genus.Size)+ # mean genus size
                 summary(mcmc.ploid)$solutions[5,2] + # annual-perennial
                 summary(mcmc.ploid)$solutions[4,2]*mean(vcv6$Hectads_shared), # for a range of hectads
               sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 
                                2*mcmc.ploid$VCV[,"T1+T2"] + 
                                Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # homoploid
  
  y2<- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
               summary(mcmc.ploid)$solutions[2,1]*x + # mean of the genetic distance
               summary(mcmc.ploid)$solutions[3,1] + # effect of heteroploidy
               summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size) + # mean genus size
               summary(mcmc.ploid)$solutions[5,1] + # annual perennial
               summary(mcmc.ploid)$solutions[4,1]*mean(vcv6$Hectads_shared), # for a range of hectads
             sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # heteroploid
  y2up<- pnorm(summary(mcmc.ploid)$solutions[1,3] + # intercept
                 summary(mcmc.ploid)$solutions[2,3]*x + # mean of the genetic distance
                 summary(mcmc.ploid)$solutions[7,3]*mean(vcv6$Genus.Size) + # mean genus size
                 summary(mcmc.ploid)$solutions[3,3] + # effect of heteroploidy
                 summary(mcmc.ploid)$solutions[5,3] + # annual-perennial
                 summary(mcmc.ploid)$solutions[4,3]*mean(vcv6$Hectads_shared), # for a range of hectads
               sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # heteroploid
  y2lo<- pnorm(summary(mcmc.ploid)$solutions[1,2] + # intercept
                 summary(mcmc.ploid)$solutions[2,2]*x + # mean of the genetic distance
                 summary(mcmc.ploid)$solutions[7,2]*mean(vcv6$Genus.Size) + # mean genus size
                 summary(mcmc.ploid)$solutions[3,2] + # effect of heteroploidy
                 summary(mcmc.ploid)$solutions[5,2] + # annual-perennial
                 summary(mcmc.ploid)$solutions[4,2]*mean(vcv6$Hectads_shared), # for a range of hectads
               sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + Bbar(tree.ploid, "Species")*mcmc.ploid$VCV[,"Sp1+Sp2"]))) # heteroploid
  
  data <- melt(data.table(x=x, Homoploid=y, Heteroploid=y2), measure.vars = c("Homoploid", "Heteroploid"), 
               variable.name = "Parental Ploidy")
  data$lCI <- c(ylo, y2lo)
  data$uCI <- c(yup, y2up)
  
  ggplot(data, aes(x = x, y = value, colour = `Parental Ploidy`)) + 
    geom_line(size = 2) +
    geom_line(aes(y = lCI), lty = 2, size =1.2) +
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = `Parental Ploidy`), alpha = 0.15) +
    geom_line(aes(y = uCI), lty = 2, size =1.2) +
    geom_vline(xintercept = vcv6[, .(mean(dist2)), by = .(Genus1)][, .(mean(V1))][1]$V1, lty =2, col = cbPalette[7], size = 2) +
    xlab(label = expression(paste("Branch length between a pair of species"))) +
    ylab(label = "Probability of forming a hybrid")+
    theme_bw()+
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.text.y =  element_text(size = 30),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_manual(values = c(cbPalette[1], cbPalette[2]))+
    scale_color_manual(values = c(cbPalette[1], cbPalette[2]))
}
ploid_plot_2()
ggsave(filename = "../figures/Supplementary/Cross_ploid_branch.pdf", plot = ploid_plot_2(), width = 12, height = 9, device = "pdf")

# phylogenetic variance ~58%
p3 <- posterior.mode(treevar.ploid/rowSums(cbind(treevar.ploid, 
                                           2*mcmc.ploid$VCV[,"T1+T2"],
                                           mcmc.ploid$VCV[,"units"])))
hpd3 <- HPDinterval(treevar.ploid/rowSums(cbind(treevar.ploid, 
                                           2*mcmc.ploid$VCV[,"T1+T2"],
                                           mcmc.ploid$VCV[,"units"])))
# species variance ~35%
p4 <- posterior.mode(2*mcmc.ploid$VCV[,"T1+T2"]/rowSums(cbind(treevar.ploid,
                                                   2*mcmc.ploid$VCV[,"T1+T2"],
                                                   mcmc.ploid$VCV[,"units"])))
hpd4 <- HPDinterval(mcmc.ploid$VCV[,"T1+T2"]/rowSums(cbind(treevar.ploid,
                                                   mcmc.ploid$VCV[,"T1+T2"],
                                                   mcmc.ploid$VCV[,"units"])))

# no different visually between the ploidy levels, so keep one figure as above.
vec2 <- matrix(nrow = length(x_test), ncol = length(y_test))
for(i in 1:length(x_test)){
  for(j in 1:length(y_test)){
    vec2[i,j] <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
                        summary(mcmc.ploid)$solutions[4,1]*y_test[j]+ # at mean hectad sharing
                        summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size) + # at mean genus size
                        summary(mcmc.ploid)$solutions[5,1] + #annual - perennial hybridisation
                        summary(mcmc.ploid)$solutions[2,1]*x_test[i], # effect of branch length
                      sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + treevar.ploid)))
    
    
  }
}
colnames(vec2)<- y_test
rownames(vec2)<- x_test

vec3 <- matrix(nrow = length(x_test), ncol = length(y_test))
for(i in 1:length(x_test)){
  for(j in 1:length(y_test)){
    vec3[i,j] <- pnorm(summary(mcmc.ploid)$solutions[1,1] + # intercept
                         summary(mcmc.ploid)$solutions[3,1] +
                         summary(mcmc.ploid)$solutions[4,1]*y_test[j]+ # at mean hectad sharing
                         summary(mcmc.ploid)$solutions[7,1]*mean(vcv6$Genus.Size) + # at mean genus size
                         summary(mcmc.ploid)$solutions[5,1] + #annual - perennial hybridisation
                         summary(mcmc.ploid)$solutions[2,1]*x_test[i], # effect of branch length
                       sd = sqrt(mean(mcmc.ploid$VCV[,"units"] + 2*mcmc.ploid$VCV[,"T1+T2"] + treevar.ploid)))
    
    
  }
}
colnames(vec3)<- y_test
rownames(vec3)<- x_test

grid.arrange(levelplot(t(vec2), col.regions=col_fun(1000), 
                                           main = "Homoploid hybrids",
                                           ylab = "Branch length", xlab = "Hectad sharing",
                                           scales = list(tck = c(1,0), x = list(rot=45))),
                        levelplot(t(vec3), col.regions=col_fun(1000), 
                                           main = "Heteroploid hybrids",
                                           ylab = "Branch length", xlab = "Hectad sharing",
                                           scales = list(tck = c(1,0), x = list(rot=45))),
                        ncol=2
                        
)



# table of variance components

write.csv(x = data.table(`Variance Component` = c("Model 1 Phylogenetic Variance", "Model 1 Species Variance",
                                    "Model 2 Phylogenetic Variance", "Model 2 Species Variance"),
           `Posterior Mode` = c(specify_decimal(p1, 4), specify_decimal(p2, 4), specify_decimal(p3, 4), specify_decimal(p4, 4)),
           `Lower Credible Interval` = c(specify_decimal(hpd1, 4)[1], specify_decimal(hpd2, 4)[1], specify_decimal(hpd3, 4)[1], specify_decimal(hpd4, 4)[1]),
           `Upper Credible Interval` = c(specify_decimal(hpd1, 4)[2], specify_decimal(hpd2, 4)[2], specify_decimal(hpd3, 4)[2], specify_decimal(hpd4, 4)[2])), 
          file = "../data/Model_outputs/variance_components.csv")
