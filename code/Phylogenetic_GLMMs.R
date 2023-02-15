# The genetic landscape of hybridisation in the British flora #
# Max R. Brown, Peter M. Hollingsworth, Laura Forrest, Michelle Hart, Laura Jones, Col Ford, Tim Rich, Natasha de Vere, Alex D. Twyford #

# All code Max Brown #
# Last updated: 05.02.23 #

# Libraries needed #
# install_github("Euphrasiologist/VCVglmm")
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

# useful little function that splits a character vector and then sorts it
# alphabetically after the split.
# the pbapply can be replaced with apply.
sort_names <- function(x, split = " ", collapse = "") {
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

# See the `Barcoding_Phylogeny_Dating` directory for more information.
dated_tree <- read.tree(
  file = "./Barcoding_Phylogeny_Dating/DNA_Barcoding.dated.treefile"
)

dated_tree$tip.label <- gsub(
  pattern = "_",
  replacement = " ", dated_tree$tip.label
)
# dated_tree has 1406 tips

# remove node labels from the tree
# otherwise the inverseA function fails.
dated_tree$node.label <- NULL

# multi membership model
# y <- vector of 0's and 1's
# as to whether a species pair gave
# rise to a successful hybrid or not
# Sp1 <- vector of species 1
# Sp2 <- vector of species 2
# T1 <- copy of Sp1
# T2 <- copy of Sp2
# dist <- variance covariance matrix of phylogeny
# cooc <- co-occurrence of parental taxa

# Sp1, Sp2 and the VCV
# just to generate the pairwise comparisons, vcv values discarded later.
vd <- vcv(dated_tree)

# create the data table object
# for the dated tree
vcv1d <- data.table(as.data.frame(as.matrix(vd)), keep.rownames = TRUE)
vcv2d <- melt.data.table(vcv1d)
setnames(
  x = vcv2d,
  old = c("rn", "variable", "value"), new = c("Sp1", "Sp2", "date_div")
)

# create column which we can sort on later (expensive call)
vcv2d[, On := sort_names(paste(Sp1, Sp2, sep = " "))]
# remove dist, updated dist2 below
vcv2d <- vcv2d[, -"date_div"]

# genus columns
vcv2d[, c("Genus1", "Genus2") := list(
  gsub(" .*", "", Sp1),
  gsub(" .*", "", Sp2)
)]

# these are all of the hybrids.
hybriddata <- fread(
  # a digital cleaned rendition of the `Hybrid Flora of the British Isles`
  "./data/Hybrid_flora_of_the_British_Isles/hybriddata_Stace4_removed.csv"
)
merg <- hybriddata[
  ,
  On := sort_names(
    paste(`Parent A`, `Parent B`, sep = " ")
  )
][, .(
  On = On,
  y = 1
)]
# merge the data sets
vcv3d <- merg[vcv2d, on = "On"]

# add in the annual perennial data.
# need to download from the bsbi website and standardise names.
# From https://database.bsbi.org/checklists.php
annper <- fread("./data/Trait_databases/annual_perennial_britf_fl.csv")
name_changes <- fread("./data/Hybrid_flora_of_the_British_Isles/name_changes.csv")
name_changes <- name_changes[-1, .(Species = V2, New.Species = V3)]
# standardise the names
# I dont know when this will be deprecated
annpernames <- Taxonstand::TPL(annper$prefTaxonName, silent = FALSE)
# set as a data table
setDT(annpernames)
fwrite(annpernames, file = "./data/Trait_databases/annpernames.csv")
# bind the updated names to the annual perennial columns
annper2 <- cbind(annpernames[, .(Species = paste(New.Genus, New.Species, sep = " "))], annper[, .(AnnPer = dataValue)])

# NA values will be introduced. Unfortunately add manually.
annper2.1 <- unique(annper2[name_changes, on = .(Species)][, .(Species = New.Species, AnnPer)])
# to see NA's
annper2.1[is.na(AnnPer)]
annual.perennial <- c("p", "p", "a", "a", "a", "p", "a", "a", "p", "a", "p", "p", "p", "p", "p", "p", "p", "p", rep("p", 9), "p", "a")
annper2.1[is.na(AnnPer), AnnPer := annual.perennial]


# merge with the names in the phylogeny
annper3 <- unique(annper2.1)[unique(vcv3[, .(Species = Sp1)]), on = "Species"]

annper3d <- unique(annper2.1)[unique(vcv3d[, .(Species = Sp1)]), on = "Species"]

# why are there 1409 in the new annper, not 1406...
# remove the table from below
annper3[Species %in% annper3[duplicated(Species)]$Species][c(1, 3, 5)]
setkey(annper3, Species, AnnPer)

annper3d[Species %in% annper3d[duplicated(Species)]$Species][c(1, 3, 5)]
setkey(annper3d, Species, AnnPer)

# and now subset.
annper4 <- annper3[!.(annper3[Species %in% annper3[duplicated(Species)]$Species][c(1, 3, 5)])]

annper4d <- annper3d[!.(annper3d[Species %in% annper3d[duplicated(Species)]$Species][c(1, 3, 5)])]

fwrite(annper4, "./data/Trait_databases/annual_perennial_britf_fl_updated.csv")
fwrite(annper4d, "./data/Trait_databases/annual_perennial_britf_fl_updated_d.csv")

# merge with vcv3
vcv3.1 <- vcv3[annper4, on = .(Sp1 = Species)]
vcv3.2 <- vcv3.1[annper4, on = .(Sp2 = Species)]

vcv3.1d <- vcv3d[annper4d, on = .(Sp1 = Species)]
vcv3.2d <- vcv3.1d[annper4d, on = .(Sp2 = Species)]

# make biennial perennial
vcv3.2[, AnnPer := ifelse(AnnPer == "b", yes = "p", no = AnnPer)]
vcv3.2[, i.AnnPer := ifelse(i.AnnPer == "b", yes = "p", no = i.AnnPer)]

vcv3.2d[, AnnPer := ifelse(AnnPer == "b", yes = "p", no = AnnPer)]
vcv3.2d[, i.AnnPer := ifelse(i.AnnPer == "b", yes = "p", no = i.AnnPer)]

# make annual perennial parental combination
vcv3.2[, Annual_Perennial := paste(AnnPer, i.AnnPer, sep = "-")]

vcv3.2d[, Annual_Perennial := paste(AnnPer, i.AnnPer, sep = "-")]
# make p-a equal to a-p
vcv3.2[, Annual_Perennial := ifelse(Annual_Perennial == "p-a", yes = "a-p", no = Annual_Perennial)]

vcv3.2d[, Annual_Perennial := ifelse(Annual_Perennial == "p-a", yes = "a-p", no = Annual_Perennial)]
# remove AnnPer and i.AnnPer
vcv3.2 <- vcv3.2[, -c("AnnPer", "i.AnnPer")]

vcv3.2d <- vcv3.2d[, -c("AnnPer", "i.AnnPer")]
# filter for congeneric
vcv4 <- vcv3.2[Genus1 == Genus2]

vcv4d <- vcv3.2d[Genus1 == Genus2]
# make zero values for y
vcv4[, y := ifelse(test = is.na(y), yes = 0, no = y)]

vcv4d[, y := ifelse(test = is.na(y), yes = 0, no = y)]
# duplicate Sp1, Sp2 for species specific variance estimates
vcv4[, c("T1", "T2") := list(Sp1, Sp2)]

vcv4d[, c("T1", "T2") := list(Sp1, Sp2)]
# hold ye horses with removing On. Might need it yet.
# vcv4[, On := NULL]

# Jarrod Hadfield's help here! (http://jarrod.bio.ed.ac.uk/)
vcv4 <- subset(vcv4, Sp1 != Sp2)

vcv4d <- subset(vcv4d, Sp1 != Sp2)
# make the relevant columns factors
as_factor <- names(vcv4)[-c(1, 2)]
for (names in as_factor) {
  set(x = vcv4, j = names, value = as.factor(vcv4[[names]]))
}

as_factord <- names(vcv4d)[-c(1, 2)]
for (names in as_factord) {
  set(x = vcv4d, j = names, value = as.factor(vcv4d[[names]]))
}

# tree manipulation
tree.4 <- tree.3
tree.4d <- dated_tree
# remove tips that are not in both Sp1 and Sp2 (can Jarrod explain this please?)
tree.5 <- ape::drop.tip(tree.4, tree.4$tip.label[
  !tree.4$tip.label %in% union(vcv4$Sp1, vcv4$Sp2)
])
tree.5$node.label <- NULL

tree.5d <- ape::drop.tip(tree.4d, tree.4d$tip.label[!tree.4d$tip.label %in% union(vcv4d$Sp1, vcv4d$Sp2)])
tree.5d$node.label <- NULL
# Distance!
D <- ape::cophenetic.phylo(x = tree.5)
Dd <- ape::cophenetic.phylo(x = tree.5d)
# cool little mapply function here
vcv4$dist2 <- mapply(as.character(vcv4$Sp1), as.character(vcv4$Sp2), FUN = function(x, y) {
  D[x, y]
})

vcv4d$dist2 <- mapply(as.character(vcv4d$Sp1), as.character(vcv4d$Sp2), FUN = function(x, y) {
  Dd[x, y]
})

# pull out the inverse matrix...
Ainv1 <- inverseA(tree.5, scale = FALSE)$Ainv

Ainv1d <- inverseA(tree.5d, scale = FALSE)$Ainv
# make factor levels the same
vcv4$T1 <- factor(vcv4$T1, levels = levels(vcv4$T2))
vcv4$Sp1 <- factor(vcv4$Sp1, levels = levels(vcv4$Sp2))

vcv4d$T1 <- factor(vcv4d$T1, levels = levels(vcv4d$T2))
vcv4d$Sp1 <- factor(vcv4d$Sp1, levels = levels(vcv4d$Sp2))

# load(file = "./data/Backups/vcv4.RData")

## NOT USED IN THE MANUSCRIPT (but here for completeness.) ##

prior.h1 <- list(
  R = list(V = diag(1), fix = 1),
  G = list(
    G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1) * 100),
    G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1) * 100)
  )
)

# gonna be a long run
# replace vcv4 with vcv4d
mcmc.h1 <- MCMCglmm(y ~ dist2 + Annual_Perennial,
  family = "threshold",
  random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
  ginverse = list(
    Sp1 = Ainv1d,
    Sp2 = Ainv1d
  ),
  data = vcv4d,
  prior = prior.h1,
  pr = TRUE,
  nitt = 13000 * 100,
  thin = 10 * 100,
  burnin = 3000 * 100
)

# 6110 pairwise zero combinations
# 558 pairwise one combinations

##### Part 2: Add in taxon co-occurrence #####

# this model and the ploidy one below are the two models described in the manuscript

taco <- fread("./data/Trait_databases/results20180326172536.csv")
taco$prefTaxonName <- gsub(x = taco$prefTaxonName, replacement = "", pattern = " [[:upper:]].*|\\(.*|=.*")
# remove hybrids
taco2 <- taco[!grepl(pattern = " x | x$", x = taco$prefTaxonName), ]
# this destroys half of the hybrid name but tha's okay because we dont care about hybrids
taco2$key <- gsub(x = taco2$prefTaxonName, pattern = " $", replacement = "")
# probably would be a good time to do this here.
# updat_sp <- TPL(splist = taco2$key)
load(file = "./data/Backups/updat_sp.RData")

taco3 <- taco2[, .(
  Species = setDT(updat_sp)[, .(Species = paste(New.Genus, New.Species))]$Species,
  Hectads = dataValue
)]
# hectads are now a list
taco3[, Hectads := strsplit(Hectads, split = ",")]
# merge duplicates
taco4 <- unique(taco3[, .(unlist(Hectads)), , by = "Species"])[, .(paste(V1, collapse = ",")), by = "Species"]
# turn hectads back to list!
taco4 <- taco4[, Hectads := strsplit(V1, split = ",")][, -"V1"]

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
all_combs <- sjoin[, .(
  Species1 = Species, Species2 = i.Species,
  Hectads_shared = pbapply::pblapply(
    1:dim(sjoin)[1],
    function(x) length(intersect(sjoin$Hectads[[x]], sjoin$i.Hectads[[x]]))
  )
)]
# make the 'On' column
### 17.10.19 something is causing NA values to appear in Hectads_shared... ###
all_combs[, On := sort_names(paste(Species1, Species2, sep = " "))]
# join to vcv4
vcv5 <- all_combs[, .(Hectads_shared, On)][vcv4, on = "On"]

vcv5d <- all_combs[, .(Hectads_shared, On)][vcv4d, on = "On"]
# turn list to vector.
vcv5[, Hectads_shared := unlist(as.numeric(as.character(vcv5$Hectads_shared)))]

vcv5d[, Hectads_shared := unlist(as.numeric(as.character(vcv5d$Hectads_shared)))]
# vcv4 now updated to vcv5 with taxon-hectad co-occurrence. vcv5 for main model!

### MAIN MODEL ###

# make the levels of each of T1 and Sp1 the same as T2 and Sp2.
vcv5$T1 <- factor(vcv5$T1, levels = levels(vcv5$T2))
vcv5$Sp1 <- factor(vcv5$Sp1, levels = levels(vcv5$Sp2))
vcv5 <- subset(vcv5, !is.na(Hectads_shared))

vcv5d$T1 <- factor(vcv5d$T1, levels = levels(vcv5d$T2))
vcv5d$Sp1 <- factor(vcv5d$Sp1, levels = levels(vcv5d$Sp2))
vcv5d <- subset(vcv5d, !is.na(Hectads_shared))

## add in genus size ##
# rowbind each of the species and genera, then find unique
gs <- unique(rbind(vcv5[, .(Species = Sp1, Genus1 = Genus1)], vcv5[, .(Species = Sp2, Genus1 = Genus2)]))

gsd <- unique(rbind(vcv5d[, .(Species = Sp1, Genus1 = Genus1)], vcv5d[, .(Species = Sp2, Genus1 = Genus2)]))
# add the genus size
gs[, Genus.Size := .N, by = .(Genus1)]
gsd[, Genus.Size := .N, by = .(Genus1)]

vcv5 <- vcv5[unique(gs[, .(Genus1, Genus.Size)]), on = .(Genus1)]

vcv5d <- vcv5d[unique(gsd[, .(Genus1, Genus.Size)]), on = .(Genus1)]
# save(vcv5d, file = "./data/Backups/vcv5d.RData")
load(file = "./data/Backups/vcv5d.RData")

### test run on ploidy dataset... ###
tree.vcv5 <- ape::drop.tip(tree.4, tree.4$tip.label[!tree.4$tip.label %in% union(vcv5$Sp1, vcv5$Sp2)])
tree.vcv5 # 1110 species
# pull out the inverse matrix
Ainv.vcv5 <- inverseA(tree.vcv5, scale = FALSE)$Ainv

tree.vcv5d <- ape::drop.tip(tree.4d, tree.4d$tip.label[!tree.4d$tip.label %in% union(vcv5d$Sp1, vcv5d$Sp2)])
tree.vcv5d # 1110 species
# pull out the inverse matrix
Ainv.vcv5d <- inverseA(tree.vcv5d, scale = FALSE)$Ainv
# save(Ainv.vcv5d, file = "./data/Backups/Ainv.vcv5d.RData")

# basic stats
# of all pairwise combinations, how many produce hybrids? 7.7%
vcv5d[duplicated(On)][y == 1, .(N = .N)] / vcv5d[duplicated(On)][, .(N = .N)]
# how many multi species genera are there?
mult.sp <- unique(rbind(vcv5d[, .(Genus = Genus1, Genus.Size, y)], vcv5d[, .(Genus = Genus2, Genus.Size, y)]))[order(Genus.Size)]

# looks like 96 in this dataset.
supp_perc_hybrids <- rbind(
  mult.sp[, .(y = sum(y)), by = .(Genus)][y > 0],
  mult.sp[, .(y = sum(y)), by = .(Genus)][y == 0]
)[, .(Genus, Hybrid_genus = y)]

fwrite(supp_perc_hybrids, file = "./data/supp_perc_hybrids.csv")

# what are the top hybridising genera in this analysis?
tops <- unique(vcv5[y == 1 & duplicated(On), .(.N, Genus.Size), by = .(Genus1)][order(-N)])
# order by genus size
ggplot(tops, aes(x = Genus.Size, y = N)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = Genus1[c(T, NA, NA, NA)]))
tops[order(-Genus.Size)][1:20]
# and how many hybrid combinations are there?
unique(vcv5[, .(y, Sp1, Sp2, On)])[!duplicated(On)][y == 1]


# of the parental species, how many are perennial? 68%
table(unique(vcv5d[y == 1, .(Sp1, AnnPer = gsub("(.){2}$", "", Annual_Perennial))])$AnnPer)
# numbers of plants forming hybrids or not broken down by parental species life history
vcv5[duplicated(On), .(n = .N), by = .(y, Annual_Perennial)][order(y)]
vcv5[Annual_Perennial == "a-p" & duplicated(On) & y == 1]
vcv5[Annual_Perennial == "a-a" & duplicated(On) & y == 1]

# are perennial only genera larger than annual only genera? Only 2 of the top twenty largest genera with unique life histories are annual.
# in this study... or of the species we looked at.
unique(vcv5[, .(Genus = Genus1, Genus.Size, Annual_Perennial)])[order(-Genus.Size)][, N := .N, by = .(Genus)][N == 1][1:20]


prior.vcv5d <- list(
  R = list(V = diag(1), fix = 1),
  G = list(
    G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1) * 1000),
    G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1) * 1000)
  )
)

# this takes 3+ hours...
#
mcmc.vcv5d <- MCMCglmm(y ~ dist2 + Hectads_shared + Annual_Perennial + Genus.Size,
  family = "threshold",
  random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
  ginverse = list(
    Sp1 = Ainv.vcv5d,
    Sp2 = Ainv.vcv5d
  ),
  data = vcv5d,
  prior = prior.vcv5d,
  pr = TRUE,
  nitt = 13000 * 100,
  thin = 10 * 100,
  burnin = 3000 * 100
)
save(mcmc.vcv5d, file = "./data/Backups/mcmc.vcv5d.RData")
load("./data/Backups/mcmc.vcv5d.RData")
summary(mcmc.vcv5d)

# relative importance of coefficients
sddist2 <- sd(vcv5d$dist2)
sdHectads_shared <- sd(vcv5d$Hectads_shared)
sdGenus <- sd(vcv5d$Genus.Size)

post_sd_vcv5d <- VCVglmm::posterior.sd(mcmc.vcv5d)
fwrite(post_sd_vcv5d, "./data/Model_outputs/post_sd_vcv5d.csv")

# seems like genetic distance mode important [not sure about the annual-perennial categorical covariates]
write.csv(
  x = data.table(
    covar = rownames(summary(mcmc.vcv5d)$solutions)[c(2, 3, 6)],
    post.mean = summary(mcmc.vcv5d)$solutions[c(2, 3, 6)],
    sd = c(sddist2, sdHectads_shared, sdGenus),
    scaled.post.mean = summary(mcmc.vcv5d)$solutions[c(2, 3, 6), 1] * c(sddist2, sdHectads_shared, sdGenus)
  ),
  file = "./data/Model_outputs/mcmc.vcv.5d.standardised.covars.csv"
)

write.csv(
  x = specify_decimal(summary(mcmc.vcv5d)$solutions, 4),
  file = "./data/Model_outputs/mcmc.vcv5d.summary.csv"
)

write.csv(
  x = cbind(c("Chi2", "df", "p"), specify_decimal(as.data.table(VCVglmm::Wald.test.auto(mcmc.vcv5d)), 4)),
  file = "./data/Model_outputs/mcmc.vcv5d.wald.tests.csv"
)

# effect of annual and perennial
# phylogenetic components here
treevar.vcv5d <- VCVglmm::Bbar(phylo = tree.vcv5d, type = "Species") * mcmc.vcv5d$VCV[, "Sp1+Sp2"]

# keeping in mind, with Annual Perennial added...
a_a <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
  summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) + # at mean distance
  summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$dist2) + # at mean genus size
  summary(mcmc.vcv5d)$solutions[3, 1] * mean(vcv5d$Hectads_shared), # at mean hectads
sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
) # accounting for phylogeny
p_p <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
  summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) + # at mean distance
  summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$dist2) + # at mean genus size
  summary(mcmc.vcv5d)$solutions[5, 1] + # perennial perennial effect
  summary(mcmc.vcv5d)$solutions[3, 1] * mean(vcv5d$Hectads_shared), # at mean hectads
sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
) # accounting for phylogeny
# 1.32
# updated is 1.30
p_p / a_a

# hectad sharing for perennial plants.
lower <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
  summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) + # at mean distance
  summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$dist2) + # at mean genus size
  summary(mcmc.vcv5d)$solutions[5, 1] + # perennial perennial effect
  summary(mcmc.vcv5d)$solutions[3, 1] * 0, # at lowest hectad sharing
sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
) # accounting for phylogeny
higher <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
  summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) + # at mean distance
  summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$dist2) + # at mean genus size
  summary(mcmc.vcv5d)$solutions[5, 1] + # perennial perennial effect
  summary(mcmc.vcv5d)$solutions[3, 1] * 3755, # at mean hectads
sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
) # accounting for phylogeny
# 2.51 times
# 2.42 updated
higher / lower

# phylogenetic variance explained/signal
p1 <- posterior.mode(treevar.vcv5d / rowSums(cbind(
  treevar.vcv5d,
  2 * mcmc.vcv5d$VCV[, "T1+T2"],
  mcmc.vcv5d$VCV[, "units"]
)))
hpd1 <- HPDinterval(treevar.vcv5d / rowSums(cbind(
  treevar.vcv5d,
  2 * mcmc.vcv5d$VCV[, "T1+T2"],
  mcmc.vcv5d$VCV[, "units"]
)))
# species variance
p2 <- posterior.mode(2 * mcmc.vcv5d$VCV[, "T1+T2"] / rowSums(cbind(
  treevar.vcv5d,
  2 * mcmc.vcv5d$VCV[, "T1+T2"],
  mcmc.vcv5d$VCV[, "units"]
)))
hpd2 <- HPDinterval(2 * mcmc.vcv5d$VCV[, "T1+T2"] / rowSums(cbind(
  treevar.vcv5d,
  2 * mcmc.vcv5d$VCV[, "T1+T2"],
  mcmc.vcv5d$VCV[, "units"]
)))

# heatmap of the joint distribution of hectad sharing and branch length.

x_test <- seq(0, 0.3, 0.01)
y_test <- seq(0, 3000, 100)

summary(vcv5d$dist2)

x_testd <- seq(0, 200, 6.5)
y_testd <- seq(0, 3000, 100)

vec <- matrix(nrow = length(x_testd), ncol = length(y_testd))
for (i in 1:length(x_testd)) {
  for (j in 1:length(y_testd)) {
    vec[i, j] <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
      summary(mcmc.vcv5d)$solutions[3, 1] * y_testd[j] + # at mean hectad sharing
      summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$Genus.Size) + # at mean genus size
      summary(mcmc.vcv5d)$solutions[4, 1] + # annual - perennial hybridisation
      summary(mcmc.vcv5d)$solutions[2, 1] * x_testd[i], # effect of branch length
    sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
    )
  }
}
colnames(vec) <- y_testd
rownames(vec) <- x_testd

col_fun <- colorRampPalette(colors = c(cbPalette[1], cbPalette[2]))

pdf(file = "./figures/Supplementary/Correlation_genetic_geographic_mod_output_vcv5d.pdf")
lattice::levelplot(t(vec),
  col.regions = col_fun(1000),
  ylab = "Branch length (Mya)", xlab = "Hectad sharing",
  scales = list(tck = c(1, 0), x = list(rot = 45))
)
dev.off()

# the effect of branch length
mcmc.vcv5.plot <- function() {
  # x <- seq(0,0.3,0.0001)
  x <- seq(0, 200, 0.1)
  y <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
    summary(mcmc.vcv5d)$solutions[3, 1] * mean(vcv5d$Hectads_shared) + # at mean hectad sharing
    summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[4, 1] + # annual - perennial hybridisation
    summary(mcmc.vcv5d)$solutions[2, 1] * x, # effect of branch length
  sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
  ) # accounting for phylogeny
  yup <- pnorm(summary(mcmc.vcv5d)$solutions[1, 3] + # intercept
    summary(mcmc.vcv5d)$solutions[3, 3] * mean(vcv5d$Hectads_shared) + # at mean hectad sharing
    summary(mcmc.vcv5d)$solutions[6, 3] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[4, 3] + # annual - perennial hybridisation
    summary(mcmc.vcv5d)$solutions[2, 3] * x, # effect of branch length
  sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
  ) # accounting for phylogeny
  ylo <- pnorm(summary(mcmc.vcv5d)$solutions[1, 2] + # intercept
    summary(mcmc.vcv5d)$solutions[3, 2] * mean(vcv5d$Hectads_shared) + # at mean hectad sharing
    summary(mcmc.vcv5d)$solutions[6, 2] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[4, 2] + # annual - perennial hybridisation
    summary(mcmc.vcv5d)$solutions[2, 2] * x, # effect of branch length
  sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
  ) # accounting for phylogeny
  ggplot(data.table(x = x, y = y, yup = yup), aes(x = x)) +
    geom_line(size = 2, aes(y = y)) +
    geom_line(size = 1, lty = 2, aes(y = yup)) +
    geom_line(size = 1, lty = 2, aes(y = ylo)) +
    geom_vline(xintercept = vcv5d[, .(mean(dist2)), by = .(Genus1)][, .(mean(V1))][1]$V1, lty = 2, col = cbPalette[7], size = 2) +
    # xlab(label = "Branch length distance between two parental species") +
    xlab(label = "Divergence time between two parental species (Mya)") +
    ylab(label = "Probability of forming a hybrid") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    xlim(0, 200)
  # xlim(0,0.3)
}
mcmc.vcv5.plot()

ggsave(filename = "./figures/Supplementary/Distance_hectads_branch_model1d.pdf", plot = mcmc.vcv5.plot(), width = 12, height = 9, device = "pdf")

# the effect of hectad sharing

x <- seq(0, 3000, 1)
y <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
  summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) + # mean branch length
  summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$Genus.Size) + # at mean genus size
  summary(mcmc.vcv5d)$solutions[3, 1] * x, # hectads shared
sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
) # accounting for phylogeny

c("Min" = min(y), "Max" = max(y), "Mean", mean(y))

mcmc.vcv5.plot2 <- function() {
  x <- seq(0, 3000, 1)
  y <- pnorm(summary(mcmc.vcv5d)$solutions[1, 1] + # intercept
    summary(mcmc.vcv5d)$solutions[2, 1] * mean(vcv5d$dist2) + # mean branch length
    summary(mcmc.vcv5d)$solutions[6, 1] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[3, 1] * x, # hectads shared
  sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
  ) # accounting for phylogeny
  yup <- pnorm(summary(mcmc.vcv5d)$solutions[1, 3] + # intercept
    summary(mcmc.vcv5d)$solutions[2, 3] * mean(vcv5d$dist2) + # mean branch length
    summary(mcmc.vcv5d)$solutions[6, 3] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[3, 3] * x, # hectads shared
  sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
  ) # accounting for phylogeny
  ylo <- pnorm(summary(mcmc.vcv5d)$solutions[1, 2] + # intercept
    summary(mcmc.vcv5d)$solutions[2, 2] * mean(vcv5d$dist2) + # mean branch length
    summary(mcmc.vcv5d)$solutions[6, 2] * mean(vcv5d$Genus.Size) + # at mean genus size
    summary(mcmc.vcv5d)$solutions[3, 2] * x, # hectads shared
  sd = sqrt(mean(mcmc.vcv5d$VCV[, "units"] + 2 * mcmc.vcv5d$VCV[, "T1+T2"] + treevar.vcv5d))
  ) # accounting for phylogeny
  ggplot(data.table(x = x, y = y), aes(x = x, y = y)) +
    geom_line(size = 2) +
    geom_line(size = 1, lty = 2, aes(y = yup)) +
    geom_line(size = 1, lty = 2, aes(y = ylo)) +
    geom_vline(xintercept = mean(vcv5d$Hectads_shared), lty = 2, col = "red", size = 2) +
    xlab(label = expression(paste("Pairwise overlap in distribution", " km"^2))) +
    ylab(label = "Probability of forming a hybrid") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30)
    )
}
mcmc.vcv5.plot2()
## not in supplementary materials ##
ggsave(filename = "./figures/Supplementary/Pairwise_overlap_model1d.pdf", plot = mcmc.vcv5.plot2(), width = 12, height = 9, device = "pdf")

# extract the posterior modes, which are the greatest in value?
tree_modesd <- VCVglmm::Solapply(mcmc.vcv5d, mean)[Group == "Sp1"][order(-Grouped_Value)]
fwrite(
  tree_modesd[, .(species_or_node = Variable, posterior_mean = Grouped_Value)],
  "./data/Model_outputs/posterior_mean_phylogenetic_blups_vcv5d.csv"
)
# can we subset the important nodes? 1053, 68 (top), 450, 225
tree_modes2d <- tree_modesd[grepl("Node", Variable)][, Variable := substring(Variable, 9)]
# tree_modes3d <- tree_modes2d[Variable %in% c(1059, 75, 207, 457)]
tree_modes3d <- tree_modes2d[Variable %in% c(749)]

pdf(file = "./figures/Supplementary/two_phyogeniesd.pdf", width = 16, height = 16)
# two clades. Many nodes will be subtrees of other trees.
par(mfrow = c(1, 2), mai = c(0.05, 0.05, 0.05, 0.05))
ape::plot.phylo(treeio::tree_subset(node = 749, tree = tree.vcv5d), font = 4, edge.width = 3)
ape::plot.phylo(treeio::tree_subset(node = 226, tree = tree.vcv5d), font = 4, edge.width = 3)
dev.off()

### PLOIDY BELOW ###

# now add into the mix the ploidy level data.
# *be careful: note to self* after doing the intial add
# manually check against the hybrid flora to check the hybridising taxa ploidy levels are correct
# this can be done by filtering on y == 1 and unique(On)!
# added a few more species from the Kew C-value data base 21.11.19
chromosome_data6 <- fread("./data/Trait_databases/chromosome_data_updated181119.csv")

# update species names
# name changes 2 from above
# so 865 species in the phylogeny for which ploidy is known with certainty
chromosome_data6[name_changes2, on = .(Species)][!is.na(Ploidy)] %>%
  ggplot(aes(x = Ploidy)) +
  geom_histogram()
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
all_combs2 <- sjoin2[, .(
  Species1 = Species, Species2 = i.Species,
  Cross_Ploid = pbapply::pblapply(
    1:dim(sjoin2)[1],
    function(x) {
      ifelse(test = sjoin2$Ploidy[[x]] == sjoin2$i.Ploidy[[x]],
        yes = "Homoploid",
        no = "Heteroploid"
      )
    }
  ),
  Mean_Ploid = pbapply::pblapply(
    1:dim(sjoin2)[1],
    function(x) {
      mean(c(sjoin2$Ploidy[[x]], sjoin2$i.Ploidy[[x]]))
    }
  )
)]
# make the 'On' column
all_combs2[, On := sort_names(paste(Species1, Species2, sep = " "))]
# join to main data
# vcv6d <- all_combs2[, .(Cross_Ploid, On)][vcv5d, on = "On"]
vcv6d <- all_combs2[, .(Cross_Ploid, Mean_Ploid, On)][vcv5d, on = "On"]
# turn list to vector.
vcv6d[, Cross_Ploid := unlist(as.factor(as.character(vcv6d$Cross_Ploid)))]
vcv6d <- vcv6d[Cross_Ploid != "NULL"]
# always check structure
str(vcv6d)
# change these four columns to factor columns (perhaps this is the problem, as.factor, factor or as.character?)
# for_factor <- c("Sp1", "Sp2", "T1", "T2")
# for(col in for_factor){
#  set(vcv6, j=col, value=as.character(vcv6[[col]]))
# }
# make the levels of each of T1 and Sp1 the same as T2 and Sp2.
vcv6d$T1 <- factor(vcv6d$T1, levels = levels(vcv6d$T2))
vcv6d$Sp1 <- factor(vcv6d$Sp1, levels = levels(vcv6d$Sp2))
vcv6d <- subset(vcv6d, !is.na(Hectads_shared))
# relevel so homoploid is baseline?
vcv6d$Mean_Ploid <- unlist(vcv6d$Mean_Ploid)
vcv6d$Cross_Ploid <- relevel(x = vcv6d$Cross_Ploid, ref = "Homoploid")
save(vcv6d, file = "./data/Backups/vcv6d.RData")

# load("./data/Backups/vcv6d.RData")

### test run on ploidy dataset... ###
# from 568 species in the old run, to 684 species.
tree.ploidd <- ape::drop.tip(tree.4d, tree.4d$tip.label[!tree.4d$tip.label %in% union(vcv6d$Sp1, vcv6d$Sp2)])
tree.ploidd
# pull out the inverse matrix
Ainv.ploid <- inverseA(tree.ploidd, scale = FALSE)$Ainv
save(Ainv.ploid, file = "./data/Backups/Ainv.ploidd.RData")

prior.h1d <- list(
  R = list(V = diag(1), fix = 1),
  G = list(
    G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1) * 1000),
    G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1) * 1000)
  )
)

# gonna be a long run (it was.)
# 165 genera, up from 141!
# 684 species, up from 568 (nearly half British Flora!)
mcmc.ploidd <- MCMCglmm(y ~ dist2 + Cross_Ploid + Hectads_shared + Annual_Perennial + Genus.Size,
  family = "threshold",
  random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
  ginverse = list(
    Sp1 = Ainv.ploid,
    Sp2 = Ainv.ploid
  ),
  data = vcv6d,
  prior = prior.h1d,
  pr = TRUE,
  nitt = 13000 * 100,
  thin = 10 * 100,
  burnin = 3000 * 100
)
summary(mcmc.ploidd)

save(mcmc.ploidd, file = "./data/Backups/mcmc.ploidd.RData")
load("./data/Backups/mcmc.ploidd.RData")

write.csv(
  x = specify_decimal(summary(mcmc.ploidd)$solutions, 4),
  file = "./data/Model_outputs/mcmc.ploidd.summary.csv"
)

write.csv(
  x = cbind(c("Chi2", "df", "p"), specify_decimal(as.data.table(VCVglmm::Wald.test.auto(mcmc.ploidd)), 4)),
  file = "./data/Model_outputs/mcmc.ploidd.wald.tests.csv"
)

VCVglmm::MCMCfixplot(mod = mcmc.ploidd)

# What is the difference between the same ploidy and cross ploidy levels?
# as tree wasnt scaled, we need to multiply the phylogenetic variance by the tree variance.
treevar.ploidd <- VCVglmm::Bbar(phylo = tree.ploidd, type = "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]

# keeping in mind, with Annual Perennial added...
same_ploidy <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
  summary(mcmc.ploidd)$solutions[2, 1] * mean(vcv6d$dist2) + # at mean distance
  summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
  summary(mcmc.ploidd)$solutions[4, 1] * mean(vcv6d$Hectads_shared), # at mean hectads
sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + treevar.ploidd))
) # accounting for phylogeny
cross_ploidy <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
  summary(mcmc.ploidd)$solutions[2, 1] * mean(vcv6d$dist2) + # at mean distance
  summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
  summary(mcmc.ploidd)$solutions[3, 1] + # cross ploid effect
  summary(mcmc.ploidd)$solutions[4, 1] * mean(vcv6d$Hectads_shared), # at mean hectads
sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + treevar.ploidd))
) # accounting for phylogeny
# 1.308
# now 1.36
same_ploidy / cross_ploidy

# what about the effect of annual perennial, at mean genetic distance
# at mean hectad sharing at same ploidy level
annual_effect <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept (annual-annual parents)
  summary(mcmc.ploidd)$solutions[2, 1] * mean(vcv6d$dist2) + # at mean distance
  summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
  summary(mcmc.ploidd)$solutions[4, 1] * mean(vcv6d$Hectads_shared), # at mean hectads
sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + treevar.ploidd))
) # accounting for phylogeny
perennial_effect <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
  summary(mcmc.ploidd)$solutions[6, 1] * mean(vcv6d$dist2) + # perennial-perennial parents
  summary(mcmc.ploidd)$solutions[2, 1] * mean(vcv6d$dist2) + # at mean distance
  summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
  summary(mcmc.ploidd)$solutions[4, 1] * mean(vcv6d$Hectads_shared), # at mean hectads
sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + treevar.ploidd))
) # accounting for phylogeny
# perennial-perennial effect is 1.02 times higher probability of creating successful hybrid
# than annual-annual effect. This is not significant.
# new model is 5.4
perennial_effect / annual_effect


# visualise the fixed effects from raw data.
mcmc.ploidvis <- vcv6d[, .(
  Prob = mean(y),
  SE.Prob = sd(y) / sqrt(.N)
), by = .(Cross_Ploid)]

ggplot(mcmc.ploidvis, aes(y = Prob, x = Cross_Ploid)) +
  geom_point() +
  geom_errorbar(aes(ymin = Prob - SE.Prob, ymax = Prob + SE.Prob), width = 0.3) +
  xlab(label = "Parental Ploidy difference") +
  ylab(label = "Probability of forming a hybrid") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# predicted model outputs
ploid_plot <- function() {
  x <- seq(0, 3000, 1)
  y <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
    summary(mcmc.ploidd)$solutions[2, 1] * mean(vcv6d$dist2) + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[5, 1] + # annual-perennial
    summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[4, 1] * x, # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] +
    2 * mcmc.ploidd$VCV[, "T1+T2"] +
    Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # homoploid

  yup <- pnorm(summary(mcmc.ploidd)$solutions[1, 3] + # intercept
    summary(mcmc.ploidd)$solutions[2, 3] * mean(vcv6d$dist2) + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[5, 3] + # annual-perennial
    summary(mcmc.ploidd)$solutions[7, 3] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[4, 3] * x, # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] +
    2 * mcmc.ploidd$VCV[, "T1+T2"] +
    Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # homoploid
  ylo <- pnorm(summary(mcmc.ploidd)$solutions[1, 2] + # intercept
    summary(mcmc.ploidd)$solutions[2, 2] * mean(vcv6d$dist2) + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[7, 2] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[5, 2] + # annual-perennial
    summary(mcmc.ploidd)$solutions[4, 2] * x, # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] +
    2 * mcmc.ploidd$VCV[, "T1+T2"] +
    Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # homoploid

  y2 <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
    summary(mcmc.ploidd)$solutions[2, 1] * mean(vcv6d$dist2) + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[3, 1] + # effect of heteroploidy
    summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[5, 1] + # annual perennial
    summary(mcmc.ploidd)$solutions[4, 1] * x, # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # heteroploid
  y2up <- pnorm(summary(mcmc.ploidd)$solutions[1, 3] + # intercept
    summary(mcmc.ploidd)$solutions[2, 3] * mean(vcv6d$dist2) + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[7, 3] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[3, 3] + # effect of heteroploidy
    summary(mcmc.ploidd)$solutions[5, 3] + # annual-perennial
    summary(mcmc.ploidd)$solutions[4, 3] * x, # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # heteroploid
  y2lo <- pnorm(summary(mcmc.ploidd)$solutions[1, 2] + # intercept
    summary(mcmc.ploidd)$solutions[2, 2] * mean(vcv6d$dist2) + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[7, 2] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[3, 2] + # effect of heteroploidy
    summary(mcmc.ploidd)$solutions[5, 2] + # annual-perennial
    summary(mcmc.ploidd)$solutions[4, 2] * x, # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # heteroploid

  data <- melt(data.table(x = x, Homoploid = y, Heteroploid = y2),
    measure.vars = c("Homoploid", "Heteroploid"),
    variable.name = "Parental Ploidy"
  )
  data$lCI <- c(ylo, y2lo)
  data$uCI <- c(yup, y2up)

  ggplot(data, aes(x = x, y = value, colour = `Parental Ploidy`)) +
    geom_line(size = 2) +
    geom_line(aes(y = lCI), lty = 2, size = 1.2) +
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = `Parental Ploidy`), alpha = 0.15) +
    geom_line(aes(y = uCI), lty = 2, size = 1.2) +
    geom_vline(xintercept = mean(vcv6d$Hectads_shared), lty = 2, col = cbPalette[7], size = 2) +
    xlab(label = expression(paste("Pairwise overlap in distribution", " km"^2))) +
    ylab(label = "Probability of forming a hybrid") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = c(cbPalette[1], cbPalette[2])) +
    scale_color_manual(values = c(cbPalette[1], cbPalette[2]))
}
ploid_plot()
ggsave(filename = "./figures/Supplementary/Cross_ploidd_hectads.pdf", plot = ploid_plot(), width = 12, height = 9, device = "pdf")

ploid_plot_2 <- function() {
  # x <- seq(0, 0.3,0.01)
  x <- seq(0, 120, 0.1)
  y <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
    summary(mcmc.ploidd)$solutions[2, 1] * x + # genetic distance
    summary(mcmc.ploidd)$solutions[5, 1] + # annual-perennial
    summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[4, 1] * mean(vcv6d$Hectads_shared), # mean of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] +
    2 * mcmc.ploidd$VCV[, "T1+T2"] +
    Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # homoploid

  yup <- pnorm(summary(mcmc.ploidd)$solutions[1, 3] + # intercept
    summary(mcmc.ploidd)$solutions[2, 3] * x + #  genetic distance
    summary(mcmc.ploidd)$solutions[5, 3] + # annual-perennial
    summary(mcmc.ploidd)$solutions[7, 3] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[4, 3] * mean(vcv6d$Hectads_shared), # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] +
    2 * mcmc.ploidd$VCV[, "T1+T2"] +
    Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # homoploid
  ylo <- pnorm(summary(mcmc.ploidd)$solutions[1, 2] + # intercept
    summary(mcmc.ploidd)$solutions[2, 2] * x + #  genetic distance
    summary(mcmc.ploidd)$solutions[7, 2] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[5, 2] + # annual-perennial
    summary(mcmc.ploidd)$solutions[4, 2] * mean(vcv6d$Hectads_shared), # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] +
    2 * mcmc.ploidd$VCV[, "T1+T2"] +
    Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # homoploid

  y2 <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
    summary(mcmc.ploidd)$solutions[2, 1] * x + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[3, 1] + # effect of heteroploidy
    summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[5, 1] + # annual perennial
    summary(mcmc.ploidd)$solutions[4, 1] * mean(vcv6d$Hectads_shared), # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # heteroploid
  y2up <- pnorm(summary(mcmc.ploidd)$solutions[1, 3] + # intercept
    summary(mcmc.ploidd)$solutions[2, 3] * x + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[7, 3] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[3, 3] + # effect of heteroploidy
    summary(mcmc.ploidd)$solutions[5, 3] + # annual-perennial
    summary(mcmc.ploidd)$solutions[4, 3] * mean(vcv6d$Hectads_shared), # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # heteroploid
  y2lo <- pnorm(summary(mcmc.ploidd)$solutions[1, 2] + # intercept
    summary(mcmc.ploidd)$solutions[2, 2] * x + # mean of the genetic distance
    summary(mcmc.ploidd)$solutions[7, 2] * mean(vcv6d$Genus.Size) + # mean genus size
    summary(mcmc.ploidd)$solutions[3, 2] + # effect of heteroploidy
    summary(mcmc.ploidd)$solutions[5, 2] + # annual-perennial
    summary(mcmc.ploidd)$solutions[4, 2] * mean(vcv6d$Hectads_shared), # for a range of hectads
  sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + Bbar(tree.ploidd, "Species") * mcmc.ploidd$VCV[, "Sp1+Sp2"]))
  ) # heteroploid

  data <- melt(data.table(x = x, Homoploid = y, Heteroploid = y2),
    measure.vars = c("Homoploid", "Heteroploid"),
    variable.name = "Parental Ploidy"
  )
  data$lCI <- c(ylo, y2lo)
  data$uCI <- c(yup, y2up)

  ggplot(data, aes(x = x, y = value, colour = `Parental Ploidy`)) +
    geom_line(size = 2) +
    geom_line(aes(y = lCI), lty = 2, size = 1.2) +
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = `Parental Ploidy`), alpha = 0.15) +
    geom_line(aes(y = uCI), lty = 2, size = 1.2) +
    geom_vline(xintercept = vcv6d[, .(mean(dist2)), by = .(Genus1)][, .(mean(V1))][1]$V1, lty = 2, col = cbPalette[7], size = 2) +
    xlab(label = expression(paste("Divergence time between a pair of species (Mya)"))) +
    ylab(label = "Probability of forming a hybrid") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = c(cbPalette[1], cbPalette[2])) +
    scale_color_manual(values = c(cbPalette[1], cbPalette[2]))
}
ploid_plot_2()
ggsave(filename = "./figures/Supplementary/Cross_ploidd_branch.pdf", plot = ploid_plot_2(), width = 12, height = 9, device = "pdf")

# phylogenetic variance ~58%
p3 <- posterior.mode(treevar.ploidd / rowSums(cbind(
  treevar.ploidd,
  2 * mcmc.ploidd$VCV[, "T1+T2"],
  mcmc.ploidd$VCV[, "units"]
)))
hpd3 <- HPDinterval(treevar.ploidd / rowSums(cbind(
  treevar.ploidd,
  2 * mcmc.ploidd$VCV[, "T1+T2"],
  mcmc.ploidd$VCV[, "units"]
)))
# species variance ~35%
p4 <- posterior.mode(2 * mcmc.ploidd$VCV[, "T1+T2"] / rowSums(cbind(
  treevar.ploidd,
  2 * mcmc.ploidd$VCV[, "T1+T2"],
  mcmc.ploidd$VCV[, "units"]
)))
hpd4 <- HPDinterval(mcmc.ploidd$VCV[, "T1+T2"] / rowSums(cbind(
  treevar.ploidd,
  mcmc.ploidd$VCV[, "T1+T2"],
  mcmc.ploidd$VCV[, "units"]
)))

# no different visually between the ploidy levels, so keep one figure as above.
vec2 <- matrix(nrow = length(x_testd), ncol = length(y_testd))
for (i in 1:length(x_testd)) {
  for (j in 1:length(y_testd)) {
    vec2[i, j] <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
      summary(mcmc.ploidd)$solutions[4, 1] * y_testd[j] + # at mean hectad sharing
      summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # at mean genus size
      summary(mcmc.ploidd)$solutions[5, 1] + # annual - perennial hybridisation
      summary(mcmc.ploidd)$solutions[2, 1] * x_testd[i], # effect of branch length
    sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + treevar.ploidd))
    )
  }
}
colnames(vec2) <- y_testd
rownames(vec2) <- x_testd

vec3 <- matrix(nrow = length(x_testd), ncol = length(y_testd))
for (i in 1:length(x_testd)) {
  for (j in 1:length(y_testd)) {
    vec3[i, j] <- pnorm(summary(mcmc.ploidd)$solutions[1, 1] + # intercept
      summary(mcmc.ploidd)$solutions[3, 1] +
      summary(mcmc.ploidd)$solutions[4, 1] * y_testd[j] + # at mean hectad sharing
      summary(mcmc.ploidd)$solutions[7, 1] * mean(vcv6d$Genus.Size) + # at mean genus size
      summary(mcmc.ploidd)$solutions[5, 1] + # annual - perennial hybridisation
      summary(mcmc.ploidd)$solutions[2, 1] * x_testd[i], # effect of branch length
    sd = sqrt(mean(mcmc.ploidd$VCV[, "units"] + 2 * mcmc.ploidd$VCV[, "T1+T2"] + treevar.ploidd))
    )
  }
}
colnames(vec3) <- y_testd
rownames(vec3) <- x_testd

grid.arrange(levelplot(t(vec2),
  col.regions = col_fun(1000),
  main = "Homoploid hybrids",
  ylab = "Branch length", xlab = "Hectad sharing",
  scales = list(tck = c(1, 0), x = list(rot = 45))
),
levelplot(t(vec3),
  col.regions = col_fun(1000),
  main = "Heteroploid hybrids",
  ylab = "Branch length", xlab = "Hectad sharing",
  scales = list(tck = c(1, 0), x = list(rot = 45))
),
ncol = 2
)



# table of variance components

write.csv(
  x = data.table(
    `Variance Component` = c(
      "Model 1 Phylogenetic Variance", "Model 1 Species Variance",
      "Model 2 Phylogenetic Variance", "Model 2 Species Variance"
    ),
    `Posterior Mode` = c(specify_decimal(p1, 4), specify_decimal(p2, 4), specify_decimal(p3, 4), specify_decimal(p4, 4)),
    `Lower Credible Interval` = c(specify_decimal(hpd1, 4)[1], specify_decimal(hpd2, 4)[1], specify_decimal(hpd3, 4)[1], specify_decimal(hpd4, 4)[1]),
    `Upper Credible Interval` = c(specify_decimal(hpd1, 4)[2], specify_decimal(hpd2, 4)[2], specify_decimal(hpd3, 4)[2], specify_decimal(hpd4, 4)[2])
  ),
  file = "./data/Model_outputs/variance_componentsd.csv"
)
