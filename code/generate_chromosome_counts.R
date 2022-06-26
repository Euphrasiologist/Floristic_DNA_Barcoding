## Chromosome number additions to the hybrid flora, but mainly for cross ploidy analysis ##

library(data.table)
library(stringr)
library(ggplot2)
library(Taxonstand)
library(taxonlookup)
library(numbers)
library(dplyr)
library(tidyr)

## Extracting ploidy data

chr_nums <- fread("./data/Trait_databases/Chromosome_numbers_update.txt", fill = TRUE)

# this extracts taxonomic names and 2N
chr_nums2<-tidyr::separate(chr_nums, V1, sep = "(?=[[:upper:]]{1}[[:lower:]]{2,}[[:blank:]][[:lower:]]+)|(?=2N)|^(?!:).*$", into = c("Data", "Extract"),
                           convert = TRUE, remove = FALSE, fill = "right")
chr_nums2[nchar(Extract) > 2]$Extract
# 4th row after 2N is the chromosome number, 7493
length(chr_nums2[which(chr_nums2$Extract == "2N")+4]$V1) # chromosome numbers
chr_nums2[which(chr_nums2$Extract == "2N")[2:7493]-9]$Extract # -5, -6, -7, -8
# make the data table, logic is that there are specific intervals anchored at the 2N flag.
chr_nums3 <- data.table(Chromosome_numbers = chr_nums2[which(chr_nums2$Extract == "2N")+4]$V1,
                        Extract1 = ifelse(nchar(chr_nums2[which(chr_nums2$Extract == "2N")-5]$Extract) == 0, NA, chr_nums2[which(chr_nums2$Extract == "2N")-5]$Extract),
                        Extract2 = ifelse(nchar(chr_nums2[which(chr_nums2$Extract == "2N")-6]$Extract) == 0, NA, chr_nums2[which(chr_nums2$Extract == "2N")-6]$Extract),
                        Extract3 = ifelse(nchar(chr_nums2[which(chr_nums2$Extract == "2N")-7]$Extract) == 0, NA, chr_nums2[which(chr_nums2$Extract == "2N")-7]$Extract),
                        Extract4 = ifelse(nchar(chr_nums2[which(chr_nums2$Extract == "2N")-8]$Extract) == 0, NA, chr_nums2[which(chr_nums2$Extract == "2N")-8]$Extract),
                        Extract5 = ifelse(nchar(c(NA, chr_nums2[which(chr_nums2$Extract == "2N")[2:7493]-9]$Extract)) == 0, NA, c(NA, chr_nums2[which(chr_nums2$Extract == "2N")[2:7493]-9]$Extract)))


# what do all of these chromosome numbers look like?
setDT(as.data.frame(table(chr_nums3$Chromosome_numbers)))[order(-Freq)][101:200]

regs <- chr_nums3$Chromosome_numbers

# remove asterisks
regs <- gsub(pattern = "([0-9]{1-3})\\*", "\\1", regs)
# remove c., c. , ca , ca, ca., ca. ,...
regs <- gsub(pattern = "c\\.([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "c\\. ([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "ca([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "ca ([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "ca\\.([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "ca\\. ([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "cf\\.([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "cf\\. ([0-9]{1-3})", "\\1", regs)
regs <- gsub(pattern = "c([0-9]{1-3})", "\\1", regs)
# remove ?
regs <- gsub(pattern = "([0-9]{1-3})\\?", "\\1", regs)
# remove bracketed expressions
regs <- gsub(pattern = "([0-9]{1-3})\\(.*", "\\1", regs)

# coerce to as.integer
regs <- as.integer(regs)

chr_nums3$Chromosome_numbers <- regs

# remove pesky NA's
chr_nums4 <- unique(chr_nums3[, .(Species = dplyr::coalesce(chr_nums3$Extract1, chr_nums3$Extract2, chr_nums3$Extract3, chr_nums3$Extract4, chr_nums3$Extract5),
                                  Chromosome_numbers = Chromosome_numbers)])

# weird character removal
chr_nums4$Species <- gsub(x = chr_nums4$Species, replacement = "x", pattern = "×")

# update taxonomy now
updated_taxonomy <- Taxonstand::TPL(chr_nums4$Species)
#save(updated_taxonomy, file = "./hybridpropensity/data/updated_taxonomy.RData")
setDT(updated_taxonomy)
# bind with chr_nums4
chr_nums4.1 <- cbind(chr_nums4, updated_taxonomy[, .(New.Species = paste(New.Genus, New.Species, sep = " "), 
                                                     New.Infraspecific,
                                                     New.Infraspecific.rank)])

# 1722 species with individual chromosome counts (this number does actually include duplicates...)
# added 168 taxa by updating taxonomy here...
# update 12.11.19 1985 taxa now.
# chr_nums5 <- chr_nums4[!is.na(Chromosome_numbers) & !grepl("forma|subsp.|var.| × ", Species)]
chr_nums5 <- chr_nums4.1[!is.na(Chromosome_numbers) & nchar(New.Infraspecific) == 0][,-c("New.Infraspecific", "New.Infraspecific.rank", "Species")]
setnames(x = chr_nums5, old = "New.Species", new = "Species")

chr_nums5[, Genus := gsub(" .*", "", Species)]
higher_orders <- setDT(taxonlookup::lookup_table(chr_nums5$Species))
setnames(higher_orders, old = c("genus", "family", "order", "group"), new = c("Genus", "Family", "Order", "Group"))
# merge with higher orders and tidy, 1838 species with 'complete cases'
chr_nums6 <- higher_orders[chr_nums5, on = "Genus"]
chr_nums6 <- chr_nums6[complete.cases(chr_nums6),]
# keep chr_nums6 to look at all chromosome numbers regardless of multiple species records of differing ploidy
chr_nums7 <- chr_nums6[, .(Chromosome_numbers_update = paste(sort(unique(Chromosome_numbers)), collapse = ",")), by = "Species"][chr_nums6, on = "Species"]

## insert ccdb estimates here ##
ccdb_brit_fl <- fread(file = "./data/Trait_databases/ccdb_brit_fl.csv")

# why did I as numeric the column here? To remove species with multiple ploidy levels.
chr_nums7[,Chromosome_numbers_update2 := as.numeric(Chromosome_numbers_update)]
# sort between aneuploidy and ploidy variation
ploid_aneu_var <- chr_nums7[grepl(pattern = ",", x = Chromosome_numbers_update)]
ploid_aneu_var2 <- lapply(X = strsplit(x = ploid_aneu_var$Chromosome_numbers_update, split = ","), FUN = function(x){
  # get all pairwise combinations of ploidies/aneuploidies
  x <- as.numeric(x)
  combs <- combn(x, 2)
  # find the largest common denominator of each pair
  vec <- apply(combs, 2, numbers::mGCD)
  # get the modes out of this
  Modes <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
  }
  # but only the ones greater than 3, (1 & 2 indicate aneuploidy)
  M <- Modes(vec[vec > 3])
  # if there are two modes, take the lowest, as that is more likely to be the base number
  if(any(duplicated(M))){
    M <- max(M)
  } else
    M <- min(M)
  
  if(any(vec %in% c(1,2)) & is.na(M)){
    x <- "Aneuploidy"
  }
  if(any(vec %in% c(1,2) & !is.na(M))){
    x <- "Aneuploidy and Polyploidy"
  }
  if(!is.na(M)){
    x <- "Polyploidy"
  }
  if(is.na(M)){
    M <- "Unknown base"
  }
  
  return(list(x,M))
  
})

### BUG HERE, URTICA ASSIGNED WRONGLY... ###

ploid_aneu_var[, `:=`(Aneuploid = lapply(ploid_aneu_var2, "[[", 1),
                      Base = lapply(ploid_aneu_var2, "[[", 2))]
# 135 species with multiple ploidies
# update 12.11.19 137 species with multiple ploidies
unique(ploid_aneu_var[,-c("Aneuploid", "Base", "Chromosome_numbers", "Chromosome_numbers_update2")])
# 56 species displaying some level of aneuploidy
# 59 species 12.11.19
unique(ploid_aneu_var[Base == "Unknown base", -c(7,9,10)])
# 114 genera (up from 99)
# 118 genera 12.11.19
length(table(ploid_aneu_var$Genus))
# chr_nums7 is all aneuploid and polyploid species in the UK! Or there are chromosome count errors
chr_nums8 <- chr_nums7[complete.cases(Chromosome_numbers_update2)]
# distribution of chromosome numbers
# overall
ggplot(chr_nums8, aes(x = Chromosome_numbers))+geom_density()
# split by orders with more than 30 species
ggplot(chr_nums8[, .(N = .N), by = "Family"][chr_nums8, on = "Family"][N > 20][,Family := paste(Family, "(", N, ")")], aes(x = Chromosome_numbers))+
  geom_density(aes(group = Family), position = position_dodge())+
  scale_x_continuous(breaks = seq(0, 400, 5))+theme_bw()+facet_wrap(~Family)

# how many species in chromosome data not in the crossplod?
setdiff(chr_nums8$Species, crossplod6$Taxa) # 373
setdiff(crossplod6$Taxa, chr_nums8$Species) # 635 in 
## merge with known data (?)
tomerge1 <- crossplod6[, .(Species = Taxa, 
                           Hybrid_flora_chromosome = PublishedChromosomeNumber, 
                           Ploidy1 = PublishedPloidyLevel,
                           Ploidy2 = Ploidy)]
tomerge2 <- data.table(Ploidy2 = c("Diploid", "Tetraploid", "Pentaploid", "Hexaploid", "Octoploid"), Ploidy = c(2,4,5,6,8))[tomerge1, on = "Ploidy2"]

tomerge3 <- tomerge2[, .(Species = Species,
                         Hybrid_flora_chromosome = Hybrid_flora_chromosome,
                         Ploidy = dplyr::coalesce(tomerge1$Ploidy1, tomerge1$Ploidy))]
# tomerge3 is from the phylogeny dataset
setdiff(chr_nums8$Species, tomerge3$Species) # 373 species are unique to the chromosome dataset

ploid_data <- chr_nums8[,.(Species = Species, Chromosome_numbers_update2)][tomerge3, on = "Species"]

ploid_data[, Hybrid_flora_chromosome := ifelse(Hybrid_flora_chromosome == 0, NA, Hybrid_flora_chromosome)]

ploid_data[, .(Species = Species,
               Chromosome = dplyr::coalesce(Chromosome_numbers_update2, Hybrid_flora_chromosome),
               Ploidy = Ploidy)][is.na(Chromosome)]

#estimate ploidy from two sources, merge.
# work out base numbers from ploid_data

ploid_data[, Base := Chromosome_numbers_update2/Ploidy]
ploid_data[, Base2 := Hybrid_flora_chromosome/Ploidy]

# add Euphrasia
ploid_data[grepl("Euphrasia", Species, TRUE)]$Base <- 11
ploid_data[grepl("Euphrasia", Species, TRUE)]$Chromosome_numbers_update2 <- c(rep(44, 9), 22,22,22,44,44,22,44,44,44)

ploid_data[is.na(Chromosome_numbers_update2) & is.na(Hybrid_flora_chromosome),][order(Species)]$Species

# base numbers from genera from ploid_aneu_var, but these are all multi-ploidied
base_from_ploid_aneu_var <- unique(ploid_aneu_var[, .(Species = Species,
                                                      Base3 = as.numeric(Base))][!is.na(Base3)])[order(Species)]
base_from_ploid_aneu_var[grepl("Urtica", Species)]
Gbase_from_ploid_aneu_var <- base_from_ploid_aneu_var[, Genus := gsub(" .*", "", Species)][, .(Base3 = unique(Base3)), by = "Genus"]
# what should really be done is for each genus try each of the base chromosome numbers
# one should divide properly

ploid_data2 <- base_from_ploid_aneu_var[ploid_data, on = "Species"]

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

ploid_data2[, Genus := gsub(" .*", "", Species)]
ploid_data3 <- ploid_data2[!is.na(Base3) | !is.na(Base2) | !is.na(Base)]

ploid_data4 <- ploid_data3[, .(Species = Species, 
                               Genus = Genus,
                               Ploidy = Ploidy,
                               Base = ifelse(is.wholenumber(Base), Base, NA),
                               Base2 = ifelse(is.wholenumber(Base2), Base2, NA),
                               Base3 = ifelse(is.wholenumber(Base3), Base3), NA)][,-"V7"]

ploid_data5 <- ploid_data4[, lapply(list(Base, Base2, Base3), function(x){
  paste(na.omit(x), collapse = ",")
}), by = "Genus"]

ploid_data6 <- cbind(ploid_data5$Genus, apply(ploid_data5[,-"Genus"], 1, function(x) paste(x, collapse = ",")))
# yes!
ploid_data7 <- data.table(Genus = ploid_data5$Genus, 
                          Base_numbers = lapply(X = strsplit(x = ploid_data6[,2], split = ","), FUN = function(x) unique(as.numeric(x[grepl("[[:digit:]]+", x)]))))
# genera with multiple basic numbers
ploid_data7[lapply(ploid_data7$Base_numbers, length) >1]

# distribution of plant base numbers in the UK
# to extract each as its own row
ploid_data7[, .(unlist(Base_numbers)), by = "Genus"] %>%
  ggplot(aes(x = V1))+geom_density()

# 285 base number estimates... 330 now (+45 estimates)
ploid_data7[, .(unlist(Base_numbers)), by = "Genus"]

# let's get to the ploidy bit now...
PLOID1 <- ploid_data[,-c("Base", "Base2")][, lapply(list(Chromosome_numbers_update2, Hybrid_flora_chromosome), function(x){
  paste(na.omit(x), collapse = ",")
}), by = "Species"]

PLOID2 <- cbind(PLOID1$Species, apply(PLOID1[,-"Species"], 1, function(x) paste(x, collapse = ",")))

PLOID3 <- data.table(Species = PLOID1$Species, 
                     Combined_Chromosome_Counts = lapply(X = strsplit(x = PLOID2[,2], split = ","), 
                                                         FUN = function(x) unique(as.numeric(x[grepl("[[:digit:]]+", x)]))))
# merge on genus
PLOID3[, Genus := gsub(" .*", "", Species)]

PLOID4 <- ploid_data7[PLOID3, on = "Genus"]

PLOID4$Base_numbers <- lapply(PLOID4$Base_numbers, as.numeric)
PLOID4$Combined_Chromosome_Counts <- lapply(PLOID4$Combined_Chromosome_Counts, as.numeric)


# unlist base numbers
unlist_base <- PLOID4[, .(Base = unlist(Base_numbers)), by = c("Genus", "Species")]
# unlist both
unlist_chrom <- PLOID4[, .(Chrom = unlist(Combined_Chromosome_Counts)), by = c("Genus", "Species")]

ploidy_1 <- unlist_base[unlist_chrom, on = c("Genus", "Species")]

# these need to be checked...
# Erophila verna, n = 26? https://www.jstor.org/stable/25063898?seq=3#metadata_info_tab_contents [jstor.org]
# Urtica dioica n = 13
# 

# reassign ploidies!
ploidy_1 <- ploidy_1[, Ploidy_level := Chrom/Base][is.wholenumber(Ploidy_level)]

# distribution of ploidy levels (apparently some haploids..? wrong... and some very highly polyploid species... also wrong...)
ploidy_prov1 <- ploidy_1[Ploidy_level < 13 & Ploidy_level > 1] %>%
  ggplot(aes(x = Ploidy_level))+geom_histogram(binwidth = 1, colour = "black")+scale_x_continuous(breaks = seq(0,12,1))+
  theme_bw()+
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_text(size = 20),
        axis.title.y.left = element_text(size = 20),
        legend.title = element_text(size = 20))+
  xlab(label = "Ploidy level") +
  ylab(label = "Count")

ggsave(filename = "./hybridpropensity/figs/prov_ploidy1.jpeg", plot = ploidy_prov1, 
       device = "jpeg", width = 5, height = 6, units = "in")

# add higher taxa levels and plot
higher_taxa <- setDT(taxonlookup::lookup_table(species_list = ploidy_1$Species))
setnames(x = higher_taxa, old = c("genus", "family", "order", "group"), new = c("Genus", "Family", "Order", "Group"))

ploidy_prov2 <- ploidy_1[higher_taxa, on = "Genus"][, N := .N, by = "Family"][N > 16 & Ploidy_level < 13 & Ploidy_level > 1] %>%
  ggplot(aes(x = Ploidy_level))+
  geom_histogram(binwidth = 1)+
  scale_x_continuous(breaks = seq(0,12,1))+
  facet_wrap(~Family)+
  theme_bw()+
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_text(size = 20),
        axis.title.y.left = element_text(size = 20),
        legend.title = element_text(size = 20))+
  xlab(label = "Ploidy level") +
  ylab(label = "Count")

ggsave(filename = "./hybridpropensity/figs/prov_ploidy2.jpeg", plot = ploidy_prov2, 
       device = "jpeg", width = 10, height = 6, units = "in")


# maybe reconcile this...
setdiff(tomerge3$Species, ploidy_1$Species)
# cakile as a subsp, raphanus as subsp

# save ploidy_1
# clean viola!!!!!
# Tilia diploid

data.table::fwrite(x = ploidy_1, file = "./hybridpropensity/data/ploidy_estimates_rbge2.csv")

######### Use data from hybriddata #############

# hybriddata is from Part_One_finaldata.R. Carex is notably missing all ploidy information. [fix]
ploid_dat_1 <- hybriddata[, .(ParentA = `Parent A`,
                              PloidyA = `Ploidy A`,
                              ChromA = `Chromosome number A`,
                              ParentB = `Parent B`,
                              PloidyB = `Ploidy B`,
                              ChromB = `Chromosome Number B`)]
# get all species in one table
ploid_dat_2 <- rbind(ploid_dat_1[,.SD,.SDcols = patterns("A$")],
                     ploid_dat_1[,.SD,.SDcols = patterns("B$")], use.names = FALSE)
ploid_dat_2 <- unique(ploid_dat_2)
setnames(ploid_dat_2, c("ParentA", "PloidyA", "ChromA"), c("Species", "Ploidy", "Chromosome_number"))
ploid_dat_2[, Base := as.numeric(Chromosome_number)/as.numeric(Ploidy)]

ploid_dat_3 <- ploid_dat_2[is.wholenumber(Base)]

# first estimate from chr_nums7, additional chromosome numbers
# a few approaches now
# use base numbers to estimate new ploidies 1. From ploid_data7[, .(unlist(Base_numbers)), by = "Genus"]
# or 2. from ploid_dat_3 in the dataset chr_nums7
# Simply concatenate ploid_dat_3 and ploidy_1 (estimates from hybriddata and rbge respectively.)

### Estimate new ploidies ###
ploid_dat_3[, Genus := gsub(" .*", "", Species)]
ploid_ests_1 <- unique(ploid_dat_3[, .(Genus, Base)])
ploid_ests_2 <- ploid_data7[, .(Base = unlist(Base_numbers)), by = "Genus"]
# all genus level base chromosome number estimates
ploid_ests_3 <- unique(rbind(ploid_ests_1,ploid_ests_2))
# merge with chromosome number data
chromosome_data <- chr_nums7[, .(Species = Species_Update, Genus, Family, Order, Chromosome_numbers)]
chromosome_data2 <- chromosome_data[ploid_ests_3, on = "Genus", allow.cartesian = TRUE]
# shows the genera for which there is no chromosome level data: chromosome_data2[is.na(Species)] 
chromosome_data3 <- chromosome_data2[, Ploidy := as.numeric(Chromosome_numbers)/as.numeric(Base)]
chromosome_data3 <- chromosome_data3[!is.na(Species)]
chromosome_data3 <- chromosome_data3[is.wholenumber(Ploidy)]
# remove species with multiple ploidies? Or if there are multiple ploidies choose the lowest?
# former for now. Much easier!
chromosome_data3[, N := .N, by = "Species"]
chromosome_data3[N == 1] # 637 in this data set
# stuff in ploid_dat_3 (hybriddata) not in ploidy estimates
setdiff(ploid_dat_3$Species, chromosome_data3[N == 1]$Species)
setdiff(ploidy_1$Species, chromosome_data3[N == 1]$Species)

# bind with ploid_dat_3
chromosome_data4 <- unique(rbind(ploid_dat_3[, .(Species, Ploidy, Chromosome_number, Base)], 
                                 ploidy_1[, .(Species, Ploidy = Ploidy_level, Chromosome_number = Chrom, Base)],
                                 chromosome_data3[N == 1][, .(Species, Ploidy, Chromosome_number = Chromosome_numbers, Base)]))
# 995 ploidy estimates for species. If it is from a species with multiple ploidies, it is from the hybrid flora.
chromosome_data5 <- chromosome_data4[, N := .N, by = .(Species)][N == 1]

# error corrections
chromosome_data5[Species == "Rorippa islandica"]$Ploidy <- 2
chromosome_data5[Species == "Erodium moschatum"]$Ploidy <- 2
chromosome_data5[Species == "Erodium maritimum"]$Ploidy <- 2
chromosome_data5[Species == "Campanula trachelium"]$Ploidy <- 2
chromosome_data5[Species == "Spergularia media"]$Ploidy <- 2

table(chromosome_data5$Ploidy)

# species in hybrid data not in the estimated numbers
add <- setdiff(ploid_dat_3$Species, chromosome_data5$Species)
# to final data, add all of the hybrid combinations from the raw data in terms of their ploidy.
# the rest can be used with this. This is because we *know* the ploidy of hybrids for some species 
# but they may be estimated wrongly.
add2 <- ploid_dat_3[Species %in% add][-c(1,2,3,4,26)][order(Species)][-c(9,10,12,13,22,23,28)][, .(Ploidy, Chromosome_number, Base, N =.N), by = "Species"]

# could be more cleaning e.g. checking through all names etc...
# add2 originates from the hybrid flora, chromosome_data5 from the cytology database.
chromosome_data6 <- unique(rbind(chromosome_data5, add2))

for_numeric <- c("Ploidy", "Chromosome_number", "Base")

for (col in for_numeric){
  set(chromosome_data6, j=col, value=as.numeric(chromosome_data6[[col]]))
}
# just a quick visualisation. Might be used for cross ploidy manuscript too...
ggplot(chromosome_data6, aes(x = Ploidy))+geom_density()
# write data, about as final as it will get.
fwrite(x = chromosome_data6, file = "./hybridpropensity/data/chromosome_data_updated151019.csv")

# add in higher taxonomic ranks
higher_orders2 <- setDT(taxonlookup::lookup_table(species_list = chromosome_data6$Species))
setnames(higher_orders2, old = c("genus", "family", "order", "group"), new = c("Genus", "Family", "Order", "Group"))
# make genus column
chromosome_data6[, Genus := gsub(" .*", "", Species)]
chromosome_data7 <- higher_orders2[chromosome_data6, on = "Genus"]
# mycelis and lamiastrum have NA's
chromosome_data7[match(chromosome_data7$Order, NA) == 1]

Mycelis <- c("Mycelis", "Asteraceae", "Asterales", "Angiosperms", "Mycelis muralis", 2, 18, 9, 1)
chromosome_data7[match(x = "Mycelis", table = chromosome_data7$Genus), names(chromosome_data7) := as.list(Mycelis)]
Lamiastrum <- c("Lamiastrum", "Lamiaceae", "Lamiales", "Angiosperms", "Lamiastrum galeobdolon", 2, 18, 9, 1)
chromosome_data7[match(x = "Lamiastrum", table = chromosome_data7$Genus), names(chromosome_data7) := as.list(Lamiastrum)]

# Chromosome numbers
chromosome_data7[, .(MeanC = mean(Chromosome_number),
                     N = .N), by = "Order"][order(-MeanC)][N >5]
# Ploidy
chromosome_data7[, .(MeanP = mean(Ploidy),
                     N = .N), by = "Family"][order(-MeanP)][N >5]

# number of species.
dim(chromosome_data7)[1]
# proportion of British Flora [how to calculate???]
length(intersect(chromosome_data7$Species, gsub("_", " ", tree$tip.label)))/length(tree$tip.label)
# those in out chromosome data but not in the phylogeny
length(setdiff(chromosome_data7$Species, tree.3$tip.label))

# add genus, family, order sizes.
chromosome_data7 <- chromosome_data7[, .(Genus_size = .N), by = .(Genus)][chromosome_data7, on = .(Genus)][, -"N"]
chromosome_data7 <- chromosome_data7[, .(Family_size = .N), by = .(Family)][chromosome_data7, on = .(Family)]
chromosome_data7 <- chromosome_data7[, .(Order_size = .N), by = .(Order)][chromosome_data7, on = .(Order)]

updated_ploidy_graph <- chromosome_data7 %>%
  ggplot(aes(x = Ploidy))+geom_histogram(binwidth = 1, colour = "black")+scale_x_continuous(breaks = seq(0,30,1))+
  theme_bw()+
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_text(size = 20),
        axis.title.y.left = element_text(size = 20),
        legend.title = element_text(size = 20))+
  xlab(label = "Ploidy level") +
  ylab(label = "Count")

updated_ploidy_graph2 <- chromosome_data7[Family_size > 25 & Genus != "Rubus"] %>%
  ggplot(aes(x = Ploidy))+geom_histogram(binwidth = 1, colour = "black")+scale_x_continuous(breaks = seq(0,30,1))+
  theme_bw()+
  facet_wrap(~Family)+
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_text(size = 20),
        axis.title.y.left = element_text(size = 20),
        legend.title = element_text(size = 20))+
  xlab(label = "Ploidy level") +
  ylab(label = "Count")

ggsave(plot = updated_ploidy_graph, filename = "./hybridpropensity/figs/Updated_ploidy.pdf", device = "pdf")
### Phylogenetic signal of ploidy ###

tree.3 <- read.tree(file = "./hybridpropensity/data/tree.3")
tree.3$tip.label <- gsub("_", " ", tree.3$tip.label)

tree.3.ploid <- ape::drop.tip(tree.3, tree.3$tip.label[!tree.3$tip.label %in% chromosome_data7$Species])

tree.3.ploid$node.label <- NULL
# pull out the inverse matrix... 
Ainv.tree.3.ploid<-inverseA(tree.3.ploid, scale=FALSE)$Ainv

prior.Ainv.tree.3.ploid <-list(R=list(V=diag(1), nu=0.02), 
                               G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*100),
                                      G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*100)))
# create phylogenetic variable
chromosome_data8 <- chromosome_data7[Species %in% tree.3.ploid$tip.label]
chromosome_data8$Species <- factor(chromosome_data8$Species)
levels(chromosome_data8$Species) <- tree.3.ploid$tip.label
chromosome_data8[, animal := Species]

# should we correct for genus size?
mcmc.tree.3.ploid <- MCMCglmm(Ploidy ~ 1, 
                              family = "poisson",
                              random = ~ Species + animal,
                              ginverse = list(animal = Ainv.tree.3.ploid),
                              data= chromosome_data8,
                              prior=prior.Ainv.tree.3.ploid,
                              pr=TRUE, 
                              nitt = 13000*50,
                              thin = 10*50,
                              burnin = 3000*50)
summary(mcmc.tree.3.ploid)

# because scale = FALSE on the inverseA call
treevar.tree.3.ploid <- mean(diag(vcv(tree.3.ploid)))*mcmc.tree.3.ploid$VCV[,"animal"]
# estimates of phylogenetic signal
posterior.mode(treevar.tree.3.ploid/rowSums(mcmc.tree.3.ploid$VCV))
HPDinterval(treevar.tree.3.ploid/rowSums(mcmc.tree.3.ploid$VCV))
