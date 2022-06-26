# The genetic landscape of hybridisation in the British flora #
# Max R. Brown, Peter M. Hollingsworth, Laura Forrest, Michelle Hart, Laura Jones, Col Ford, Tim Rich, Natasha de Vere, Alex D. Twyford #

# All code Max Brown #
# Last updated: 04.06.20 #

# Libraries needed #
library(ape)
library(seqinr)
library(data.table)
library(Taxonstand)
library(ggplot2)

`%!in%` <- function(x,y)!('%in%'(x,y))

# need to make names correct.
# in making the names correct, we lose some data, but we still have more than
# enough to actually answer the questions we want.
stace4 <- fread("./data/Trait_databases/DTOL_updated_list_070220.csv")
# make the double barrels separate by a dash
stace4[, Stace_4_Species := gsub("^([[:alnum:]]* [[:alnum:]]*) ", "\\1-", Stace_4_Species)]
stace4[, Plantlist := gsub("^([[:alnum:]]* [[:alnum:]]*) ", "\\1-", Plantlist)]
stace4[, `taxon name` := gsub("^([[:alnum:]]* [[:alnum:]]*) ", "\\1-", `taxon name`)]

#### CPDNA SEQUENCES ####
sequ <- read.alignment("./data/Alignment_fastas/final_merged_for_pairwise_dist.fasta", format = "fasta")
# get sequence names and remove underscore
sequ_names <- as.data.table(sequ[[2]])
sequ_names[, Species := gsub("_", " ", V1)][, V1 := NULL]
# hit against the plant list
#sequ_names2 <- Taxonstand::TPL(sequ_names$Species, silent = FALSE)
#save(sequ_names2, file = "../data/Backups/sequ_names2.RData")
load("./data/Backups/sequ_names2.RData")
sequ_names[, Plantlist := paste(sequ_names2$New.Genus, sequ_names2$New.Species, sep = " ")]
# deal with double barrel names
sequ_names[, Plantlist := gsub("^([[:alnum:]]* [[:alnum:]]*) ", "\\1-", Plantlist)]
# add ID column
sequ_names[, ID := 1:nrow(sequ_names)]


# merge names. This merges stace 4 with the names from the alignment
sequ_names3 <- sequ_names[stace4[,.(Stace_4_Species, Plantlist, `taxon name`)], on = .(Plantlist)]
sequ_names4 <- sequ_names3[!is.na(ID)]
# I have lost some data here but not much.
sequ_names4[ID %!in% setdiff(1:1421, sequ_names4$ID)]
# 1345 species in this dataset
sequ_names5 <- sequ_names4[!duplicated(Stace_4_Species)]

# all that remains now is to subset the sequence data, and replace the names with Stace!
str(sequ)
sequ[[1]] <- 1345 # dist.alignment freaks out if you dont change this...
sequ[[2]] <- sequ_names5$Stace_4_Species 
sequ[[3]] <- sequ[[3]][sequ_names5$ID]

dist.sequ <- dist.alignment(sequ)

#### ITS SEQUENCES ####
sequ2 <- read.alignment("./data/Alignment_fastas/final_merged_ITS_extraction.fasta", format = "fasta")

sequ_names.2 <- as.data.table(sequ2[[2]])
sequ_names.2[, Species := gsub("_", " ", V1)][, V1 := NULL]
# hit against the plant list
sequ_names2.2 <- Taxonstand::TPL(sequ_names.2$Species, silent = FALSE)
sequ_names.2[, Plantlist := paste(sequ_names2.2$New.Genus, sequ_names2.2$New.Species, sep = " ")]
# deal with double barrel names
sequ_names.2[, Plantlist := gsub("^([[:alnum:]]* [[:alnum:]]*) ", "\\1-", Plantlist)]
# add ID column
sequ_names.2[, ID := 1:nrow(sequ_names.2)]


# merge names. This merges stace 4 with the names from the alignment
sequ_names3.2 <- sequ_names.2[stace4[,.(Stace_4_Species, Plantlist, `taxon name`)], on = .(Plantlist)]
sequ_names4.2 <- sequ_names3.2[!is.na(ID)]
# I have lost some data here but not much.
sequ_names4.2[ID %!in% setdiff(1:1354, sequ_names4.2$ID)]
# 1282 species in this dataset
sequ_names5.2 <- sequ_names4.2[!duplicated(Stace_4_Species)]

# all that remains now is to subset the sequence data, and replace the names with Stace!
sequ2[[1]] <- 1282 # dist.alignment freaks out if you dont change this...
sequ2[[2]] <- sequ_names5.2$Stace_4_Species 
sequ2[[3]] <- sequ2[[3]][sequ_names5.2$ID]

dist.sequ2 <- dist.alignment(sequ2)

# I think NaNs mean that the data is missing in the alignment. So it should be 
# removed.

##### Part 1: Plastid distances #####

##### Part 1.1: Genetic distance for hybridising/non hybridising taxa #####

# now format the data so that it is in the same format as in Part_Four_Models.R
# put into long data format and set as data table
dist.sequ.cp <- setDT(as.data.frame(as.matrix(dist.sequ)), keep.rownames = TRUE)
# melt as before
dist.sequ.cp2 <- melt(dist.sequ.cp)
# remove missing data (NaN's)
dist.sequ.cp2 <- dist.sequ.cp2[!is.nan(value),]

# change names to remove '_' 's
dist.sequ.cp2$rn <- gsub(pattern = "_", replacement = " ", x = dist.sequ.cp2$rn)
dist.sequ.cp2$variable <- gsub(pattern = "_", replacement = " ", x = dist.sequ.cp2$variable)
# check that the spaces worked
nchar(dist.sequ.cp2[1,2])

# import hybriddata
hybriddata <- fread("./data/Hybrid_flora_of_the_British_Isles/hybriddata_Stace4_removed.csv")
hybriddata$BOTH <- paste(hybriddata$`Parent A`, hybriddata$`Parent B`, sep = " ")
dist.sequ.cp2$BOTH <- paste(dist.sequ.cp2$rn, dist.sequ.cp2$variable, sep = " ")

# remove reciprocal hybrids. Long winded but works?
# strsplit BOTH on space. Takes a while!
dist.sequ.cp2$test2 <- unlist(lapply(lapply(strsplit(x = dist.sequ.cp2$BOTH, split = " "), sort), function(x) paste(x, collapse = "")))
dist.sequ.cp3 <- unique(dist.sequ.cp2, by = "test2")

# which hybrid combinations are in the distance matrix?
# create a logical vector
hybriddata$test2<- unlist(lapply(lapply(strsplit(x = hybriddata$BOTH, split = " "), sort), function(x) paste(x, collapse = "")))
dist.sequ.cp3[,test := dist.sequ.cp3$test2 %in% hybriddata$test2]

# filter data when taxa is not congeneric
# create genus columns for parent 1 and 2
dist.sequ.cp3$Genus1 <- sub(" .*", "", dist.sequ.cp3$rn) 
dist.sequ.cp3$Genus2 <- sub(" .*", "", dist.sequ.cp3$variable) 

# filter where genus of both parents are the same!
dist.sequ.cp4 <- dist.sequ.cp3[Genus1 == Genus2 & rn != variable]

# plot
ggplot(dist.sequ.cp4, aes(x=test, y=value))+ geom_jitter(col = "tomato", alpha=0.1)+
  geom_boxplot(alpha = 0)+
  theme_bw(base_line_size = 0)+
  ylab(label = "Genetic distance (Plastid)")+
  xlab(label = "Hybrid formed between a pair of taxa?")

ggplot(dist.sequ.cp4, aes(x = value))+
  geom_histogram(aes(fill = test), position = position_dodge())

# top/lowest 10 genetic diverged
# top/lowest 10 genetic diverged
dist.sequ.cp4[test == TRUE][order(-value)][1:20,]
dist.sequ.cp4[test == TRUE][order(value)][1:20,]

##### Part 1.2: Ranking genera for mean genetic distance #####
sum.dist.sequ.cp3 <- dist.sequ.cp4[test == TRUE, .(Mean_dist = mean(value),
                                          SE = sd(value)/sqrt(.N)), by = "Genus1"][order(-Mean_dist)]

# then plot the results

ggplot(sum.dist.sequ.cp3, aes(x = reorder(Genus1, Mean_dist), y=Mean_dist))+
  geom_point()+
  geom_errorbar(aes(ymin = Mean_dist - SE, ymax = Mean_dist + SE))+
  #geom_vline(xintercept = 82, col="red", lty=2)+
  theme_bw(base_line_size = 0)+
  theme(axis.text.x.bottom = element_text(hjust = 1, angle = 60))+
  ylab(label = "Mean Genetic Distance (Plastid)")+
  xlab(label = "Hybridising Genus Rank Order")

##### Part 1.3: Regressing mean genetic distance vs size of genus #####

# for all genera in the British Flora
sum.dist.sequ.cp4_all <- dist.sequ.cp4[, .(Mean_dist = mean(value),
                              SE = sd(value)/sqrt(.N)), by = "Genus1"][order(-Mean_dist)]
# change the genus name
colnames(sum.dist.sequ.cp4_all)[1] <- "Genus"
# merge again
resizCP <- data3[sum.dist.sequ.cp4_all, on = "Genus"]
# but we want to add a hybridising vs non hybridising genera...
resizCP$Hybrid_Genus <- resizCP[[1]] %in% regenCP[[1]]

# now plot
ggplot(resizCP, mapping = aes(y = Mean_dist, x = Genus_Size, group = Hybrid_Genus))+
  #geom_errorbar(mapping = aes(ymax = Mean_co_dist +SE, ymin = Mean_co_dist-SE), width = 0.5)+
  geom_point(aes(colour = Hybrid_Genus), size = 2)+
  ggrepel::geom_text_repel(mapping = aes(label= ifelse(test = resizCP$Genus %in% c("Carex",
                                                                                  "Euphrasia",
                                                                                  "Rosa",
                                                                                  "Rumex",
                                                                                  "Trifolium",
                                                                                  "Sedum",
                                                                                  "Salix",
                                                                                  "Ranunculus",
                                                                                  "Sorbus",
                                                                                  "Juncus"), yes = as.character(resizCP$Genus), no = NA) ), point.padding = 0.5, box.padding = 0.5)+
  theme_bw()+
  xlab(label = "Genus size")+
  ylab(label = "Mean Genetic Distance (Plastid)")





##### Part 2: ITS distances #####

##### Part 2.1: Genetic distance for hybridising/non hybridising taxa #####

# now format the data so that it is in the same format as in Part_Four_Models.R
# put into long data format and set as data table
dist.sequ.its <- setDT(as.data.frame(as.matrix(dist.sequ2)), keep.rownames = TRUE)
# melt as before
dist.sequ.its2 <- melt(dist.sequ.its)
# remove missing data (NaN's)
dist.sequ.its2 <- dist.sequ.its2[!is.nan(value),]

# change names to remove '_' 's
dist.sequ.its2$rn <- gsub(pattern = "_", replacement = " ", x = dist.sequ.its2$rn)
dist.sequ.its2$variable <- gsub(pattern = "_", replacement = " ", x = dist.sequ.its2$variable)
# check that the spaces worked
nchar(dist.sequ.its2[1,2])

# hybriddata from Part_One_finaldata.R
hybriddata$BOTH <- paste(hybriddata$`Parent A`, hybriddata$`Parent B`, sep = " ")
dist.sequ.its2$BOTH <- paste(dist.sequ.its2$rn, dist.sequ.its2$variable, sep = " ")

# remove reciprocal hybrids. Long winded but works?
# strsplit BOTH on space. Takes a while! look at that beautiful double lapply!!
dist.sequ.its2$test2 <- unlist(lapply(lapply(strsplit(x = dist.sequ.its2$BOTH, split = " "), sort), function(x) paste(x, collapse = "")))
dist.sequ.its3 <- unique(dist.sequ.its2, by = "test2")

# which hybrid combinations are in the distance matrix?
# create a logical vector
hybriddata$test2<- unlist(lapply(lapply(strsplit(x = hybriddata$BOTH, split = " "), sort), function(x) paste(x, collapse = "")))
dist.sequ.its3[,test := dist.sequ.its3$test2 %in% hybriddata$test2]

# filter data when taxa is not congeneric
# create genus columns for parent 1 and 2
dist.sequ.its3$Genus1 <- sub(" .*", "", dist.sequ.its3$rn) 
dist.sequ.its3$Genus2 <- sub(" .*", "", dist.sequ.its3$variable) 

# filter where genus of both parents are the same!
dist.sequ.its4 <- dist.sequ.its3[Genus1 == Genus2 & rn != variable]

# plot
ggplot(dist.sequ.its4, aes(x=test, y=value))+ geom_jitter(col = "tomato", alpha=0.1)+
  geom_boxplot(alpha = 0)+
  theme_bw(base_line_size = 0)+
  ylab(label = "Genetic distance (ITS)")+
  xlab(label = "Hybrid formed between a pair of taxa?")

# top/lowest 10 genetic diverged
dist.sequ.its4[test == TRUE][order(-value)][1:20,]
dist.sequ.its4[test == TRUE][order(value)][1:20,]


# have a look at the genus Potamogeton.
# Is hybridisation restricted to closely related taxa?
# Or if hybridisation has homogenised species differences by introgression or chloroplast capture...
# Potamegeton produces sterile hybrids so it cannot be explanation two.
library(dplyr)

dist.sequ.its4[Genus1 == "Potamogeton"] %>%
  ggplot(aes(x = test, y = value))+
  geom_jitter(col = "tomato", alpha=0.5)+
  geom_boxplot(alpha = 0)+
  theme_bw(base_line_size = 0)+
  ylab(label = "Genetic distance (ITS)")+
  xlab(label = "Hybrid formed between a pair of taxa?")
  # 5 hybrid combinations with sequences
dim(dist.sequ.cp3[Genus1 == "Potamogeton" & Genus2 == "Potamogeton"])[1] # out of 231 possible combinations

# import the tree
TREE <- read.tree(file = "./Barcoding_Phylogeny_Dating/DNA_Barcoding.dated.treefile")
D <- cophenetic.phylo(TREE)
dist.tree <- setDT(as.data.frame(as.matrix(D)), keep.rownames = TRUE)
dist.tree2 <- melt(dist.tree)
dist.tree2 <- dist.tree2[!is.nan(value),]
dist.tree2$rn <- gsub(pattern = "_", replacement = " ", x = dist.tree2$rn)
dist.tree2$variable <- gsub(pattern = "_", replacement = " ", x = dist.tree2$variable)

dist.tree2$BOTH <- paste(dist.tree2$rn, dist.tree2$variable, sep = " ")

dist.tree2$test2 <- unlist(lapply(lapply(strsplit(x = dist.tree2$BOTH, split = " "), sort), function(x) paste(x, collapse = "")))
dist.tree3 <- unique(dist.tree2, by = "test2")

dist.tree3[,test := dist.tree3$test2 %in% hybriddata$test2]

# filter data when taxa is not congeneric
# create genus columns for parent 1 and 2
dist.tree3$Genus1 <- sub(" .*", "", dist.tree3$rn) 
dist.tree3$Genus2 <- sub(" .*", "", dist.tree3$variable) 

# filter where genus of both parents are the same!
dist.tree4 <- dist.tree3[Genus1 == Genus2 & rn != variable]

# both boxplots on the same plot

dist.sequ.its4[, Marker := "ITS"]
dist.sequ.cp4[, Marker := "Plastid"]
dist.tree4[, Marker := "Tree"]

bothbox <- rbind(dist.sequ.its4, dist.sequ.cp4, dist.tree4)

library(cowplot)

(bothboxplot <- ggplot(bothbox[Marker %in% c("ITS", "Plastid")], aes(x=test, y=value))+ geom_jitter(col = "tomato", alpha=0.1)+
  geom_boxplot(alpha = 0, size = 1.5)+
  facet_wrap(~ Marker)+
  theme_bw(base_line_size = 0)+
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y =  element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
    scale_x_discrete(limits = rev(c("TRUE", "FALSE")), labels = rev(c("Yes", "No")))+
  ylab(label = "Genetic distance") + xlab(label = ""))

breaks <- c(-2.302585, 0, 2.302585, 3.912023, 5.298317) # 1, 50, 100, 150, 200
labels <- c(0.1, round(c(exp(0), 10, exp(3.912023), exp(5.298317)), 0))

(treeboxplot <- ggplot(bothbox[Marker %in% c("Tree")], aes(x=test, y=log(value)))+ geom_jitter(col = "tomato", alpha=0.1)+
    geom_boxplot(alpha = 0, size = 1.5)+
    facet_wrap(~ Marker)+
    theme_bw(base_line_size = 0)+
    theme(strip.background = element_rect(fill = "white", colour = "white"), 
          strip.text = element_text(size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.text.y =  element_text(size = 30),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))+
    scale_x_discrete(limits = rev(c("TRUE", "FALSE")), labels = rev(c("Yes", "No")))+
    scale_y_continuous(
      breaks = breaks,
      labels = labels) +
    ylab(label = "Divergence time (Mya)") + 
    xlab(label = ""))

# merge plots
library(grid)
library(gridExtra)
x.grob <- textGrob("Hybrid formed between two taxa?", 
                   gp=gpar(col="black", fontsize=25, fontfamily = "Helvetica"))
both_distances <- cowplot::plot_grid(bothboxplot, treeboxplot, rel_widths = c(1.5, 1))
save_both_distances <- grid.arrange(arrangeGrob(both_distances, bottom = x.grob))
# 
ggsave(filename = "./figures/Main/Genetic_distance_boxplotsd.pdf", plot = bothboxplot, 
       device = "pdf", width = 11, height = 11, units = "in")
ggsave(filename = "./figures/Main/All_distances_boxplots.pdf", plot = save_both_distances, 
       device = "pdf", width = 11, height = 11, units = "in")

##### Part 2.2: Ranking genera for mean genetic distance #####

sum.dist.sequ.its <- dist.sequ.its4[test == TRUE, .(Mean_dist = mean(value),
                                                   SE = sd(value)/sqrt(.N),
                                                   N =.N), by = "Genus1"][order(-Mean_dist)]

# then plot the results

ggplot(sum.dist.sequ.its, aes(x = reorder(Genus1, Mean_dist), y=Mean_dist))+
  geom_point()+
  geom_errorbar(aes(ymin = Mean_dist - SE, ymax = Mean_dist + SE))+
  #geom_vline(xintercept = 82, col="red", lty=2)+
  theme_bw()+
  theme(axis.text.x.bottom = element_text(hjust = 1, angle = 60))+
  ylab(label = "Mean Genetic Distance (ITS)")+
  xlab(label = "Hybridising Genus Rank Order")

##### Part 2.3: Regressing mean genetic distance vs size of genus #####

# for all genera in the British Flora
sum.dist.sequ.its4_all <- dist.sequ.its4[, .(Mean_dist = mean(value),
                                           SE = sd(value)/sqrt(.N)), by = "Genus1"][order(-Mean_dist)]
# change the genus name
colnames(sum.dist.sequ.its4_all)[1] <- "Genus"
# merge again
resizITS <- data3[sum.dist.sequ.its4_all, on = "Genus"]
regenITS <- data3[sum.dist.sequ.its, on = .(Genus = Genus1)]
# but we want to add a hybridising vs non hybridising genera...
resizITS$Hybrid_Genus <- resizITS[[1]] %in% regenITS[[1]]

# now plot
ggplot(resizITS, mapping = aes(y = Mean_dist, x = Genus_Size, group = Hybrid_Genus))+
  #geom_errorbar(mapping = aes(ymax = Mean_co_dist +SE, ymin = Mean_co_dist-SE), width = 0.5)+
  geom_point(aes(colour = Hybrid_Genus), size = 2)+
  ggrepel::geom_text_repel(mapping = aes(label= ifelse(test = resizITS$Genus %in% c("Carex",
                                                                                    "Euphrasia",
                                                                                    "Rosa",
                                                                                    "Rumex",
                                                                                    "Trifolium",
                                                                                    "Sedum",
                                                                                    "Salix",
                                                                                    "Ranunculus",
                                                                                    "Sorbus",
                                                                                    "Juncus"), yes = as.character(resizITS$Genus), no = NA) ), box.padding = 0.5, point.padding = 0.5)+
  theme_bw()+
  xlab(label = "Genus size")+
  ylab(label = "Mean Genetic Distance (ITS)")







##### Part 3: Correlation in genetic distance between the two loci #####

# they are (maybe quite obviously..? very highly correlated)
cor.dist <- dist.sequ.its4[,c(3,5,6)][dist.sequ.cp4[,c(3,5,6)], on = "test2"]

# but more interestingly we can then colour by hybrid combinations
ggplot(cor.dist[!is.na(value) & !is.na(i.value),], mapping = aes(x = value, y = i.value, colour = test))+
  geom_point(aes(alpha = test), size =3)+
  theme_bw()+
  xlab(label = "ITS distances")+
  ylab(label = "Plastid distances")+
  guides(colour=guide_legend(title="Hybrid formed?"),
         alpha = guide_legend(title="Hybrid formed?"))



