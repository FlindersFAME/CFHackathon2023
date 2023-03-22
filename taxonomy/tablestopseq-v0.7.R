# CF Kraken output -> Pavian -> Phyloseq
# Taxonomy Team
## data were seperated into a "read counts" table and a "taxonomy table"
## metadata was added up and merged to produce phyloseq object
# Read data and set libraries ----------

#setwd() #change to match folder hierarchy in wd

reads <- read.csv("reads.csv")
taxonomy <- read.csv("taxonomy.csv")
metadata <- read.csv("metadata2.csv")

# library(metagenomeSeq)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(dplyr)

### convert table to tabular split version
taxtable <- taxonomy %>%
  as_tibble() %>%
  separate(taxon, sep=">", c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Subspecies"))

taxtable_m <- as.matrix(taxtable)

# Create phyloseq obj ---------
otu <- otu_table(reads, taxa_are_rows = T)
otu <- otu[,-1]
tax <- tax_table(taxtable_m)

# add zeros
otu[is.na(otu)] <- 0

ps <-  phyloseq(otu, tax)
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

rownames(metadata) <- metadata$sample
meta <- sample_data(metadata)

setdiff(rownames(metadata), (colnames(otu)))

ps.1 <- merge_phyloseq(ps, meta)
ps.1
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

## Ranks and levels ---------------------------
rank_names(ps.1) # [1] "taxID"      "Kingdom"    "Phylum"     "Class"      "Order"      "Family"     "Genus"      "Species"    "Subspecies"
# sort( as.character( unique( tax_table(ps)[, "Phylum"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Class"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Order"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Family"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Genus"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Species"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Subspecies"] ) ))

## Rarefaction ------------------------------
rare.object <- ps.1

min(sample_sums(ps.1) ) # 18196
min(taxa_sums(ps.1)) # 1

max(sample_sums(ps.1) ) # 521423
max(taxa_sums(ps.1)) # 430405

sort(sample_sums(ps.1))
# X698917_20180128_S_Run4barcode03 X698917_20190119_S_Run7barcode05 X698917_20171207_S_Run7barcode10 X788707_20180313_S_Run2barcode10 X788707_20180301_S_Run2barcode07
# 18196                            23475                            47061                            57090                            110912
# X788707_20171213_S_Run9barcode04 X658355_20180321_S_Run1Barcode01 X658355_20171204_S_Run6barcode11 X788707_20181116_S_Run8barcode08 X658355_20180122_S_Run8barcode05
# 155532                           156437                           191382                         205677                           521423  


#RAREFY PACKAGE
# library(Rarefy)

library(vegan)
# # Rarefaction plot
# par(xpd=T, mar=par()$mar+c(0,0,0,10))
dev.off()

# ps.1@sam_data$patient.id <- as.character(ps.1@sam_data$patient.id)

# Plot for rarefaction curve
# rare.plot <- rarecurve(t(otu_table(ps.1)), label = F, ylab = "Number of features", xlab = "Number of reads")
# legend(locator(1),legend=c("),pt.bg=cols,pch=21,bty="n",ncol=1,cex = 0.75 ,pt.cex = 0.75, title = "Treatment", title.adj = 0)
# abline(v = 11976, col = "red", lty = 4)
# dev.off()

#Rarefaction step
seed <- 123
ps.rare <- rarefy_even_depth(ps.1, sample.size = 18196,
                             rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

# `set.seed(123)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(123); .Random.seed` for the full vector
# ...
# 1966OTUs were removed because they are no longer 
# present in any sample after random subsampling
# 
# ...

# Diversity: Calculate observed ASVs and the Shannon diversity index--------
# obs.ps.raw <- plot_richness(ps.1, measures=c("Observed"))
# obs.ps.raw
# str(obs.ps)
# 
# chao.ps.raw <- plot_richness(ps.1, measures=c("Chao1"))
# chao.ps.raw
# str(chao.ps)
# 
# shan.ps.raw <- plot_richness(ps.1, measures=c("Shannon"))
# shan.ps.raw
# colnames(ps.1@sam_data)

# Ouput file with rarefied data
obs.ps.rare <- plot_richness(ps.rare, measures=c("Observed"))
chao.ps.rare <- plot_richness(ps.rare, measures=c("Chao1"))
shan.ps.rare <- plot_richness(ps.rare, measures=c("Shannon"))

out.r1 <- data.frame(
  sample=shan.ps.rare$data$samples,
  observed=obs.ps.rare$data$value,
  shannon=shan.ps.rare$data$value,
  chao1=chao.ps.rare$data$value,
  patient.id=chao.ps.rare$data$patient.id,
  date=chao.ps.rare$data$date
)

out.r1$patient.id <- as.factor(out.r1$patient.id)
out.r1$date <- as.Date(out.r1$date, format = "%d/%m/%Y")
str(out.r1)
out.r1
#                              sample observed     shannon     chao patient.id       date
# 1  X658355_20171204_S_Run6barcode11      781    2.161164 405.0909     658355 2017-12-04
# 2  X658355_20180122_S_Run8barcode05     1901    2.327209 687.1000     658355 2018-01-22
# 3  X658355_20180321_S_Run1Barcode01      322    1.837654 195.5556     658355 2018-03-21
# 4  X698917_20171207_S_Run7barcode10      612    2.780461 700.3036     698917 2017-12-09
# 5  X698917_20180128_S_Run4barcode03      324    2.776247 697.5789     698917 2018-01-28
# 6  X698917_20190119_S_Run7barcode05      297    2.510563 456.1765     698917 2019-01-19
# 7  X788707_20171213_S_Run9barcode04      586    2.246014 348.8966     788707 2017-12-13
# 8  X788707_20180301_S_Run2barcode07      522    2.099383 816.1538     788707 2018-03-01
# 9  X788707_20180313_S_Run2barcode10      338    2.098574 473.0000     788707 2018-03-12
# 10 X788707_20181116_S_Run8barcode08      970    2.126763 595.5161     788707 2018-11-16

shanon.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=shannon, colour = patient.id))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Shannon's diversity index (rarefied)")+
  xlab("Patient ID")
shanon.div.rare.plot

chao.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=chao1, colour = patient.id))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Chao1 diversity index (rarefied)")+  
  xlab("Patient ID")
chao.div.rare.plot

library(ggpubr)
shan.chao.div.rare <- ggarrange(shanon.div.rare.plot,chao.div.rare.plot)
# ggsave("Diversity-shannon-chao-rare.pdf", plot =shan.chao.div.rare, width = 17)

# Chao over time
chao.time.plot <- ggplot(out.r1, aes(x= date, y=chao1, colour = patient.id))+
  geom_point()+facet_wrap(.~patient.id, scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()

# Ordination --------
ps.rare
ps.rare@sam_data$patient.id <- as.character(ps.rare@sam_data$patient.id)
ps.rare@sam_data$date <- as.character(ps.rare@sam_data$date)

set.seed(123)
PCoA_ord <- ordinate(ps.rare, method="PCoA", distance="bray", weighted=TRUE)
PCoA_all.plot <- plot_ordination(ps.rare, PCoA_ord, label = "date")

PCoA_all <- PCoA_all.plot + 
  geom_point(aes(fill = as.character(patient.id)), size=2.5, shape = 21, colour = "black")+
  theme_bw()
PCoA_all
# ggsave("Ordination_PCoA-bray-curtis_patients-rarefied-data.pdf", plot =PCoA_all, height = 5.5, width = 6)

PCoA_ord
PCoA_ord$vectors
dev.off()

# Relative abundance plots -------------------------
ps.1
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

ps.rare
# otu_table()   OTU Table:         [ 1231 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 1231 taxa by 9 taxonomic ranks ]

rank_names(ps.1) # "taxID"      "Kingdom"    "Phylum"     "Class"      "Order"      "Family"     "Genus"      "Species"    "Subspecies"

(all.K.psq <- ps.1)
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

## remove taxa not assigned as Bacteria
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","cellular organisms", "Eukaryota", "root", "unclassified", "Viruses") )
bacteria.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
bacteria.psq
# otu_table()   OTU Table:         [ 2620 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 2620 taxa by 9 taxonomic ranks ]

## remove taxa not assigned as Cellular organims
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" " "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","Eukaryota", "Bacteria", "root", "unclassified", "Viruses") )
cellular.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
cellular.psq
# otu_table()   OTU Table:         [ 3 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3 taxa by 9 taxonomic ranks ]

## remove taxa not assigned as Eukaryota
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","cellular organisms", "Bacteria", "root", "unclassified", "Viruses") )
eukaryota.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
eukaryota.psq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 428 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 428 taxa by 9 taxonomic ranks ]

## Relative abundance plots All divisions in one: -------------------------
rel_abun.Phylum.phy_obj <- transform_sample_counts(tax_glom(ps.1, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel_abun.Phylum.phy_obj
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3175 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3175 taxa by 9 taxonomic ranks ]

rel_abun.Genus.phy_obj <- transform_sample_counts(tax_glom(ps.1, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel_abun.Genus.phy_obj
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2558 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 2558 taxa by 9 taxonomic ranks ]

# Identify significantly trending phyla and those with substantial proportion
phylum.plot <- plot_bar(rel_abun.Phylum.phy_obj , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# phylum.plot <- phylum.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.pdf", plot = phylum.plot, height = 20, width = 10)

genus.plot <- plot_bar(rel_abun.Genus.phy_obj , x = "sample", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.pdf", plot = genus.plot, height = 49, width = 10)

## Relative abundance plots top taxa: -------------------------

### Bacteria -------------
#### Phylum -------------
rel.abun.bact.phy <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel.abun.bact.phy

bact.phylum.plot <- plot_bar(rel.abun.bact.phy , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")#+ 
  theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.pdf", plot = phylum.plot, height = 20, width = 10)

# Break down phyla into more manageable numbers
p <- bact.phylum.plot
# str(p$data)

## Identify rare Phyla
hist(p$data$Abundance)
out <- p$data
# str(out)

## Set threshold for major phyla
summary(out$Abundance)

## how many Phyla at thresholds > X% in any sample ?
#  modify X until number of major phyla is manageable, e.g. ~10 or so
X <- 0.15

major_phyla.bact <- levels(as.factor(as.character( out$Phylum[which(out$Abun >= X )]  )))
major_phyla.bact
#   [1] " Bacteroidota"                "Actinomycetota"               "Bacillota"                    "Bacteria candidate phyla"    
#   [5] "Bacteroidetes/Chlorobi group" "Fusobacteria"                 "Pseudomonadota"               "Tenericutes"                 
#   [9] "Terrabacteria"                "Terrabacteria group"     

sel.row <- which(out$Phylum %in% major_phyla.bact) # qty 
length(sel.row) #24960

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 99.84085 % of relative abundance is covered
100 - 99.84085 # 0.15915 % left representing rare phyla
rel.abun.cov.bact.phyla <- 99.84085
perc.rare.bact.phyla <- 0.15915

str(out) # do this below if $Phylum is a factor
out$Phylum <- as.character(out$Phylum)
out$Phylum[-sel.row] <- "Other minor phyla"
out$Phylum <- as.factor(out$Phylum)

str(p$data)

sel.row.plot <- which(p$data$Phylum %in% major_phyla.bact) 
str(sel.row.plot)
p$data$Phylum <- as.character(p$data$Phylum)
p$data$Phylum[ -sel.row.plot ] <- "Other minor phyla"
p$data$Phylum <- as.factor(p$data$Phylum)

levels(p$data$Phylum)
#   [1] " Bacteroidota"                "Actinomycetota"               "Bacillota"                    "Bacteria candidate phyla"    
#   [5] "Bacteroidetes/Chlorobi group" "Fusobacteria"                 "Other minor phyla"            "Pseudomonadota"              
#   [9] "Tenericutes"                  "Terrabacteria"                "Terrabacteria group"     

p$labels$fill <- "Phyla"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor phyla to last position, retitle if required
old_phyla<- c(
   " Bacteroidota"                ,"Actinomycetota"    ,           "Bacillota"          ,          "Bacteria candidate phyla"  ,  
   "Bacteroidetes/Chlorobi group" ,"Fusobacteria"      ,           "Other minor phyla"  ,          "Pseudomonadota"            ,  
   "Tenericutes"                  ,"Terrabacteria"     ,           "Terrabacteria group"    )      

## re-label phyla for plot
new_phyla<- c(
  "Bacteroidota"                ,"Actinomycetota"    ,           "Bacillota"          ,          "Bacteria candidate phyla"  ,  
  "Bacteroidetes/Chlorobi group" ,"Fusobacteria"      ,           "Other minor phyla"  ,          "Pseudomonadota"            ,  
  "Tenericutes"                  ,"Terrabacteria"     ,           "Terrabacteria group"    )     
#change to character, then revert to factor afterwards
p$data$Phylum <- as.character(p$data$Phylum)

for (i in 1:length(old_phyla)) {
  #i<-1
  p$data$Phylum[which(p$data$Phylum == old_phyla[i])] <- new_phyla[i]
}

bact.phyla.top.plot <- p + theme_bw() +
  theme(
    #axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1) 
    axis.text.x  = element_text(angle=60, vjust=1, hjust = 1) 
  )+
  labs(x = NULL, y = "Relative abundance (%)" )

bact.phyla.top.plot
# ggsave("Bacterial-relative-abundance-phyla.pdf", plot = bact.phyla.top.plot)

#### Genus -------------
rel.abun.bact.gen <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.bact.gen

bact.genus.plot <- plot_bar(rel.abun.bact.gen , x = "sample", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.pdf", plot = genus.plot, height = 49, width = 10)

# Break down genera into more manageable numbers
p <- bact.genus.plot
# str(p$data)

## Identify rare Genera
hist(p$data$Abundance)
out <- p$data
# str(out)

## Set threshold for major genera
summary(out$Abundance)

## how many Genera at thresholds > X% in any sample ?
#  modify X until number of major genera is manageable, e.g. ~30 or so
X <- 0.7

major_genera.bact <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.bact
#   [1] "Abiotrophia"                               "Acinetobacter"                            
#   [3] "Actinomyces"                               "Alloprevotella"                           
#   [5] "Bacillales Family XI. Incertae Sedis"      "Bacillus"                                 
#   [7] "Burkholderia"                              "Corynebacterium"                          
#   [9] "Cronobacter"                               "Cutibacterium"                            
#   [11] "Dolosicoccus"                              "Eubacteriales Family XIII. Incertae Sedis"
#   [13] "Granulicatella"                            "Klebsiella/Raoultella group"              
#   [15] "Lancefieldella"                            "Lawsonella"                               
#   [17] "Leptotrichia"                              "Levilactobacillus"                        
#   [19] "Listeria"                                  "Mycoplasmopsis"                           
#   [21] "Oribacterium"                              "Porphyromonas"                            
#   [23] "Prevotella"                                "Psychrobacter"                            
#   [25] "Rothia"                                    "Solobacterium"                            
#   [27] "Staphylococcus"                            "Streptococcus"                            
#   [29] "unclassified Lachnospiraceae"              "Veillonella"   

sel.row <- which(out$Genus %in% major_genera.bact) # qty 
length(sel.row) #7390

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 93.29846 % of relative abundance is covered
100 - 93.29846 # 6.70154 % left representing rare genera
rel.abun.cov.bact.genera <- 93.29846
perc.rare.bact.genera <- 6.70154

str(out) # do this below if $Genus is a factor
out$Genus <- as.character(out$Genus)
out$Genus[-sel.row] <- "Other minor genera"
out$Genus <- as.factor(out$Genus)

str(p$data)

sel.row.plot <- which(p$data$Genus %in% major_genera.bact) 
str(sel.row.plot)
p$data$Genus <- as.character(p$data$Genus)
p$data$Genus[ -sel.row.plot ] <- "Other minor genera"
p$data$Genus <- as.factor(p$data$Genus)

levels(p$data$Genus)
#  [1] "Abiotrophia"                               "Acinetobacter"                             "Actinomyces"                              
#  [4] "Alloprevotella"                            "Bacillales Family XI. Incertae Sedis"      "Bacillus"                                 
#  [7] "Burkholderia"                              "Corynebacterium"                           "Cronobacter"                              
# [10] "Cutibacterium"                             "Dolosicoccus"                              "Eubacteriales Family XIII. Incertae Sedis"
# [13] "Granulicatella"                            "Klebsiella/Raoultella group"               "Lancefieldella"                           
# [16] "Lawsonella"                                "Leptotrichia"                              "Levilactobacillus"                        
# [19] "Listeria"                                  "Mycoplasmopsis"                            "Oribacterium"                             
# [22] "Other minor genera"                        "Porphyromonas"                             "Prevotella"                               
# [25] "Psychrobacter"                             "Rothia"                                    "Solobacterium"                            
# [28] "Staphylococcus"                            "Streptococcus"                             "unclassified Lachnospiraceae"             
# [31] "Veillonella"                              

p$labels$fill <- "Genera"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor genera to last position, retitle if required
old_genera<- c(
        "Abiotrophia"     , "Acinetobacter"                       ,  "Actinomyces"                              ,
        "Alloprevotella"  , "Bacillales Family XI. Incertae Sedis",  "Bacillus"                                 ,
        "Burkholderia"    , "Corynebacterium"                     ,  "Cronobacter"                              ,
        "Cutibacterium"   , "Dolosicoccus"                        ,  "Eubacteriales Family XIII. Incertae Sedis",
        "Granulicatella"  , "Klebsiella/Raoultella group"         ,  "Lancefieldella"                           ,
        "Lawsonella"      , "Leptotrichia"                        ,  "Levilactobacillus"                        ,
        "Listeria"        , "Mycoplasmopsis"                      ,  "Oribacterium"                             ,
        "Porphyromonas"   , "Prevotella"                          , 
        "Psychrobacter"   , "Rothia"                              ,  "Solobacterium"                            ,
        "Staphylococcus"  , "Streptococcus"                       ,  "unclassified Lachnospiraceae"             ,
        "Veillonella"     , "Other minor genera"
        )      

## re-label genera for plot
new_genera<- c(
  "Abiotrophia"     , "Acinetobacter"                       ,  "Actinomyces"                              ,
  "Alloprevotella"  , "Bacillales Family XI. Incertae Sedis",  "Bacillus"                                 ,
  "Burkholderia"    , "Corynebacterium"                     ,  "Cronobacter"                              ,
  "Cutibacterium"   , "Dolosicoccus"                        ,  "Eubacteriales Family XIII. Incertae Sedis",
  "Granulicatella"  , "Klebsiella/Raoultella group"         ,  "Lancefieldella"                           ,
  "Lawsonella"      , "Leptotrichia"                        ,  "Levilactobacillus"                        ,
  "Listeria"        , "Mycoplasmopsis"                      ,  "Oribacterium"                             ,
  "Porphyromonas"   , "Prevotella"                          , 
  "Psychrobacter"   , "Rothia"                              ,  "Solobacterium"                            ,
  "Staphylococcus"  , "Streptococcus"                       ,  "unclassified Lachnospiraceae"             ,
  "Veillonella"     , "Other minor genera"
)      

#change to character, then revert to factor afterwards
p$data$Genus <- as.character(p$data$Genus)

for (i in 1:length(old_genera)) {
  #i<-1
  p$data$Genus[which(p$data$Genus == old_genera[i])] <- new_genera[i]
}

bact.genera.top.plot <- p + theme_bw() +
  theme(
    #axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1) 
    axis.text.x  = element_text(angle=60, vjust=1, hjust = 1) 
  )+
  labs(x = NULL, y = "Relative abundance (%)" )

bact.genera.top.plot
# ggsave("Bacterial-relative-abundance-genera.pdf", plot = bact.genera.top.plot)

### Eukaryota -------------
#### Phylum -------------
rel.abun.euka.phy <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel.abun.euka.phy

euka.phylum.plot <- plot_bar(rel.abun.euka.phy , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")#+ 
# theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.pdf", plot = phylum.plot, height = 20, width = 10)

p <- euka.phylum.plot
# str(p$data)

euka.phyla.top.plot <- p + theme_bw() +
  theme(
    #axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1) 
    axis.text.x  = element_text(angle=60, vjust=1, hjust = 1) 
  )+
  labs(x = NULL, y = "Relative abundance (%)" )

euka.phyla.top.plot
# ggsave("Eukaryote-relative-abundance-phyla.pdf", plot = euka.phyla.top.plot)

#### Genus -------------
rel.abun.euka.gen <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.euka.gen

euka.genus.plot <- plot_bar(rel.abun.euka.gen , x = "sample", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free") #+ 
  # theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.pdf", plot = genus.plot, height = 49, width = 10)

# Break down genera into more manageable numbers
p <- euka.genus.plot
# str(p$data)

## Identify rare Genera
hist(p$data$Abundance)
out <- p$data
# str(out)

## Set threshold for major genera
summary(out$Abundance)

## how many Genera at thresholds > X% in any sample ?
#  modify X until number of major genera is manageable, e.g. ~30 or so
X <- 0.31

major_genera.euka <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.euka
#  [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"              "Bacillariophyceae"             
#  [5] "Bathycoccaceae"                 "Chlorellales"                   "Chlorodendrales"                "Chytridiomycota incertae sedis"
#  [9] "Coccidia"                       "CS clade"                       "Deuterostomia"                  "Dictyosteliales"               
# [13] "Haemosporida"                   "Hexamitinae"                    "Kickxellomycotina"              "Mamiellaceae"                  
# [17] "Monothalamids"                  "Mucoromycotina"                 "Oligohymenophorea"              "Pelagomonadales"               
# [21] "Phaeocystaceae"                 "Piroplasmida"                   "Protostomia"                    "saccharomyceta"                
# [25] "Spirotrichea"                   "Symbiodiniaceae"                "Tracheophyta"                   "Trebouxiales"                  
# [29] "Trypanosomatida"                "Ustilaginomycotina"    

sel.row <- which(out$Genus %in% major_genera.euka) # qty 
length(sel.row) #3770

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 99.49871 % of relative abundance is covered
100 - 99.49871 # 6.70154 % left representing rare genera
rel.abun.cov.euka.genera <- 99.49871
perc.rare.euka.genera <- 0.50129

str(out) # do this below if $Genus is a factor
out$Genus <- as.character(out$Genus)
out$Genus[-sel.row] <- "Other minor genera"
out$Genus <- as.factor(out$Genus)

str(p$data)

sel.row.plot <- which(p$data$Genus %in% major_genera.euka) 
str(sel.row.plot)
p$data$Genus <- as.character(p$data$Genus)
p$data$Genus[ -sel.row.plot ] <- "Other minor genera"
p$data$Genus <- as.factor(p$data$Genus)

levels(p$data$Genus)
#  [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"              "Bacillariophyceae"             
#  [5] "Bathycoccaceae"                 "Chlorellales"                   "Chlorodendrales"                "Chytridiomycota incertae sedis"
#  [9] "Coccidia"                       "CS clade"                       "Deuterostomia"                  "Dictyosteliales"               
# [13] "Haemosporida"                   "Hexamitinae"                    "Kickxellomycotina"              "Mamiellaceae"                  
# [17] "Monothalamids"                  "Mucoromycotina"                 "Oligohymenophorea"              "Other minor genera"            
# [21] "Pelagomonadales"                "Phaeocystaceae"                 "Piroplasmida"                   "Protostomia"                   
# [25] "saccharomyceta"                 "Spirotrichea"                   "Symbiodiniaceae"                "Tracheophyta"                  
# [29] "Trebouxiales"                   "Trypanosomatida"                "Ustilaginomycotina"   

p$labels$fill <- "Genera"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor genera to last position, retitle if required
old_genera<- c(
   "Agaricomycotina"     ,           "Anthozoa"             ,          "Apansporoblastina"  ,  "Bacillariophyceae"             ,
   "Bathycoccaceae"      ,           "Chlorellales"         ,          "Chlorodendrales"    ,  "Chytridiomycota incertae sedis",
   "Coccidia"            ,           "CS clade"             ,          "Deuterostomia"      ,  "Dictyosteliales"               ,
   "Haemosporida"        ,           "Hexamitinae"          ,          "Kickxellomycotina"  ,  "Mamiellaceae"                  ,
   "Monothalamids"       ,           "Mucoromycotina"       ,          "Oligohymenophorea"  ,
   "Pelagomonadales"     ,           "Phaeocystaceae"       ,          "Piroplasmida"       ,  "Protostomia"                   ,
   "saccharomyceta"      ,           "Spirotrichea"         ,          "Symbiodiniaceae"    ,  "Tracheophyta"                  ,
   "Trebouxiales"        ,           "Trypanosomatida"      ,          "Ustilaginomycotina" ,  "Other minor genera"
)      

## re-label genera for plot
new_genera <- c(
  "Agaricomycotina"     ,           "Anthozoa"             ,          "Apansporoblastina"  ,  "Bacillariophyceae"             ,
  "Bathycoccaceae"      ,           "Chlorellales"         ,          "Chlorodendrales"    ,  "Chytridiomycota incertae sedis",
  "Coccidia"            ,           "CS clade"             ,          "Deuterostomia"      ,  "Dictyosteliales"               ,
  "Haemosporida"        ,           "Hexamitinae"          ,          "Kickxellomycotina"  ,  "Mamiellaceae"                  ,
  "Monothalamids"       ,           "Mucoromycotina"       ,          "Oligohymenophorea"  ,
  "Pelagomonadales"     ,           "Phaeocystaceae"       ,          "Piroplasmida"       ,  "Protostomia"                   ,
  "saccharomyceta"      ,           "Spirotrichea"         ,          "Symbiodiniaceae"    ,  "Tracheophyta"                  ,
  "Trebouxiales"        ,           "Trypanosomatida"      ,          "Ustilaginomycotina" ,  "Other minor genera"
  )         

#change to character, then revert to factor afterwards
p$data$Genus <- as.character(p$data$Genus)

for (i in 1:length(old_genera)) {
  #i<-1
  p$data$Genus[which(p$data$Genus == old_genera[i])] <- new_genera[i]
}

euka.genera.top.plot <- p + theme_bw() +
  theme(
    #axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1) 
    axis.text.x  = element_text(angle=60, vjust=1, hjust = 1) 
  )+
  labs(x = NULL, y = "Relative abundance (%)" )

euka.genera.top.plot
# ggsave("Eukaryote-relative-abundance-genera.pdf", plot = euka.genera.top.plot)

getwd()
