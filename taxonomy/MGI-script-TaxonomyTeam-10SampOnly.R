#  MGI script
## Taxonomy team
.libPaths()
R.Version() # "R version 4.2.2 (2022-10-31)"
citation()  # R Core Team (2022)

# Read data and libraries ----------
setwd("...") #change to match folder hierarchy in wd
datadir <- "..." # set path for a place to send saved files
library(dplyr);packageVersion("dplyr") # '1.1.0'
library(tidyr);packageVersion("tidyr") # ‘1.3.0’
library(phyloseq);packageVersion("phyloseq") # '1.42.0’
library(ggplot2);packageVersion("ggplot2") # '3.4.1’
library(ggpubr);packageVersion("ggpubr") # '0.6.0’
library(vegan);packageVersion("vegan") # '2.6.4'

# Requires an abundance table, taxonomy table, and metadata table
reads <- read.csv("reads_finalMGI_Sample.csv")
taxonomy <- read.csv("taxonomy_finalMGI_Sample.csv")
metadata <- read.csv("metadata_finalMGI_Sample.csv")

### convert table to tabular split version
taxtable <- taxonomy %>%
  as_tibble() %>%
  separate(taxon, sep=">", c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Subspecies"))

taxtable_m <- as.matrix(taxtable)

# Create phyloseq object ---------
otu <- otu_table(reads, taxa_are_rows = T)
otu <- otu[,-1]
tax <- tax_table(taxtable_m)

otu[is.na(otu)] <- 0

ps <-  phyloseq(otu, tax)
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8176 taxa and 10 samples ]
# tax_table()   Taxonomy Table:    [ 8176 taxa by 9 taxonomic ranks ]

rownames(metadata) <- metadata$sample
meta <- sample_data(metadata)

# sanity check
setdiff(rownames(metadata), (colnames(otu)))

ps.0 <- merge_phyloseq(ps, meta)
ps.0
# otu_table()   OTU Table:         [ 8176 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 8176 taxa by 9 taxonomic ranks ]

## Filter 10% of reads  --------
library(microViz);packageVersion("microViz") # ‘0.10.6’
ps.1 <- tax_filter(ps.0, min_prevalence=0.1)
# Proportional min_prevalence given: 0.1 --> min 1/10 samples.

# ## Normalise: Log  --------------
# Remove zeros and log transform
ps.1.2 <- ps.1
ps.1.2@otu_table[ps.1.2@otu_table == 0] <- NA
ps.1.2@otu_table <- log(ps.1.2@otu_table)+1

# read zeros
ps.1.2@otu_table[is.na(ps.1.2@otu_table)] <- 0

## Ranks and levels ---------------------------
rank_names(ps.1) # [1] "taxID"      "Kingdom"    "Phylum"     "Class"      "Order"      "Family"     "Genus"      "Species"    "Subspecies"
# sort( as.character( unique( tax_table(ps)[, "Phylum"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Class"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Order"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Family"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Genus"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Species"] ) ))
# sort( as.character( unique( tax_table(ps)[, "Subspecies"] ) ))

# ## IF: Rarefaction - ALL HASHED ------------------------------
# rare.object <- ps.1
# 
# min(sample_sums(ps.1) ) # 101145
# min(taxa_sums(ps.1)) # 0
# 
# max(sample_sums(ps.1) ) # 2415533
# max(taxa_sums(ps.1)) # 1485523
# 
# sort(sample_sums(ps.1))
# # X698917_20180128 X698917_20190119 X698917_20171207 X788707_20180313 X788707_20180301 X658355_20171204 X788707_20171213 
# # 101145           121234           226621           251540           522834           777230           809289 
# # X658355_20180321 X788707_20181116 X658355_20180122 
# # 830567           933179          2415533 
# 
# library(vegan)
# # # Rarefaction plot
# # par(xpd=T, mar=par()$mar+c(0,0,0,10))
# # dev.off()
# 
# # ps.1@sam_data$patient.id <- as.character(ps.1@sam_data$patient.id)
# 
# # Plot for rarefaction curve
# rare.plot <- rarecurve(t(otu_table(ps.1)), label = F, ylab = "Number of features", xlab = "Number of reads")
# # legend(locator(1),legend=c("),pt.bg=cols,pch=21,bty="n",ncol=1,cex = 0.75 ,pt.cex = 0.75, title = "Treatment", title.adj = 0)
# # abline(v = 11976, col = "red", lty = 4)
# # dev.off()
# 
# #Rarefaction step
# seed <- 123
# ps.rare <- rarefy_even_depth(ps.1, sample.size = 101145,
#                              rngseed = seed, replace = FALSE, 
#                              trimOTUs = TRUE, verbose = TRUE)
# 
# # `set.seed(123)` was used to initialize repeatable random subsampling.
# # Please record this for your records so others can reproduce.
# # Try `set.seed(123); .Random.seed` for the full vector
# # ...
# # 6424OTUs were removed because they are no longer 
# # present in any sample after random subsampling
# # 
# # ...

# If rarefying - make sure to remove normalisation steps above.
# ps.2 <- ps.rare

# If not rarefying
ps.2 <- ps.1.2 


# Diversity: Shannon diversity indexーChao1--------
ps.2 # gives log normalised data
ps.1 # gives unnormalised data

# obs.ps.raw <- plot_richness(ps.1, measures=c("Observed"))
# obs.ps.raw１
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
obs.ps.1 <- plot_richness(ps.1, measures=c("Observed"))
chao.ps.1 <- plot_richness(ps.1, measures=c("Chao1"))
shan.ps.2 <- plot_richness(ps.2, measures=c("Shannon"))

out.r1 <- data.frame(
  sample=shan.ps.2$data$samples,           # sample ID
  observed=obs.ps.1$data$value,            # On log normalised data
  shannon=shan.ps.2$data$value,            # On log normalised data
  chao1=chao.ps.1$data$value,              # On unnormalised data
  patient.id=shan.ps.2$data$patient.id,    # patient.id
  date=shan.ps.2$data$date                 # sampling date
)

out.r1$patient.id <- as.factor(out.r1$patient.id)
out.r1$eff_no_spp <- exp(out.r1$shannon) # calculate effective no of species
out.r1$date <- as.Date(out.r1$date, format = "%d/%m/%Y")
str(out.r1)
out.r1
#              sample observed  shannon     chao1 patient.id       date eff_no_spp
# 1  X658355_20171204     1045 6.687466 1715.5658     658355 2017-12-04   802.2871
# 2  X658355_20180122     2347 7.519050 3951.1385     658355 2018-01-22  1842.8160
# 3  X658355_20180321      508 5.958327  761.0147     658355 2018-03-21   386.9622
# 4  X698917_20171207      863 6.489932 1440.5231     698917 2017-12-09   658.4784
# 5  X698917_20180128      527 5.991813  913.2881     698917 2018-01-28   400.1394
# 6  X698917_20190119      493 5.936671  768.5185     698917 2019-01-19   378.6723
# 7  X788707_20171213      868 6.521157 1442.4545     788707 2017-12-13   679.3638
# 8  X788707_20180301      767 6.391243 1260.8462     788707 2018-03-01   596.5975
# 9  X788707_20180313      538 6.032534  922.2462     788707 2018-03-12   416.7700
# 10 X788707_20181116     1305 6.921119 2360.3987     788707 2018-11-16  1013.4538

shanon.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=shannon, colour = patient.id))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Shannon's diversity index")+
  xlab("Patient ID")
shanon.div.rare.plot

chao.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=chao1, colour = patient.id))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Chao1 diversity index")+  
  xlab("Patient ID")
chao.div.rare.plot

exp.shanon.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=eff_no_spp, colour = patient.id))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("exp(Shannon's diversity index)")+
  xlab("Patient ID")
exp.shanon.div.rare.plot

library(ggpubr)
shan.chao.div.rare <- ggarrange(shanon.div.rare.plot,chao.div.rare.plot, common.legend = TRUE)
shan.chao.div.rare
# ggsave("Diversity-shannon-chao-rare.MGI.pdf", plot =shan.chao.div.rare, path = datadir)

# Chao1 over time
chao1.time.plot <- ggplot(out.r1, aes(x= date, y=chao1, colour = patient.id))+
  geom_point(size=2.5)+
  facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chao1.time.plot 

chao1.time.plot2 <- ggplot(out.r1, aes(x= date, y=chao1, colour = patient.id))+
  geom_point(size=2.5)+
  # facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chao1.time.plot2

shan.time.plot <- ggplot(out.r1, aes(x= date, y=shannon, colour = patient.id))+
  geom_point(size=2.5)+
  facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
shan.time.plot 

shan.time.plot2 <- ggplot(out.r1, aes(x= date, y=shannon, colour = patient.id))+
  geom_point(size=2.5)+
  # facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
shan.time.plot2 

shan.chao.time.rare <- ggarrange(shan.time.plot,chao1.time.plot, common.legend = TRUE, ncol = 1, align = "hv")
shan.chao.time.rare
# ggsave("Diversity-time-shannon-chao-rare.MGI.pdf", plot =shan.chao.time.rare, path = datadir)

# Ordination on normalised data -----------
# All taxa in samples
ps.2

set.seed(123)
PCoA_ord <- ordinate(ps.2, method="PCoA", distance="bray", weighted=TRUE)
PCoA_all.plot <- plot_ordination(ps.2, PCoA_ord, label = "date")

PCoA_all <- PCoA_all.plot + 
  geom_point(aes(fill = as.character(patient.id)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all
# ggsave("Ordination_PCoA-bray-curtis_patients-all-data.MinION.pdf", plot =PCoA_all, height = 3.5, width = 5.2, path = datadir)

PCoA_ord
PCoA_ord$vectors

# Bacteria only
levels(factor(tax_table(ps.2)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(ps.2)[, "Kingdom"] %in%  c("Archaea","cellular organisms", "Eukaryota", "root", "unclassified", "Viruses") )
bacteria.psq <- prune_taxa(ps.2, taxa = row.names(tax_table(ps.2)[-rem_taxa, ]) )
bacteria.psq

bacteria.psq@sam_data$patient.id <- as.character(bacteria.psq@sam_data$patient.id)
bacteria.psq@sam_data$date <- as.character(bacteria.psq@sam_data$date)

set.seed(123)
PCoA_ord.bact <- ordinate(bacteria.psq, method="PCoA", distance="bray", weighted=TRUE)
PCoA_bact.plot <- plot_ordination(bacteria.psq, PCoA_ord.bact, label = "date")

PCoA_bact <- PCoA_bact.plot + 
  geom_point(aes(fill = as.character(patient.id)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_bact
# ggsave("Ordination_PCoA-bray-curtis_patients-bacterialonly-data.MinION.pdf", plot =PCoA_bact, height = 3.5, width = 5.2, path = datadir)

PCoA_ord.bact
PCoA_ord.bact$vectors
dev.off()

# Eukaryotes only
levels(factor(tax_table(ps.2)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(ps.2)[, "Kingdom"] %in%  c("Archaea","cellular organisms", "Bacteria", "root", "unclassified", "Viruses") )
eukaryota.psq <- prune_taxa(ps.2, taxa = row.names(tax_table(ps.2)[-rem_taxa, ]) )
eukaryota.psq

eukaryota.psq@sam_data$patient.id <- as.character(eukaryota.psq@sam_data$patient.id)
eukaryota.psq@sam_data$date <- as.character(eukaryota.psq@sam_data$date)

set.seed(123)
PCoA_ord.euka <- ordinate(eukaryota.psq, method="PCoA", distance="bray", weighted=TRUE)
PCoA_euka.plot <- plot_ordination(eukaryota.psq, PCoA_ord.euka, label = "date")

PCoA_euka <- PCoA_euka.plot + 
  geom_point(aes(fill = as.character(patient.id)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_euka
# ggsave("Ordination_PCoA-bray-curtis_patients-eukaryotesonly-data.MinION.pdf", plot =PCoA_bact, height = 3.5, width = 5.2, path = datadir)

PCoA_ord.euka
PCoA_ord.euka$vectors
dev.off()

## Statistics ----------
library(car);packageVersion("car") #‘3.1.1’
par(mfrow=c(2,2))

#Shannon's
aov.patient.sh <- aov(shannon~patient.id, out.r1)

plot(aov.patient.sh)
shapiro.test(resid(aov.patient.sh)) #W = 0.98192, p-value = 0.9746 Non-normal
leveneTest(shannon~patient.id, out.r1) # F= 0.9534, Pr(>F) = 0.4303
summary(aov.patient.sh)
#             Df Sum Sq Mean Sq F value Pr(>F)
# patient.id   2 0.5114  0.2557   0.989  0.418
# Residuals    7 1.8091  0.2584  

#Chao1
aov.patient.ch <- aov(chao1~patient.id, out.r1)

plot(aov.patient.ch)
shapiro.test(resid(aov.patient.ch)) #W = 0.93506, p-value = 0.4994
leveneTest(chao1~patient.id, out.r1) # F= 1.2944, Pr(>F) = 0.3324
summary(aov.patient.ch)
#              Df Sum Sq Mean Sq F value Pr(>F)
# patient.id   2 1842679  921340   0.956  0.429
# Residuals    7 6746665  963809      

# Test hypothesis that microbiota vary (with different centroids) by patient
# Calculate bray-curtis distance matrix - ALL
set.seed(123)
bray.rare <- phyloseq::distance(ps.2, method = "bray")
sam.df <- data.frame(sample_data(ps.2))

str(sam.df)
dev.off()

set.seed(123)
adonis2(bray.rare ~ patient.id, data = sam.df)  
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = bray.rare ~ patient.id, data = sam.df)
#            Df SumOfSqs     R2      F Pr(>F)   
# patient.id  1  0.13942 0.10206 0.9093  0.492
# Residual    8  1.22660 0.89794              
# Total       9  1.36602 1.00000                  

#Test for beta dispersion
beta <- betadisper(bray.rare, sam.df$patient.id)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#            Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.005967 0.0029833 0.2579    999  0.802
# Residuals  7 0.080983 0.0115690          

# Calculate bray-curtis distance matrix - Bacteria only
set.seed(123)
bray.rare <- phyloseq::distance(bacteria.psq, method = "bray")
sam.df <- data.frame(sample_data(bacteria.psq))

str(sam.df)
# dev.off()

set.seed(123)
adonis2(bray.rare ~ patient.id, data = sam.df)  
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = bray.rare ~ patient.id, data = sam.df)
#            Df SumOfSqs     R2      F Pr(>F)   
# patient.id  2  0.54088 0.35241 1.9046  0.026 *
# Residual    7  0.99394 0.64759                
# Total       9  1.53482 1.00000         

#Test for beta dispersion
beta <- betadisper(bray.rare, sam.df$patient.id)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     2 0.006596 0.0032978 0.2555    999  0.805
# Residuals  7 0.090343 0.0129062       

# Calculate bray-curtis distance matrix - Eukaryotes only
set.seed(123)
bray.rare <- phyloseq::distance(eukaryota.psq, method = "bray")
sam.df <- data.frame(sample_data(eukaryota.psq))

str(sam.df)
# dev.off()

set.seed(123)
adonis2(bray.rare ~ patient.id, data = sam.df)  
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = bray.rare ~ patient.id, data = sam.df)
#            Df SumOfSqs      R2      F Pr(>F)
# patient.id  2  0.30865 0.33628 1.7733   0.11
# Residual    7  0.60919 0.66372              
# Total       9  0.91785 1.00000  

#Test for beta dispersion
beta <- betadisper(bray.rare, sam.df$patient.id)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     2 0.014273 0.0071367 0.7838    999  0.493
# Residuals  7 0.063739 0.0091055   

# Relative abundance plots -------------------------
ps.1 # unnormalised
# otu_table()   OTU Table:         [ 3647 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3647 taxa by 9 taxonomic ranks ]

ps.2 # normalised
# otu_table()   OTU Table:         [ 3647 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3647 taxa by 9 taxonomic ranks ]

rank_names(ps.2) # "taxID"      "Kingdom"    "Phylum"     "Class"      "Order"      "Family"     "Genus"      "Species"    "Subspecies"

(all.K.psq <- ps.2)
# otu_table()   OTU Table:         [ 3647 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3647 taxa by 9 taxonomic ranks ]

## remove taxa not assigned as Bacteria
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","cellular organisms", "Eukaryota", "root", "unclassified", "Viruses") )
bacteria.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
bacteria.psq
# otu_table()   OTU Table:         [ 2948 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 2948 taxa by 9 taxonomic ranks ]

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
# otu_table()   OTU Table:         [ 536 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 536 taxa by 9 taxonomic ranks ]

## Relative abundance plots All divisions in one: -------------------------
rel_abun.Phylum.phy_obj <- transform_sample_counts(tax_glom(ps.2, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel_abun.Phylum.phy_obj
# otu_table()   OTU Table:         [ 3614 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3614 taxa by 9 taxonomic ranks ]

rel_abun.Genus.phy_obj <- transform_sample_counts(tax_glom(ps.2, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel_abun.Genus.phy_obj
# otu_table()   OTU Table:         [ 3052 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3052 taxa by 9 taxonomic ranks ]

# Identify significantly trending phyla and those with substantial proportion
phylum.plot <- plot_bar(rel_abun.Phylum.phy_obj , x = "date", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")
# phylum.plot <- phylum.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.MGI.pdf", plot = phylum.plot, height = 20, width = 10, path = datadir)

genus.plot <- plot_bar(rel_abun.Genus.phy_obj , x = "date", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.MGI.pdf", plot = genus.plot, height = 49, width = 10, path = datadir)

## Relative abundance plots top taxa: -------------------------

### Bacteria -------------
#### Phylum -------------
rel.abun.bact.phy <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel.abun.bact.phy

bact.phylum.plot <- plot_bar(rel.abun.bact.phy , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")#+ 
# theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.MGI.pdf", plot = phylum.plot, height = 20, width = 10)

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
X <- 0.099

major_phyla.bact <- levels(as.factor(as.character( out$Phylum[which(out$Abun >= X )]  )))
major_phyla.bact
# [1]  [1] "Acidobacteria"           "Bacteria_incertae_sedis" "Elusimicrobia"           "environmental_samples"  
# [5] "FCB_group"               "Fusobacteria"            "Pseudomonadota"          "PVC_group"              
# [9] "Spirochaetes"            "Synergistetes"           "Terrabacteria_group"     "unclassified_Bacteria"  

sel.row <- which(out$Phylum %in% major_phyla.bact) # qty 
length(sel.row) #29060

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 99.87422 % of relative abundance is covered
100 - 99.87422 # 0.12578 % left representing rare phyla
rel.abun.cov.bact.phyla <- 99.87422
perc.rare.bact.phyla <- 0.12578

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
# [1] "Acidobacteria"           "Bacteria_incertae_sedis" "Elusimicrobia"           "environmental_samples"  
# [5] "FCB_group"               "Fusobacteria"            "Other minor phyla"       "Pseudomonadota"         
# [9] "PVC_group"               "Spirochaetes"            "Synergistetes"           "Terrabacteria_group"    
# [13] "unclassified_Bacteria"  

p$labels$fill <- "Phyla"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor phyla to last position, retitle if required
old_phyla<- c(
   "Acidobacteria"         ,  "Bacteria_incertae_sedis", "Elusimicrobia"      ,     "environmental_samples" , 
   "FCB_group"             ,  "Fusobacteria"           ,     "Pseudomonadota"        , 
   "PVC_group"             ,  "Spirochaetes"           , "Synergistetes"      ,     "Terrabacteria_group"   , 
   "unclassified_Bacteria" ,  "Other minor phyla"
   )

## re-label phyla for plot
new_phyla<- c(
  "Acidobacteria"         ,  "Bacteria_incertae_sedis", "Elusimicrobia"      ,     "environmental_samples" , 
  "FCB_group"             ,  "Fusobacteria"           ,     "Pseudomonadota"        , 
  "PVC_group"             ,  "Spirochaetes"           , "Synergistetes"      ,     "Terrabacteria_group"   , 
  "unclassified_Bacteria" ,  "Other minor phyla"
)
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
# ggsave("Bacterial-relative-abundance-phyla.MGI.pdf", plot = bact.phyla.top.plot, path = datadir)

#### Genus -------------
rel.abun.bact.gen <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.bact.gen

bact.genus.plot <- plot_bar(rel.abun.bact.gen , x = "date", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.MGI.pdf", plot = genus.plot, height = 49, width = 10, path = datadir)

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
X <- 0.8

major_genera.bact <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.bact
#  [1] "Acinetobacter"                "Aerococcaceae"                "Bacillaceae"                 
#  [4] "Bacillales_incertae_sedis"    "Bacteroidales"                "Carnobacteriaceae"           
#  [7] "Corynebacteriaceae"           "Cronobacter"                  "Curvibacter"                 
# [10] "Enterococcaceae"              "Escherichia"                  "Eubacteriales_incertae_sedis"
# [13] "Flavobacteriales"             "Klebsiella/Raoultella_group"  "Lachnospiraceae"             
# [16] "Lactobacillaceae"             "Lawsonellaceae"               "Listeriaceae"                
# [19] "Massilia_group"               "Micrococcaceae"               "Mycoplasmataceae"            
# [22] "Peptostreptococcaceae"        "Propionibacteriaceae"         "Pseudomonas"                 
# [25] "Psychrobacter"                "Salmonella"                   "Staphylococcaceae"           
# [28] "Streptococcaceae"             "unclassified_Eubacteriales"   "unclassified_Sinobacteraceae"
# [31] "Veillonellaceae"   

sel.row <- which(out$Genus %in% major_genera.bact) # qty 
length(sel.row) #14230

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 78.29535 % of relative abundance is covered
100 - 78.29535 # 21.70465 % left representing rare genera
rel.abun.cov.bact.genera <- 78.29535
perc.rare.bact.genera <- 21.70465

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
#  [1] "Acinetobacter"                "Aerococcaceae"                "Bacillaceae"                 
#  [4] "Bacillales_incertae_sedis"    "Bacteroidales"                "Carnobacteriaceae"           
#  [7] "Corynebacteriaceae"           "Cronobacter"                  "Curvibacter"                 
# [10] "Enterococcaceae"              "Escherichia"                  "Eubacteriales_incertae_sedis"
# [13] "Flavobacteriales"             "Klebsiella/Raoultella_group"  "Lachnospiraceae"             
# [16] "Lactobacillaceae"             "Lawsonellaceae"               "Listeriaceae"                
# [19] "Massilia_group"               "Micrococcaceae"               "Mycoplasmataceae"            
# [22] "Other minor genera"           "Peptostreptococcaceae"        "Propionibacteriaceae"        
# [25] "Pseudomonas"                  "Psychrobacter"                "Salmonella"                  
# [28] "Staphylococcaceae"            "Streptococcaceae"             "unclassified_Eubacteriales"  
# [31] "unclassified_Sinobacteraceae" "Veillonellaceae"                          

p$labels$fill <- "Genera"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor genera to last position, retitle if required
old_genera<- c(
     "Acinetobacter"               , "Aerococcaceae"               , "Bacillaceae"                 ,
     "Bacillales_incertae_sedis"   , "Bacteroidales"               , "Carnobacteriaceae"           ,
     "Corynebacteriaceae"          , "Cronobacter"                 , "Curvibacter"                 ,
     "Enterococcaceae"             , "Escherichia"                 , "Eubacteriales_incertae_sedis",
     "Flavobacteriales"            , "Klebsiella/Raoultella_group" , "Lachnospiraceae"             ,
     "Lactobacillaceae"            , "Lawsonellaceae"              , "Listeriaceae"                ,
     "Massilia_group"              , "Micrococcaceae"              , "Mycoplasmataceae"            ,
     "Peptostreptococcaceae"       , "Propionibacteriaceae"        ,
     "Pseudomonas"                 , "Psychrobacter"               , "Salmonella"                  ,
     "Staphylococcaceae"           , "Streptococcaceae"            , "unclassified_Eubacteriales"  ,
     "unclassified_Sinobacteraceae", "Veillonellaceae"  ,"Other minor genera")

 
 
 ## re-label genera for plot
new_genera<- c(
  "Acinetobacter"               , "Aerococcaceae"               , "Bacillaceae"                 ,
  "Bacillales_incertae_sedis"   , "Bacteroidales"               , "Carnobacteriaceae"           ,
  "Corynebacteriaceae"          , "Cronobacter"                 , "Curvibacter"                 ,
  "Enterococcaceae"             , "Escherichia"                 , "Eubacteriales_incertae_sedis",
  "Flavobacteriales"            , "Klebsiella/Raoultella_group" , "Lachnospiraceae"             ,
  "Lactobacillaceae"            , "Lawsonellaceae"              , "Listeriaceae"                ,
  "Massilia_group"              , "Micrococcaceae"              , "Mycoplasmataceae"            ,
  "Peptostreptococcaceae"       , "Propionibacteriaceae"        ,
  "Pseudomonas"                 , "Psychrobacter"               , "Salmonella"                  ,
  "Staphylococcaceae"           , "Streptococcaceae"            , "unclassified_Eubacteriales"  ,
  "unclassified_Sinobacteraceae", "Veillonellaceae"  ,"Other minor genera")

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
# ggsave("Bacterial-relative-abundance-genera.MGI.pdf", plot = bact.genera.top.plot, path = datadir)

### Eukaryota -------------
#### Phylum -------------
rel.abun.euka.phy <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel.abun.euka.phy

euka.phylum.plot <- plot_bar(rel.abun.euka.phy , x = "date", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")#+ 
# theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.MGI.pdf", plot = phylum.plot, height = 20, width = 10, path = datadir)

p <- euka.phylum.plot
# str(p$data)

euka.phyla.top.plot <- p + theme_bw() +
  theme(
    #axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1) 
    axis.text.x  = element_text(angle=60, vjust=1, hjust = 1) 
  )+
  labs(x = NULL, y = "Relative abundance (%)" )

euka.phyla.top.plot
# ggsave("Eukaryote-relative-abundance-phyla.MGI.pdf", plot = euka.phyla.top.plot, path = datadir)

#### Genus -------------
rel.abun.euka.gen <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.euka.gen

euka.genus.plot <- plot_bar(rel.abun.euka.gen , x = "date", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free") #+ 
  # theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.MGI.pdf", plot = genus.plot, height = 49, width = 10)

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
#  modify X until number of major genera is manageable, e.g. ~10 or so
X <- 0.3

major_genera.euka <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.euka
#   [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"             
# [4] "Bacillariophyceae"              "Bathycoccaceae"                 "Charales"                      
# [7] "Chlorellales"                   "Chlorodendrales"                "Chytridiomycota_incertae_sedis"
# [10] "Coccidia"                       "CS_clade"                       "Deuterostomia"                 
# [13] "Dictyosteliales"                "Haemosporida"                   "Hexamitinae"                   
# [16] "Mamiellaceae"                   "Monothalamids"                  "Mucoromycotina"                
# [19] "Oligohymenophorea"              "Pelagomonadales"                "Peronosporaceae"               
# [22] "Phaeocystaceae"                 "Piroplasmida"                   "Protostomia"                   
# [25] "saccharomyceta"                 "Spirotrichea"                   "Symbiodiniaceae"               
# [28] "Tracheophyta"                   "Trypanosomatida"                "Ustilaginomycotina"            

sel.row <- which(out$Genus %in% major_genera.euka) # qty 
length(sel.row) #4630

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 99.25094 % of relative abundance is covered
100 - 99.25094 # 0.74906 % left representing rare genera
rel.abun.cov.euka.genera <- 99.25094
perc.rare.euka.genera <- 0.74906

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
#  [5] "Bathycoccaceae"                 "Chlorellales"                   "Chlorodendrales"                "Chytridiomycota_incertae_sedis"
#  [9] "Coccidia"                       "CS_clade"                       "Deuterostomia"                  "Dictyosteliales"               
# [13] "Haemosporida"                   "Hexamitinae"                    "Mamiellaceae"                   "Monothalamids"                 
# [17] "Mucoromycotina"                 "Oligohymenophorea"              "Other minor genera"             "Pelagomonadales"               
# [21] "Peronosporaceae"                "Phaeocystaceae"                 "Piroplasmida"                   "Protostomia"                   
# [25] "saccharomyceta"                 "Spirotrichea"                   "Symbiodiniaceae"                "Tracheophyta"                  
# [29] "Trypanosomatida"                "Ustilaginomycotina"   

p$labels$fill <- "Genera"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor genera to last position, retitle if required
old_genera<- c(
  "Agaricomycotina"   ,    "Anthozoa"          ,  "Apansporoblastina"             ,
  "Bacillariophyceae" ,    "Bathycoccaceae"    ,  "Charales"                      ,
  "Chlorellales"      ,    "Chlorodendrales"   ,  "Chytridiomycota_incertae_sedis",
  "Coccidia"          ,    "CS_clade"          ,  "Deuterostomia"                 ,
  "Dictyosteliales"   ,    "Haemosporida"      ,  "Hexamitinae"                   ,
  "Mamiellaceae"      ,    "Monothalamids"     ,  "Mucoromycotina"                ,
  "Oligohymenophorea" ,    "Pelagomonadales"   ,
  "Peronosporaceae"   ,    "Phaeocystaceae"    ,  "Piroplasmida"                  ,
  "Protostomia"       ,    "saccharomyceta"    ,  "Spirotrichea"                  ,
  "Symbiodiniaceae"   ,    "Tracheophyta"      ,  "Trypanosomatida"               ,
  "Ustilaginomycotina",    "Other minor genera"
  )

## re-label genera for plot
new_genera <- c(
  "Agaricomycotina"   ,    "Anthozoa"          ,  "Apansporoblastina"             ,
  "Bacillariophyceae" ,    "Bathycoccaceae"    ,  "Charales"                      ,
  "Chlorellales"      ,    "Chlorodendrales"   ,  "Chytridiomycota_incertae_sedis",
  "Coccidia"          ,    "CS_clade"          ,  "Deuterostomia"                 ,
  "Dictyosteliales"   ,    "Haemosporida"      ,  "Hexamitinae"                   ,
  "Mamiellaceae"      ,    "Monothalamids"     ,  "Mucoromycotina"                ,
  "Oligohymenophorea" ,    "Pelagomonadales"   ,
  "Peronosporaceae"   ,    "Phaeocystaceae"    ,  "Piroplasmida"                  ,
  "Protostomia"       ,    "saccharomyceta"    ,  "Spirotrichea"                  ,
  "Symbiodiniaceae"   ,    "Tracheophyta"      ,  "Trypanosomatida"               ,
  "Ustilaginomycotina",    "Other minor genera"
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
# ggsave("Eukaryote-relative-abundance-genera.MGI.pdf", plot = euka.genera.top.plot, path = datadir)

# getwd()
