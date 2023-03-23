# CF Kraken output -> Pavian -> Phyloseq
# MinION data (10 samples)
## Taxonomy team

# Provides diversity of microbial communities, microbial community composition, and relativa abundance plots to show major phyla and genera

.libPaths()
R.Version() # "R version 4.2.2 (2022-10-31)"
citation()  # R Core Team (2022)

# Read data and libraries ----------
setwd("~/Documents/PhD/Collaborations/CF Hackathon/PseqFiles") #change to match folder hierarchy in wd

library(dplyr);packageVersion("dplyr") # '1.1.0'
library(tidyr);packageVersion("tidyr") # ‘1.3.0’
library(phyloseq);packageVersion("phyloseq") # '1.42.0’
library(ggplot2);packageVersion("ggplot2") # '3.4.1’
library(ggpubr);packageVersion("ggpubr") # '0.6.0’

# Requires an abundance table, taxonomy table, and metadata table
reads <- read.csv("reads.csv")
taxonomy <- read.csv("taxonomy.csv")
metadata <- read.csv("metadata2.csv")

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
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

rownames(metadata) <- metadata$sample
meta <- sample_data(metadata)

# sanity check
setdiff(rownames(metadata), (colnames(otu)))

ps.0 <- merge_phyloseq(ps, meta)
ps.0
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

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
ps.2 # gives log rarefied data
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
#                              sample observed  shannon     chao1 patient.id       date eff_no_spp
# 1  X658355_20171204_S_Run6barcode11      781 6.408909 1404.8938     658355 2017-12-04   607.2310
# 2  X658355_20180122_S_Run8barcode05     1901 7.315678 3638.5197     658355 2018-01-22  1503.6917
# 3  X658355_20180321_S_Run1Barcode01      322 5.507915  547.2128     658355 2018-03-21   246.6363
# 4  X698917_20171207_S_Run7barcode10      612 6.169140 1208.2619     698917 2017-12-09   477.7748
# 5  X698917_20180128_S_Run4barcode03      324 5.528412  697.5789     698917 2018-01-28   251.7439
# 6  X698917_20190119_S_Run7barcode05      297 5.457093  528.0000     698917 2019-01-19   234.4150
# 7  X788707_20171213_S_Run9barcode04      586 6.143228 1038.9765     788707 2017-12-13   465.5539
# 8  X788707_20180301_S_Run2barcode07      522 6.013915 1044.1587     788707 2018-03-01   409.0817
# 9  X788707_20180313_S_Run2barcode10      338 5.561981  680.3750     788707 2018-03-12   260.3380
# 10 X788707_20181116_S_Run8barcode08      970 6.640759 2028.8783     788707 2018-11-16   765.6762

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
# ggsave("Diversity-shannon-chao-rare.MinION.pdf", plot =shan.chao.div.rare)

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
# ggsave("Diversity-time-shannon-chao-rare.MinION.pdf", plot =shan.chao.time.rare)

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
ggsave("Ordination_PCoA-bray-curtis_patients-all-data.MinION.pdf", plot =PCoA_all, height = 3.5, width = 5.2)

PCoA_ord
PCoA_ord$vectors

# Bacteria only
levels(factor(tax_table(ps.2)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(ps.2)[, "Kingdom"] %in%  c("Archaea","cellular organisms", "Eukaryota", "root", "unclassified", "Viruses") )
bacteria.psq <- prune_taxa(ps.2, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
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
# ggsave("Ordination_PCoA-bray-curtis_patients-bacterialonly-data.MinION.pdf", plot =PCoA_bact, height = 3.5, width = 5.2)

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
# ggsave("Ordination_PCoA-bray-curtis_patients-eukaryotesonly-data.MinION.pdf", plot =PCoA_bact, height = 3.5, width = 5.2)

PCoA_ord.euka
PCoA_ord.euka$vectors
dev.off()

## Statistics ----------
library(car);packageVersion("car") #‘3.1.1’
par(mfrow=c(2,2))

#Shannon's
aov.patient.sh <- aov(shannon~patient.id, out.r1)

plot(aov.patient.sh)
shapiro.test(resid(aov.patient.sh)) #W = 0.98266, p-value = 0.9778 Non-normal
leveneTest(shannon~patient.id, out.r1) # F= 0.803, Pr(>F) = 0.4853
summary(aov.patient.sh)
#             Df Sum Sq Mean Sq F value Pr(>F)
# patient.id   2 0.7211  0.3606   0.997  0.416
# Residuals    7 2.5323  0.3618 

#Chao1
aov.patient.ch <- aov(chao1~patient.id, out.r1)

plot(aov.patient.ch)
shapiro.test(resid(aov.patient.ch)) #W = 0.92323, p-value = 0.3847
leveneTest(chao1~patient.id, out.r1) # F= 1.2338, Pr(>F) = 0.3475
summary(aov.patient.ch)
#              Df Sum Sq Mean Sq F value Pr(>F)
# patient.id   2 1707463  853731   0.941  0.435
# Residuals    7 6351647  907378     

# Test hypothesis that microbiota vary (with different centroids) by patient
# Calculate bray-curtis distance matrix
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
# patient.id  1  0.23161 0.11514 1.041  0.401
# Residual    8  1.77988 0.88486             
# Total       9  2.01148 1.00000       

#Test for beta dispersion
beta <- betadisper(bray.rare, sam.df$patient.id)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#            Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.005901 0.0029503 0.2148    999  0.825
# Residuals  7 0.096158 0.0137369     

# Relative abundance plots -------------------------
ps.1 # unnormalised
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

ps.2 # normalised
# otu_table()   OTU Table:         [ 3197 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3197 taxa by 9 taxonomic ranks ]

rank_names(ps.2) # "taxID"      "Kingdom"    "Phylum"     "Class"      "Order"      "Family"     "Genus"      "Species"    "Subspecies"

(all.K.psq <- ps.2)
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
rel_abun.Phylum.phy_obj <- transform_sample_counts(tax_glom(ps.2, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel_abun.Phylum.phy_obj
# otu_table()   OTU Table:         [ 3175 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 3175 taxa by 9 taxonomic ranks ]

rel_abun.Genus.phy_obj <- transform_sample_counts(tax_glom(ps.2, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel_abun.Genus.phy_obj
# otu_table()   OTU Table:         [ 2558 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 2558 taxa by 9 taxonomic ranks ]

# Identify significantly trending phyla and those with substantial proportion
phylum.plot <- plot_bar(rel_abun.Phylum.phy_obj , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# phylum.plot <- phylum.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.MinION.pdf", plot = phylum.plot, height = 20, width = 10)

genus.plot <- plot_bar(rel_abun.Genus.phy_obj , x = "sample", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.MinION.pdf", plot = genus.plot, height = 49, width = 10)

## Relative abundance plots top taxa: -------------------------

### Bacteria -------------
#### Phylum -------------
rel.abun.bact.phy <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel.abun.bact.phy

bact.phylum.plot <- plot_bar(rel.abun.bact.phy , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")#+ 
# theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.MinION.pdf", plot = phylum.plot, height = 20, width = 10)

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
X <- 0.8

major_phyla.bact <- levels(as.factor(as.character( out$Phylum[which(out$Abun >= X )]  )))
major_phyla.bact
# [1] " Bacteroidota"                "Actinomycetota"               "Bacillota"                   
# [4] "Bacteroidetes/Chlorobi group" "Fusobacteria"                 "Pseudomonadota"              
# [7] "Tenericutes"                  "Terrabacteria"                "Terrabacteria group"

sel.row <- which(out$Phylum %in% major_phyla.bact) # qty 
length(sel.row) #64080

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 96.64789 % of relative abundance is covered
100 - 96.64789 # 3.35211 % left representing rare phyla
rel.abun.cov.bact.phyla <- 96.64789
perc.rare.bact.phyla <- 3.35211

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
#  [1] " Bacteroidota"                "Actinomycetota"               "Bacillota"                   
#  [4] "Bacteroidetes/Chlorobi group" "Fusobacteria"                 "Other minor phyla"           
#  [7] "Pseudomonadota"               "Tenericutes"                  "Terrabacteria"               
# [10] "Terrabacteria group"   

p$labels$fill <- "Phyla"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor phyla to last position, retitle if required
old_phyla<- c(
  " Bacteroidota"               , "Actinomycetota", "Bacillota"       ,            
  "Bacteroidetes/Chlorobi group", "Fusobacteria",          
  "Pseudomonadota"              , "Tenericutes" , "Terrabacteria"     ,          
  "Terrabacteria group" ,  "Other minor phyla" 
)

## re-label phyla for plot
new_phyla<- c(
  "Bacteroidota"               , "Actinomycetota", "Bacillota"       ,            
  "Bacteroidetes/Chlorobi group", "Fusobacteria",          
  "Pseudomonadota"              , "Tenericutes" , "Terrabacteria"     ,          
  "Terrabacteria group" ,  "Other minor phyla" 
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
# ggsave("Bacterial-relative-abundance-phyla-lognorm.MinION.pdf", plot = bact.phyla.top.plot)

#### Genus -------------
rel.abun.bact.gen <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.bact.gen

bact.genus.plot <- plot_bar(rel.abun.bact.gen , x = "sample", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")+ 
  theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.MinION.pdf", plot = genus.plot, height = 49, width = 10)

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
X <- 1.19

major_genera.bact <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.bact
# [1] "Abiotrophia"                          "Acinetobacter"                        "Bacillales Family XI. Incertae Sedis"
# [4] "Bacillus"                             "Clostridioides"                       "Corynebacterium"                     
# [7] "Cronobacter"                          "Curvibacter"                          "Cutibacterium"                       
# [10] "Enterococcus"                         "Escherichia"                          "Klebsiella/Raoultella group"         
# [13] "Lawsonella"                           "Levilactobacillus"                    "Listeria"                            
# [16] "Massilia group"                       "Mycobacterium"                        "Mycoplasmopsis"                      
# [19] "Oribacterium"                         "Porphyromonas"                        "Prevotella"                          
# [22] "Pseudoleptotrichia"                   "Pseudomonas"                          "Psychrobacter"                       
# [25] "Salmonella"                           "Staphylococcus"                       "Streptococcus"                       
# [28] "unclassified Sinobacteraceae"         "unclassified Spongiibacteraceae"  

sel.row <- which(out$Genus %in% major_genera.bact) # qty 
length(sel.row) #28340

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 58.32988 % of relative abundance is covered
100 - 58.32988 # 41.67012 % left representing rare genera
rel.abun.cov.bact.genera <- 58.32988
perc.rare.bact.genera <- 41.67012

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
#  [1] "Abiotrophia"                          "Acinetobacter"                        "Bacillales Family XI. Incertae Sedis"
#  [4] "Bacillus"                             "Clostridioides"                       "Corynebacterium"                     
#  [7] "Cronobacter"                          "Curvibacter"                          "Cutibacterium"                       
# [10] "Enterococcus"                         "Escherichia"                          "Klebsiella/Raoultella group"         
# [13] "Lawsonella"                           "Levilactobacillus"                    "Listeria"                            
# [16] "Massilia group"                       "Mycobacterium"                        "Mycoplasmopsis"                      
# [19] "Oribacterium"                         "Other minor genera"                   "Porphyromonas"                       
# [22] "Prevotella"                           "Pseudoleptotrichia"                   "Pseudomonas"                         
# [25] "Psychrobacter"                        "Salmonella"                           "Staphylococcus"                      
# [28] "Streptococcus"                        "unclassified Sinobacteraceae"         "unclassified Spongiibacteraceae"                               

p$labels$fill <- "Genera"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor genera to last position, retitle if required
old_genera<- c(
  "Abiotrophia"    , "Acinetobacter"                 ,     "Bacillales Family XI. Incertae Sedis",
  "Bacillus"       , "Clostridioides"                ,     "Corynebacterium"                     ,
  "Cronobacter"    , "Curvibacter"                   ,     "Cutibacterium"                       ,
  "Enterococcus"  ,  "Escherichia"                  ,      "Klebsiella/Raoultella group"         ,
  "Lawsonella"    ,  "Levilactobacillus"            ,      "Listeria"                            ,
  "Massilia group",  "Mycobacterium"                ,      "Mycoplasmopsis"                      ,
  "Oribacterium"  ,  "Other minor genera"           ,      "Porphyromonas"                       ,
  "Prevotella"    ,  "Pseudoleptotrichia"           ,      "Pseudomonas"                         ,
  "Psychrobacter" ,  "Salmonella"                   ,      "Staphylococcus"                      ,
  "Streptococcus" ,  "unclassified Sinobacteraceae" ,      "unclassified Spongiibacteraceae"    , "Other minor genera"
)
## re-label genera for plot
new_genera<- c(
  "Abiotrophia"    , "Acinetobacter"                 ,     "Bacillales Family XI. Incertae Sedis",
  "Bacillus"       , "Clostridioides"                ,     "Corynebacterium"                     ,
  "Cronobacter"    , "Curvibacter"                   ,     "Cutibacterium"                       ,
  "Enterococcus"  ,  "Escherichia"                  ,      "Klebsiella/Raoultella group"         ,
  "Lawsonella"    ,  "Levilactobacillus"            ,      "Listeria"                            ,
  "Massilia group",  "Mycobacterium"                ,      "Mycoplasmopsis"                      ,
  "Oribacterium"  ,  "Other minor genera"           ,      "Porphyromonas"                       ,
  "Prevotella"    ,  "Pseudoleptotrichia"           ,      "Pseudomonas"                         ,
  "Psychrobacter" ,  "Salmonella"                   ,      "Staphylococcus"                      ,
  "Streptococcus" ,  "unclassified Sinobacteraceae" ,      "unclassified Spongiibacteraceae"    , "Other minor genera"
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
ggsave("Bacterial-relative-abundance-genera.MinION.pdf", plot = bact.genera.top.plot, width = 10)

### Eukaryota -------------
#### Phylum -------------
rel.abun.euka.phy <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Phylum" ), function(x) 100*x / sum(x))
rel.abun.euka.phy

euka.phylum.plot <- plot_bar(rel.abun.euka.phy , x = "sample", fill = "Phylum")+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free")#+ 
# theme(legend.position = "bottom")
# ggsave("relabund.all.phylum.plot.MinION.pdf", plot = phylum.plot, height = 20, width = 10)

p <- euka.phylum.plot
# str(p$data)

euka.phyla.top.plot <- p + theme_bw() +
  theme(
    #axis.text.x  = element_text(angle=90, vjust=0.5, hjust = 1) 
    axis.text.x  = element_text(angle=60, vjust=1, hjust = 1) 
  )+
  labs(x = NULL, y = "Relative abundance (%)" )

euka.phyla.top.plot
# ggsave("Eukaryote-relative-abundance-phyla.MinION.pdf", plot = euka.phyla.top.plot)

#### Genus -------------
rel.abun.euka.gen <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.euka.gen

euka.genus.plot <- plot_bar(rel.abun.euka.gen , x = "sample", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free") #+ 
# theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.all.genus.plot.MinION.pdf", plot = genus.plot, height = 49, width = 10)

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
X <- 0.7689

major_genera.euka <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.euka#  [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"              "Bacillariophyceae"             
#  [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"             
#  [4] "Bacillariophyceae"              "Bathycoccaceae"                 "Charales"                      
#  [7] "Chlorellales"                   "Chlorodendrales"                "Chytridiomycota incertae sedis"
# [10] "Coccidia"                       "CS clade"                       "Deuterostomia"                 
# [13] "Dictyosteliales"                "Haemosporida"                   "Hexamitinae"                   
# [16] "Kickxellomycotina"              "Mamiellaceae"                   "Monothalamids"                 
# [19] "Mucoromycotina"                 "Oligohymenophorea"              "Pelagomonadales"               
# [22] "Phaeocystaceae"                 "Piroplasmida"                   "Protostomia"                   
# [25] "saccharomyceta"                 "Spirotrichea"                   "Symbiodiniaceae"               
# [28] "Tracheophyta"                   "Trebouxiales"                   "Trypanosomatida"               
# [31] "Ustilaginomycotina" 

sel.row <- which(out$Genus %in% major_genera.euka) # qty 
length(sel.row) #3780

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 98.59056 % of relative abundance is covered
100 - 98.59056 # 1.40944 % left representing rare genera
rel.abun.cov.euka.genera <- 98.59056
perc.rare.euka.genera <- 1.40944

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
#  [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"             
#  [4] "Bacillariophyceae"              "Bathycoccaceae"                 "Charales"                      
#  [7] "Chlorellales"                   "Chlorodendrales"                "Chytridiomycota incertae sedis"
# [10] "Coccidia"                       "CS clade"                       "Deuterostomia"                 
# [13] "Dictyosteliales"                "Haemosporida"                   "Hexamitinae"                   
# [16] "Kickxellomycotina"              "Mamiellaceae"                   "Monothalamids"                 
# [19] "Mucoromycotina"                 "Oligohymenophorea"              "Other minor genera"            
# [22] "Pelagomonadales"                "Phaeocystaceae"                 "Piroplasmida"                  
# [25] "Protostomia"                    "saccharomyceta"                 "Spirotrichea"                  
# [28] "Symbiodiniaceae"                "Tracheophyta"                   "Trebouxiales"                  
# [31] "Trypanosomatida"                "Ustilaginomycotina"  

p$labels$fill <- "Genera"
p$labels$x <- NULL

#change to character, then revert to factor afterwards
# copy & paste & modify as required below. Shift Other minor genera to last position, retitle if required
old_genera<- c(
  "Agaricomycotina"   ,     "Anthozoa"           ,       "Apansporoblastina"             ,
  "Bacillariophyceae" ,     "Bathycoccaceae"     ,       "Charales"                      ,
  "Chlorellales"      ,     "Chlorodendrales"    ,       "Chytridiomycota incertae sedis",
  "Coccidia"          ,     "CS clade"           ,       "Deuterostomia"                 ,
  "Dictyosteliales"   ,     "Haemosporida"       ,       "Hexamitinae"                   ,
  "Kickxellomycotina" ,     "Mamiellaceae"       ,       "Monothalamids"                 ,
  "Mucoromycotina"    ,     "Oligohymenophorea"  ,       
  "Pelagomonadales"   ,     "Phaeocystaceae"     ,       "Piroplasmida"                  ,
  "Protostomia"       ,     "saccharomyceta"     ,       "Spirotrichea"                  ,
  "Symbiodiniaceae"   ,     "Tracheophyta"       ,       "Trebouxiales"                  ,
  "Trypanosomatida"   ,     "Ustilaginomycotina" ,       "Other minor genera"
)

## re-label genera for plot
new_genera <- c(
  "Agaricomycotina"   ,     "Anthozoa"           ,       "Apansporoblastina"             ,
  "Bacillariophyceae" ,     "Bathycoccaceae"     ,       "Charales"                      ,
  "Chlorellales"      ,     "Chlorodendrales"    ,       "Chytridiomycota incertae sedis",
  "Coccidia"          ,     "CS clade"           ,       "Deuterostomia"                 ,
  "Dictyosteliales"   ,     "Haemosporida"       ,       "Hexamitinae"                   ,
  "Kickxellomycotina" ,     "Mamiellaceae"       ,       "Monothalamids"                 ,
  "Mucoromycotina"    ,     "Oligohymenophorea"  ,       
  "Pelagomonadales"   ,     "Phaeocystaceae"     ,       "Piroplasmida"                  ,
  "Protostomia"       ,     "saccharomyceta"     ,       "Spirotrichea"                  ,
  "Symbiodiniaceae"   ,     "Tracheophyta"       ,       "Trebouxiales"                  ,
  "Trypanosomatida"   ,     "Ustilaginomycotina" ,       "Other minor genera"
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
# ggsave("Eukaryote-relative-abundance-genera.MinION.pdf", plot = euka.genera.top.plot)

getwd()

