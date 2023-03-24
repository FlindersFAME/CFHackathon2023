#  MinION script
## Taxonomy team
# Provides diversity analyais (alpha and beta), plus relative anundance plots, with a few statistical tests 

.libPaths()
R.Version() # "R version 4.2.2 (2022-10-31)"
citation()  # R Core Team (2022)

# Read data and libraries ----------
setwd("~/Documents/PhD/Collaborations/CF Hackathon/PseqFiles/MinION full") #change to match folder hierarchy in wd
datadir <- "~/Documents/PhD/Collaborations/CF Hackathon/PseqFiles/MinION full/datadir"

library(dplyr);packageVersion("dplyr") # '1.1.0'
library(tidyr);packageVersion("tidyr") # ‘1.3.0’
library(phyloseq);packageVersion("phyloseq") # '1.42.0’
library(ggplot2);packageVersion("ggplot2") # '3.4.1’
library(ggpubr);packageVersion("ggpubr") # '0.6.0’
library(readxl);packageVersion("readxl") # '1.4.2’
library(vegan);packageVersion("vegan") # '2.6.4’

# Requires an abundance table, taxonomy table, and metadata table
reads <- read.csv("reads_finalMinION_All_corrected.csv")
metadata <- readRDS("CF_Metadata_Table-[JCJ-OD-v-2023-03-23--1352].RDS")
# taxonomy <- read.csv("taxonomy_finalMinION_All.csv")

genus.list <- read.table(file = 'Genus_CladeReads.tsv', sep = '\t', header = TRUE)
colnames(genus.list)[colnames(genus.list)=="name"] <- "Genus"
genus.list2 <- genus.list %>% select(Genus, taxID, lineage)

phylum.list <- read.table(file = 'Phylum_CladeReads.tsv', sep = '\t', header = TRUE)
colnames(phylum.list)[colnames(phylum.list)=="name"] <- "Phylum"
phylum.list2 <- phylum.list %>% select(Phylum, taxID, lineage)

taxonomy <- genus.list2

metadata <- subset(metadata, minion_ID != "NA")
nrow(metadata)

### convert table to tabular split version
taxtable <- taxonomy %>%
  as_tibble() %>%
  separate(lineage, sep=">", c("Group", "Kingdom"))
taxtable.2 <- taxtable %>% select(taxID, Group, Kingdom, Genus)

taxtable_m <- as.matrix(taxtable.2)

unique(taxtable$Kingdom)
unique(taxtable$Group)

# Create phyloseq object ---------
otu <- otu_table(reads, taxa_are_rows = T)
otu <- otu[,-1]
tax <- tax_table(taxtable_m)

otu[is.na(otu)] <- 0

colnames(otu) <- gsub('X','',colnames(otu))
colnames(otu) <- strtrim(colnames(otu),18)
colnames(otu)[21:59] <- strtrim(colnames(otu)[21:59],17)
colnames(otu)

# sanity check
setdiff(metadata$unique_ID, (colnames(otu)))
# 770590_20170925_S and 1068841_20180306_S do not exist, remove from metadata

metadata <- as.data.frame(metadata)
metadata <- subset(metadata, unique_ID != "770590_20170925_S") 
metadata <- subset(metadata, unique_ID != "1068841_20180306_S") 
nrow(metadata) #59

rownames(metadata) <- metadata$unique_ID
meta <- sample_data(metadata)

ps <-  phyloseq(otu, tax)
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2168 taxa and 59 samples ]
# tax_table()   Taxonomy Table:    [ 2168 taxa by 4 taxonomic ranks ]

rownames(meta)

ps.0 <- merge_phyloseq(ps, meta)
ps.0
# otu_table()   OTU Table:         [ 2168 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 2168 taxa by 4 taxonomic ranks ]

## Filter 10% of reads  --------
library(microViz);packageVersion("microViz") # ‘0.10.6’
ps.1 <- tax_filter(ps.0, min_prevalence=0.1)
# Proportional min_prevalence given: 0.1 --> min 6/59 samples.

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
sort(sample_sums(ps.1))
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
  sample=obs.ps.1$data$unique_ID,          # sample ID
  observed=obs.ps.1$data$value,            # On log normalised data
  shannon=shan.ps.2$data$value,            # On log normalised data
  chao1=chao.ps.1$data$value,              # On unnormalised data
  patient.id=shan.ps.2$data$Patient.x,     # patient.id
  date=shan.ps.2$data$Date.x,              # sampling date
  antibioticsYN=shan.ps.2$data$Antibiotics_YN.x,              # Yes No antibiotics
  age=shan.ps.2$data$Age.groups,           # Age
  IP.OP=shan.ps.2$data$IP.vs.OP              # In patient vs out patient
)

glimpse(out.r1)
# Rows: 59
# Columns: 9
# $ sample        <chr> 
# $ observed      <dbl> 
# $ shannon       <dbl> 
# $ chao1         <dbl> 
# $ patient.id    <dbl> 
# $ date          <dttm>
# $ antibioticsYN <chr> 
# $ age           <dbl> 
# $ IP.OP         <chr> 

out.r1$patient.id <- as.factor(out.r1$patient.id)
out.r1$antibioticsYN <- as.integer(out.r1$antibioticsYN)
out.r1$eff_no_spp <- exp(out.r1$shannon) # calculate effective no of species
out.r1$date <- as.Date(out.r1$date, format = "%d/%m/%Y")
str(out.r1)
head(out.r1)
#               sample observed  shannon     chao1 patient.id       date antibioticsYN age IP.OP eff_no_spp
# 1 1447437_20171212_S      520 5.993719  778.1481    1447437 2017-12-12             0   5    OP   400.9028
# 2 1128691_20171218_S      151 4.813398  241.0000    1128691 2017-12-18             1   7    OP   123.1494
# 3 1128691_20180116_S      231 5.214521  328.5000    1128691 2018-01-16             1   7    OP   183.9238
# 4 1590009_20171212_S     1200 6.842575 1712.7219    1590009 2017-12-12             0   5    OP   936.8986
# 5 1282052_20180206_S      804 6.447013 1145.4519    1282052 2018-02-06             1   5    OP   630.8154
# 6 1316935_20180417_S      951 6.604602 1376.6842    1316935 2018-04-17            NA   5    OP   738.4859

shanon.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=shannon))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Shannon's diversity index")+
  xlab("Patient ID")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
shanon.div.rare.plot

chao.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=chao1))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Chao1 diversity index")+  
  xlab("Patient ID")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
chao.div.rare.plot

exp.shanon.div.rare.plot <- ggplot(out.r1, aes(x= patient.id, y=eff_no_spp))+
  geom_boxplot(colour = "black")+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("exp(Shannon's diversity index)")+
  xlab("Patient ID")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
exp.shanon.div.rare.plot

library(ggpubr)
shan.chao.div.rare <- ggarrange(shanon.div.rare.plot,chao.div.rare.plot, common.legend = TRUE, ncol = 1)
shan.chao.div.rare
# ggsave("Diversity-shannon-chao-59samp.MinION-all.pdf", plot =shan.chao.div.rare, path = datadir)

# Identify longitudinal patients
out.r1.single <- out.r1 %>% group_by(patient.id) %>% filter(n() < 2)

out.r1.longitude <- out.r1 %>% group_by(patient.id) %>% filter(n() > 1)
long.patient.vec <- unique(out.r1.longitude$patient.id) 
long.patient.vec #1447437 1128691 1651490 1565754 1593973 658355  673895  698917  748160  748699  785991  788707  825012 

# Patient diversity over time
chao1.time.plot <- ggplot(out.r1.longitude, aes(x= date, y=chao1, colour = patient.id))+
  geom_point(size=2.5)+
  facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chao1.time.plot 

chao1.time.plot2 <- ggplot(out.r1.longitude, aes(x= date, y=chao1, colour = patient.id))+
  geom_point(size=2.5)+
  # facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chao1.time.plot2

shan.time.plot <- ggplot(out.r1.longitude, aes(x= date, y=shannon, colour = patient.id))+
  geom_point(size=2.5)+
  facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
shan.time.plot 

shan.time.plot2 <- ggplot(out.r1.longitude, aes(x= date, y=shannon, colour = patient.id))+
  geom_point(size=2.5)+
  # facet_wrap(.~patient.id, scales = "free")+ 
  theme_bw()+
  geom_line()+
  scale_x_date(date_labels =  "%m-%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
shan.time.plot2 

# Patient diversity by age
shanon.age.plot <- ggplot(na.omit(out.r1.single), aes(x= age, y=shannon, colour = age))+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Shannon's diversity index")+
  xlab("Age")
shanon.age.plot

chao1.age.plot <- ggplot(na.omit(out.r1.single), aes(x= age, y=chao1, colour = age))+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("Chao1 diversity index")+
  xlab("Age")
chao1.age.plot

exp.shanon.age.plot <- ggplot(na.omit(out.r1.single), aes(x= age, y=eff_no_spp, colour = age))+
  geom_point(size=2.5)+
  theme_bw()+
  ylab("exp(Shannon's diversity index)")+
  xlab("Age")
exp.shanon.age.plot

# Antibiotics
shanon.anti.plot <- ggplot(na.omit(out.r1.single), aes(x= as.factor(antibioticsYN), y=shannon, colour = as.factor(antibioticsYN)))+
  geom_boxplot(colour = "black",outlier.shape = NA)+
  geom_jitter(size=2.5, width = 0.2)+
  theme_bw()+
  ylab("Shannon's diversity index")+
  xlab("Antibiotics")
shanon.anti.plot

chao1.anti.plot <- ggplot(na.omit(out.r1.single), aes(x= as.factor(antibioticsYN), y=chao1, colour = as.factor(antibioticsYN)))+
  geom_boxplot(colour = "black",outlier.shape = NA)+
  geom_jitter(size=2.5, width = 0.2)+
  theme_bw()+
  ylab("Chao1 diversity index")+
  xlab("Antibiotics")
chao1.anti.plot

ggarrange(shanon.anti.plot,chao1.anti.plot, common.legend = TRUE)

exp.shanon.anti.plot <- ggplot(na.omit(out.r1.single), aes(x= as.factor(antibioticsYN), y=eff_no_spp, colour = as.factor(antibioticsYN)))+
  geom_boxplot(colour = "black")+
  geom_jitter(size=2.5, width = 0.2)+
  theme_bw()+
  ylab("exp(Shannon's diversity index)")+
  xlab("Antibiotics")
exp.shanon.anti.plot

# In patient - Out patient
shanon.ipop.plot <- ggplot(na.omit(out.r1.single), aes(x= as.factor(IP.OP), y=shannon, colour = as.factor(IP.OP)))+
  geom_boxplot(colour = "black",outlier.shape = NA)+
  geom_jitter(size=2.5, width = 0.2)+
  theme_bw()+
  ylab("Shannon's diversity index")+
  xlab("In-patient vs Out-patient")
shanon.ipop.plot

chao1.ipop.plot <- ggplot(na.omit(out.r1.single), aes(x= as.factor(IP.OP), y=chao1, colour = as.factor(IP.OP)))+
  geom_boxplot(colour = "black",outlier.shape = NA)+
  geom_jitter(size=2.5, width = 0.2)+
  theme_bw()+
  ylab("Chao1 diversity index")+
  xlab("In-patient vs Out-patient")
chao1.ipop.plot

exp.shanon.ipop.plot <- ggplot(na.omit(out.r1.single), aes(x= as.factor(IP.OP), y=eff_no_spp, colour = as.factor(IP.OP)))+
  geom_boxplot(colour = "black")+
  geom_jitter(size=2.5, width = 0.2)+
  theme_bw()+
  ylab("exp(Shannon's diversity index)")+
  xlab("In-patient vs Out-patient")
exp.shanon.ipop.plot

# Ordination on normalised data -----------
# All taxa in samples
ps.2
# otu_table()   OTU Table:         [ 2288 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 2288 taxa by 9 taxonomic ranks ]

set.seed(123)
PCoA_ord <- ordinate(ps.2, method="PCoA", distance="bray", weighted=TRUE)
PCoA_all.plot <- plot_ordination(ps.2, PCoA_ord)

# Patient ID ordination
PCoA_all <- PCoA_all.plot + 
  geom_point(aes(fill = as.character(Patient.x)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all
# ggsave("Ordination_PCoA-bray-curtis_patients-all-data.MinION.pdf", plot =PCoA_all, height = 3.5, width = 5.2)

# ps.2.long.only <- prune_samples(sample_names(ps.2) == long.samp.list, ps.2)

#list of the samples to keep
long.samp.list <- out.r1.longitude$sample

#prune to new phyloseq object
ps.2.long.only <- subset_samples(ps.2, unique_ID %in% long.samp.list)
ps.2.long.only@sam_data$Date.x <- as.Date(ps.2.long.only@sam_data$Date.x)

set.seed(123)
PCoA_ord.long <- ordinate(ps.2.long.only, method="PCoA", distance="bray", weighted=TRUE)
PCoA_all.long.plot <- plot_ordination(ps.2.long.only, PCoA_ord.long, label = "Date.x")

# Patient ID ordination - longitudinal data only
# Define the number of colors you want
library(randomcoloR)
n <- 13
set.seed(98)
palette <- distinctColorPalette(n)

PCoA_all.long <- PCoA_all.long.plot + 
  # geom_line(aes(colour = Patient.x, group = Patient.x))+
  geom_path(aes(group = Patient.x, colour = as.character(Patient.x)), size = 1, arrow = arrow(length = unit(0.3, "cm"), type = "closed"))+
  geom_point(aes(fill = as.character(Patient.x)), size=3, shape = 21, colour = "black")+
  scale_fill_manual(values = palette, name = "Patient ID")+
  scale_colour_manual(values = palette)+
  theme_bw()+
  guides(colour="none")
PCoA_all.long

## Other comparison plots -----------
PCoA_all.ab <- PCoA_all.plot + 
  geom_point(aes(fill = as.factor(Antibiotics_YN.x)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all.ab
# ggsave("PCoA_all.minion.antibiotics.pdf", plot = PCoA_all.ab)

PCoA_all.ipop <- PCoA_all.plot + 
  geom_point(aes(fill = as.factor(IP.vs.OP)), size=3, shape = 21, colour = "black")+
  scale_fill_manual(name = NULL,labels=c("In-patient", "Out-patient"), values=c("#F8766D", "#00BFC4"))+
  theme_bw()
PCoA_all.ipop

PCoA_all.age <- PCoA_all.plot + 
  geom_point(aes(fill = Age), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all.age

PCoA_all.room <- PCoA_all.plot + 
  geom_point(aes(fill = Room), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all.room

PCoA_all.hosp <- PCoA_all.plot + 
  geom_point(aes(fill = Hospital), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all.hosp

PCoA_all.PA <- PCoA_all.plot + 
  geom_point(aes(fill = Paediatric.vs.Adult), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_all.PA

PCoA_all.Sex <- PCoA_all.plot + 
  geom_point(aes(fill = Gender), size=3, shape = 21, colour = "black")+
  scale_fill_manual(name = "Sex",labels=c("Female", "Male"), values=c("#F8766D", "#00BFC4"))+
  theme_bw()
PCoA_all.Sex

PCoA_ord
PCoA_ord$vectors

# Currently can not run below script due to issues with taxonomy table issues--------
(all.K.psq <- ps.2)
# Bacteria only
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","Duplodnaviria", "Eukaryota", "unclassified&nbsp;viruses", "Varidnaviria") )
bacteria.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
bacteria.psq
# otu_table()   OTU Table:         [ 1046 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1046 taxa by 4 taxonomic ranks ]

bacteria.psq@sam_data$Patient.x <- as.character(bacteria.psq@sam_data$Patient.x)
bacteria.psq@sam_data$Date.x <- as.character(bacteria.psq@sam_data$Date.x)

set.seed(123)
PCoA_ord.bact <- ordinate(bacteria.psq, method="PCoA", distance="bray", weighted=TRUE)
PCoA_bact.plot <- plot_ordination(bacteria.psq, PCoA_ord.bact, label = "Date.x")

PCoA_bact <- PCoA_bact.plot + 
  geom_point(aes(fill = as.character(Patient.x)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_bact
# ggsave("Ordination_PCoA-bray-curtis_patients-bacterialonly-data.MinION.pdf", plot =PCoA_bact, height = 3.5, width = 5.2)

PCoA_bact.ab <- PCoA_bact.plot + 
  geom_point(aes(fill = as.factor(Antibiotics_YN.x)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_bact.ab

PCoA_ord.bact
PCoA_ord.bact$vectors
dev.off()

# Eukaryotes only
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","Duplodnaviria", "Bacteria", "unclassified&nbsp;viruses", "Varidnaviria") )
eukaryota.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
eukaryota.psq
# otu_table()   OTU Table:         [ 499 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 499 taxa by 4 taxonomic ranks ]

eukaryota.psq@sam_data$Patient.x <- as.character(eukaryota.psq@sam_data$Patient.x)
eukaryota.psq@sam_data$Date.x <- as.character(eukaryota.psq@sam_data$Date.x)

set.seed(123)
PCoA_ord.euka <- ordinate(eukaryota.psq, method="PCoA", distance="bray", weighted=TRUE)
PCoA_euka.plot <- plot_ordination(eukaryota.psq, PCoA_ord.euka, label = "Date.x")

PCoA_euka <- PCoA_euka.plot + 
  geom_point(aes(fill = as.character(Patient.x)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_euka
# ggsave("Ordination_PCoA-bray-curtis_patients-eukaryotesonly-data.MinION.pdf", plot =PCoA_bact, height = 3.5, width = 5.2)

PCoA_euka.ab <- PCoA_euka.plot + 
  geom_point(aes(fill = as.factor(Antibiotics_YN.x)), size=3, shape = 21, colour = "black")+
  theme_bw()
PCoA_euka.ab

PCoA_ord.euka
PCoA_ord.euka$vectors
dev.off()

## Statistics ----------
library(car);packageVersion("car") #‘3.1.1’
par(mfrow=c(2,2))

#Shannon's
aov.patient.sh <- aov(shannon~patient.id, out.r1)

plot(aov.patient.sh)
shapiro.test(resid(aov.patient.sh)) #W = 0.83623, p-value = 1.442e-06 Non-normal
leveneTest(shannon~patient.id, out.r1) # F= 0.7621, Pr(>F) = 0.7691
summary(aov.patient.sh)
#             Df Sum Sq Mean Sq F value Pr(>F)
# patient.id  39 10.979 0.28151   2.979 0.00652 **
# Residuals   19  1.796 0.09451 

#Chao1
aov.patient.ch <- aov(chao1~patient.id, out.r1)

plot(aov.patient.ch)
shapiro.test(resid(aov.patient.ch)) #W = 0.81453, p-value = 3.77e-07 # non-normal
leveneTest(chao1~patient.id, out.r1) # F= 0.6918, Pr(>F) = 0.838
summary(aov.patient.ch)
#              Df Sum Sq Mean Sq F value Pr(>F)
#  patient.id  39 4595487  117833   2.569 0.0151 *
#  Residuals   19  871444   45865   

# Test hypothesis that microbiota vary (with different centroids) by patient
# Calculate bray-curtis distance matrix
set.seed(123)
bray.rare <- phyloseq::distance(ps.2, method = "bray")
sam.df <- data.frame(sample_data(ps.2))

str(sam.df)
dev.off()

set.seed(123)
adonis2(bray.rare ~ Patient.x, data = sam.df)  
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = bray.rare ~ patient.id, data = sam.df)
#           Df SumOfSqs      R2      F Pr(>F)
# Patient.x  1   0.1454 0.02295 1.3387  0.178
# Residual  57   6.1897 0.97705              
# Total     58   6.3350 1.00000   

#Test for beta dispersion
beta <- betadisper(bray.rare, sam.df$Patient.x)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#            Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups    39 0.60830 0.0155975 2.719    999  0.093 .
# Residuals 19 0.10899 0.0057365  

# Relative abundance plots -------------------------
ps.1 # unnormalised
# otu_table()   OTU Table:         [ 1603 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1603 taxa by 4 taxonomic ranks ]

ps.2 # normalised
# otu_table()   OTU Table:         [ 1603 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1603 taxa by 4 taxonomic ranks ]

rank_names(ps.2) # "taxID"   "Group"   "Kingdom" "Genus"  

(all.K.psq <- ps.2)
# otu_table()   OTU Table:         [ 1603 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1603 taxa by 4 taxonomic ranks ]

## remove taxa not assigned as Bacteria
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","Duplodnaviria", "Eukaryota", "unclassified&nbsp;viruses", "Varidnaviria") )
bacteria.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
bacteria.psq
# otu_table()   OTU Table:         [ 1046 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1046 taxa by 4 taxonomic ranks ]

## remove taxa not assigned as Eukaryota
levels(factor(tax_table(all.K.psq)[, "Kingdom"])) # "[1] "Archaea" "Bacteria" "cellular organisms" "Eukaryota" "root" "unclassified" "Viruses"   
rem_taxa <- which(tax_table(all.K.psq)[, "Kingdom"] %in%  c("Archaea","Duplodnaviria", "Bacteria", "unclassified&nbsp;viruses", "Varidnaviria") )
eukaryota.psq <- prune_taxa(all.K.psq, taxa = row.names(tax_table(all.K.psq)[-rem_taxa, ]) )
eukaryota.psq
# otu_table()   OTU Table:         [ 499 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 499 taxa by 4 taxonomic ranks ]

## Genera across all kingdoms: -------------------------
rel_abun.Genus.phy_obj <- transform_sample_counts(tax_glom(all.K.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel_abun.Genus.phy_obj
# otu_table()   OTU Table:         [ 1603 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1603 taxa by 4 taxonomic ranks ]

all.genus.plot <- plot_bar(rel_abun.Genus.phy_obj , x = "unique_ID", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~Patient.x, scales = "free")+
  theme(legend.position = "none")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
ggsave("relabund.all.genus.plot.MinION.pdf", plot = all.genus.plot, height = 20, width = 10, path=datadir)

## Relative abundance plots top taxa: -------------------------
### Bacteria -------------
#### Genus -------------
rel.abun.bact.gen <- transform_sample_counts(tax_glom(bacteria.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.bact.gen

bact.genus.plot <- plot_bar(rel.abun.bact.gen , x = "unique_ID", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~Patient.x, scales = "free")+ 
  theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("relabund.bacteria.genus.plot.MinION.pdf", plot = bact.genus.plot, height = 49, width = 10)

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
#  modify X until number of major genera is manageable enough to explore in the plot, e.g. ~30 or so
X <- 1.46

major_genera.bact <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.bact
#  [1] "Abiotrophia"              "Aerococcus"               "Aliarcobacter"           
#  [4] "Anaeromyxobacter"         "Candidatus Mycosynbacter" "Candidatus Symbiothrix"  
#  [7] "Cupriavidus"              "Dysgonomonas"             "Ellagibacter"            
# [10] "Francisella"              "Globicatella"             "Granulicatella"          
# [13] "Gulosibacter"             "Krasilnikovia"            "Lactobacillus"           
# [16] "Lactococcus"              "Mobiluncus"               "Olivibacter"             
# [19] "Parolsenella"             "Pasteurella"              "Pedobacter"              
# [22] "Periweissella"            "Pluralibacter"            "Prevotellamassilia"      
# [25] "Pseudoleptotrichia"       "Pseudoprevotella"         "Psychrobacter"           
# [28] "Slackia"                  "Spongiactinospora"        "Tetrasphaera"    "Slackia"                  "Spongiactinospora"        "Tetrasphaera"    
length(major_genera.bact) # 30

sel.row <- which(out$Genus %in% major_genera.bact) # qty 
length(sel.row) #1770

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 17.16459 % of relative abundance is covered
100 - 17.16459 # 82.83541 % left representing rare genera
rel.abun.cov.bact.genera <- 17.16459
perc.rare.bact.genera <- 82.83541

major_genera.bact # 30 genera

######## major genera - absolute abundance
rem_taxa.aa <- which(tax_table(bacteria.psq)[, "Genus"] %in%  major_genera.bact )
bacteria.maj.gen.psq.aa <- prune_taxa(bacteria.psq, taxa = row.names(tax_table(all.K.psq)[rem_taxa.aa, ]) )
bacteria.maj.gen.psq.aa
# otu_table()   OTU Table:         [ 30 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 30 taxa by 4 taxonomic ranks ]

sample_sums(bacteria.maj.gen.psq.aa@otu_table)
sum(bacteria.maj.gen.psq.aa@otu_table)

#### major genera - Relative abundance
rel.abun.bact.gen # is gives relative proportion of each feature
# otu_table()   OTU Table:         [ 1046 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 1046 taxa by 4 taxonomic ranks ]

rem_taxa.ra <- which(tax_table(rel.abun.bact.gen)[, "Genus"] %in%  major_genera.bact )
bacteria.maj.gen.psq.ra <- prune_taxa(rel.abun.bact.gen, taxa = row.names(tax_table(all.K.psq)[rem_taxa.ra, ]) )
bacteria.maj.gen.psq.ra
# otu_table()   OTU Table:         [ 30 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 30 taxa by 4 taxonomic ranks ]

# Generate OTU tables
rel.gen.otu <- bacteria.maj.gen.psq.ra@otu_table
aa.gen.otu <- bacteria.maj.gen.psq.aa@otu_table

# bacteria.maj.gen.psq.ra.2 <- tax_glom(bacteria.maj.gen.psq.ra,)
bacteria.maj.gen.psq.ra

sample_sums(rel.gen.otu)
sample_sums(aa.gen.otu)

ra.genera.bact <- psmelt(rel.abun.bact.gen)

ra.agglom.gen.bact <- ra.genera.bact %>% group_by (unique_ID, Genus) %>% 
  mutate (genus_sum = sum(Abundance)) %>% 
  distinct (Sample, Genus, genus_sum) %>%
  filter (Genus %in% major_genera.bact)

library(RColorBrewer)
# Define the number of colors you want
library(randomcoloR)
n <- 30
set.seed(123)
palette <- distinctColorPalette(n)

here.palette <- c("#61416B", "#DE6F90", "#DF5F46", "#8DB95B", "#DCDE82", "#DCC089", "#76ECAA", "#E3AACC", "#C869DB", "#B9C6E1",
                  "#E540A9", "#7F51E1", "#E8D5DA", "#DAE552", "#519D81", "#E8AC4E", "#AF9CD8", "#E19B8A", "#B4E9E3", "#E2E3C9",
                  "#D33DE9", "#85EC44", "#C6EAB2", "#D98BD5", "#998D83", "#6AA3DB", "#6EC5D8", "#677DDF", "#69E7D4", "#75E371"
                  )

#library(ggplot2)
MinION.all.sampl.bact.ra <- ggplot(ra.agglom.gen.bact, aes(x=unique_ID, y = genus_sum, fill=Genus))+
  scale_fill_manual(values = palette)+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  # facet_grid(.~Patient.x)+ # No patient variable
  ylab("Relative abundance (%)")
# ggsave("Relative-abundance-MajGen-bact.MinION.allsamp30.pdf", plot = MinION.all.sampl.bact.ra, width = 20, path = datadir)

### Eukaryota -------------
#### Genus -------------
rel.abun.euka.gen <- transform_sample_counts(tax_glom(eukaryota.psq, taxrank = "Genus" ), function(x) 100*x / sum(x))
rel.abun.euka.gen

euka.genus.plot <- plot_bar(rel.abun.euka.gen , x = "unique_ID", fill = "Genus" )+
  geom_bar(stat="identity")+
  facet_grid(.~patient.id, scales = "free") #+ 
# theme(legend.position = "bottom")
# genus.plot <- genus.plot + theme(legend.position = "bottom")
# ggsave("Relative-abundance-all-Euka.MinION.allsamp30.bact.pdf", plot = euka.genus.plot, height = 49, width = 10)

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
X <- 3.2

major_genera.euka <- levels(as.factor(as.character( out$Genus[which(out$Abun >= X )]  )))
major_genera.euka#  [1] "Agaricomycotina"                "Anthozoa"                       "Apansporoblastina"              "Bacillariophyceae"             
#  [1] "Artibeus"        "Blastocystis"    "Bradysia"        "Calocera"        "Canis"          
#  [6] "Castor"          "Corethrella"     "Ditylum"         "Fonticula"       "Gelidium"       
# [11] "Genlisea"        "Hylobates"       "Lamprigera"      "Leishmania"      "Marmota"        
# [16] "Meleagris"       "Onchocerca"      "Ophiophagus"     "Ovis"            "Pan"            
# [21] "Plutella"        "Portunus"        "Pteropus"        "Sciurus"         "Serinus"        
# [26] "Spodoptera"      "Tetraselmis"     "Thyridium"       "Trichinella"     "Wickerhamomyces"

length(major_genera.euka)

sel.row <- which(out$Genus %in% major_genera.euka) # qty 
length(sel.row) #1770

100*sum(out$Abun[sel.row] )/sum(out$Abun) # 33.63406 % of relative abundance is covered
100 - 33.63406 # 66.36594 % left representing rare genera
rel.abun.cov.euka.genera <- 33.63406
perc.rare.euka.genera <- 66.36594

major_genera.euka # 30 genera

######## major genera - absolute abundance
rem_taxa.aa <- which(tax_table(eukaryota.psq)[, "Genus"] %in%  major_genera.euka )
eukayotes.maj.gen.psq.aa <- prune_taxa(eukaryota.psq, taxa = row.names(tax_table(eukaryota.psq)[rem_taxa.aa, ]) )
eukayotes.maj.gen.psq.aa
# otu_table()   OTU Table:         [ 30 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 30 taxa by 4 taxonomic ranks ]

sample_sums(eukayotes.maj.gen.psq.aa@otu_table)
sum(eukayotes.maj.gen.psq.aa@otu_table)

#### major genera - Relative abundance
rel.abun.euka.gen # is gives relative proportion of each feature
# otu_table()   OTU Table:         [ 499 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 499 taxa by 4 taxonomic ranks ]

rem_taxa.ra <- which(tax_table(rel.abun.euka.gen)[, "Genus"] %in%  major_genera.euka )
eukayotes.maj.gen.psq.ra <- prune_taxa(rel.abun.euka.gen, taxa = row.names(tax_table(rel.abun.euka.gen)[rem_taxa.ra, ]) )
eukayotes.maj.gen.psq.ra
# otu_table()   OTU Table:         [ 30 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 152 sample variables ]
# tax_table()   Taxonomy Table:    [ 30 taxa by 4 taxonomic ranks ]

# Generate OTU tables
rel.gen.otu.euka <- eukayotes.maj.gen.psq.ra@otu_table
aa.gen.otu.euka <- eukayotes.maj.gen.psq.aa@otu_table

eukayotes.maj.gen.psq.ra

sample_sums(rel.gen.otu.euka)
sample_sums(aa.gen.otu.euka)

ra.genera.euka <- psmelt(rel.abun.euka.gen)

ra.agglom.gen.euka <- ra.genera.euka %>% group_by (unique_ID, Genus) %>% 
  mutate (genus_sum = sum(Abundance)) %>% 
  distinct (Sample, Genus, genus_sum) %>%
  filter (Genus %in% major_genera.euka)

# Define the number of colors you want
library(randomcoloR)
n <- 30
set.seed(123)
palette <- distinctColorPalette(n)

here.palette <- c("#61416B", "#DE6F90", "#DF5F46", "#8DB95B", "#DCDE82", "#DCC089", "#76ECAA", "#E3AACC", "#C869DB", "#B9C6E1",
                  "#E540A9", "#7F51E1", "#E8D5DA", "#DAE552", "#519D81", "#E8AC4E", "#AF9CD8", "#E19B8A", "#B4E9E3", "#E2E3C9",
                  "#D33DE9", "#85EC44", "#C6EAB2", "#D98BD5", "#998D83", "#6AA3DB", "#6EC5D8", "#677DDF", "#69E7D4", "#75E371"
)

#library(ggplot2)
MinION.all.sampl.euka.ra <- ggplot(ra.agglom.gen.euka, aes(x=unique_ID, y = genus_sum, fill=Genus))+
  scale_fill_manual(values = palette)+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  # facet_grid(.~Patient.x)+ # No patient variable
  ylab("Relative abundance (%)")
# ggsave("Relative-abundance-MajGen-Euka.MinION.allsamp30.pdf", plot = MinION.all.sampl.euka.ra, width = 20, path = datadir)
