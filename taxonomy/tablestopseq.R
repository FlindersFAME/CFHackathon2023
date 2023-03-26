# CF Kraken output -> Pavian -> Phyloseq
# 20/03/23 - Flinders 10.03 - Chris Keneally
## data were seperated into a "read counts" table and a "taxonomy table
setwd("~/Documents/Postgrad/CF Hackathon/PavianOut/Final") #change to match folder hierarchy in wd
reads <- read.csv("reads.csv")
taxonomy <- read.csv("taxonomy.csv")

### convert table to tabular split version
library(tidyr)
taxtable <- taxonomy %>%
  as_tibble() %>%
  separate(taxon, sep=">", c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Subspecies"))

taxtable_m <- as.matrix(taxtable)

#create phyloseq obj
library(phyloseq)
otu <- otu_table(reads, taxa_are_rows = T)
otu <- otu[,-1]
tax <- tax_table(taxtable_m)

ps <-  phyloseq(otu, tax)
ps
#3197 tax, 11 samples

