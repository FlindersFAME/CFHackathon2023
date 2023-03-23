# Alluvial plots - taxa by sample interval
# 21/03/23 - Chris Keneally - CF Hackathon 2023
# Input file format found inside "PavianOut Folder"
# You can format Kraken outputs to this using pavian - check to make sure you are using clade counts!

## Required packages
require('ggalluvial')
require('ggplot2')
# install easyalluvial
# devtools::install_github("erblast/easyalluvial")
require(easyalluvial)

## setwd
setwd("~/Documents/Postgrad/CF Hackathon/PavianOut") #change to match folder hierarchy in wd

## Read in counts data
data <- read.table('genus1.tsv', header = T)

## Count data (genus level) are split into seperate tables (per patient)
## Select 'names', and 'counts' columns per patient, duplicate names columns * intervals
pt1 <- data %>% select(c(1,5,6,7))
pt2 <- data %>% select(c(1,8,9,10))
pt3 <- data %>% select(c(1,11,12,13,14))

## rename columns to sample intervals
pt1_1 <- pt1 %>% 
  rename(first = 2, second = 3, third = 4)
pt2_1 <- pt2 %>%
  rename(first = 2, second = 3, third = 4)
pt3_1 <- pt3 %>%
  rename(first = 2, second = 3, third = 4, fourth = 5)

pt1_1[is.na(pt1_1)] <- 0
pt2_1[is.na(pt2_1)] <- 0
pt3_1[is.na(pt3_1)] <- 0

# calculate the total number of reads per sample
count_cols <- c("first", "second", "third")
total_reads <- colSums(pt1_1[, count_cols], na.rm=TRUE)

# convert the counts to relative abundance
pt1_1[, count_cols] <- pt1_1[, count_cols] / total_reads

# calculate the relative abundance of each taxon across all samples
taxon_abundance <- rowSums(pt1_1[, count_cols], na.rm=TRUE)

# select the taxa with relative abundances above 1%
pt1_filtered <- pt1_1[taxon_abundance >= 0.01, ]

## Then, data are converted to long format for ggalluvial (requires long format)
#pt1_1$name <- factor(pt1_1$name)
keycol <- "interval"
valuecol <- "count"
gathercols <- c("first", "second", "third")
pt1_long <- gather(pt1_filtered, keycol, valuecol, gathercols)

#figure out the bin labels (this is also passed to the plot function in next step to count bin numbers)
bin_labels = pt1_filtered$name

pt1_long$valuecol <- pt1_long$valuecol*100

ordered_by_abund = pt1_long %>%
  group_by(name) %>%
  count() %>%
  arrange( n ) %>%
  .[['name']]

#Plot1
pt1plot <- alluvial_long(pt1_long, key = keycol, value = valuecol, id = name, fill_by = 'value',
              bin_labels = bin_labels, bins = length(bin_labels), col_vector_value = palette_filter( greys = F),
              order_levels_fill = ordered_by_abund) + 
  theme_void()

# calculate the total number of reads per sample for pt2
count_cols <- c("first", "second", "third")
total_reads <- colSums(pt2_1[, count_cols], na.rm=TRUE)

# convert the counts to relative abundance
pt2_1[, count_cols] <- pt2_1[, count_cols] / total_reads

# calculate the relative abundance of each taxon across all samples
taxon_abundance <- rowSums(pt2_1[, count_cols], na.rm=TRUE)

# select the taxa with relative abundances above 1%
pt2_filtered <- pt2_1[taxon_abundance >= 0.01, ]

## Then, data are converted to long format for ggalluvial (requires long format)
#pt1_1$name <- factor(pt1_1$name)
keycol <- "interval"
valuecol <- "count"
gathercols <- c("first", "second", "third")

pt2_long <- gather(pt2_filtered, keycol, valuecol, gathercols)

bin_labels = pt2_filtered$name

alluvial_long(pt2_long, key = keycol, value = valuecol, id = name, fill_by = 'value',
                         bin_labels = bin_labels, bins = length(bin_labels), 
                         col_vector_value = palette_filter( greys = F)) + theme_void()

#Issue with alluvial_long() - removes levels seemingly indiscriminately needs troubleshooting

