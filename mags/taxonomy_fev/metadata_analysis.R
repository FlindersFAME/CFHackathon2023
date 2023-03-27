library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(ggrepel)

### Loading in data, setting column formats, etc. 
metadata_orig = read.csv('~/Documents/metadata_fame_hackathon.csv')
metadata = data.frame(cbind(cbind(metadata_orig[, 1:6], metadata_orig[, 9:10], metadata_orig[, 19:20], metadata_orig[, 22:23], metadata_orig[, 34:36]), apply(cbind(metadata_orig[, 7:8], metadata_orig[, 11:18], metadata_orig[, 21], metadata_orig[, 24:33], metadata_orig[, 37]), 2, as.numeric))); colnames(metadata)[26] = 'NTM'; colnames(metadata)[37] = 'CFLD'
metadata$Date = as.Date(metadata$Date, '%d/%m/%Y')
samples = read.csv('~/Documents/sequenced_samples.csv'); samples = samples[, 1:3]; samples = subset(samples, MGI_ID_from_fastq.gz != "")
metadata = subset(metadata, unique_ID %in% samples$unique_ID)

### Basic histograms of variables
hist(metadata$Age, main = "Age distribution of sequenced samples", xlab = "Age bracket")
barplot(table(table(metadata$Patient)), main = "Sampling distribution", xlab = "Number of samples per individual patient", ylab = "Frequency")


### Plots of longitudinal data
meta_2plus = subset(metadata, Patient %in% (rownames(data.frame(which(table(metadata$Patient) >= 2)))))
meta_2plus = meta_2plus[order(meta_2plus$Patient, meta_2plus$Date),]
meta_2plus$first_Patient = !duplicated(meta_2plus$Patient)
days = c()
for (i in 1:nrow(meta_2plus)) {
  if (meta_2plus$first_Patient[i]) {
    days = append(days, 0)
    i -> patient_start
  } else {
      days = append(days, meta_2plus$Date[i] - meta_2plus$Date[patient_start])
    }
}
meta_2plus$first_Patient = NULL; meta_2plus$Days = days
meta_2plus$log2days = log2(meta_2plus$Days)

plot = na.omit(cbind(meta_2plus[, 1:17], meta_2plus[, 21:23], meta_2plus$log2days)); plot$log2days = plot$`meta_2plus$log2days`; plot$`meta_2plus$log2days` = NULL
plot1 = ggplot(plot, aes(x = log2days, y = CH4_Corrected, group = Patient, colour = Gender)) + geom_line() + geom_point(aes(shape = IP.vs.OP), size = 2.5) + ylab("CH4 (corrected)") + theme_classic() + theme(legend.position = "none") + xlab("log2(days)")
plot2 = ggplot(plot, aes(x = log2days, y = H2_Corrected, group = Patient, colour = Gender)) + geom_line() + geom_point(aes(shape = IP.vs.OP), size = 2.5) + ylab("H2 (corrected)") + theme_classic() + theme(legend.position = "none") + xlab("log2(days)")
plot = na.omit(cbind(meta_2plus[, 1:17], meta_2plus[, 31:34], meta_2plus$log2days)); plot$log2days = plot$`meta_2plus$log2days`; plot$`meta_2plus$log2days` = NULL
plot3 = ggplot(plot, aes(x = log2days, y = FEV1, group = Patient, colour = Gender)) + geom_line() + geom_point(aes(shape = IP.vs.OP), size = 2.5) + ylab("FEV1") + theme_classic() + theme(legend.position = "none") + xlab("log2(days)")
plot = ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
annotate_figure(plot, top = text_grob("Longitudinal clinical data", face = "bold"))
