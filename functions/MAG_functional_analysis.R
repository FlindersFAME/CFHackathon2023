#R for the MAGs datasets

#load in the file

workingDir <- "C:\\Users\\belm9\\OneDrive - Flinders\\Documents\\CF Hackathon\\MAGs\\"
CF_counts <- read.csv(file=paste(workingDir,"CF_MAGs_subsystems_counts.csv",sep=""),head=TRUE,sep=",")

#create subset of data where qCPR results are present
#CF_63 <- subset(CF_counts, CF_counts$CF_bin=="63")

CF_counts$log10gene_count <- log10(CF_counts$Gene.Count)
CF_counts$Superclass <- ifelse(CF_counts$Superclass == "", NA, CF_counts$Superclass)
CF_counts <- subset(CF_counts, !is.na(CF_counts$Superclass))
CF_counts$CF_bin <- factor(CF_counts$CF_bin)

###violin plot with separated mag subset (9 mags) THIS WORKS YAY!!!!###

Fig2 <- ggplot(data=CF_counts, aes(x=CF_counts$Superclass, y=CF_counts$log10gene_count, fill=CF_counts$CF_bin))
Fig2 <- Fig2 + geom_violin(fill="grey77")
Fig2 <- Fig2 + xlab("superclass")
Fig2 <- Fig2 + ylab("gene count")
Fig2 <- Fig2 + ylim(0, 2)
Fig2 <- Fig2 + stat_boxplot(geom='errorbar',width=0.5)
Fig2 <- Fig2 + geom_boxplot(position=position_dodge(0.5), width=0.2,outlier.colour = "black", outlier.shape=8)
#Fig2 <- Fig2 + scale_fill_brewer(palette="Pastel2")
#Fig2 <- Fig2 + geom_hline(yintercept=0, linetype="dashed", color="red")
#Fig2 <- Fig2 + geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5))
Fig2 <- Fig2 + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Fig2 <- Fig2 + theme(text=element_text(size=10))
Fig2

## stacked bar chart of antimicrobial genes - this one works - problem was text size! :) need to change colours and collapse the rows
.
tiff(file = "mx.tiff", width = 10, height = 10, units = "in", res = 300)
mx = ggplot(CF_resistance, aes(x=CF_resistance$CF_bin, fill = CF_resistance$Subsystem.Name, y = CF_resistance$log10gene_count)) +
  geom_bar(stat= "identity", colour = "black") +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face = "bold"),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 5, face = "bold"),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "MAG", y= "log(10)Count", fill = "SubSystem") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_viridis_d(option = "viridis")



### NEW DATA



##counts for ALL MAGs

all_mags <- read.csv(file=paste(workingDir,"all_MAGs_subsystems_counts_with_category_withclassification.csv",sep=""),head=TRUE,sep=",")
names(all_mags) <- trimws(names(all_mags))

resallmags <- all_mags[all_mags$gene %in% c("Aminoglycosidemodifyingenzymes", "Beta-lactamasesAmblerclassA", "Antibiotictargets", "Bacitracinresistance", "Bicyclomycinresistancecluster", "Fosfomycinresistance", "Fusidicacidresistance", "MLSKOresistance", "Mupirocinresistance", "ResistancetoDaptomycin", "Resistancetothefluoroquinolonesnorfloxacinandciprofloxacin", "ResistancetoTriclosan", "TeicoplaninresistanceinStaphylococci", "Tetracyclineresistance,all.mechanisms", "VraTSRandLiaFSRthreecomponentregulatorysystems", "BetalactamasesAmberclassA"), ]
resallmags$logcount <- log10(resallmags$count)
resallmags$mag <- factor(resallmags$mag)

## trying to change the names in the legend - THIS WORKS!!!

mx = ggplot(resallmags, aes(x=resallmags$mag, fill = resallmags$gene, y = resallmags$logcount)) +
  geom_bar(stat= "identity") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1,),
        legend.text = element_text(size = 7, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 10, face = "bold"), legend.title = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold")) +
 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "MAG", y= "log(10)Count", fill = "SubSystem") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_viridis_d(option = "viridis", labels = c("Aminoglycoside Modifying Enzymes", "Antibiotic Targets", "Bacitracin Resistance", "Beta-Lactamase", "Bicyclomycin Resistance", "Fosfomycin Resistance", "Fusidicacid Resistance", "Mupirocin Resistance", "Daptomycin Resistance", "Fluoroquinolones, Norfloxacin and Ciprofloxacin Resistance", "Triclosan Resistance", "Teicoplanin Resistance", "Tetracycline Resistance", "VraTSR and LiaFSR three-component reg. systems"))
#scale_fill_viridis_d(option = "viridis")
mx

