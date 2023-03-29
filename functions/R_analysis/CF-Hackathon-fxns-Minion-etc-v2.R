#########################
#
# CF Hackathon - Functions - MinION data
# Craig Liddicoat Mar-2023
# DRAFT R code for data analysis
# Work in progress
#########################

# record library and version info
.libPaths() # "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"

R.Version()
# "R version 4.2.2 (2022-10-31)"
citation()
# R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL
# https://www.R-project.org/.


# not all libraries are used yet ...

library(readxl); packageVersion("readxl") # '1.4.1'
library(plyr); packageVersion("plyr") # '1.8.8'
library(dplyr); packageVersion("dplyr") # '1.0.10'
library(vegan);packageVersion("vegan") # '2.6.4'
library(phyloseq); packageVersion("phyloseq") # '1.42.0'
library(ggplot2); packageVersion("ggplot2") # '3.4.0'
library(grid); packageVersion("grid") #  '4.2.2'
library(reshape2); packageVersion("reshape2") # '1.4.4'
library(tidyr); packageVersion("tidyr") # '1.2.1'
library(ggforce); packageVersion("ggforce") # '0.4.1'
library(ggrepel); packageVersion("ggrepel") # '0.9.2'

library(stringdist); packageVersion("stringdist") # ‘0.9.10’
library(stringr); packageVersion("stringr") # ‘1.5.0’
library(doParallel); packageVersion("doParallel") # '1.0.17'

#library(zCompositions); packageVersion("zCompositions") # '1.4.0.1'
#library(propr); packageVersion("propr") # '4.2.6'

library(RColorBrewer); packageVersion("RColorBrewer") # '1.1.3'
library(ggpp); packageVersion("ggpp") # ‘0.5.0’ # https://cran.r-project.org/web/packages/ggpp/vignettes/grammar-extensions.html

library(corrplot)                  ;packageVersion("corrplot") #  '0.92'
library(caret)                     ;packageVersion("caret") # '6.0.93'
library(MASS)                     ;packageVersion("MASS") # ‘7.3.58.1’
library(ggsignif); packageVersion("ggsignif") # '0.6.4'
library(moments)                  ;packageVersion("moments") # ‘0.14.1’
library(ANCOMBC); packageVersion("ANCOMBC") # ‘2.0.1’
library(grDevices); packageVersion("grDevices") #  '4.2.2'
library(ggbiplot); packageVersion("ggbiplot") #  ‘0.55’
library(viridis); packageVersion("viridis") #  ‘0.6.2’

library(FSA); packageVersion("FSA") # '0.9.3'
library(rcompanion); packageVersion("rcompanion") # '2.4.18'


# up to line ...


#########################
##save.image("/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion/WORKSPACE-v2.RData")
##      load("/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion/WORKSPACE-v2.RData")
#########################

workdir <- "/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion"
setwd(workdir)
getwd()


par.default <- par()



#### Patient Metadata
#-------------------------

#### read in patient metadata

meta <- read_excel("/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/CF_Metadata_Table.xlsx",
                  sheet=1, range="A1:AK164", col_names = TRUE)
meta <- as.data.frame(meta)
str(meta)

row.names(meta) <- meta$unique_ID

names(meta)
# [1] "unique_ID"              "Patient"                "Date"                  
# [4] "IP vs OP"               "Hospital"               "Room"                  
# [7] "Age"                    "Age groups"             "Paediatric vs Adult"   
# [10] "Gender"                 "H2_Uncorrected"         "CH4_Uncorrected"       
# [13] "CO2"                    "H2_Corrected"           "CH4_Corrected"         
# [16] "CH4/H2 ratio_corrected" "Corr."                  "Antibiotics_YN"        
# [19] "Antibiotics (duration)" "Culture Result"         "NTM"                   
# [22] "Previous 12 months"     "Others"                 "IgE"                   
# [25] "Spec IgE"               "Spec IgG"               "Precipitins"           
# [28] "FVC"                    "FEV1"                   "Best FEV1"             
# [31] "FEV1/best FEV1"         "CFRD"                   "PI"                    
# [34] "CF gene 1"              "CF gene 2"              "Notes"                 
# [37] "CFLD" 


# change these to numeric

vars_to_numeric <- c(
  "H2_Uncorrected"   ,      "CH4_Uncorrected"   ,    
  "CO2"              ,      "H2_Corrected"      ,     "CH4_Corrected"  ,       
  "CH4/H2 ratio_corrected", "Corr."      ,           
  "IgE"   ,                
  "Spec IgE"    ,           "Spec IgG"      ,         "Precipitins"    ,       
  "FVC"           ,         "FEV1"    
  
)

meta[ ,vars_to_numeric] <- lapply(X = meta[ ,vars_to_numeric], FUN = as.numeric)

str(meta)
# 'data.frame':	163 obs. of  37 variables:
# $ unique_ID             : chr  "623361_20180123_S" "634207_20180510_S" "634207_20180517_S" "639354_20171206_S" ...
# $ Patient               : num  623361 634207 634207 639354 639354 ...
# $ Date                  : POSIXct, format: "2018-01-23" "2018-05-10" ...
# $ IP vs OP              : chr  "OP" "IP" "IP" "IP" ...
# $ Hospital              : chr  "RAH" "WCH" "WCH" "WCH" ...
# $ Room                  : chr  "Chest Clinic 9" "Adol Rm9" "Adol Rm9" "Adolescent 10" ...
# $ Age                   : num  18 17 17 17 17 17 17 17 16 16 ...
# $ Age groups            : num  4 3 3 3 3 3 3 3 3 3 ...
# $ Paediatric vs Adult   : chr  "Adult" "Paediatric" "Paediatric" "Paediatric" ...
# $ Gender                : chr  "M" "F" "F" "F" ...
# $ H2_Uncorrected        : num  6 39 26 10 9 22 31 3 11 14 ...
# $ CH4_Uncorrected       : num  2 15 10 2 2 8 5 1 2 2 ...
# $ CO2                   : num  4.1 4.2 4 3.6 3.5 3.8 4.2 4 3.6 3.8 ...
# $ H2_Corrected          : num  8 50 36 15 14 32 40 4 17 20 ...
# $ CH4_Corrected         : num  3 19 14 3 3 12 7 1 3 3 ...
# $ CH4/H2 ratio_corrected: num  0.375 0.38 0.389 0.2 0.214 ...
# $ Corr.                 : num  1.34 1.27 1.37 1.52 1.57 1.44 1.3 1.48 1.52 1.44 ...
# $ Antibiotics_YN        : chr  "0" "1" "1" "1" ...
# $ Antibiotics (duration): chr  "No antibiotics for >2/12" "Yes (IV tobramycin/? For 1/7)" "Yes (IV tobramycin/? For 8/7)" "IV tob (x1 dose)" ...
# $ Culture Result        : chr  "P. aeruginosa (mucoid), S. maltophilia, A. fumigatus, oral flora" "Oral flora; Aspergillus flavus" "Oral flora; Aspergillus flavus" "P. aeruginosa (mucoid)" ...
# $ NTM                   : chr  "0" "0" "0" "0" ...
# $ Previous 12 months    : chr  "2; 6; 10" "1; 6; 11" "1; 6; 11" "2" ...
# $ Others                : chr  NA "Penicillium; Enterobacter cloacae" "Penicillium; Enterobacter cloacae" NA ...
# $ IgE                   : num  NA 1750 NA 107 NA NA 11 NA NA 175 ...
# $ Spec IgE              : num  NA 31 NA 4.7 NA NA 0 NA NA 1.9 ...
# $ Spec IgG              : num  NA 122 NA NA NA NA NA NA NA NA ...
# $ Precipitins           : num  NA 1 NA 2 NA NA 1 NA NA 1 ...
# $ FVC                   : num  92 87 84 85 94 92 83 87 89 90 ...
# $ FEV1                  : num  84 79 79 53 69 91 86 89 86 87 ...
# $ Best FEV1             : num  89 94 94 77 77 91 87 89 98 96 ...
# $ FEV1/best FEV1        : num  0.944 0.84 0.84 0.688 0.896 ...
# $ CFRD                  : num  0 2 2 0 0 0 1 1 2 2 ...
# $ PI                    : num  1 0 0 0 0 1 1 1 1 1 ...
# $ CF gene 1             : chr  "F508" "F508" "F508" "F508" ...
# $ CF gene 2             : chr  "F508" "A455" "A455" "p.Arg851" ...
# $ Notes                 : chr  NA "2 week admission for acute on chronic bronchitis" "2 week admission for acute on chronic bronchitis" "2 week admission for 25% decline in PFT" ...
# $ CFLD                  : chr  "0" "0" "0" "0" ...

plot(meta$`Age groups`, meta$FEV1)
plot(meta$`Age groups`, meta$`FEV1/best FEV1`)
plot(meta$`Age groups`, meta$`CH4/H2 ratio_corrected`)



## later need to subselect only those sample IDs that occur in our MinION data
## e.g. select on names(fxns) that match up to row.names(meta)
## some trimming will be needed
## e.g. this format: "X658355_20171204" ... vs this format: "1845116_20180403_S" 

#sel <- which(names(fxns) %in% row.names(meta))


names(meta)

# keyvars <- c(
#   "IP vs OP"           , #   "Hospital"             ,  "Room"        ,          
#   #"Age"                 ,   
#   #"Paediatric vs Adult"   ,
#   "Gender"               ,    
#   "H2_Corrected"        ,   "CH4_Corrected"   ,      
#   "CH4/H2 ratio_corrected", # "Corr."      ,           
#   "Antibiotics_YN"   ,     
#   # "Culture Result"     ,    
#   "NTM"            ,       
#   "IgE"                  , 
#   "Spec IgE"          ,     "Spec IgG"     ,          "Precipitins"    ,       
#   "FEV1"              ,
#   "FEV1/best FEV1"    # ,    "CFRD"        ,           "PI"    ,                
#   #"CF gene 1"         ,     "CF gene 2"    ,            
#   #"CFLD" 
# )
# keyvars
# # [1] "IP vs OP"               "Gender"                
# # [3] "H2_Corrected"           "CH4_Corrected"         
# # [5] "CH4/H2 ratio_corrected" "Antibiotics_YN"        
# # [7] "NTM"                    "IgE"                   
# # [9] "Spec IgE"               "Spec IgG"              
# # [11] "Precipitins"            "FEV1"                  
# # [13] "FEV1/best FEV1"  
# 
# length(keyvars) # 13
# 
# meta.melt <- melt(meta, id.vars = c("unique_ID", "Patient", "Age groups"))
# # Using sample, age_category, site_name, city, year_planted, restoration_status as id variables
# 
# names(meta.melt)
# # "unique_ID"  "Patient"    "Age groups" "variable"   "value"
# 
# sel.plot <- which(meta.melt$variable %in% keyvars)
# 
# p <- ggplot( meta.melt[sel.plot, ], aes(x = `Age groups`, y = value )) +
#   #geom_boxplot()+
#   geom_point() +
#   #geom_boxplot(outlier.shape = NA)+
#   #geom_jitter(size=1.5,width = 0.15, alpha=0.25) +
#   facet_wrap(facets = vars(variable), nrow = 2, scales = "free")+
#   theme_bw()+
#   theme(
#     #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1),
#     
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     
#     strip.text = element_text(size = rel(0.7)),
#     #strip.text.y = element_text(size = rel(1.1)),
#     strip.background = element_blank()
#     )
# 
# p
# 
# #grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=15, fontface="bold") )
# 
# dev.print(tiff, file = paste0(workdir,"/plots/","metadata-values.tiff"), width = 22, height = 22, units = "cm", res=600, compression="lzw", type = "cairo")


sel <- which(meta$sample == "788707_20180313")
meta$`IP vs OP`[sel] # NA
meta$`IP vs OP`[sel] <- "IP"

shapes.IPvsOP <- c("IP" = 16, "OP" = 1)

no_patients <- length(unique(meta$Patient)) # 76
col.XPatients <- sample( rainbow(no_patients) )

# include 'X' on named patient list to play nicely with ggplot scale_color_manual()
names(col.XPatients) <- paste0("X",unique(meta$Patient))
meta$XPatient <- paste0("X",meta$Patient)
meta$XPatient

meta$Hospital <- factor(meta$Hospital, levels = c("WCH", "RAH"), ordered = TRUE)

plot(meta$Age, meta$FEV1)
plot(meta$Age, meta$`IP vs OP`)

names(meta)
meta2 <- meta[ , c( "unique_ID"        ,      "Patient"           ,     "Date"       ,           
                    "IP vs OP"          ,     "Hospital"          ,     "Room"      ,            
                    "Age"                ,    "Age groups"         ,    "Paediatric vs Adult"   ,
                    "Gender"         ,
                    "H2_Corrected"    ,       "CH4_Corrected"     ,    
                    "CH4/H2 ratio_corrected" ,
                    "Antibiotics_YN"        ,
                    
                    "FEV1"                ,   "Best FEV1"      ,       
                    "FEV1/best FEV1"      ,   "CFRD"        ,                  
                    "sample"              ,   "XPatient")]
ok <- complete.cases(meta2)
sel.ok <- which(ok==TRUE)
meta2 <- meta2[sel.ok, ]

p <- ggplot(data = meta2, aes(x = Age, y = FEV1))+
  geom_point(aes(shape = `IP vs OP`, color = XPatient))+
  geom_path(aes(group = XPatient, color = XPatient), alpha = 0.4, arrow = arrow(angle = 30, length = unit(0.2, "cm")))+
  scale_shape_manual(values = shapes.IPvsOP)+
  scale_color_manual(values = col.XPatients)+
  geom_smooth(method = "loess", na.rm = TRUE, alpha = 0.4)+
  guides(color = "none")+
  theme_bw()+
  facet_wrap(facets = vars(Hospital), nrow = 1, scales = "fixed", drop = TRUE)+
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.8,"line"),
    legend.margin = margin(t = 1, l = 0, r = 0, b = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(jpeg, file = paste0(workdir,"/plots/","Metadata-FEV1-vs-Age-with-IP-OP.jpeg"), width = 16, height = 9, units = "cm", res=450, type="cairo")

plot(meta2$Age, meta2$`FEV1/best FEV1`)

p <- ggplot(data = meta2, aes(x = Age, y = `FEV1/best FEV1`))+
  geom_point(aes(shape = `IP vs OP`, color = XPatient))+
  geom_path(aes(group = XPatient, color = XPatient), alpha = 0.4, arrow = arrow(angle = 30, length = unit(0.2, "cm")))+
  scale_shape_manual(values = shapes.IPvsOP)+
  scale_color_manual(values = col.XPatients)+
  geom_smooth(method = "loess", na.rm = TRUE, alpha = 0.4)+
  guides(color = "none")+
  theme_bw()+
  geom_hline(yintercept=1, linetype = "dashed", color = "red")+
  facet_wrap(facets = vars(Hospital), nrow = 1, scales = "fixed", drop = TRUE)+
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.8,"line"),
    legend.margin = margin(t = 1, l = 0, r = 0, b = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(jpeg, file = paste0(workdir,"/plots/","Metadata-FEV1toBestFEV1-vs-Age-with-IP-OP.jpeg"), width = 16, height = 9, units = "cm", res=450, type="cairo")




#-------------------------

#### Fxn Minion - Read in Functions
#-------------------------

# need to prepare:
# - OTU table
# - TAX table
# - metadata sample table



# raw <- read.csv(file = "/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion/10sample_uniref50counts_with_SEED.tsv", 
#                      header = TRUE, sep = "\t")

raw <- read.csv(file = "/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion/all_minion_sample_uniref50counts_with_SEED.tsv", 
                header = TRUE, sep = "\t")

dim(raw) # 395195     66
length(unique(raw$X)) # 395195
row.names(raw) <- raw$X
raw <- raw[ ,-1]

names(raw)
names(raw) <- gsub(pattern = "_S_(.+)nooverlap", replacement = "", x = names(raw))
names(raw)
names(raw) <- gsub(pattern = "X", replacement = "", x = names(raw))
names(raw)
# [1] "1068841_20180306" "1112926_20171212" "1128691_20171218" "1128691_20180116" "1255498_20171212"
# [6] "1282052_20180206" "1316935_20180417" "1316979_20171215" "1447437_20171006" "1470026_20180502"
# [11] "1565754_20171128" "1565754_20180403" "1586713_20180309" "1593967_20180424" "1593973_20180427"
# [16] "1593973_20180504" "1598281_20180508" "1651490_20171010" "1651490_20180206" "1834617_20180501"
# [21] "1845116_20180403" "623361_20180123"  "639354_20171206"  "650003_20180207"  "658355_20171204" 
# [26] "658355_20180122"  "658355_20180321"  "673895_20180122"  "673895_20180205"  "676138_20180130" 
# [31] "698917_20171207"  "698917_20180128"  "698917_20190119"  "715927_20180205"  "748160_20180321" 
# [36] "748160_20180329"  "748699_20180329"  "748699_20180410"  "752797_20170927"  "753522_20180606" 
# [41] "756934_20181218"  "763742_20180129"  "768745_20171123"  "770590_20170925"  "770590_20180115" 
# [46] "778851_20171204"  "785991_20171129"  "785991_20171206"  "785991_20180321"  "788707_20171213" 
# [51] "788707_20180301"  "788707_20180313"  "788707_20181116"  "802971_20180605"  "825012_20181120" 
# [56] "825012_20181126"  "892355_20180123"  "895293_20180502"  "983493_20180123"  "FIG_function"    
# [61] "Uniref.function"  "Subsystem"        "Subclass"         "Class"            "Superclass"

( nreads0 <- length(unique(row.names(raw))) ) # 395195
length(unique(raw$FIG_function)) # 10666
length(unique(raw$Uniref.function)) # 53685
length(unique(raw$Subsystem)) # 578
length(unique(raw$Subclass)) # 111
length(unique(raw$Class)) # 32
length(unique(raw$Superclass)) # 15



fxns <- raw[ ,sort( c("1068841_20180306", "1112926_20171212", "1128691_20171218", "1128691_20180116" ,"1255498_20171212",
                      "1282052_20180206", "1316935_20180417", "1316979_20171215", "1447437_20171006" ,"1470026_20180502",
                      "1565754_20171128", "1565754_20180403", "1586713_20180309", "1593967_20180424" ,"1593973_20180427",
                      "1593973_20180504", "1598281_20180508", "1651490_20171010", "1651490_20180206" ,"1834617_20180501",
                      "1845116_20180403", "623361_20180123" , "639354_20171206" , "650003_20180207"  ,"658355_20171204" ,
                      "658355_20180122" , "658355_20180321" , "673895_20180122" , "673895_20180205"  ,"676138_20180130" ,
                      "698917_20171207" , "698917_20180128" , "698917_20190119" , "715927_20180205"  ,"748160_20180321" ,
                      "748160_20180329" , "748699_20180329" , "748699_20180410" , "752797_20170927"  ,"753522_20180606" ,
                      "756934_20181218" , "763742_20180129" , "768745_20171123" , "770590_20170925"  ,"770590_20180115" ,
                      "778851_20171204" , "785991_20171129" , "785991_20171206" , "785991_20180321"  ,"788707_20171213" ,
                      "788707_20180301" , "788707_20180313" , "788707_20181116" , "802971_20180605"  ,"825012_20181120" ,
                      "825012_20181126" , "892355_20180123" , "895293_20180502" , "983493_20180123" 
                
                )) ]

names(fxns)
str(fxns)
# 'data.frame':	395195 obs. of  59 variables:
#   $ 1068841_20180306: num  181 25 13 3 2 1 323 81 225 40 ...
# $ 1112926_20171212: num  43 0 3 1 1 1 30 11 43 12 ...
# $ 1128691_20171218: num  0 0 0 0 0 0 0 0 0 0 ...
# $ 1128691_20180116: num  0 0 0 0 0 0 0 0 0 0 ...
# $ 1255498_20171212: num  868 30 17 6 5 0 287 225 662 92 ...
# $ 1282052_20180206: num  113 2 0 0 0 0 73 18 68 37 ...
# $ 1316935_20180417: num  138 3 1 1 1 0 106 42 93 16 ...
# $ 1316979_20171215: num  83 40 0 0 0 0 65 36 124 15 ...
# $ 1447437_20171006: num  859 5 3 4 0 0 543 195 669 112 ...
# $ 1470026_20180502: num  4 0 0 2 0 0 14 5 7 0 ...
# $ 1565754_20171128: num  6 2 0 0 0 0 5 0 4 1 ...
# $ 1565754_20180403: num  12 0 0 0 0 0 1 22 54 10 ...
# $ 1586713_20180309: num  734 49 0 0 0 0 434 107 335 53 ...
# $ 1593967_20180424: num  161 23 1 0 0 0 262 102 274 49 ...
# $ 1593973_20180427: num  253 7 0 1 0 0 7 40 118 17 ...
# $ 1593973_20180504: num  492 2 0 0 0 0 211 122 353 61 ...
# $ 1598281_20180508: num  196 43 0 0 0 0 143 64 185 21 ...
# $ 1651490_20171010: num  28 11 2 0 0 0 18 5 18 8 ...
# $ 1651490_20180206: num  13 1 0 0 0 0 4 3 9 14 ...
# $ 1834617_20180501: num  262 8 0 0 0 0 50 96 196 53 ...
# $ 1845116_20180403: num  86 10 0 0 0 0 109 41 95 22 ...
# $ 623361_20180123 : num  486 75 26 4 0 0 124 131 348 59 ...
# $ 639354_20171206 : num  38 14 0 0 0 0 30 19 44 7 ...
# $ 650003_20180207 : num  0 0 0 0 0 0 0 0 0 0 ...
# $ 658355_20171204 : num  528 26 0 0 0 0 51 134 466 62 ...
# $ 658355_20180122 : num  883 61 0 0 0 0 497 144 448 71 ...
# $ 658355_20180321 : num  625 7 0 0 0 0 218 75 271 43 ...
# $ 673895_20180122 : num  39 7 0 0 0 0 0 8 11 1 ...
# $ 673895_20180205 : num  0 2 0 0 0 0 0 1 5 0 ...
# $ 676138_20180130 : num  60 1 0 0 0 0 48 32 149 20 ...
# $ 698917_20171207 : num  119 0 0 2 0 0 50 17 42 5 ...
# $ 698917_20180128 : num  0 0 0 0 0 0 1 0 0 1 ...
# $ 698917_20190119 : num  3 0 0 0 0 0 0 2 0 0 ...
# $ 715927_20180205 : num  79 16 4 0 0 0 25 22 83 11 ...
# $ 748160_20180321 : num  0 4 0 0 0 0 26 26 61 13 ...
# $ 748160_20180329 : num  0 0 0 0 0 0 6 1 4 2 ...
# $ 748699_20180329 : num  294 43 0 0 0 0 283 102 239 50 ...
# $ 748699_20180410 : num  224 0 0 0 0 0 109 57 139 20 ...
# $ 752797_20170927 : num  0 0 0 0 0 0 0 2 0 0 ...
# $ 753522_20180606 : num  109 34 0 0 0 0 131 59 94 12 ...
# $ 756934_20181218 : num  198 38 0 0 0 0 256 117 436 76 ...
# $ 763742_20180129 : num  155 2 0 0 0 0 98 57 132 18 ...
# $ 768745_20171123 : num  138 0 0 0 0 0 31 26 54 5 ...
# $ 770590_20170925 : num  0 26 35 10 5 0 24 15 48 9 ...
# $ 770590_20180115 : num  1066 60 8 1 0 ...
# $ 778851_20171204 : num  51 0 1 0 0 0 106 44 131 15 ...
# $ 785991_20171129 : num  185 4 0 0 0 0 53 5 37 7 ...
# $ 785991_20171206 : num  144 9 0 0 0 0 35 21 78 10 ...
# $ 785991_20180321 : num  427 1 66 24 12 3 163 71 201 33 ...
# $ 788707_20171213 : num  177 1 0 0 0 0 145 104 152 24 ...
# $ 788707_20180301 : num  284 0 0 0 0 0 262 110 195 31 ...
# $ 788707_20180313 : num  139 0 0 0 0 0 117 61 105 27 ...
# $ 788707_20181116 : num  348 10 0 0 0 0 199 88 269 38 ...
# $ 802971_20180605 : num  34 32 0 0 0 0 55 19 84 12 ...
# $ 825012_20181120 : num  48 0 0 0 0 0 59 38 135 25 ...
# $ 825012_20181126 : num  17 0 0 0 0 0 669 178 648 104 ...
# $ 892355_20180123 : num  10 0 0 0 0 0 6 2 11 12 ...
# $ 895293_20180502 : num  399 168 44 13 5 0 31 266 405 53 ...
# $ 983493_20180123 : num  136 19 0 0 0 0 54 63 113 18 ...


length(unique(raw$FIG_function)) # 10666
length(unique(raw$Uniref.function)) # 53685
length(unique(raw$Subsystem)) # 578
length(unique(raw$Subclass)) # 111
length(unique(raw$Class)) # 32
length(unique(raw$Superclass)) # 15

tax <- raw[ ,c("Superclass" , # 15
               "Class", # 32
               "Subclass", # 111
               "Subsystem", # 578
               "FIG_function", # 10666
               "Uniref.function"  )] # 53685

## cells with missing data across all functional annotation 
## do not play nicely with phyloseq

sel.missing <- which(is.na(tax[]), arr.ind = TRUE)
length(sel.missing) # 0

sel.missing <- which(tax[]=="", arr.ind = TRUE)
length(sel.missing) # 723036

tax[sel.missing] <- NA

ok <- complete.cases(tax)
sel.ok <- which(ok==TRUE)
length(sel.ok) # 334942

100*length(sel.ok)/nreads0 # 84.7536 % retained
100 - 100*length(sel.ok)/nreads0 # 15.2464 % of reads have no functional annotation & were excluded

tax.clean <- tax[sel.ok, ]

# also clean fxns table

identical(row.names(tax), row.names(fxns)) # TRUE

fxns.clean <- fxns[sel.ok, ]



dim(tax.clean) # 334942      6

nreads1 <- dim(tax.clean)[1]



unique( tax.clean$Superclass )
# [1] "unknown"                             "Protein processing"                 
# [3] "RNA processing"                      "Energy"                             
# [5] "Dna processing"                      "Stress response, defense, virulence"
# [7] "Metabolism"                          "Membrane transport"                 
# [9] "Cellular processes"                  "Regulation and cell signaling"      
# [11] "DNA processing"                      "Rna processing"                     
# [13] "Cell envelope"                       "Miscellaneous" 

unique( tax.clean$Class )
# [1] "unknown"                                                      "Protein Synthesis"                                           
# [3] "RNA Processing"                                               "Energy and Precursor Metabolites Generation"                 
# [5] "DNA Processing"                                               "Stress Response, Defense and Virulence"                      
# [7] "Fatty Acids, Lipids, and Isoprenoids"                         "Cofactors, Vitamins, Prosthetic Groups"                      
# [9] "Carbohydrates"                                                "Secondary Metabolism"                                        
# [11] "Membrane Transport"                                           "Iron acquisition and metabolism"                             
# [13] "Cell Cycle, Cell Division and Death"                          "Regulation and Cell signaling"                               
# [15] "Nucleosides and Nucleotides"                                  "Amino Acids and Derivatives"                                 
# [17] "Sulfur Metabolism"                                            "Respiration"                                                 
# [19] "Clustering-based subsystems"                                  "Protein Fate (folding, modification, targeting, degradation)"
# [21] "Prokaryotic cell type differentiation"                        "Nitrogen Metabolism"                                         
# [23] "Cell Envelope, Capsule and Slime layer"                       "Phosphate Metabolism"                                        
# [25] "Metabolite damage and its repair or mitigation"               "Prophages, Transposable elements, Plasmids"                  
# [27] "Miscellaneous"                                                "Microbial communities"                                       
# [29] "Experimental Subsystems"                                      "Photosynthesis"                                              
# [31] "Motility and Chemotaxis" 

unique( tax.clean$Subclass )


## samples

head( meta$unique_ID )
# [1] "623361_20180123_S" "634207_20180510_S" "634207_20180517_S" "639354_20171206_S" "639354_20171213_S"
# [6] "642660_20180601_S"

meta$sample <- gsub(pattern = "_S", replacement = "", x = meta$unique_ID)

dim(fxns) # 395195     59

length(which(names(fxns.clean) %in% meta$sample )) # 54

## corrected sample names ???

sample_names_v0 <- names(fxns.clean)

crx <- read_excel("/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/CF_Metadata_Table.xlsx",
                      sheet="Corrected patient samples list", range="A1:C17", col_names = TRUE)
crx <- as.data.frame(crx)
str(crx)
# 'data.frame':	16 obs. of  3 variables:
#   $ Sequencer type : chr  "MGI" "MGI" "MGI" "MinION" ...
# $ Old sample name: chr  "649354_20170206_S" "652927_20180215_S" "658355_20180301_S" "698917_20190119_S" ...
# $ New sample name: chr  "639354_20171206_S" "715927_20180226_S" "658355_20180327_S" "698917_20180119_S" ...

crx <- crx[which(crx$`Sequencer type` == "MinION"), ]

crx$old <- gsub(pattern = "_S", replacement = "", x = crx$`Old sample name` )
crx$new <- gsub(pattern = "_S", replacement = "", x = crx$`New sample name` )

length(which(meta$sample %in% crx$old)) # 0
length(which(meta$sample %in% crx$new)) # 5 - i.e. meta already contains the new sample name

length(which(names(fxns) %in% crx$old)) # 5
length(which(names(fxns) %in% crx$new)) # 0

length(crx$new) # 5

meta.fixed <- meta

fxns.fixed <- fxns.clean

# replace $old for $new sample name
for (i in 1:length(crx$old)) {
  #i<-1
  this_old_name <- crx$old[i]
  sel_fxns_name <- which(names(fxns) == this_old_name)
  names(fxns.fixed)[sel_fxns_name] <- crx$new[i]
  print(paste0("completed ",i))
}

row.names(meta.fixed) <- gsub(pattern = "_S", replacement = "", x = row.names(meta.fixed))

length(which(row.names(meta.fixed) %in% names(fxns.fixed))) # 59

identical(row.names(meta.fixed), names(fxns.fixed)) # FALSE

samps <- meta.fixed

samps <- samps[ names(fxns.fixed), ]
dim(samps) # 59 38


#-------------------------

#### Fxn Minion - get into Phyloseq object
#-------------------------

# fxns.fixed - is equiv to OTU table

# tax.clean - is equiv to TAX table

# samps - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.clean )
dim(tax.m) # 334942      6

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( fxns.fixed )
dim(otu.m)
# 334942     59

identical(row.names(otu.m), row.names(tax.m)) # TRUE

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334942 taxa and 59 samples ]
# tax_table()   Taxonomy Table:    [ 334942 taxa by 6 taxonomic ranks ]


## metadata ? need dataframe with row.names as sample_names

SAMP <- sample_data(samps)

phy = merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334942 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 334942 taxa by 6 taxonomic ranks ]


head(taxa_names(phy))
# [1] "UniRef50_V8BD69"     "UniRef50_T0V120"     "UniRef50_E7MQS7"     "UniRef50_A0A7X7V2F8"
# [5] "UniRef50_A0A7S3ID01" "UniRef50_Q99XH4"

head(phy@tax_table)

getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion"

saveRDS(object = phy, file = "phy-phyloseq-object-All-Minion-fxns-longreads.RDS")

#phy <- readRDS("phy-phyloseq-object-All-Minion-fxns-longreads.RDS")


#-------------------------

#### Fxn Minion - Alpha & Beta diversity
#-------------------------

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334942 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 334942 taxa by 6 taxonomic ranks ]

sort( sample_sums(phy))
# 1128691_20171218  752797_20170927 1565754_20171128  698917_20180128  673895_20180205 
# 1477            10350            19108            22735            23606 
# 1128691_20180116  673895_20180122  698917_20180119 1565754_20180403  768745_20171123 
# 24948            64378            72391            73361            95739 
# 715927_20180205  650003_20180207 1651490_20171010  639354_20171206  748160_20180329 
# 118141           127817           135676           155588           164784 
# 676138_20180130  825012_20181120  753522_20180606  698917_20171207  802971_20180605 
# 165596           181287           182176           186625           188754 
# 1447437_20171212 1593973_20180427  788707_20180313  763742_20180129  748160_20180321 
# 204863           211216           255520           261702           265226 
# 892355_20180123  983493_20180123 1651490_20180206  778851_20171204 1316935_20180417 
# 269082           272348           278404           285296           308226 
# 1651490_20171215  785991_20171206  785991_20171129  788707_20181116  748699_20180410 
# 338206           339469           350971           394198           415731 
# 1845116_20180403  788707_20180301 1068841_20180306 1586713_20180309  658355_20180321 
# 430344           439492           443192           445698           458201 
# 748699_20180329 1588281_20180508 1593973_20180504  623361_20180123 1593967_20180424 
# 468837           505293           522735           531627           536559 
# 1470026_20180502  788707_20171213  770590_20170925  658355_20180122 1834617_20180501 
# 551294           657472           699018           718436           721901 
# 895293_20180502  658355_20171204  825012_20181126  756934_20181218  770590_20180115 
# 729768           757781           925008           925650           960896 
# 785991_20180321 1447437_20171006 1590009_20171212 1282052_20180206 
# 1050528          1053432          1136483          1196208 


summary(sample_sums(phy))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1477  165190  308226  395438  534093 1196208

sort(table(phy@sam_data$Patient), decreasing = TRUE)
# 788707  658355  698917  785991 1651490  673895  748160  748699  770590  825012 1128691 1447437 
# 4       3       3       3       3       2       2       2       2       2       2       2 
# 1565754 1593973  623361  639354  650003  676138  715927  752797  753522  756934  763742  768745 
# 2       2       1       1       1       1       1       1       1       1       1       1 
# 778851  802971  892355  895293  983493 1068841 1282052 1316935 1470026 1586713 1588281 1590009 
# 1       1       1       1       1       1       1       1       1       1       1       1 
# 1593967 1834617 1845116 
# 1       1       1 

patients_with_3ormore_samps <- c("788707", "658355", "698917", "785991", "1651490")


## Normalize data by rarefying

# remove sample with low library size: "1128691_20171218"

keep_samps <- sample_names(phy)[ -which(sample_names(phy) == "1128691_20171218")]

phy2 <- prune_samples(keep_samps ,phy)

min(taxa_sums(phy2)) # 0
# prune taxa that have zero sequence reads
phy2 <- prune_taxa(taxa = taxa_sums(phy2) > 0, x = phy2)
phy2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334927 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 334927 taxa by 6 taxonomic ranks ]

min(sample_sums(phy2)) # 10350
min(taxa_sums(phy2)) # 1





# rarefy #1
seed <- 123
r1.ps <- rarefy_even_depth(phy2, sample.size = min(sample_sums(phy2)),
                            rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)


min(taxa_sums(r1.ps)) # 1
sample_sums(r1.ps) # all 10350
ntaxa(r1.ps) # 95357

shan.r1.ps <- plot_richness(r1.ps, measures=c("Shannon"))
shan.r1.ps
str(shan.r1.ps$data)

out <- data.frame(sample=shan.r1.ps$data$sample,
                  shannon=shan.r1.ps$data$value,
                  patient=shan.r1.ps$data$Patient,
                  in_vs_out_patient = shan.r1.ps$data$IP.vs.OP
                  )

# calculate effective no of species
out$eff_no_fxn <- exp(out$shannon)
unique(out$patient)
# [1] 1068841 1447437 1128691 1590009 1282052 1316935 1651490 1470026 1565754 1586713 1593967 1593973 1588281 1834617
# [15] 1845116  623361  639354  650003  658355  673895  676138  698917  715927  748160  748699  752797  753522  756934
# [29]  763742  768745  770590  778851  785991  788707  802971  825012  892355  895293  983493

str(out)
# 'data.frame':	58 obs. of  5 variables:
# $ sample           : chr  "1068841_20180306" "1447437_20171212" "1128691_20180116" "1590009_20171212" ...
# $ shannon          : num  8.75 8.72 4.49 8.67 8.39 ...
# $ patient          : num  1068841 1447437 1128691 1590009 1282052 ...
# $ in_vs_out_patient: chr  "OP" "OP" "OP" "OP" ...
# $ eff_no_fxn       : num  6334 6147.5 89.5 5812.1 4383.5 ...

# out$sample <- gsub(pattern = "X", replacement = "", x = out$sample)
out$sample <- factor(out$sample, levels = sort(out$sample), ordered=TRUE)
# out$patient <- gsub(pattern = "X", replacement = "", x = out$patient)

out

shapes <- c( "OP"= 1, "IP"=16) # OP = outpatient = open circle; IP = inpatient = filled circle


table(out$in_vs_out_patient, useNA = "ifany")
# IP   OP <NA> 
# 25   32    1 

# fix missing IP/OP data ??
sel <- which(is.na(out$in_vs_out_patient))
out[sel, ]
#             sample  shannon patient in_vs_out_patient eff_no_fxn Xpatient
# 51 788707_20180313 8.511191  788707              <NA>    4970.08  X788707
out$in_vs_out_patient[sel] # NA
out$in_vs_out_patient[sel] <- "OP"

sel <- which(phy@sam_data$sample == "788707_20180313")
phy@sam_data$`IP vs OP`[sel] # NA
phy@sam_data$`IP vs OP`[sel] <- "OP"

sel <- which(phy2@sam_data$sample == "788707_20180313")
phy2@sam_data$IP.vs.OP[sel] # NA
phy2@sam_data$IP.vs.OP[sel] <- "OP"

table(out$in_vs_out_patient, useNA = "ifany")
# IP OP 
# 25 33

# colours?
# http://www.sthda.com/english/wiki/the-elements-of-choosing-colors-for-great-data-visualization-in-r
#library(colortools)
#col1 <- wheel("steelblue", num = (dim(out)[1])/2)

# https://stackoverflow.com/questions/13665551/generate-pairs-of-bright-and-dark-colours-for-ggplot2
library(grDevices)

dim(out)[1] # 58
# no of patients
no_patients <- length(unique(out$patient)) # 39

col1 <- rainbow(no_patients)

# include 'X' on named patient list to play nicely with ggplot scale_color_manual()
names(col1) <- paste0("X",unique(out$patient))
out$Xpatient <- paste0("X",out$patient)

names(out)

head(col1)
# X1068841  X1447437  X1128691  X1590009  X1282052  X1316935 
# "#FF0000" "#FF2700" "#FF4E00" "#FF7600" "#FF9D00" "#FFC400" 

p <- ggplot(data=out, aes(x=sample, y=eff_no_fxn, color=Xpatient, shape = in_vs_out_patient)) +
  geom_point() +
  geom_line(aes(group = Xpatient), alpha=0.5)+
  ggtitle("Alpha diversity of Uniref Functions\n(rarefied to minimum library size)")+
  #geom_boxplot() +
  #geom_jitter(size=1.5, width = 0.15) +
  #theme(axis.text.x  = element_text(angle=90, hjust=1, vjust = 0.5) ) +
  labs(x = NULL, y = "Effective no of functions\nbased on exp(Shannon's diversity)")+
  theme_bw()+
  scale_shape_manual(values = shapes, name = "In-patient or out-patient")+
  scale_color_manual(values = col1)+
  guides(color = "none")+
  theme(
    #legend.position = c(0.15,0.18),
    legend.position = "bottom",
    legend.key.height = unit(0.8,"line"),
    legend.title = element_text(size = rel(0.9)),
    axis.text.x  = element_text(angle=60, hjust=1, vjust = 1, size = rel(0.8)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p

#dev.print(tiff, file = paste0(workdir,"/plots/","Minion-Uniref-Function-Alpha-diversity.tiff"), width = 12, height = 9, units = "cm", res=450, compression="lzw",type="cairo")
dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-ALL-Uniref-Function-Alpha-diversity.jpeg"), width = 20, height = 12, units = "cm", res=450, type="cairo")


## consider FEV1

out2 <- out

out2$FEV1 <- shan.r1.ps$data$FEV1
out2$FEV1.best.FEV1 <- shan.r1.ps$data$FEV1.best.FEV1
out2$Hospital <- shan.r1.ps$data$Hospital
out2$XPatient <- out2$Xpatient
out2$IP.vs.OP <- shan.r1.ps$data$IP.vs.OP

sel <- which(is.na(out2$Hospital)==TRUE)
out2[sel, ]
# sample  shannon patient in_vs_out_patient eff_no_fxn Xpatient FEV1 FEV1.best.FEV1
# 51 788707_20180313 8.511191  788707                IP    4970.08  X788707   NA             NA
# Hospital XPatient IP.vs.OP
# 51     <NA>  X788707     <NA>

# from inspection of other patient data for 788707, set as WCH IP
out2$Hospital[sel] <- "WCH"
out2$IP.vs.OP[sel] <- "IP"

plot(out2$eff_no_fxn, out2$FEV1)
plot(out2$eff_no_fxn, out2$FEV1.best.FEV1)

p <- ggplot(data=out2, aes(x=eff_no_fxn, y=FEV1))+
  geom_point(aes(shape = IP.vs.OP, color = XPatient))+
  geom_path(aes(group = Xpatient, color = XPatient), alpha = 0.4, arrow = arrow(angle = 30, length = unit(0.2, "cm")))+
  scale_shape_manual(values = shapes.IPvsOP)+
  scale_color_manual(values = col.XPatients)+
  geom_smooth(method = "loess", na.rm = TRUE, alpha = 0.4)+
  guides(color = "none")+
  theme_bw()+
  xlab("Effective no of Uniref50 functional genes")+
  facet_wrap(facets = vars(Hospital), nrow = 1, scales = "fixed", drop = TRUE)+
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.8,"line"),
    legend.margin = margin(t = 1, l = 0, r = 0, b = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(jpeg, file = paste0(workdir,"/plots/","MinION-Functions-Alpha-diversity-vs-FEV1-with-Hospitals-IP-OP.jpeg"), width = 16, height = 9, units = "cm", res=450, type="cairo")



### compare with taxonomic diversity???

divtax.shan <- read.csv(file = "/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/Fxns_minion/alphadiversity_shannon_bacteria_genus.csv", header = TRUE)
class(divtax.shan) # "data.frame"
divtax.shan$sample_id <- gsub(pattern = "_S", replacement = "", x = divtax.shan$sample_id)


dim(out2) # 58 11
dim(divtax.shan) # 59 10

out2$divtaxShan <- NA

for (i in 1:dim(out2)[1]) {
  #i<-1
  this_samp <- as.character( out2$sample[i] ) # NOTE SOME ARE MISSING DUE TO MISMATCH sample names???!!!!
  sel.div <- which(as.character(divtax.shan$sample_id) == this_samp)
  if (length(sel.div)==1) {
    out2$divtaxShan[i] <- divtax.shan$Shannon_H[sel.div]
  }
  print(paste0("completed ",i))
}

out2$eff_no_taxa <- exp(out2$divtaxShan)

p <- ggplot(data=out2, aes(x=eff_no_taxa, y=eff_no_fxn))+
  geom_point(aes(shape = IP.vs.OP, color = XPatient))+
  geom_path(aes(group = Xpatient, color = XPatient), alpha = 0.4, arrow = arrow(angle = 30, length = unit(0.2, "cm")))+
  scale_shape_manual(values = shapes.IPvsOP)+
  scale_color_manual(values = col.XPatients)+
  geom_smooth(method = "loess", na.rm = TRUE, alpha = 0.4)+
  guides(color = "none")+
  theme_bw()+
  xlab("Effective no of genera")+
  ylab("Effective no of Uniref50\nfunctional genes")+
  
  facet_wrap(facets = vars(Hospital), nrow = 1, scales = "fixed", drop = TRUE)+
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.8,"line"),
    legend.margin = margin(t = 1, l = 0, r = 0, b = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(jpeg, file = paste0(workdir,"/plots/","MinION-Taxonomic-Diversity-vs-Functions-Alpha-diversity--with-Hospitals-IP-OP.jpeg"), width = 16, height = 9, units = "cm", res=450, type="cairo")


## FEV1 vs Taxonomic diversity??
p <- ggplot(data=out2, aes(x=eff_no_taxa, y=FEV1))+
  geom_point(aes(shape = IP.vs.OP, color = XPatient))+
  geom_path(aes(group = Xpatient, color = XPatient), alpha = 0.4, arrow = arrow(angle = 30, length = unit(0.2, "cm")))+
  scale_shape_manual(values = shapes.IPvsOP)+
  scale_color_manual(values = col.XPatients)+
  geom_smooth(method = "loess", na.rm = TRUE, alpha = 0.4)+
  guides(color = "none")+
  theme_bw()+
  xlab("Effective no of genera")+
  facet_wrap(facets = vars(Hospital), nrow = 1, scales = "fixed", drop = TRUE)+
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.8,"line"),
    legend.margin = margin(t = 1, l = 0, r = 0, b = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(jpeg, file = paste0(workdir,"/plots/","MinION-Taxonomic-Alpha-diversity-vs-FEV1-with-Hospitals-IP-OP.jpeg"), width = 16, height = 9, units = "cm", res=450, type="cairo")






## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(r1.ps, "PCoA", "bray")

ord


str(r1.ps@sam_data)
names(sample_data(r1.ps))


sel <- which(r1.ps@sam_data$sample == "788707_20180313")
r1.ps@sam_data$IP.vs.OP[sel] # NA
r1.ps@sam_data$IP.vs.OP[sel] <- "OP"


#p <- plot_ordination(r1.ps, ord, type="samples", color="age_category", shape = "city")
p <- plot_ordination(r1.ps, ord, type="samples", color="patient")
p

str(p$data)

p$labels$x # "Axis.1   [12.8%]"
x_lab <- "PCo1 (12.8%)"

p$labels$y # "Axis.2   [7.2%]"
y_lab <- "PCo2 (7.2%)"


p_df <- p$data

p_df$XPatient <- paste0("X",p_df$Patient)
p_df$sample <- factor( p_df$sample, levels = sort(p_df$sample), ordered=TRUE )
p_df$sample_select_display <- NA
sel <- which(p_df$Patient %in% patients_with_3ormore_samps)
p_df$sample_select_display[sel] <- as.character( p_df$sample[sel] )

p_df$IP.vs.OP

  
p <- ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = XPatient))+
  ggtitle("Beta diversity of Uniref Functions\n(PCoA on Bray-Curtis distances\nrarefied to minimum library size)")+
  
  theme_bw()+
  xlim(-0.5, 0.4) + ylim(-0.4, 0.4)+
  #geom_polygon(aes(group=patient),alpha=0.1)+ # ,  inherit.aes = FALSE fill=NA,  alpha=0.2,
  
  #geom_line(aes(group=XPatient),alpha=0.6)+ # ,  inherit.aes = FALSE fill=NA,  alpha=0.2,
  
  #geom_path(aes(group=XPatient),alpha=0.6, arrow = arrow(type = "open", angle = 30, length = unit(0.2, "cm")))+ # ,  inherit.aes = FALSE fill=NA,  alpha=0.2,
  geom_path(aes(group=XPatient),alpha=0.6, linetype = "dotted", arrow = arrow(type = "open", angle = 30, length = unit(0.2, "cm")))+ # ,  inherit.aes = FALSE fill=NA,  alpha=0.2,
  
  geom_point(aes(shape = IP.vs.OP))+
  #geom_text(aes(x = Axis.1, y = Axis.2,label=label_display), nudge_x = 0.005, nudge_y = 0.005, size=2.5, inherit.aes = FALSE)+
  
  #geom_text(aes(x = Axis.1, y = Axis.2,label=sample_select_display), hjust=0, vjust=0, size=3, inherit.aes = FALSE)+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_shape_manual(values = shapes, name = "In-patient or out-patient")+
  scale_color_manual(values = col1)+
  guides(color = "none")+
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
p


#dev.print(tiff, file = paste0(workdir,"/plots/","Minion-Uniref-Functions-Beta-diversity-PCoA-Bray.tiff"), width = 16, height = 16, units = "cm", res=450, compression="lzw", type = "cairo")
dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-Uniref-ALL-Functions-Beta-diversity-PCoA-Bray.jpeg"), width = 18, height = 20, units = "cm", res=450,  type = "cairo")

dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-Uniref-ALL-Functions-Beta-diversity-PCoA-Bray--no-labels.jpeg"), width = 18, height = 20, units = "cm", res=450,  type = "cairo")

#-------------------------


#### Fxn Minion - explore functional themes - Heatmaps ??
#-------------------------

phy_in <- phy2

phy2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334927 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 334927 taxa by 6 taxonomic ranks ]

rank_names(phy_in)
# [1] "Superclass"      "Class"           "Subclass"        "Subsystem"       "FIG_function"   
# [6] "Uniref.function"

sum(sample_sums(phy_in)) # 23329371 = 23,329,371

sort( table(phy_in@tax_table[ ,"Superclass"]) , decreasing = TRUE)
# unknown                          Metabolism                  Protein processing 
# 222094                               34382                               18620 
# Energy                      Dna processing Stress response, defense, virulence 
# 11229                               10791                                7957 
# Cellular processes                  Membrane transport                      Rna processing 
# 7631                                6675                                5286 
# Cell envelope                      RNA processing       Regulation and cell signaling 
# 4243                                1863                                1449 
# Miscellaneous                      DNA processing 
# 1366                                1341 

# remove reads with unknown Superclass
sel.rem_taxa <- which(phy_in@tax_table[ ,"Superclass"] == "unknown")
head(phy_in@tax_table[sel.rem_taxa])
keep_taxa <- row.names(phy_in@tax_table)[-sel.rem_taxa]
head(phy_in@tax_table[ keep_taxa ,"Superclass"])

phy3 <- prune_taxa(keep_taxa, phy_in)
sum(sample_sums(phy3)) # 8759816 = 8,759,816

100*8759816/23329371 # 37.54844 % of reads remaining

100*(sum(sample_sums(phy2)) - sum(sample_sums(phy3)))/sum(sample_sums(phy2))
# 62.45156 % of reads removed due to unknown Superclass

# # use phy_in = phy3 from here !!!!!!!!!!

phy_in <- phy3
phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 112833 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 112833 taxa by 6 taxonomic ranks ]

phy.relabun.Superclass <- transform_sample_counts( tax_glom( phy_in, taxrank = "Superclass" ), function(x) 100*x / sum(x) )
phy.relabun.Superclass
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 13 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 13 taxa by 6 taxonomic ranks ]

phy.relabun.Class <- transform_sample_counts( tax_glom( phy_in, taxrank = "Class" ), function(x) 100*x / sum(x) )
phy.relabun.Class
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 31 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 31 taxa by 6 taxonomic ranks ]

phy.relabun.Subclass <- transform_sample_counts( tax_glom( phy_in, taxrank = "Subclass" ), function(x) 100*x / sum(x) )
phy.relabun.Subclass
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 133 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 133 taxa by 6 taxonomic ranks ]

phy.relabun.Subsystem <- transform_sample_counts( tax_glom( phy_in, taxrank = "Subsystem" ), function(x) 100*x / sum(x) )
phy.relabun.Subsystem
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 586 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 586 taxa by 6 taxonomic ranks ]

phy.relabun.FIG_function <- transform_sample_counts( tax_glom( phy_in, taxrank = "FIG_function" ), function(x) 100*x / sum(x) )
phy.relabun.FIG_function
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2928 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 38 sample variables ]
# tax_table()   Taxonomy Table:    [ 2928 taxa by 6 taxonomic ranks ]

phy.relabun.Uniref.function <- transform_sample_counts( tax_glom( phy_in, taxrank = "Uniref.function" ), function(x) 100*x / sum(x) )
phy.relabun.Uniref.function



phy_in <- phy.relabun.Superclass


#same order as factor for group
samp_order <- sample_names(phy_in)[order(phy_in@sam_data$sample)]
samp_order

## heatmap for Superclass

p <- phyloseq::plot_heatmap(phy_in, sample.label = "sample", taxa.label = "Superclass" ,
                            sample.order = samp_order )
pp <- p + guides(fill=guide_legend(title= "Relative\nabundance (%)" ))

pp

dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-Uniref-ALL-Functions-Heatmap-Superclass.jpeg"), width = 24, height = 20, units = "cm", res=450,  type = "cairo")



phy_in <- phy.relabun.Class

rank_names(phy_in) 
# "Superclass"      "Class"           "Subclass"        "Subsystem"       "FIG_function"    "Uniref.function"



## heatmap for Class
p <- phyloseq::plot_heatmap(phy_in, sample.label = "sample", taxa.label = "Class" ,
                            sample.order = samp_order )
pp <- p + guides(fill=guide_legend(title= "Relative\nabundance (%)" ))

pp

dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-Uniref-ALL-Functions-Heatmap-Class.jpeg"), width = 30, height = 30, units = "cm", res=450,  type = "cairo")




phy_in <- phy.relabun.Subclass

## heatmap for Subclass
p <- phyloseq::plot_heatmap(phy_in, sample.label = "sample", taxa.label = "Subclass" ,
                            sample.order = samp_order )
pp <- p + guides(fill=guide_legend(title= "Relative\nabundance (%)" ))

pp

dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-Uniref-ALL-Functions-Heatmap-Subclass.jpeg"), width = 30, height = 40, units = "cm", res=450,  type = "cairo")





## still to do: 
#  - explore themes (Classes, Subclasses, etc...) of interest??
#  - determine associations b/w health outcomes (e.g. FEV1, etc) with functional components??



#-------------------------


#### Relabun bar plots
#-------------------------
rank_names(phy) # "Superclass"      "Class"           "Subclass"        "Subsystem"       "FIG.function"    "Uniref.function"



phy.relabun <- transform_sample_counts( tax_glom( phy, taxrank = "Superclass" ), function(x) 100*x / sum(x) )
phy.relabun
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 14 taxa and 10 samples ]
# sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 14 taxa by 6 taxonomic ranks ]

data <- psmelt(phy.relabun) # create dataframe from phyloseq object
str(data)
# 'data.frame':	140 obs. of  6 variables:
# $ OTU       : chr  "UniRef50_A0A1X1L0A4" "UniRef50_A0A1X1L0A4" "UniRef50_A0A1X1L0A4" "UniRef50_A0A1X1L0A4" ...
# $ Sample    : chr  "X698917_20190119" "X658355_20180321" "X788707_20171213" "X698917_20180128" ...
# $ Abundance : num  89.6 71.3 63.9 63.6 62.2 ...
# $ sample    : chr  "X698917_20190119" "X658355_20180321" "X788707_20171213" "X698917_20180128" ...

#simple way to rename phyla with < 1% abundance
data$Superclass[data$Abundance < 0.01] <- "< 1% abund."

#plot with condensed phyla into "< 1% abund" category
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Superclass))
p + geom_bar(aes(), stat="identity", position="stack")+
  theme_bw()+
  theme(
    
    legend.key.height = unit(0.8,"line"),
    legend.title = element_text(size = rel(0.9)),
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(0.8)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)

dev.print(jpeg, file = paste0(workdir,"/plots/","Minion-Superclass-Relabunbarplotp.jpeg"), width = 18, height = 14, units = "cm", res = 400,type="cairo")


#-------------------------


#### AMR profiles - draft
#-------------------------

amr <- read_excel("/Users/lidd0026/WORKSPACE/PROJ/CFHackathon2023-FAME/AMR/presence_absence_abricate_minionReads_amrfinderplust all data.xlsx",
                  sheet=1, range="A1:CR62", col_names = TRUE)
amr <- as.data.frame(amr)
str(amr)

sel.zero <- which(amr[]==".", arr.ind = TRUE)
head( amr[sel.zero] ) # "." "." "." "." "." "."
amr[sel.zero] <- "0"

# any NA values?
sel.na <- which(is.na(amr[]), arr.ind = TRUE)
amr[sel.na] # NA
amr[sel.na] <- "0"


names(amr)
# [1] "#FILE"          "aac(6')-IIc"    "aac(6')-Im"     "aacA-ENT1"      "ant(4')-Ia"     "ant(6)-Ia"     
# [7] "aph(2'')-IIa"   "aph(2'')-If"    "aph(2'')-Ih"    "aph(3')-IIIa"   "blaCSP-1"       "blaI_of_Z"     
# etc...

# replace FileID - write as row.name
# https://stackoverflow.com/questions/12677178/regular-expression-with-wildcards-to-match-any-character
head(gsub(pattern = "_S_(.+)txt", replacement = "", x = amr$`#FILE`))
# [1] "amrfinder/1068841_20180306" "amrfinder/1112926_20171212" "amrfinder/1128691_20171218"
# [4] "amrfinder/1128691_20180116" "amrfinder/1255498_20171212" "amrfinder/1282052_20180206"

# make row.names as sample ids
row.names(amr) <- gsub(pattern = "_S_(.+)txt", replacement = "", x = amr$`#FILE`)
head(row.names(amr))
# same as above

row.names(amr) <- gsub(pattern = "amrfinder/", replacement = "", x = row.names(amr))

#remove 1st column
head( amr[ ,1] )
amr <- amr[ ,-1]

# now replace non-zero strings with counts of cells with entries like this: 99.26;97.77;98.32;98.70;98.88;95.90;97.21;97.21;

amr.tempcopy <- amr

sel <- which(!amr[]=="0", arr.ind = TRUE)
length(sel) # 1648
head(amr[sel])
# [1] "96.91"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# [2] "99.26;97.77;98.32;98.70;98.88;95.90;97.21;97.21;98.32;95.72;99.26;98.88;97.58;98.70;98.70;99.07;97.58;97.39;97.58;99.81;96.46;98.51;97.95;99.44;98.14;94.60"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
# [3] "98.36"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# [4] "97.54;98.18;98.44;98.57;99.35;98.96;99.61;97.67;86.64;98.96;96.50;97.92;95.72;98.18;98.83;96.50;97.51;97.02;83.07;98.31;97.92;97.0
# etc ... some very long ones

for (i in 1:length(sel)) {
  #i<-3
  this_string <- amr[sel][i]
  this_count <- length(unlist(strsplit(x = this_string, split = ";")))
  amr[sel][i] <- this_count
  print(paste0("completed ",i))
}

str(amr)

# need to convert character fields into numeric counts
amr.num <- amr
amr.num[ ,names(amr.num)] <- lapply(FUN = as.numeric, X = amr[ ,names(amr)])

# transpose so that columns are samples, rows are fxn/taxa
amr.num <- t(amr.num)

str(amr.num)
# num [1:95, 1:61] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:95] "aac(6')-IIc" "aac(6')-Im" "aacA-ENT1" "ant(4')-Ia" ...
# ..$ : chr [1:61] "1068841_20180306" "1112926_20171212" "1128691_20171218" "1128691_20180116" ...


## amr taxonomy ?

amr.tax <- data.frame(AMRgene = row.names(amr.num))
row.names(amr.tax) <- amr.tax$AMRgene

## later join AMR 'taxonomy' info ???


## read into Phyloseq object

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#tax.m <- as.matrix( tax.fxn )
tax.m <- as.matrix( amr.tax )
dim(tax.m) # 95 1

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
#otu.m <- as.matrix( sfx.wide )
otu.m <- as.matrix( amr.num )
dim(otu.m)
# 95 61

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 95 taxa and 61 samples ]
# tax_table()   Taxonomy Table:    [ 95 taxa by 1 taxonomic ranks ]


## metadata ? need dataframe with row.names as sample_names
SAMP = data.frame(sample = sample_names(phy))
row.names(SAMP) <- SAMP$sample
SAMP$Patient <- gsub(pattern = "_(.+)", replacement = "", x = SAMP$sample)
SAMP$label_3patients <- NA
sel <- grep(pattern = "788707|698917|658355",x = SAMP$sample)
SAMP$label_3patients[sel] <- SAMP$sample[sel]

div <- plot_richness(phy,measures = c("Observed","Shannon"))
div
head( div$data )
#            samples variable value se
# 1 1068841_20180306 Observed    12 NA
# 2 1112926_20171212 Observed    14 NA
# 3 1128691_20171218 Observed     4 NA
# 4 1128691_20180116 Observed     0 NA
# 5 1255498_20171212 Observed    16 NA
# 6 1282052_20180206 Observed    23 NA
class(div$data) # "data.frame"

div.wide <- dcast(div$data, formula = samples ~ variable, value.var = "value")
identical(row.names(SAMP),div.wide$samples) # TRUE
names(div.wide) # "samples"  "Observed" "Shannon" 
names(div.wide) <- c("samples", "No. observed", "Shannon")

SAMP <- cbind(SAMP, div.wide[ ,c("No. observed", "Shannon")])
SAMP <- sample_data(SAMP)


phy = merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 95 taxa and 61 samples ]
# sample_data() Sample Data:       [ 61 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 95 taxa by 1 taxonomic ranks ]


head(phy@sam_data)
# Sample Data:        [6 samples by 5 sample variables]:
#                            sample Patient label_3patients No..observed  Shannon
# 1068841_20180306 1068841_20180306 1068841            <NA>           12 1.779070
# 1112926_20171212 1112926_20171212 1112926            <NA>           14 1.744328
# 1128691_20171218 1128691_20171218 1128691            <NA>            4 1.332179
# 1128691_20180116 1128691_20180116 1128691            <NA>            0 0.000000
# 1255498_20171212 1255498_20171212 1255498            <NA>           16 2.027752
# 1282052_20180206 1282052_20180206 1282052            <NA>           23 1.924043


## Jaccard & Bray-Curtis plots ??

sort( sample_sums(phy) )
# 1128691_20180116  752797_20170927  673895_20180205 1128691_20171218  698917_20190119  673895_20180122 
# 0                3                4                5               11               15 
# 650003_20180207 1565754_20171128  748160_20180329  802971_20180605  748160_20180321 1565754_20180403 
# 21               23               33               57               60               62 
# 698917_20171207  763742_20180129  778851_20171204  715927_20180205  698917_20180128  983493_20180123 
# 64               95               99              102              104              115 
# 788707_20171213 1112926_20171212  768745_20171123  753522_20180606  748699_20180410  639354_20171206 
# 121              124              130              132              141              143 
# 1068841_20180306 1593973_20180427  785991_20171129  825012_20181120 1316935_20180417 1651490_20171010 
# 187              201              204              204              216              226 
# 788707_20180301  892355_20180123 1470026_20180502  788707_20180313 1845116_20180403  676138_20180130 
# 233              234              253              253              259              271 
# 785991_20171206  748699_20180329  770590_20170925  788707_20181116  623361_20180123 1586713_20180309 
# 284              290              332              337              364              370 
# 756934_20181218 1593967_20180424  895293_20180502 1593973_20180504 1834617_20180501  658355_20180321 
# 376              381              390              435              452              469 
# 1255498_20171212  658355_20180122  825012_20181126  875028_20180115 1316979_20171215  658355_20171204 
# 586              617              638              712              742              772 
# 785991_20180321  770590_20180115 1447437_20171006 1651490_20180206  642660_20180601 1598281_20180508 
# 824              830              863              873             1085             1419 
# 1282052_20180206 
# 2435 


## ordination plot
## PCoA + Bray-Curtis

r1.ps <- phy

set.seed(123)
# ord <- ordinate(r1.ps, "PCoA", "bray") # missing values are not allowed with argument 'na.rm = FALSE'
ord <- ordinate(r1.ps, "PCoA", "bray", na.rm = TRUE)



p <- plot_ordination(r1.ps, ord, type="samples", color="Shannon")
p

p$labels$x # "Axis.1   [22%]"
x_lab <- "PCo1 (22%)"

p$labels$y # "Axis.2   [13.7%]"
y_lab <- "PCo2 (13.7%)"

#temp <- r1.ps
p_df <- p$data



p <- ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = Shannon))+
  theme_bw()+
  geom_point()+
  xlab(x_lab) + ylab(y_lab)+
  geom_text(aes(label= label_3patients), size = 3)+
  guides(colour = guide_colorbar(title = "Shannon\ndiversity"))+
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

#dev.print(tiff, file = paste0(workdir,"/plots/","Alcoa-ASV-Bray-Curtis-NMDS-ordination.tiff"), width = 12.4, height = 10, units = "cm", res=600, compression="lzw", type = "cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","AMR-Bray-Curtis-PCoA-ordination.tiff"), width = 12, height = 14, units = "cm", res=450, compression="lzw", type = "cairo")



## ordination plot
## PCoA + Jaccard

r1.ps <- phy

set.seed(123)
# ord <- ordinate(r1.ps, "PCoA", "bray") # missing values are not allowed with argument 'na.rm = FALSE'
#ord <- ordinate(r1.ps, "PCoA", "bray", na.rm = TRUE)
ord <- ordinate(r1.ps, "PCoA", "jaccard", na.rm = TRUE)

p <- plot_ordination(r1.ps, ord, type="samples", color="Shannon")
p

p$labels$x # "Axis.1   [15.9%]"
x_lab <- "PCo1 (15.9%)"

p$labels$y # "Axis.2   [11.5%]"
y_lab <- "PCo2 (11.5%)"

#temp <- r1.ps
p_df <- p$data



p <- ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = Shannon))+
  theme_bw()+
  geom_point()+
  xlab(x_lab) + ylab(y_lab)+
  geom_text(aes(label= label_3patients), size = 3)+
  guides(colour = guide_colorbar(title = "Shannon\ndiversity"))+
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","AMR-Jaccard-PCoA-ordination.tiff"), width = 12, height = 14, units = "cm", res=450, compression="lzw", type = "cairo")




# amr.num

amr.num[1:5, 1:5]
#             1068841_20180306 1112926_20171212 1128691_20171218 1128691_20180116 1255498_20171212
# aac(6')-IIc                0                0                0                0                0
# aac(6')-Im                 0                0                0                0                0
# aacA-ENT1                  0                0                0                0                0
# ant(4')-Ia                 0                0                0                0                0
# ant(6)-Ia                  0                0                0                0                0


pco.bray<-cmdscale(vegdist(amr.num,method='bray'))
plot(pco.bray)
pco.bray.binary<-cmdscale(vegdist(amr.num,method='bray',binary=T))
plot(pco.bray.binary)
pco.jac<-cmdscale(vegdist(amr.num,method='jaccard')) # , eig = TRUE
plot(pco.jac)
pco.jac.binary<-cmdscale(vegdist(amr.num,method='jaccard',binary=T)) # , eig = TRUE
plot(pco.jac.binary)




### HEatmaps??
data <- as.matrix( t(amr.num) )

heatmap(data, scale="row", cexRow=0.15, cexCol=0.1, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
# cexRow=1.5,
# labRow=rownames(data), 

# Plot a legend in bottom right part of heatmap
legend(x = "topleft" , legend = c("low", "medium", "high"),
       cex = 0.8, fill = colorRampPalette(brewer.pal(8, "Blues"))(3))

dev.print(tiff, file = paste0(workdir,"/plots/","AMR-Heatmap-sampleRows-amrColumns-scale-rows.tiff"), width = 22, height = 16, units = "cm", res=450, compression="lzw", type = "cairo")



#-------------------------


