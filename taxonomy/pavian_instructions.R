# Pavian alternate install - may work where it did not previously.
library(devtools)
install_github("fbreitwieser/pavian")
library(pavian)

#Run Pavian GUI
pavian::runApp(port=5000)

#Can run as R Shiny app or be accessed at http://127.0.0.1:5000 
## Instructions for generating tables
## 1) Run Pavian (above)
## 2) upload .taxa_report files for each sample
## 3) open 'comparison' tab
## 4) Select desired tax rank (Green buttons, top left)
## 5) Make sure 'reads' are selected (Blue buttons)
## 6) Choose clade or taxon reads (you can include both - Yellow buttons on top right)
## 7) When you've got everything click download button on bottom left