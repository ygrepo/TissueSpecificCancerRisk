#### Process TCGA data to retain only primary tumors and filter out WGD samples ###
library(tidyverse)
library(readxl)
## read data
cnv_arm <- read.table("../../../Huang_lab_data/TCGA_PanCanAtlas_2018/CNV_Taylor_CancerCell2018/PANCAN_ArmCallsAndAneuploidyScore_092817.txt",h=T,sep="\t",stringsAsFactors=FALSE)
str(cnv_arm)
primary_only_cnv <- cnv_arm[grep("-01$", cnv_arm$Sample), ] #keep only primary tumors
str(primary_only_cnv)
cnv_aoi <- cnv_arm [,c(1,2, 28,34)] #keep only chr13q and chr17q
head(cnv_aoi)
colnames(cnv_aoi) [3:4] <- c ("Chr13q","Chr17q")
sum(duplicated(cnv_aoi$Sample)) #double check to ensure no duplicated samples

tn="../data/PANCAN_ArmCallsAndAneuploidyScore_092817_filtered.txt"
write.table(cnv_aoi, quote=F, sep="\t", file = tn, row.names = F)

# read Genome Doubling data
gn_db <- read_excel("../../CNV_immune/data/1-s2.0-S1535610818301119-mmc2.xlsx")[,1:6]
sample_no_db <- gn_db$Sample[gn_db$Genome_doublings == 0]

# aol_noWGD
cnv_aoi_noWGD <- cnv_aoi %>% filter(Sample %in% sample_no_db)
str(cnv_aoi)
str(cnv_aoi_noWGD)
tn="../data/PANCAN_ArmCallsAndAneuploidyScore_092817_filtered_noWGD.txt"
write.table(cnv_aoi_noWGD, quote=F, sep="\t", file = tn, row.names = F)