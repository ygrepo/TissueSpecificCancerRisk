#### Process ICGC data to retain only primary tumors and filter out WGD samples ###
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(readxl)
system("mkdir out")

### Read data ###
# arm-level
cnv_arm <- read.table("../data/broad_values_by_arm.rmcnv.pt_170207.txt",h=T,sep="\t",stringsAsFactors=FALSE,check.names = FALSE)

# Patient data
clin_f = "../data/pcawg_sample_sheet.tsv"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)
# we can probably get cancer type through project code
table(clin$dcc_project_code)
clin$cancer_type = gsub("-.*","",clin$dcc_project_code)

# get tcga filtering info
clin2_f = "../data/pcawg_donor_clinical_August2016_v9.xlsx"
tcga_info <- read_excel (clin2_f)
colnames(tcga_info)[1] <- "donor_unique_id"
clin <- merge(clin, tcga_info [,c(1,5)], by = "donor_unique_id")
table(clin$dcc_specimen_type)
clin_primary = clin[grep("Primary tumour", clin$dcc_specimen_type), ] #filter to primary tumor only
sum(is.na(clin_primary$tcga_donor_uuid))
clin_primary_noTCGA = clin_primary[is.na(clin_primary$tcga_donor_uuid),]

# wgd data
wgd = read.table("../data/consensus.20170217.purity.ploidy.txt",h=T,sep="\t",stringsAsFactors=FALSE,check.names = FALSE)
colnames(wgd)[1] = "aliquot_id"

### merge CNV data to patient/cancer type data ###
cnv_arm_m = melt(cnv_arm)
colnames(cnv_arm_m) = c("chromosome_arm","aliquot_id","gistic_value")

# match based on aliquot_id
cnv_arm_m_clin = merge(cnv_arm_m,clin_primary_noTCGA[,c("aliquot_id","cancer_type")],by="aliquot_id")
filtered_data <- subset(cnv_arm_m_clin, chromosome_arm %in% c("13q", "17q"))
cnv_arm_m_clin_t <- as.data.frame(dcast(filtered_data, aliquot_id + cancer_type ~ chromosome_arm, value.var = "gistic_value"))
tn="../data/broad_values_by_arm.rmcnv.pt_170207_filtered.txt"
write.table(cnv_arm_m_clin_t, quote=F, sep="\t", file = tn, row.names = F)

# remove WGD samples
cnv_arm_m_clin_wgd = merge(cnv_arm_m_clin_t,wgd[,c("aliquot_id","wgd_status")],by="aliquot_id") 
cnv_arm_m_clin_rmwgd = cnv_arm_m_clin_wgd[cnv_arm_m_clin_wgd$wgd_status == "no_wgd",]
cnv_arm_m_clin_rmwgd = cnv_arm_m_clin_rmwgd[,colnames(cnv_arm_m_clin_rmwgd)!="wgd_status"]
tn="../data/broad_values_by_arm.rmcnv.pt_170207_filtered_noWGD.txt"
write.table(cnv_arm_m_clin_rmwgd, quote=F, sep="\t", file = tn, row.names = F)