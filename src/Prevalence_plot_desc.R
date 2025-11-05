########## Figure 1. Freq of 13q/17q in various cancer types
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
source("./global_aes_out.R")
library("patchwork")

# Function to summarize chromosome arms data
summarize_chr_arms <- function(data) {
  # Identify chromosome arm columns
  chrom_arms <- names(data)[3:ncol(data)]
  # Initialize an empty list to store results
  results <- list()
  
  # Loop over each cancer type and chromosome arm
  for (cancer_type in unique(data$Type)) {
    data_cancer <- data %>% filter(Type == cancer_type)
    for (chrom_arm in chrom_arms) {
      # Total non-NA samples for this chromosome arm
      non_na_samples <- sum(!is.na(data_cancer[[chrom_arm]]))
      
      # Number of samples deleted (<= -1) for this chromosome arm
      deleted_samples <- sum(data_cancer[[chrom_arm]] <= -1, na.rm = TRUE)
      
      # Calculate the proportion of deleted samples
      proportion_deleted <- ifelse(non_na_samples > 0, deleted_samples / non_na_samples, 0)
      
      description = paste(cancer_type, " (",  deleted_samples, "/", non_na_samples,")", sep="")
      
      # Add the results to the list
      results <- append(results, list(data.frame(Cancer_Type = cancer_type, 
                                                 Chromosome_Arm = chrom_arm, 
                                                 Total_Non_NA_Samples = non_na_samples, 
                                                 Samples_Deleted = deleted_samples, 
                                                 Proportion_Deleted = proportion_deleted,
                                                 Description = description)))
    }
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}

## Function to plot the prevalence
plot_prevalence <- function(data) {
  data_name <- deparse(substitute(data))
  
  chr13q <- data %>% filter(Total_Non_NA_Samples >= 10, Chromosome_Arm %in% c("Chr13q","13q"))
  chr17q <- data %>% filter(Total_Non_NA_Samples >= 10, Chromosome_Arm %in% c("Chr17q","17q"))
  
  p1 <- ggplot(chr13q, aes(x = reorder(Description, desc(Proportion_Deleted)), y=Proportion_Deleted*100)) +
    geom_bar(stat = 'identity',fill = "#606060", width =0.7, alpha =0.7) + 
    theme_nogrid() +
    theme(axis.text.x = element_text(size = 9, angle = 90),
          axis.title.x = element_blank()) +
    labs(y="Prevalence (%)")
  fn = paste0('out/Prevalence_Figures/',data_name,'_chr13q.pdf')
  ggsave(fn,w = 6, h = 3)
  
  p2 <- ggplot(chr17q, aes(x = reorder(Description, desc(Proportion_Deleted)), y=Proportion_Deleted*100)) +
    geom_bar(stat = 'identity',fill = "#606060",width =0.7, alpha =0.7) + 
    theme_nogrid() +
    theme(axis.text.x = element_text(size = 9, angle = 90),
          axis.title.x = element_blank()) +
    labs(y="Prevalence (%)")
  fn = paste0('out/Prevalence_Figures/',data_name,'_chr17q.pdf')
  ggsave(fn,w = 6, h = 3)
}

## plot each data
TCGA = read.table("../data/PANCAN_ArmCallsAndAneuploidyScore_092817_filtered.txt",h=T,sep="\t",stringsAsFactors=FALSE)
cat("\nTCGA full data\n")
str(TCGA)
TCGA_filtered = summarize_chr_arms(TCGA)
write.table(TCGA_filtered, quote=F, sep = "\t", file ="out/TCGA_prevalence.txt", row.names = F)
plot_prevalence(TCGA_filtered)
print(TCGA_filtered[order(TCGA_filtered$Chromosome_Arm, -TCGA_filtered$Proportion_Deleted), c(2,5,6) ])

TCGA_noWGD = read.table("../data/PANCAN_ArmCallsAndAneuploidyScore_092817_filtered_noWGD.txt",h=T,sep="\t",stringsAsFactors=FALSE)
cat("\nTCGA noWGD data\n")
str(TCGA_noWGD)
TCGA_filtered_noWGD = summarize_chr_arms(TCGA_noWGD)
write.table(TCGA_filtered_noWGD, quote=F, sep = "\t", file ="out/TCGA_noWGD_prevalence.txt", row.names = F)
plot_prevalence(TCGA_filtered_noWGD)
print(TCGA_filtered_noWGD[order(TCGA_filtered_noWGD$Chromosome_Arm, -TCGA_filtered_noWGD$Proportion_Deleted), c(2,5,6) ])

ICGC = read.table("../data/broad_values_by_arm.rmcnv.pt_170207_filtered.txt",h=T,sep="\t",stringsAsFactors=FALSE)
cat("\nICGC full data\n")
str(ICGC)
colnames(ICGC) = c("aliquot_id","Type","13q","17q")
ICGC_filtered = summarize_chr_arms(ICGC)
write.table(ICGC_filtered, quote=F, sep = "\t", file ="out/ICGC_prevalence.txt", row.names = F)
plot_prevalence(ICGC_filtered)
print(ICGC_filtered[order(ICGC_filtered$Chromosome_Arm, -ICGC_filtered$Proportion_Deleted), c(2,5,6)])

ICGC_noWGD = read.table("../data/broad_values_by_arm.rmcnv.pt_170207_filtered_noWGD.txt",h=T,sep="\t",stringsAsFactors=FALSE)
cat("\nICGC noWGD data\n")
str(ICGC_noWGD)
colnames(ICGC_noWGD) = c("aliquot_id","Type","13q","17q")
ICGC_filtered_noWGD = summarize_chr_arms(ICGC_noWGD)
write.table(ICGC_filtered_noWGD, quote=F, sep = "\t", file ="out/ICGC_noWGD_prevalence.txt", row.names = F)
plot_prevalence(ICGC_filtered_noWGD)
print(ICGC_filtered_noWGD[order(ICGC_filtered_noWGD$Chromosome_Arm, -ICGC_filtered_noWGD$Proportion_Deleted), c(2,5,6)])