##########################################################
## Project: RBMX mutations
## Author: J Clark
## Date: June 21st, 2024
##########################################################

##########################################################
## Set working directory
##########################################################
setwd("C:/Users/jclark/iCloudDrive/The Hub/Projets/P5_RBMX Cancer database/data/")
output_path <- "C:/Users/jclark/iCloudDrive/The Hub/Projets/P5_RBMX Cancer database/analysis/21062024/output"
set.seed(326)

#################################################################################
## Load Packages
#################################################################################
library(data.table)
library(tidyverse)
library(readxl)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(lawstat)
library(table1)

#################################################################################
## Load Data
#################################################################################
# RBMX cases
rbmx_manifest <- fread("RBMX/gdc_manifest.2024-06-07.txt")
head(rbmx_manifest)

# RBMX clinical data
rbmx_clinical <- fread("RBMX/clinical.tsv")
head(rbmx_clinical)

# Define columns to retain from MAF files
cols1 <- c("case_id", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand", 
          "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", 
          "Tumor_Seq_Allele2", "dbSNP_RS", "Mutation_Status", "Gene", "One_Consequence", 
          "SIFT", "PolyPhen", "IMPACT")

# Define columns to retain from clinical file
cols2 <- c("case_id", "case_submitter_id", "project_id", "gender", "race",
           "age_at_diagnosis", "icd_10_code", "primary_diagnosis")

# Define columns in an MAF file to make sure we select the correct files in loop
test_file_path <- paste0("RBMX/", rbmx_manifest$id[1], "/", rbmx_manifest$filename[rbmx_manifest$id == rbmx_manifest$id[1]])
maf_cols <- colnames(fread(test_file_path))

# Load gene panels (from Kyle)
panels <- read_excel("Gene List Mission Bio Panels_formatted for R.xlsx")
cancer_genes <- c(panels[[3]], panels[[8]], panels[[13]]) %>% na.omit()


#################################################################################
## Organize analytical dataset
#################################################################################
analytical_dataset <- data.frame(matrix(ncol = 27, nrow = 0))

# Create additional cols & add placeholders for data from MAF and clinical files
colnames(analytical_dataset) <- c("manifest_id", "CT_mutation", "GA_mutation", "Type of mutation",
                         "Num_cancer_gene_mutations", "CTs_genome", "CTs_cancer", 
                         "GAs_genome", "GAs_cancer", cols1)

process_file <- function(file, manifest_id) {
  if (!file.exists(file_path)) {
    warning(paste("File does not exist:", file_path))
    return(NULL)
  }
  
  file <- fread(file_path, showProgress = TRUE)
  if (!all(colnames(file) %in% maf_cols)) {
    warning(paste("File does not contain all required columns:", file_path))
    return(NULL)
  }
  
  file <- file %>% select(all_of(cols1))
  file$manifest_id <- manifest_id
  
  file <- file %>%
    mutate(
      CT_mutation = ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "T", 1, 0),
      GA_mutation = ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "A", 1, 0),
      `Type of mutation` = ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "T", "C>T",
                           ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "A", "C>A",
                           ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "G", "C>G",
                           ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "T", "G>T",
                           ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "C", "G>C",
                           ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "A", "G>A", "Other")))))),
      CTs_genome = sum(CT_mutation[Hugo_Symbol != "RBMX"]),
      CTs_cancer = sum(CT_mutation[Hugo_Symbol %in% cancer_genes]),
      GAs_genome = sum(GA_mutation[Hugo_Symbol != "RBMX"]),
      GAs_cancer = sum(GA_mutation[Hugo_Symbol %in% cancer_genes]),
      Num_cancer_gene_mutations = sum(Hugo_Symbol %in% cancer_genes)
    ) %>%
    filter(Hugo_Symbol == "RBMX") %>%
    select(manifest_id, Num_cancer_gene_mutations, CT_mutation, GA_mutation, `Type of mutation`,
           CTs_genome, CTs_cancer, GAs_genome, GAs_cancer, all_of(cols1))
  
  return(file)
}

for (i in rbmx_manifest$id) {
  file_path <- paste0("RBMX/", i, "/", rbmx_manifest$filename[rbmx_manifest$id == i])
  file <- fread(file_path)
  processed_file <- process_file(file, i)
  analytical_dataset <- rbind(analytical_dataset, processed_file)
}


# Identify any cases with multiple RBMX mutations -- to avoid counting as separate observations 
dups <- analytical_dataset %>% 
  select("case_id") %>%
  group_by(case_id) %>%
  summarize(n = n()) # 116 unique IDs
dups

length(which(dups$n == 1)) # 97 cases with 1 mutation 

# Randomly select unique cases (removing 2nd, 3rd, ... nth iteration of RBMX mutations)
#analytical_dataset <- distinct(analytical_dataset[sample(1:nrow(analytical_dataset)), ], 
#           case_id, 
#           .keep_all = TRUE) # N = 116

# Select same ID and SIFT pairings from last analysis -- where cases were previously selected randomly -- to ensure
# the same cases are included
orig_ids <- fread("../analysis/archive/07062024/output/analytical_dataset.csv") %>%
  select(c("case_id", "manifest_id", "SIFT", "PolyPhen", "IMPACT", "Variant_Classification"))
analytical_dataset <- left_join(orig_ids,
             analytical_dataset,
             by = c("case_id", "manifest_id", "SIFT", "PolyPhen", "IMPACT", "Variant_Classification"))

rbmx_clinical <- distinct(rbmx_clinical[sample(1:nrow(rbmx_clinical)), ], 
                          case_id, 
                          .keep_all = TRUE) # N = 116 as well


# Finally, merge analytical_dataset (maf file) with clinical data
analytical_dataset <- left_join(analytical_dataset,
            rbmx_clinical %>% select(all_of(cols2)),
            by = "case_id")

#################################################################################
## Create variables for pathogenicity -- using both initial and new coding definitions
#################################################################################
analytical_dataset <- analytical_dataset %>%
  mutate(
    `PolyPhen (0=benign; 1=possibly or probably damaging)` = 
      ifelse(str_detect(PolyPhen, c("benign")), 0, 
             ifelse(str_detect(PolyPhen, "probably_damaging|possibly_damaging"), 1, NA)),
    
    `PolyPhen (0=benign/possibly damaging; 1=probably damaging)` = 
      ifelse(str_detect(PolyPhen, c("benign|possibly_damaging")), 0, 
             ifelse(str_detect(PolyPhen, "probably_damaging"), 1, NA)),
    
    `PolyPhen [0-1]` = as.numeric(gsub("[a-z()]", "", PolyPhen)),
    
    `SIFT (0=tolerated; 1=deleterious)` = 
      ifelse(str_detect(SIFT, "tolerated"), 0, 
             ifelse(str_detect(SIFT, "deleterious"), 1, NA)),
    
    `SIFT [0-1]` = as.numeric(gsub("[a-z()]", "", SIFT)),
    
    `IMPACT (0=low/modifier; 1=moderate/high)` = 
      ifelse(str_detect(IMPACT, "LOW|MODIFIER"), 0, 
             ifelse(str_detect(IMPACT, "HIGH|MODERATE"), 1, NA)),
    
    `IMPACT (0=low/moderate/modifier; 1=high)` = 
      ifelse(str_detect(IMPACT, "LOW|MODERATE|MODIFIER"), 0, 
             ifelse(str_detect(IMPACT, "HIGH"), 1, NA))
  )

write.csv(analytical_dataset, file.path(output_path, "/analytical_dataset.csv"), row.names = FALSE)

#################################################################################
## Compile entire dataset -- not filtering for a single row @ RBMX (for downstream data viz)
#################################################################################
fulldataset <- data.frame(matrix(ncol = 27, nrow = 0))

# Create additional cols & add placeholders for data from MAF and clinical files
colnames(fulldataset) <- c("manifest_id", "CT_mutation", "GA_mutation", "Type of mutation",
                         "Num_cancer_gene_mutations", "CTs_genome", "CTs_cancer", 
                         "GAs_genome", "GAs_cancer", cols1)

process_entire_file <- function(file, manifest_id) {
  if (!file.exists(file_path)) {
    warning(paste("File does not exist:", file_path))
    return(NULL)
  }
  
  file <- fread(file_path, showProgress = TRUE)
  if (!all(colnames(file) %in% maf_cols)) {
    warning(paste("File does not contain all required columns:", file_path))
    return(NULL)
  }
  
  file <- file %>% select(all_of(cols1))
  file$manifest_id <- manifest_id
  
  file <- file %>%
    mutate(
      CT_mutation = ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "T", 1, 0),
      GA_mutation = ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "A", 1, 0),
      `Type of mutation` = ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "T", "C>T",
                           ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "A", "C>A",
                           ifelse(Reference_Allele == "C" & Tumor_Seq_Allele2 == "G", "C>G",
                           ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "T", "G>T",
                           ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "C", "G>C",
                           ifelse(Reference_Allele == "G" & Tumor_Seq_Allele2 == "A", "G>A", "Other")))))),
      CTs_genome = sum(CT_mutation[Hugo_Symbol != "RBMX"]),
      CTs_cancer = sum(CT_mutation[Hugo_Symbol %in% cancer_genes]),
      GAs_genome = sum(GA_mutation[Hugo_Symbol != "RBMX"]),
      GAs_cancer = sum(GA_mutation[Hugo_Symbol %in% cancer_genes]),
      Num_cancer_gene_mutations = sum(Hugo_Symbol %in% cancer_genes)
    ) %>%
    #filter(Hugo_Symbol == "RBMX") %>%
    select(manifest_id, Num_cancer_gene_mutations, CT_mutation, GA_mutation, `Type of mutation`,
           CTs_genome, CTs_cancer, GAs_genome, GAs_cancer, all_of(cols1))
  
  return(file)
}

for (i in rbmx_manifest$id) {
  file_path <- paste0("RBMX/", i, "/", rbmx_manifest$filename[rbmx_manifest$id == i])
  file <- fread(file_path)
  entire_processed_file <- process_entire_file(file, i)
  fulldataset <- rbind(fulldataset, entire_processed_file)
}


## Merge with select columns from analytical_dataset, to include pathogenicity information
fulldataset <- left_join(fulldataset,
                         analytical_dataset %>% 
                           select("case_id", `PolyPhen (0=benign; 1=possibly or probably damaging)`,
                                  `PolyPhen (0=benign/possibly damaging; 1=probably damaging)`,
                                  `PolyPhen [0-1]`,
                                  `SIFT (0=tolerated; 1=deleterious)`,
                                  `SIFT [0-1]`,
                                  `IMPACT (0=low/modifier; 1=moderate/high)`, 
                                  `IMPACT (0=low/moderate/modifier; 1=high)`),
                         by="case_id")
                                  

write.csv(fulldataset, file.path(output_path, "/full_dataset.csv"), row.names = FALSE)
