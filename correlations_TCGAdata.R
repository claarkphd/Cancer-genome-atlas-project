##########################################################
## Project: RBMX mutations
## Author: J Clark
## Date: June 21st, 2024
##########################################################

##########################################################
## Set working directory
##########################################################
setwd("C:/Users/jclark/iCloudDrive/The Hub/Projets/P5_RBMX Cancer database/analysis/21062024/output/")
output_path <- "C:/Users/jclark/iCloudDrive/The Hub/Projets/P5_RBMX Cancer database/analysis/21062024/output/"
set.seed(326)

#################################################################################
#################################################################################
##                    L O A D   P A C K A G E S   A N D
##                                  D A T A 
#################################################################################
#################################################################################
library(data.table)
library(tidyverse)
library(readxl)
library(ggplot2)
#install.packages("hrbrthemes")
#install.packages("viridis")
#install.packages("lawstat")
#install.packages("table1")
library(hrbrthemes)
library(viridis)
library(lawstat)
library(table1)


# LOAD RBMX DATA
rbmx <- fread("analytical_dataset.csv") 
which(duplicated(rbmx$manifest_id) == TRUE) # duplicates have been removed

#################################################################################
#################################################################################
##                  S P E A R M A N   C O R R E L A T I O N S
##                            R B M X    A C T I V I T Y
#################################################################################
#################################################################################
variable <- c("PolyPhen (0=benign; 1=possibly or probably damaging)",
              "PolyPhen (0=benign/possibly damaging; 1=probably damaging)",
              "PolyPhen [0-1]",
              "SIFT (0=tolerated; 1=deleterious)",
              "SIFT [0-1]",
              "IMPACT (0=low/modifier; 1=moderate/high)", 
              "IMPACT (0=low/moderate/modifier; 1=high)")
mutation <- c("CTs_genome", "GAs_genome", "CTs_cancer", "GAs_cancer")
results <- setNames(data.frame(matrix(ncol = 4, nrow = length(variable)*length(mutation))), c("Mutation Type", "Variable", 
                                                              "Spearman correlation (ρ)", "p-value")) 
row <- 1
# Correlations - CT
for (v in variable){
  for (m in mutation){
    print(paste(v,m))
    results$`Mutation Type`[row] <- m
    results$Variable[row] <- v
    results$`Spearman correlation (ρ)`[row] <- eval(parse(text = paste0("cor.test(x = rbmx$`",v,"`, y = rbmx$",m,", method = 'spearman')$estimate")))
    results$`p-value`[row] <- eval(parse(text = paste0("cor.test(x = rbmx$`",v,"`, y = rbmx$",m,", method = 'spearman')$p.value")))
    row <- row + 1
  }
}
write.csv(results, paste0(output_path, "/spearman_correlations.csv"), row.names = FALSE)


#################################################################################
#################################################################################
##                              T A B L E   O N E 
##                           
#################################################################################
#################################################################################
rbmx$age <- as.numeric(rbmx$age_at_diagnosis)/365
colnames(rbmx)
# Construct table 1
#General formula: table1(~ var1 + var2 | marginalvar or nestvar*nestvar, data = rbmx)
Table1 <- table1(~ gender + race + age + primary_diagnosis +
         CTs_genome + GAs_genome + CTs_cancer + GAs_cancer + 
         Variant_Classification + factor(`PolyPhen (0=benign; 1=possibly or probably damaging)`)+ 
         factor(`PolyPhen (0=benign/possibly damaging; 1=probably damaging)`)+ 
         `PolyPhen [0-1]`+ 
         factor(`SIFT (0=tolerated; 1=deleterious)`)+ 
         `SIFT [0-1]`+ 
         factor(`IMPACT (0=low/modifier; 1=moderate/high)`)+  
         factor(`IMPACT (0=low/moderate/modifier; 1=high)`), data = rbmx) %>%
  as.data.frame()
write.csv(Table1, paste0(output_path,"/Table1.csv"))
