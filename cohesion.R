

rm(list = ls())
setwd("D:/phD/16s/16s_files/cohesion/WATER/")


library(psych) 
library(vegan)  
library(dplyr) 



# Read the OTU table (Behavior OTU/ASV, listed as sample)


cat("loading data...\n")
otu_raw <- read.table("water.txt", header = TRUE, row.names = 1, sep="\t")  
otu_raw

cat("??? Data reading completed! Original data size:", dim(otu_raw)[1], "species ×", dim(otu_raw)[2], "samples\n")  
# Convert to relative abundance
otu_rel <- sweep(otu_raw, 2, colSums(otu_raw), FUN = "/")  
otu_rel[is.na(otu_rel)] <- 0 

otu_rel

persistence_cutoff <- 0.1 
min_rel_abundance <- 0
max_otus <- 500

#====================================  
# 2. Species screening
#====================================  
cat("Species are being screened....\n")  
# Calculate the occurrence frequency of each species (persistence)  
persistence <- rowSums(otu_rel > 0) / ncol(otu_rel)  
# Select the species that meet the persistence condition  
otu_filtered <- otu_rel[persistence >= persistence_cutoff, ]  
cat("??? Filter by frequency of occurrence:", nrow(otu_filtered), "The preservation of a species (present to≥",   
    persistence_cutoff*100, "%samples)\n")  


otu_filtered <- otu_filtered[rowMeans(otu_filtered) >= min_rel_abundance, ]  
cat("??? Filter by average abundance:", nrow(otu_filtered), "Species retention (average abundance ≥",   
    min_rel_abundance, ")\n")  

if(nrow(otu_filtered) > max_otus) {  
  total_abundance <- rowSums(otu_filtered)  
  otu_filtered <- otu_filtered[order(total_abundance, decreasing=TRUE), ][1:max_otus, ]  
  cat("??? Only retain the ones with the highest abundance.", max_otus, "species\n")  
}  

cat("\nFinal screening results:\n")  
cat("??? Original number of OTUs:", nrow(otu_raw),   
    "\n??? Number of OTUs after screening:", nrow(otu_filtered),  
    "\n??? samples:", ncol(otu_filtered), "\n\n")  
#====================================  
# 3. Calculate the observation correlation matrix  
#====================================  
cat("Calculate the correlations between species...\n")  
# Calculate the Spearman correlation matrix between OTUs
obs_cor_result <- corr.test(t(otu_filtered),   
                            method = "spearman",   
                            adjust = "fdr",   # FDR校正P值  
                            alpha = 0.05)    
# Extract the correlation coefficient matrix

obs_cor <- obs_cor_result$r  
diag(obs_cor) <- 0  
cat("??? The relevant matrix calculation is complete!\n")  
#====================================  
# 4. Null model correction
#====================================  

set.seed(123)      # Random number seed to ensure repeatable results
num_perm <- 200    # Number of Replacement Operations (The paper suggests at least 100 times)
null_cors <- array(0, dim=c(nrow(obs_cor), ncol(obs_cor), num_perm))  

pb <- txtProgressBar(min=0, max=num_perm, style=3)  

# Execute the "row-shuffle" null model
for(i in 1:num_perm) {  
  shuffled_data <- apply(otu_filtered, 2, sample)  
  rownames(shuffled_data) <- rownames(otu_filtered)  
  perm_cor <- corr.test(t(shuffled_data), method="spearman", adjust="none")$r  
  diag(perm_cor) <- 0  
  null_cors[,,i] <- perm_cor  
  setTxtProgressBar(pb, i)
  }


close(pb)

# Calculate the average null correlation matrix (expected correlation) 
exp_cor <- apply(null_cors, c(1,2), mean)  
# The corrected correlation matrix = Observed correlation - Expected correlation (null model)
adj_cor <- obs_cor - exp_cor  
rownames(adj_cor) <- rownames(obs_cor)  
colnames(adj_cor) <- colnames(obs_cor)  
cat("\n??? Null model correction completed!\n")  


#====================================  
# 5. Calculate species connectivity (connectedness)
#====================================  
cat("Calculate species connectivity...\n")  
# Calculate the positive connectivity (r_pos) and negative connectivity (r_neg) of each species.
r_pos <- numeric(nrow(adj_cor))  
r_neg <- numeric(nrow(adj_cor))  
taxa_names <- rownames(adj_cor)  
for(i in 1:nrow(adj_cor)) {  
  taxon_correlations <- adj_cor[i, ]  
  pos_cors <- taxon_correlations[taxon_correlations > 0]  
  r_pos[i] <- ifelse(length(pos_cors) > 0, mean(pos_cors), 0)  
  neg_cors <- taxon_correlations[taxon_correlations < 0]  
  r_neg[i] <- ifelse(length(neg_cors) > 0, mean(neg_cors), 0)  
}  

connectedness_df <- data.frame(  
  Taxon = taxa_names,  
  r_pos = r_pos,  
  r_neg = r_neg  
)  
cat("??? Species connectivity calculation completed!\n")  

#====================================  
# 6. Calculate the sample cohesion
#====================================  
cat("Calculate the cohesion of the sample...\n")  

samples <- colnames(otu_filtered)  
n_samples <- length(samples)  
C_pos <- numeric(n_samples)  
C_neg <- numeric(n_samples)  
# Calculate the cohesion of each sample 
for(j in 1:n_samples) {  
  abund_vec <- otu_filtered[, j]  
  C_pos[j] <- sum(abund_vec * r_pos)  
  C_neg[j] <- sum(abund_vec * r_neg)  
}  
# Integration of Sample Cohesion Results
cohesion_df <- data.frame(  
  Sample = samples,  
  C_pos = C_pos,  
  C_neg = C_neg  
)  

#====================================  
# 7. Output
#====================================  
# Save the results of species connectivity
out_file1 <- "taxa_connectedness.txt"  
write.table(connectedness_df, out_file1,   
            sep="\t", row.names=FALSE, quote=FALSE)  
# Save the results of sample cohesion 
out_file2 <- "sample_cohesion.txt"  
write.table(cohesion_df, out_file2,   
            sep="\t", row.names=FALSE, quote=FALSE)  
# Save the corrected correlation matrix (optional)
out_file3 <- "adjusted_correlation_matrix.csv"  
write.csv(adj_cor, out_file3)  
# Summary of Results
cat("\n==== Summary of Species Connectivity Statistics ====\n")  
print(summary(connectedness_df[, c("r_pos", "r_neg")]))  
cat("\n==== Summary of Sample Cohesion Statistics ====\n")  
print(summary(cohesion_df[, c("C_pos", "C_neg")]))  
cat("\nThe result has been saved to the following file:\n")  
cat("???", out_file1, "- Species connectivity data\n")  
cat("???", out_file2, "- Sample cohesion data\n")  
cat("???", out_file3, "- Corrected correlation matrix\n") 

