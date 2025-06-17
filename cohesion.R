

rm(list = ls())
setwd("D:/phD/16s/16s_files/cohesion/gut/")

# Calculate the number of zeros in the vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}
# Calculate the average value of the negative values in the vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

# Calculate the average value of the positive values in the vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}


############## Workflow options ########################
##Set the classification persistence threshold (minimum proportion of existence) to retain the classification in the analysis

pers.cutoff <- 0.10

## Set the number of iterations of the virtual model (it is recommended to be >=200)
iter <- 200

##Choose category shuffling (tax.shuffle= T) or row shuffling (tax.shuffle= F£©
tax.shuffle <- T

#Choose whether to use a custom correlation matrix
# Note: Your correlation table must have the same number of classifications as the abundance table. There should be no empty (all zero) classification vectors in the abundance table.
#Even if you input a custom correlation table, the persistence threshold will still be applied

use.custom.cors <- F

########### Calculate cohesion########################

b <-
  read.csv("gut.csv", header = T, row.names = 1)

# If a custom correlation matrix is used, please read and check the dimensions

if(use.custom.cors == T) {
  custom.cor.mat <- read.csv("content.csv", header = T,
                             row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  print(dim(b)[2] == dim(custom.cor.mat)[2])
  }

# Reformat the data to remove empty samples and classifications
c <- as.matrix(b)
c <- c[rowSums(c) > 0,
       colSums(c) > 0]
# The total number of individuals that preserve the original sample
rowsums.orig <- rowSums(c)

# Determine the quantitative threshold of classification zeros based on the persistence threshold
zero.cutoff <- ceiling(pers.cutoff
                       * dim(c)[1])

# Delete the classifications that are below the persistence threshold
d <- c[ , apply(c, 2, zero) <
          (dim(c)[1]-zero.cutoff) ]

# Delete the samples without individuals
d <- d[rowSums(d) > 0, ]

# If a custom correlation matrix is used, update the correlation matrix
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < 
                                         (dim(c)[1]-zero.cutoff), 
                                       apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
  }


# Create the relative abundance matrix
rel.d <- d / rowsums.orig

hist(rowSums(rel.d))
# Calculate the observed correlation matrix
cor.mat.true <- cor(rel.d)
# Save the median correlation of each category
med.tax.cors <- vector()

# Calculate the expected correlation of the virtual model
# If a custom correlation matrix is used, the virtual model is skipped

if(use.custom.cors == F) {
  if(tax.shuffle) {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      # Save the correlation of each permutation
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        # Create an empty matrix with the same dimension as rel.d
        perm.rel.d <- matrix(numeric(0),
                             dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <-
          rownames(rel.d)
        colnames(perm.rel.d) <-
          colnames(rel.d)
        
        # Shuffle each category
        for(j in 1:dim(rel.d)[2]){
          perm.rel.d[, j ] <- sample(rel.d[
            ,j ])
        }
        # Keep the focus column unchanged
        perm.rel.d[, which.taxon] <- rel.d[
          , which.taxon]
        
        # Calculate the correlation matrix of the permutation matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # Save the correlation of the focus classification
        perm.cor.vec.mat <-
          cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlation
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1,
                                                median))
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
    
# Calculate the observed minus the expected correlation
if(use.custom.cors == T) {
  obs.exp.cors.mat <- custom.cor.mat.sub
  } else {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }


diag(obs.exp.cors.mat) <- 0


#### Generate the vectors of connectivity and cohesion
# Calculate the connectivity (the average value of positive and negative correlations)
connectedness.pos <-
  apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <-
  apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate the degree of cohesion (the product of relative abundance and connectivity)
cohesion.pos <- rel.d %*%
  connectedness.pos
cohesion.neg <- rel.d %*%
  connectedness.neg


####results
output <- list(connectedness.neg,
               connectedness.pos, cohesion.neg, cohesion.pos)

names(output) <- c("Negative Connectedness",
                   "Positive Connectedness", "Negative Cohesion",
                   "Positive Cohesion")
print(output)

#####Connectedness############

connectedness_data <- list(
  "Negative Connectedness" = output[["Negative
                                     Connectedness"]],
  "Positive Connectedness" = output[["Positive
                                     Connectedness"]]
  )

# Convert to a data frame and ensure that the number of rows is the same
connectedness_df <-
  as.data.frame(do.call(cbind, connectedness_data))

# Set names for each column
colnames(connectedness_df) <-names(connectedness_data)
# CSV files
write.csv(connectedness_df, file ="connectedness_data.csv", row.names = TRUE)

###############Output Cohesion data ############ extract the matrices or vectors to be merged
cohesion_data <- list(
  "Negative Cohesion" = output[["Negative Cohesion"]],
  "Positive Cohesion" = output[["Positive Cohesion"]]
  )
# Convert to a data frame and ensure that the number of rows is the same
cohesion_df <-
  as.data.frame(do.call(cbind, cohesion_data))

colnames(cohesion_df) <-names(cohesion_data)

write.csv(cohesion_df, file ="cohesion.csv", row.names = TRUE)

