
# Clear all objects
rm(list = ls())
wd = "D:/phD/16s·ÖÎöÄ£¿é/16s_files/icamp/water/"
# Define the CRAN and Bioconductor packages to be installed
cran_packages <- c("ggplot2", "ggstream", "dplyr", "tidyr")   # CRAN packages
bioc_packages <- c("phyloseq", "ggtree", "ggtreeExtra")       # Bioconductor packages

# Check and install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg)
  }
}

# Install Bioconductor packages if BiocManager is not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}

# Check and install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg)
  }
}

# Load required libraries
library(ggplot2)
library(ggstream)
library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(dplyr)
library(tidyr)

# Create a physeq object; below is an example with dummy data
otu <- read.table("otus.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
taxa <- read.table("classification.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
metadata <- read.table("treat2col.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
colnames(metadata) <- "site" # Commonly used terms: group or treatment
tree <- read.tree("tree.tre")

# Data processing
otu <- as.matrix(otu)
otu_table_obj <- otu_table(otu, taxa_are_rows = TRUE)
taxa <- as.matrix(taxa)
tax_table_obj <- tax_table(taxa)
sample_data_obj <- sample_data(metadata)

physeq <- phyloseq(otu_table_obj, tax_table_obj, sample_data_obj, tree)



# Check if the created physeq object is correct
summary(physeq)
nsamples(physeq)  # Number of samples
ntaxa(physeq)     # Number of ASVs

# Read bins information; data source from iCAMP analysis
# For details, refer to the link: 
# https://mp.weixin.qq.com/s/mnSfdFca24XLZOg_E04byQ
bins <- read.csv("Testoutput2/Test.Bin_TopTaxon.csv", row.names = 1, header = TRUE)
bins_assembly <- read.csv("Testoutput2/Test.ProcessImportance_EachBin_EachGroup.csv", header = TRUE)
bins_assembly <- bins_assembly[1:5, 4:ncol(bins_assembly)]

# Extract the "core" ASVs of bins; this will not affect the overall tree structure
# The resulting tree will generally have hundreds of branches
bin_otu_ids <- bins[, "TopTaxonID"]

# Filter the main ASVs of each bin, significantly reducing the size of the physeq object
physeq <- prune_taxa(bin_otu_ids, physeq) # Retain the main ASVs of bins

# Extract main ASVs from bins
bins_map <- bins[, "TopTaxonID", drop = FALSE]

# Sort by ASVs in the physeq object
bins_map <- bins_map[match(taxa_names(physeq), bins_map$TopTaxonID), , drop = FALSE]

# Ensure that ASVs in physeq and bins_map are consistent
identical(bins_map$TopTaxonID, taxa_names(physeq))

# Rename ASVs to bin names
taxa_names(physeq) <- rownames(bins_map)

# Validate that there are no duplicates in taxa_names or tip.labels
length(unique(taxa_names(physeq))) == length(taxa_names(physeq))
length(unique(phy_tree(physeq)$tip.label)) == length(phy_tree(physeq)$tip.label)
validObject(physeq)
head(taxa_names(physeq))

# Prepare data for assembly
bins_assembly <- as.data.frame(t(bins_assembly))
colnames(bins_assembly) <- as.character(unlist(bins_assembly[1, ]))

# Convert case: Bin -> bin. Case mismatch often causes errors and can be easily overlooked
rownames(bins_map) <- tolower(rownames(bins_map))

# Reorder bins_assembly to match the order in bins_map
bins_assembly <- bins_assembly[match(rownames(bins_map), rownames(bins_assembly)), , drop = FALSE]

# Check if row names of bins_map and bins_assembly are consistent
setdiff(rownames(bins_map), rownames(bins_assembly))  # If the result is 0, the row names match
setdiff(rownames(bins_assembly), rownames(bins_map))  # If the result is 0, the row names match

# Create a new assembly data frame with row names added as a new column
bins_assembly_new <- bins_assembly %>%
  tibble::rownames_to_column(var = "bin")  
bins_assembly_new$bin <- sub("^b", "B", bins_assembly_new$bin)

# Convert columns to numeric format
bins_assembly_new$HeS <- as.numeric(bins_assembly_new$HeS)
bins_assembly_new$HoS <- as.numeric(bins_assembly_new$HoS)
bins_assembly_new$DL <- as.numeric(bins_assembly_new$DL)
bins_assembly_new$HD <- as.numeric(bins_assembly_new$HD)
bins_assembly_new$DR <- as.numeric(bins_assembly_new$DR)

# Import relative abundance data for bins
bins_RA <- read.csv("Testoutput2/Test.Bin_TopTaxon.csv", row.names = 1, header = TRUE) 

# Convert row names to lowercase
rownames(bins_RA) <- tolower(rownames(bins_RA))

# Reorder bins_RA to match the order in bins_map
bins_RA <- bins_RA[rownames(bins_map), , drop = FALSE]
rownames(bins_RA) <- sub("^b", "B", rownames(bins_RA))

# Check if row names of bins_RA and taxa_names in physeq are identical
identical(rownames(bins_RA), taxa_names(physeq))

# Create a new bins_RA data frame with row names added as a new column
bins_RA_new <- bins_RA %>%
  tibble::rownames_to_column(var = "bin")

# Filter the top 10 phyla for visualization preparation
summarized_data <- bins_RA_new %>%
  group_by(TopTaxon.Phylum) %>%
  summarise(total_RA = sum(BinRA)) %>%
  ungroup() %>%
  arrange(desc(total_RA))  # Arrange in descending order of total relative abundance

# Select the top 10 phyla (example data has only two phyla)
top_phyla <- summarized_data$TopTaxon.Phylum[1:10]

# Create a new column to indicate whether each bin's phylum is in the top 10
bins_RA_new <- bins_RA_new %>%
  mutate(Phylum = ifelse(TopTaxon.Phylum %in% top_phyla, TopTaxon.Phylum, "others"))

# Extract the tax_table, ensuring it's in data.frame format
tax_table_df <- as.data.frame(physeq@tax_table)
original_rownames <- rownames(tax_table_df)
original_colnames <- colnames(tax_table_df)

# Update the Phylum column, setting low-abundance groups to "others"
tax_table_df$Phylum <- ifelse(tax_table_df$Phylum %in% top_phyla, tax_table_df$Phylum, "others")

# Create a new physeq object
new_physeq <- physeq

# Assign the modified tax_table back to the phyloseq object
new_physeq@tax_table <- tax_table(tax_table_df)

# Restore the row and column names of the tax_table
rownames(new_physeq@tax_table) <- original_rownames
colnames(new_physeq@tax_table) <- original_colnames

# Organize metadata for bins
bins_info <- cbind(bins_assembly, bins_RA_new$BinRA)
bins_info <- cbind(bins_info, bins_RA_new$Phylum)
rownames(bins_info) <- bins_assembly_new$bin
colnames(bins_info) <- c("HeS", "HoS", "DL", "HD", "DR", "binsRA", "Phylum")

# Extract the tree from the physeq object
new_tree <- phy_tree(physeq)
# Initialize an isolated environment using the renv package
renv::init()
# Use the picante package for the match.phylo.data function
# Note: picante may have compatibility issues with newer R versions, so switching to an isolated environment with R-4.2.2 is necessary.
library(picante)
# Match the phylogenetic tree with the metadata
metaplot <- match.phylo.data(new_tree, bins_info$Phylum)
# Exit the isolated environment
renv::deactivate()

library(phyloseq)
library(dplyr)
library(tidyr)


# Create the tree with additional information
treeWithInfo <- metaplot$phy 

###Save the tree file##
library(ape)
write.tree(treeWithInfo, file = "tree_file.nwk")

