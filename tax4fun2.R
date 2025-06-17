
install.packages(pkgs = 'Tax4Fun2_1.1.5.tar.gz', repos = NULL, source = TRUE)

library(Tax4Fun2)

###samples

# Installation Dependencies
buildDependencies(path_to_reference_data = './Tax4Fun2_ReferenceData_v2', use_force = T, install_suggested_packages = TRUE)

# Species Annotation
runRefBlast(path_to_otus = 'rep-sequences.fasta', path_to_reference_data = './Tax4Fun2_ReferenceData_v2', path_to_temp_folder = 'bac_Ref99NR', database_mode = 'Ref99NR', use_force = TRUE, num_threads = 12)



# Function Prediction ko
makeFunctionalPrediction(path_to_otu_table = 'out_table_Hem.tsv', path_to_reference_data = './Tax4Fun2_ReferenceData_v2', path_to_temp_folder = 'bac_Ref99NR-hem',database_mode = 'Ref99NR', normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)


# Calculate the Functional Redundancy Index FRI
calculateFunctionalRedundancy(path_to_otu_table = "out_table_Hem.tsv", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "bac_Ref99NR-hem", database_mode = "Ref99NR", min_identity_to_reference = 0.97)
