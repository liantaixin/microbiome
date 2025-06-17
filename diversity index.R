setwd("D:/phD/16s/16s_files")

library(vegan)
library(picante)

otu <- read.delim('feature-table-otus.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
View(otu)

 
#observed species/Richness  index
observed_species <- estimateR(otu)[1, ]

observed_species

View(observed_species)

write.csv(observed_species, file = "observed_species.csv", row.names = TRUE)

#Chao 1 index

Chao1  <- estimateR(otu)[2, ]

Chao1
write.csv(Chao1, file = "chao1.csv", row.names = TRUE)

#Shannon index
Shannon <- diversity(otu, index = 'shannon', base = exp(1))
Shannon
write.csv(Shannon, file = "Shannon.csv", row.names = TRUE)

#Simpson index
Gini_simpson  <- diversity(otu, index = 'simpson')

Gini_simpson
write.csv(Gini_simpson, file = "Gini_simpson.csv", row.names = TRUE)


#coverage rate
goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu)
goods_coverage
write.csv(goods_coverage, file = "goods_coverage.csv", row.names = TRUE)
