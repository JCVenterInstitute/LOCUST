suppressPackageStartupMessages(library(mclust))
library(ape)


args = commandArgs(trailingOnly = TRUE)

input_file = args[1]
in_aligned_fasta <- read.FASTA(input_file)
distance_matrix <- dist.dna(in_aligned_fasta, as.matrix = T)
write.table(distance_matrix, file = "all_vs_all.dist", sep = "\t", quote = F, row.names = TRUE, col.names = NA)
mds_results <- cmdscale(distance_matrix, k = 2)
mc <- Mclust(mds_results)

column_names <- list()
for (i in 1:ncol(mc$z)){
  outName = paste("Cluster ", i, sep = "")
  column_names <- c(column_names, outName)
}
colnames(mc$z) <- column_names
outMatrix <- cbind(distance_matrix, Cluster = mc$classification)
outMatrix <- cbind(outMatrix, uncertainty = mc$z)

write.table(outMatrix, file = "mclust_results.txt", sep = "\t", quote = F, row.names = TRUE, col.names = NA)
