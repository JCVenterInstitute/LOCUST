suppressPackageStartupMessages(library(mclust))
library(ape)
library(ggplot2)


args = commandArgs(trailingOnly = TRUE)

input_file = args[1]
in_aligned_fasta <- read.FASTA(input_file)
distance_matrix <- dist.dna(in_aligned_fasta, as.matrix = T)
write.table(distance_matrix, file = "all_vs_all.dist", sep = "\t", quote = F, row.names = TRUE, col.names = NA)
mds_results <- cmdscale(distance_matrix, k = 5)
colnames(mds_results) <- c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5")
mc <- Mclust(mds_results)

write.table(mds_results, file = "mds_results.txt", sep = "\t", quote = F, row.names = TRUE, col.names = NA)

png("mclust.png")
plot.Mclust(mc, what = 'density')

column_names <- list()
for (i in 1:ncol(mc$z)){
  outName = paste("Cluster ", i, sep = "")
  column_names <- c(column_names, outName)
}
colnames(mc$z) <- column_names
outMatrix <- cbind(distance_matrix, Cluster = mc$classification)
outMatrix <- cbind(outMatrix, uncertainty = mc$z)

write.table(outMatrix, file = "mclust_results.txt", sep = "\t", quote = F, row.names = TRUE, col.names = NA)
