library(Ckmeans.1d.dp)
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  stop("Must supply a sorted distances file.n", call. = FALSE)
}

input_file = args[1]
input = read.table(input_file, header = 1, sep = "\t", check.names = F)

#input_vector <- apply(as.vector(input), MARGIN = 1, function(x) x * 100)
input_vector <- as.vector(input)
genome_number <- ncol(input_vector)
results <- Ckmeans.1d.dp(input_vector, k = c(1, genome_number))
#results <- Ckmedian.1d.dp(input_vector, k = c(1, genome_number))
cluster_centers <- results$centers

#Generate Output File Names
num_of_splits <- length(unlist(strsplit(input_file, "/")))

in_path <- paste(unlist(strsplit(input_file, "/"))[-num_of_splits], collapse = "/")
sample_name <- unlist(strsplit(unlist(strsplit(input_file, "/"))[num_of_splits], "_"))[1]
out_path <- paste(in_path, sample_name, sep = "/")

#Output Cluster Cutoffs File
cutoffs_file_name = paste(out_path, "_cluster_cutoffs.txt", sep = "")
write.table(t(as.table(cluster_centers)), file = cutoffs_file_name, sep = '\t', col.names = F, row.names = F)

#Output Cluster Identity File
cluster_identity_file = paste(out_path, "_cluster_identity.txt", sep = "")

out_results <- results$cluster

input[2,] <- as.integer(out_results)
write.table(input, file = cluster_identity_file, sep = '\t', quote = F, row.names = F, col.names = T)
