library(ape)

args = commandArgs(trailingOnly = TRUE)

input_file = args[1]
in_aligned_fasta <- read.FASTA(input_file)
distance_matrix <- dist.dna(in_aligned_fasta, model = "raw", as.matrix = T)


output_prefix <- unlist(strsplit(input_file, "[.]"))[1]
output_file = paste(output_prefix, ".dist", sep = "")

write.csv(distance_matrix, file = output_file, quote = F)