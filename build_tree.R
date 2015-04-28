library(ape)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.tree <- args[2]
output.plot <- args[3]

dna <- read.dna(input.file, format="fasta")

tree <- nj(dist.dna(dna))
write.tree(tree, output.tree)

pdf(output.plot, width=8, height=6)
plot(tree)
dev.off()
