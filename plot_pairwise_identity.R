# Load required libraries.
library(PopGenome)
library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly=TRUE)
alignment.dir <- args[1]
output.file <- args[2]

# Set window sizes for analysis.
window.size <- 500
window.slide <- 100

# Read in multiple sequence alignment.
GENOME.class <- readData(alignment.dir)

# Load all "individuals" from the multiple sequence alignment.
individuals <- unlist(get.individuals(GENOME.class))

# Create all pairs of individuals (e.g., 1 and 1, 1 and 2, etc.) and
# convert the resulting data frame into a list of vectors for PopGenome.
pairs <- expand.grid(individuals, individuals, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
pairs.list <- as.list(data.frame(t(as.matrix(pairs)), stringsAsFactors=FALSE))

# Set all pairs as the "populations" to analyze.
GENOME.class <- set.populations(GENOME.class, pairs.list)

# Calculate total segregating sites per population across windows.
GENOME.class.slide <- sliding.window.transform(GENOME.class, width=window.size, jump=window.slide, type=2, whole.data=FALSE)
GENOME.class.slide <- neutrality.stats(GENOME.class.slide)

# Use segregating sites as a proxy for mismatches in a normal pairwise alignment.
# Calculate % identity from total bases and mismatches.
mismatches <- GENOME.class.slide@n.segregating.sites
total.bases <- GENOME.class.slide@n.sites
identity <- (total.bases - mismatches) / total.bases

# Create a data frame with one row per window start and population.
position <- seq(1, length.out=length(total.bases), by=window.slide)
df <- data.frame(identity=identity, position=position)
df.melt <- melt(df, id=c("position"), variable_name="population")

# Plot one line for % identity per population.
pdf(output.file, width=8, height=6)
ggplot(df.melt, aes(y=value, x=position, colour=population)) + geom_line() + theme_bw() + ylab("% Identity") + xlab("Alignment position (bp)") + guides(colour=FALSE)
dev.off()
