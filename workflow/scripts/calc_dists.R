
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(phytools))

input  <- snakemake@input
output <- snakemake@output
params <- snakemake@params

tree.files <- sort(params$trees)
alig.files <- sort(input$mafft)
dist.file  <- unlist(output)

alig.names <- basename(dirname(alig.files))
alig.lens <- lapply(alig.files, read.fasta) %>%
	lapply(`[`, 1) %>%
	lapply(unlist) %>%
	lapply(length) %>%
	setNames(alig.names)

tree.dists <- lapply(tree.files, read.tree) %>%
	do.call(c,.) %>%
	lapply(cophenetic) %>%
	lapply(as.table) %>%
	setNames(basename(dirname(tree.files)))

# Write distance matrices for ERABLE
lapply(names(alig.lens), function(gene.name) c(
	"",
	paste(
		ncol(tree.dists[[gene.name]]), alig.lens[[gene.name]]
	),
	capture.output(
		write.table(tree.dists[[gene.name]], sep = "\t", col.names = F, quote = F)
	)
)) %>% {c(length(.), unlist(.))} %>%
	cat(file = dist.file, sep = "\n")
