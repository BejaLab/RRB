
library(ape)
library(dplyr)
library(tidyr)

input <- snakemake@input
output <- snakemake@output

astral_files <- input$astral
erable_files <- input$erable
output_file  <- unlist(output)

trees <- lapply(1:length(astral_files), function(taxon) {
	astral <- read.tree(astral_files[[taxon]]) %>%
		collapse.singles
	erable <- read.tree(erable_files[[taxon]]) %>%
		collapse.singles
	astral.labels <- comparePhylo(erable, astral)$NODES %>%
		data.frame %>%
		extract(erable, into = c("erable.label", "erable.node"), regex = "(.*?) ?\\((.+)\\)", convert = F) %>%
		extract(astral, into = c("astral.label", "astral.node"), regex = "(.*?) ?\\((.+)\\)", convert = F)
	mismatches <- filter(astral.labels, astral.label != erable.label, erable.label != "")
	if (nrow(mismatches) > 1) {
		write("Labels mismatch between astral and erable", stderr())
		q(status = 1)
	}
	erable$node.label <- astral.labels$astral.label
	return(erable)
}) %>%
	lapply(chronos) %>%
	lapply(`class<-`, "phylo")

write.tree(do.call(c, trees), output_file)
