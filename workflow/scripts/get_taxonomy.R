
library(dplyr)
library(tidyr)
library(tools)
library(seqinr)
library(taxize)

input_file  <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)

read.fasta(input_file) %>%
	sapply(attr, "Annot") %>%
	data.frame %>%
	extract(1, into = "TaxID", regex = "TaxID=(\\d+)", convert = T) %>%
	distinct(TaxID) %>%
	pull %>%
	classification(db = "ncbi") %>%
	`[`(!is.na(.)) %>%
	bind_rows(.id = "Organism.ID") %>%
	write.table(output_file, sep = "\t", row.names = F, col.names = T, quote = F)
