library(dplyr)
library(tidyr)
library(tools)
library(bio3d)
library(seqinr)

with(snakemake@input, {
	pdb_files      <<- pdbs
	opm_file       <<- opm
	sequences_file <<- sequences
})
output_file <- unlist(snakemake@output)
pdbs <- file_path_sans_ext(basename(pdb_files))
data <- lapply(pdb_files, read.pdb)

aliases <- read.fasta(sequences_file) %>%
	lapply(attr, "Annot") %>%
	unlist %>%
	data.frame %>%
	extract(1, into = c("Alias", "ID", "chain", "pos"), regex = ">([^ ]+).*\\[(.+?)\\.([A-Z]):(\\d+)\\]", convert = T) %>%
	mutate(chain = ifelse(chain == "", "A", chain))

sheets <- lapply(data, `[[`, "sheet") %>%
	lapply(data.frame) %>%
	setNames(pdbs) %>%
	bind_rows(.id = "ID") %>%
	mutate(sense = ifelse("sense" %in% names(.), as.numeric(sense), NA), sense = ifelse(sense < 0, "-", "+"))
helices <- lapply(data, `[[`, "helix") %>%
	lapply(data.frame) %>%
	setNames(pdbs) %>%
	bind_rows(.id = "ID")
tms <- read.table(opm_file, sep = "\t", header = T) %>%
	separate_rows(segments, sep = ",") %>%
	extract(segments, into = c("Name", "start", "end"), regex = "(\\d+)\\s*\\(\\s*(\\d+)\\s*-\\s*(\\d+)\\s*\\)", convert = T) %>%
	group_by(ID, chain) %>%
	mutate(Name = Name - n() + 7, Name = paste0("TM", Name))

list(sheet = sheets, helix = helices, tm = tms) %>%
	bind_rows(.id = "ss") %>%
	left_join(aliases, by = c("ID", "chain")) %>%
	filter(!is.na(Alias)) %>%
	replace_na(list(sense = ".")) %>%
	mutate(start = start - pos + 1, end = end - pos + 1) %>%
	mutate(source = "SS", score = ".", frame = ".", attrib = ifelse(is.na(Name), ".", paste0("Name=", Name))) %>%
	select(Alias, source, ss, start, end, score, sense, frame, attrib) %>%
	write.table(output_file, quote = F, sep = "\t", col.names = F, row.names = F)
