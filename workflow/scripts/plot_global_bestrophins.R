
library(ggtree)
library(dplyr)
library(tidyr)
library(treeio)
library(tools)
library(readxl)
library(seqinr)
library(ggplot2)
library(bioformatr)
library(pals)

input  <- snakemake@input
output <- snakemake@output

tree_file     <- input$tree
tab_file      <- input$tab
fasta_file    <- input$fasta
clstr_file    <- input$clstr
genes_file    <- input$genes
targetp_file  <- input$targetp
asafind_file  <- input$asafind
blast_file    <- input$blast
fasta_file    <- input$fasta
taxa_file     <- input$taxa
taxonomy_file <- input$taxonomy
output_file   <- unlist(output)

extract_uniref_cluster <- function(.data, col.name) {
	extract(.data, !!col.name, into = col.name, regex = "UniRef50_(.+)")
}
extract_pfam_cluster <- function(.data, col.name) {
	extract(.data, !!col.name, into = col.name, regex = "(.+)/")
}

uniref50 <- read.table(tab_file, sep = "\t", header = T, quote = "") %>%
	separate_rows(Cluster.members, sep = "; ") %>%
	extract_uniref_cluster("Cluster.ID") %>%
	select(Entry = Cluster.members, UniRef50_Cluster = Cluster.ID)

cdhit <- read.cdhit.clstr(clstr_file) %>%
	extract_pfam_cluster("Seq.Name") %>%
	extract_pfam_cluster("Representative") %>%
	select(UniRef50_Cluster = Seq.Name, Cluster = Representative)

refs <- read.table(genes_file, header = T, fill = T, sep = "\t") %>%
	left_join(uniref50, by = "Entry") %>%
	group_by(UniRef50_Cluster) %>%
	summarize(Gene = paste(Gene, collapse = "; "))

targetp <- read.table(targetp_file, skip = 1, comment.char = "", sep = "\t", header = T) %>%
	rename(UniRef50_Cluster = 1) %>%
	extract_uniref_cluster("UniRef50_Cluster") %>%
	filter(Prediction == "cTP") %>%
	pull(UniRef50_Cluster)

asafind <- read.table(asafind_file, header = T, sep = "\t") %>%
	separate(ASAFind.Prediction, into = c("Prediction", "Confidence"), sep = ", ") %>%
	filter(Prediction == "Plastid", Confidence == "high confidence") %>%
	extract_uniref_cluster("Identifier") %>%
	pull(Identifier)

RRBs <- read.outfmt6(blast_file) %>%
	distinct(qseqid) %>%
	extract_uniref_cluster("qseqid") %>%
	pull

bestrophins <- read.fasta(fasta_file) %>%
	sapply(attr, "Annot") %>%
	data.frame(UniRef50_Cluster = names(.), desc = .) %>%
	extract_uniref_cluster("UniRef50_Cluster") %>%
	left_join(refs, by = "UniRef50_Cluster") %>%
	mutate(cTP = case_when(UniRef50_Cluster %in% targetp ~ "primary", UniRef50_Cluster %in% asafind ~ "secondary")) %>%
	mutate(is.R_R_B = UniRef50_Cluster %in% RRBs) %>%
	extract(desc, into = c("description", "n", "dummy", "Tax", "TaxID", "RepID"), regex = ">[^ ]+ (.+) n=(\\d+)( Tax=(.+) TaxID=(\\d+))? RepID=(.+)", convert = T) %>%
	mutate(is.fragment = grepl("\\(Fragment\\)", description)) %>%
	select(-dummy) %>%
	left_join(cdhit, by = "UniRef50_Cluster") %>%
	filter(!is.na(Cluster))

taxa <- read.table(taxa_file, header = T, sep = "\t", fill = T) %>%
	mutate(no.chloroplasts = !is.na(no.chloroplasts) & no.chloroplasts == 1)
aliases <- list(
	Amoebozoa = "Other unikonts",
	Apusozoa  = "Other unikonts",
	Choanoflagellata = "Other unikonts",
	Ichthyosporea = "Other unikonts",
	Glaucocystophyceae = "Archaeplastida",
	# Cryptophyceae = "Archaeplastida",
	Viridiplantae = "Archaeplastida",
	Rhodophyta = "Archaeplastida"
)

ncbi.tax.all <- read.table(taxonomy_file, header = T, sep = "\t", comment.char = "")
ncbi.tax <- left_join(ncbi.tax.all, taxa, by = c("name", "rank")) %>%
	group_by(Organism.ID) %>%
	mutate(level_id = last(which(!is.na(level))), Group = name[level_id + level[level_id]], Group = ifelse(is.na(Group), last(name), Group)) %>%
	mutate(no.chloroplasts = any(no.chloroplasts)) %>%
	mutate(Superkingdom = first(name[rank == "superkingdom"])) %>%
	distinct(TaxID = Organism.ID, Group, Superkingdom, no.chloroplasts) %>%
	mutate(Group = recode(Group, !!!aliases))

proteins <- left_join(bestrophins, ncbi.tax, by = "TaxID") %>%
	mutate(cTP = ifelse(!is.na(no.chloroplasts) & no.chloroplasts, NA, cTP)) %>%
	group_by(Cluster, Group) %>%
	mutate(n.group = n()) %>%
	group_by(Cluster, cTP) %>%
	mutate(n.cTP = n()) %>%
	group_by(Cluster) %>%
	summarize(Gene = na.omit(Gene) %>% paste(collapse = "; ") %>% first, Group = Group[which.max(n.group)], cTP = cTP[which.max(n.cTP)], is.R_R_B = first(is.R_R_B))

tree <- read.iqtree(tree_file) %>%
	as_tibble %>%
	mutate(Cluster = label) %>%
	extract_pfam_cluster("Cluster") %>%
	left_join(proteins, by = "Cluster") %>%
	# mutate(Group_common = ifelse(n > 3, Group, NA)) %>%
	`class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))

p <- ggtree(as.treedata(tree), aes(color = Group), layout = "ape") +
	# scale_linetype_manual(values = c(Bacteria = "dotted", Eukaryota = "dashed"), na.value = "solid") +
	geom_nodepoint(aes(x = branch.x, y = branch.y, subset = !is.na(UFboot) & UFboot >= 95)) +
	geom_tippoint(aes(subset = is.R_R_B), color = "red", shape = "✱", size = 5) +
	geom_tippoint(aes(subset = !is.na(cTP), shape = cTP), color = "green", size = 5) +
	geom_tiplab(aes(subset = !is.na(Gene) & Gene != "", label = Gene), size = 5) +
	scale_shape_manual(values = list(primary = -0x25CF, secondary = -0x25C9)) +
	scale_color_manual(values = cols25(unique(tree$Group) %>% length), na.value = "gray40") +
	geom_treescale(width = 1)

ggsave(output_file, p, device = cairo_pdf, width = 20, height = 10)
