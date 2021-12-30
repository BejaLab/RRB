
library(ggtree)
library(dplyr)
library(tidyr)
library(treeio)
library(tools)
library(readxl)
library(seqinr)
library(ggnewscale)
library(pals)
library(taxize)

input  <- snakemake@input
output <- snakemake@output

tree_file     <- input$tree
selected_8TMs_file <- input$selected_8TMs
taxa_file     <- input$taxa
fasta_file    <- input$fasta
interpro_file <- input$interpro
domains_file  <- input$domains
taxonomy_file <- input$taxonomy

output_file   <- unlist(output)

clades.ignore <- c(NA, "chromerid-rhodopsins", "PK-fusions", "RB-like", "9TM", "bestrophins-complete")
selected_8TMs <- read.fasta(selected_8TMs_file) %>%
	sapply(attr, "Annot") %>%
	data.frame(label = names(.)) %>%
	extract(1, into = "TaxID", regex = "TaxID=(\\d+)", convert = T)
taxa <- read.table(taxa_file, header = T, sep = "\t", fill = T)
aliases <- list()

ncbi.tax.all <- read.table(taxonomy_file, header = T, sep = "\t", comment.char = "")

ncbi.tax <- left_join(ncbi.tax.all, taxa, by = c("name", "rank")) %>%
	group_by(Organism.ID) %>%
	mutate(level_id = last(which(!is.na(level))), Group = name[level_id + level[level_id]], Group = ifelse(is.na(Group), last(name), Group)) %>%
	mutate(Superkingdom = first(name[rank == "superkingdom"])) %>%
	ungroup %>%
	distinct(TaxID = Organism.ID, Group, Superkingdom)

genera <- list(
	Chromera        = "Colpodellida",
	Vitrella        = "Colpodellida",
	Karlodinium     = "Dinophyceae",
	Dinophysis      = "Dinophyceae",
	Paragymnodinium = "Dinophyceae",
	Gymnodinium     = "Dinophyceae",
	dinoflagellate  = "Dinophyceae",
	Cymbomonas      = "Viridiplantae",
	Phaeocystis     = "Haptista"
)

rhodopsins <- read.fasta(fasta_file) %>%
	sapply(attr, "Annot") %>%
	data.frame(label = names(.)) %>%
	extract(1, into = "Clade",   regex = "\\{(.+)\\}", remove = F) %>%
	extract(1, into = "Organism", regex = "\\[(.+)\\]") %>%
	mutate(Organism = sub(" .+", "", sub("Uncultured ", "", Organism))) %>%
	mutate(Group = recode(Organism, !!!genera))

bestrophin_domains <- filter(rhodopsins, grepl("bestrhodopsins", Clade)) %>%
	mutate(Domain = "Best") %>%
	select(Protein.accession = label, Domain)

#pos <- c(57, 61, 169)
#triad <- read.fasta("rhodopsins.trimal.K_plus_10.fasta", seqtype = "AA") %>%
#	lapply(`[`, pos) %>%
#	lapply(`%in%`, c("E","D")) %>%
#	lapply(all) %>%
#	unlist %>%
#	`[`(. == T) %>%
#	names

tree <- read.iqtree(tree_file) %>%
	as_tibble %>%
	left_join(rhodopsins, by = "label") %>%
	left_join(selected_8TMs, by = "label") %>%
	left_join(ncbi.tax, by = "TaxID") %>%
	unite(Group, Group.x, Group.y, na.rm = T) %>%
	mutate(Group = ifelse(Group == "", Clade, Group)) %>%
	#mutate(Triad = label %in% triad) %>%
	`class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))

col.names <- c(
	"Protein.accession",
	"Sequence.MD5.digest",
	"Sequence.length",
	"Analysis",
	"Signature.accession",
	"Signature.description",
	"Start.location",
	"Stop.location",
	"Score",
	"Status",
	"Date",
	"InterPro.annotations.accession",
	"InterPro.annotations.description"
)

domains <- read.table(domains_file, header = T, sep = "\t")
ipr <- read.table(interpro_file, sep = "\t", quote = "", col.names = col.names) %>%
	left_join(domains, by = "InterPro.annotations.accession") %>%
	mutate(Protein.accession = sub("UniRef50_", "", Protein.accession)) %>%
	filter(!is.na(Domain), Domain != "Rhodopsin") %>%
	bind_rows(bestrophin_domains) %>%
	distinct(label = Protein.accession, Domain)

add_data <- function(p, data, nudge) {
	p$data <- left_join(p$data, data, by = "label") %>%
		group_by(node) %>%
		mutate(num = 1:n(), angle.rad = angle * pi / 180, hjust = num * nudge * cos(angle.rad), vjust = num * nudge * sin(angle.rad))
	return(p)
}

p <- ggtree(as.treedata(tree), layout = "ape") %>% add_data(ipr, 0.05) +
	# scale_linetype_manual(values = c(Bacteria = "dotted", Eukaryota = "dashed"), na.value = "solid") +
	# scale_color_manual(values = cols25(unique(tree$Group) %>% length), na.value = "gray40") +
	geom_nodepoint(aes(x = branch.x, y = branch.y, subset = !is.na(UFboot) & UFboot >= 95), color = "#4d4d4dff") +
	geom_nodepoint(aes(x = branch.x, y = branch.y, subset = !is.na(UFboot) & UFboot >= 90 & UFboot < 95), color = "#b3b3b3ff") +
	new_scale_color() +
	geom_point2(aes(subset = !is.na(Domain), x = x + hjust, y = y + vjust, color = Domain, shape = Domain)) +
	# geom_tippoint(aes(subset = !is.na(Clade), color = Clade)) +
	geom_treescale(width = 1) # +
	#geom_tippoint(aes(subset = is.R_R_B), color = "red", shape = "✱", size = 5) +
	#geom_tippoint(aes(subset = !is.na(cTP), shape = cTP), color = "green", size = 5) +
	#geom_tiplab(aes(subset = !is.na(Gene) & Gene != "", label = Gene), size = 5) +
	#scale_shape_manual(values = list(primary = -0x25CF, secondary = -0x25C9)) +

ggsave(output_file, p, width = 20, height = 10)
