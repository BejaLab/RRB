
library(ape)
library(dplyr)
library(tidyr)
library(adephylo)
library(phytools)
library(treeio)
library(ggtree)
library(ggstance)
library(ggplot2)
library(ggpubr)

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

busco_files  <- input$busco
species_file <- input$species

astral_file  <- input$astral
erable_file  <- input$erable

newick_file <- output$newick
svg_file    <- output$svg

root.ratio <- params$root_ratio
width  <- params$width
height <- params$height
output_file <- unlist(output)

busco_names <- busco_files %>% # {path}/{name}/run_{lineage}/full_table.tsv
		dirname %>%    # {path}/{name}/run_{lineage}
		dirname %>%    # {path}/{name}
		basename       # {name}

species <- read.table(species_file, header = T, sep = "\t", fill = T) %>%
	rename(label = Name)

busco <- lapply(busco_files, read.table, sep = "\t", fill = T, quote = "", col.names = c("Busco_id", "Status", "Sequence", "Score", "Length", "url", "Description")) %>%
	setNames(busco_names) %>%
	bind_rows(.id = "label") %>%
	distinct(label, Busco_id, Status) %>%
	mutate(Status = recode(Status, Duplicated = "Complete")) %>%
	group_by(label, Status) %>%
	summarize(n = n(), .groups = "drop_last") %>%
	mutate(pct = n / sum(n) * 100) %>%
	ungroup

astral <- read.astral(astral_file)
astral@phylo$node.label <- unite(astral@data, "united", sep = ";") %>% pull
erable <- read.tree(erable_file) %>%
	collapse.singles
astral.labels <- comparePhylo(erable, astral@phylo)$NODES %>%
	data.frame %>%
	extract(erable, into = c("erable.label", "node"), regex = "(.*?) ?\\((.+)\\)", convert = T) %>%
	extract(astral.phylo, into = c("astral.label", "astral.node"), regex = "(.*?) ?\\((.+)\\)") %>%
	separate(astral.label, into = c("V1", "EN", "f1", "f2", "f3", "pp1", "pp2", "pp3", "q1", "q2", "q3", "QC", "astral.node"), sep = ";", convert = T)
mismatches <- filter(astral.labels, abs(pp1 - erable.label) > 0.01)
if (nrow(mismatches) > 1) {
	write("Labels mismatch between astral and erable", stderr())
	q(status = 1)
}
tree <- as_tibble(erable) %>%
	left_join(astral.labels, by = "node") %>%
	mutate(branch.length = case_when(
		parent == rootnode(erable) & label == "" ~ branch.length * 0.1,
		parent == rootnode(erable) & label != "" ~ branch.length * (1 - root.ratio),
		T ~ branch.length)
	) %>%
	left_join(species, by = "label") %>%
	as.treedata
busco.data <- filter(busco, label %in% tree@phylo$tip.label) %>%
	select(label, Status, pct)
busco.colors <- c(
	Complete   = "#00ba38",
	Fragmented = "#619cff",
	Missing    = "#f8766d"
)

tree.p <- ggtree(tree) +
	geom_tiplab(aes(label = paste(Species, Strain))) +
	geom_text(aes(x = branch, label = round(pp1, 2)), size = 2.5, nudge_y =  0.3) +
	geom_text(aes(x = branch, label = round(EN,  2)), size = 2.5, nudge_y = -0.3) +
	geom_treescale(width = 0.2)

p <- facet_plot(tree.p, data = busco.data, panel = "Busco", geom = geom_barh, mapping = aes(x = pct, fill = Status), stat = "identity") +
	scale_fill_manual(values = busco.colors)
p <- facet_widths(p, c(Busco = .2))
ggsave(output_file, p, width = width, height = height)
