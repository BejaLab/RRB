library(phytools)
library(treeio)
library(dplyr)
library(tidyr)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(tibble)

input <- snakemake@input
output <- snakemake@output

clades_file <- input$clades

best_core_file                   <- input$best_core
best_core_raxml_nwk_file         <- input$best_core_raxml_nwk
best_core_raxml_nhx_file         <- input$best_core_raxml_nhx
best_subset_file                 <- input$best_subset
best_subset_raxml_nwk_file       <- input$best_subset_raxml_nwk
best_subset_raxml_nhx_file       <- input$best_subset_raxml_nhx
best_subset_nt2aa_file           <- input$best_subset_nt2aa
best_w_outgroups_nt2_rooted_file <- input$best_w_outgroups_nt2_rooted
rhod_w_outgroups_rooted_file     <- input$rhod_w_outgroups_rooted
rhod_w_outgroups_nt2_rooted_file <- input$rhod_w_outgroups_nt2_rooted

output_file <- unlist(output)

clades <- read.table(clades_file, sep = "\t", header = T, comment.char = "")
colors <- distinct(clades, clade, color) %>%
	with(setNames(color, clade))

add.clade <- function(tree) {
	as_tibble(tree) %>%
		separate(label, into = c("gene", "domain"),  fill = "right", sep = "__", remove = F) %>%
		separate(label, into = c("gene2", "family"), fill = "right", extra = "drop", sep = "_", remove = F) %>%
		left_join(clades, by = "gene") %>%
		mutate(label_show = case_when(
			family == "BRL" ~ gene2,
			family %in% c("9TM", "PK") ~ sprintf("%s {%s}", gene2, family),
			domain %in% c(NA, "R", "B", "RX") ~ gene,
			T ~ sprintf("%s (%s)", gene, domain)
		)) %>%
		`class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame")) %>%
		as.treedata
}

scale_ufboot <- function() scale_radius(range = c(1,3), limits = c(70, 100))
rdigger.tree <- function(iqtree, nwk, nhx) {
	rooted.tree <- read.newick(nwk)
	outgroup <- rooted.tree %>%
		{Descendants(rooted.tree, rootnode(.), type = "children")} %>%
		{Descendants(rooted.tree, first(.),    type = "tips")} %>%
		{`[`(rooted.tree$tip.label, .[[1]])}
	lwr <- read.nhx(nhx) %>%
		as_tibble %>%
		pull(LWR) %>%
		max(na.rm = T)

	iqtree.rooted <- read.newick(iqtree) %>%
		root(outgroup, edgelabel = T, resolve.root = T) %>%
		as_tibble
	root.node <- rootnode(iqtree.rooted)$node
	root.lens <- filter(iqtree.rooted, parent == root.node, node != root.node) %>%
		pull(branch.length) %>%
		mean
	tree <- mutate(iqtree.rooted, branch.length = case_when(node == root.node ~ 0, parent == root.node ~ root.lens, T ~ branch.length)) %>%
		mutate(LWR = ifelse(node == parent, lwr, NA)) %>%
		mutate(UFboot = ifelse(node %in% parent & node != parent, label, NA), UFboot = as.numeric(UFboot)) %>%
		add.clade

	ggtree(tree, layout = "rectangular") +
		geom_tiplab(aes(subset = isTip, label = label_show, color = clade), offset = 0.01) +
		geom_text2(aes(subset = !is.na(LWR), label = LWR)) +
		geom_point2(aes(subset = !isTip & UFboot >= 95,              x = branch, size = UFboot), color = "red") +
		geom_point2(aes(subset = !isTip & UFboot >= 70 & UFboot < 95, x = branch, size = UFboot), color = "black") +
		scale_ufboot() +
		scale_color_manual(values = colors, na.value = "black") +
		geom_treescale() +
		theme(legend.position = "None")
}
epa.tree <- function(fname) {
	#epa <- read_jplace(fname)
	#epa.edges <- t(epa$edge_key) %>%
	#	data.frame %>%
	#	setNames(c("parent", "node")) %>%
	#	rownames_to_column("edge_num") %>%
	#	mutate(edge_num = as.numeric(edge_num)) %>%
	#	left_join(epa$placement_positions, by = "edge_num") %>%
	#	filter(!is.na(likelihood))
	#epa.tree <- as_tibble(epa$arbre) %>%
	#	left_join(epa.edges, by = c("node", "parent")) %>%
	#	add.clade
	epa <- read.jplace(fname) %>% add.clade
	ggtree(epa, aes(color = nplace), layout = "ape") +
		scale_colour_gradient2(high = "red", mid = "black") +
		new_scale_color() +
		geom_tiplab(aes(subset = isTip, label = label_show, color = clade), offset = 0.01) +
		scale_color_manual(values = colors, na.value = "black") +
		geom_treescale() +
		theme(legend.position = "None")
}
iqtree.tree.ape <- function(fname) {
	iqtree <- read.iqtree(fname) %>%
		add.clade
	ggtree(iqtree, layout = "ape") +
		scale_colour_gradient2(high = "red", mid = "black") +
		new_scale_color() +
		geom_point2(aes(subset = !isTip & UFboot >= 95,              x = branch.x, y = branch.y, size = UFboot), color = "red") +
		geom_point2(aes(subset = !isTip & UFboot >= 70 & UFboot < 95, x = branch.x, y = branch.y, size = UFboot), color = "black") +
		scale_ufboot() +
		new_scale_color() +
		geom_tiplab2(aes(label = label_show, color = clade), offset = 0.01) +
		scale_color_manual(values = colors, na.value = "black") +
		geom_treescale()
}

iqtree.tree.rect <- function(fname) {
	iqtree <- read.iqtree(fname) %>%
		add.clade
	ggtree(iqtree, layout = "rectangular") +
		scale_colour_gradient2(high = "red", mid = "black") +
		new_scale_color() +
		geom_point2(aes(subset = !isTip & UFboot >= 95,               x = branch, size = UFboot), color = "red") +
		geom_point2(aes(subset = !isTip & UFboot >= 70 & UFboot < 95, x = branch, size = UFboot), color = "black") +
		scale_ufboot() +
		new_scale_color() +
		geom_tiplab(aes(label = label_show, color = clade), offset = 0.01) +
		scale_color_manual(values = colors, na.value = "black") +
		geom_treescale()
}

blank <- ggplot() + theme_void()

best.codon.core.rd <- rdigger.tree(best_core_file,           best_core_raxml_nwk_file, best_core_raxml_nhx_file)
best.codon.subset.rd <- rdigger.tree(best_subset_file,       best_subset_raxml_nwk_file, best_subset_raxml_nhx_file)
best.nt2aa.subset.rd <- rdigger.tree(best_subset_nt2aa_file, best_subset_raxml_nwk_file, best_subset_raxml_nhx_file)
best.nt2aa.og      <- iqtree.tree.ape(best_w_outgroups_nt2_rooted_file)

best <- ggarrange(best.codon.core.rd, best.codon.subset.rd, best.nt2aa.subset.rd, best.nt2aa.og, nrow = 1, common.legend = T, legend = "right")

rhod.codon <- iqtree.tree.rect(rhod_w_outgroups_rooted_file)
rhod.nt2aa <- iqtree.tree.rect(rhod_w_outgroups_nt2_rooted_file)

rhod <- ggarrange(rhod.codon, rhod.nt2aa, blank, widths = c(1,1,0.5), nrow = 1, common.legend = T, legend = "right")

trees <- ggarrange(rhod, best, nrow = 2, heights = c(1.8, 1), common.legend = T, legend = "right")

ggsave(output_file, trees, height = 10, width = 16)
