
library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)
library(readxl)
library(ggtree)
library(treeio)
library(seqinr)
library(ggpubr)
library(ggrepel)
library(zoo)

with(snakemake@input, {
	fasta_file     <<- fasta
	geneconv_file  <<- geneconv
	gard_file      <<- gard
	metadata_file  <<- metadata
	RB_groups_file <<- RB_groups
	TM_file        <<- TMs
})
output_file   <- unlist(snakemake@output)
codons <- read.fasta(fasta_file, forceDNAtolower = F)
last_res <- length(codons[[1]])

metadata <- read_excel(metadata_file) %>%
	filter(!is.na(Gene)) %>%
	mutate(Gene_dash = gsub("-", "_", Gene))

TM_df <- read.table(TM_file, col.names = c("ref_seq", "name","start","end"))
ref_seq <- first(TM_df$ref_seq)

residues <- data.frame(res = as.character(codons[[ref_seq]]), aln_num = 1:last_res) %>%
	filter(res != "-") %>%
	mutate(res_num = 1:n())

TMs <- TM_df %>%
	left_join(residues, by = c(start = "res_num")) %>%
	left_join(residues, by = c(end = "res_num"))

TM_rect <- geom_rect(data = TMs, mapping = aes(xmin = aln_num.x, xmax = aln_num.y), ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.4)
TM_text <- geom_text(data = TMs, mapping = aes(x = (aln_num.x + aln_num.y) / 2, label = paste0(name, "\n")), y = -Inf, vjust = "inward", color = "gray")

scale_x <- scale_x_continuous(limits = c(1, last_res), expand = c(0, 0))
my_theme <- theme_bw() + theme(axis.title.x = element_blank(), legend.position = "none")
my_colors <- read.table(RB_groups_file, col.names = c("Group", "Color"), comm = "") %>%
	add_row(Group = "NA", Color = "darkgray") %>%
	with(scale_colour_manual(values = setNames(Color, Group)))

json <- fromJSON(gard_file)
bps <- with(json, data.frame(
	pos = names(siteBreakPointSupport) %>% as.numeric,
	support = unlist(siteBreakPointSupport)
))

roll_ident <- function(x, y, window = 15) {
	gaps <- c("x","X","n","N","-")
	data.frame(x = c(x), y = c(y)) %>%
		mutate(pos = 1:n()) %>%
		filter(! x %in% gaps, ! y %in% gaps) %>%
		mutate(ident = rollapply(x == y, window, mean, fill = NA)) %>%
		select(pos, ident) %>%
		filter(!is.na(ident))
}

format_label <- function(Gene, Domain, Completeness) {
	compl_suffix <- ifelse(Completeness > 3, "*", "")
	domain_suffix <- ifelse(Domain %in% c("R1","R2"), sprintf(" (%s)", Domain), "")
	paste0(Gene, compl_suffix, domain_suffix)
}

col.names <- c("Type", "Seq.Names", "Sim.Pvalue", "BC.KA.Pvalue", "Begin", "End", "Len", "Num.Poly", "Num.Dif", "Tot.Difs", "MisM.Pen")
frags <- read.table(geneconv_file, col.names = col.names) %>%
	filter(Type == "GI") %>%
	separate(Seq.Names, into = c("A","B"), sep = ";") %>%
	separate(A, into = c("A.Gene", "A.Domain"), sep = "__", remove = F) %>%
	separate(B, into = c("B.Gene", "B.Domain"), sep = "__", remove = F) %>%
	left_join(select(metadata, A.Gene = Gene, A.Group = Group, A.Completeness = Completeness), by = "A.Gene") %>%
	left_join(select(metadata, B.Gene = Gene, B.Group = Group, B.Completeness = Completeness), by = "B.Gene") %>%
	mutate(A.Label = format_label(A.Gene, A.Domain, A.Completeness), B.Label = format_label(B.Gene, B.Domain, B.Completeness)) %>%
	mutate(A.Label = ifelse(A.Gene == B.Gene, sprintf("%s (%s vs %s)", A.Gene, A.Domain, B.Domain), A.Label), B.Label = ifelse(A.Gene == B.Gene, "", B.Label))


p3 <- ggplot(bps) + TM_rect + TM_text +
	geom_segment(aes(x = pos, xend = pos, y = 0, yend = support)) +
	ylab("GARD model-averaged support") +
	scale_x + my_colors + my_theme

#protein <- lapply(codons, translate)

if (nrow(frags) > 0) {

	data <- list(
		`Tara-RRB` = with(codons, roll_ident(`Tara-RRB__R1`, `Tara-RRB__R2`)) %>% mutate(Group = "P"),
		`Ct-RRB` = with(codons, roll_ident(`M54302063-RB__R`, `Ct-RB__R`)) %>% mutate(Group = "C")
	) %>% bind_rows(.id = "Prot") %>%
		mutate(pos2 = lead(pos), ident2 = lead(ident)) %>%
		filter(pos2 - pos == 1)

	p1 <- ggplot(frags) + TM_rect + TM_text +
		geom_rect(aes(xmin = Begin, xmax = End, ymin = BC.KA.Pvalue, ymax = Sim.Pvalue, linetype = A.Label), color = "black", alpha = 0.1) +
		geom_text(aes(x = (Begin + End) / 2, y = (BC.KA.Pvalue + Sim.Pvalue) / 2, label = A.Label, color = A.Group), vjust = "outward") +
		geom_text(aes(x = (Begin + End) / 2, y = (BC.KA.Pvalue + Sim.Pvalue) / 2, label = paste0("\n", B.Label), color = B.Group)) +
		#scale_y_continuous(trans = "log10") +
		ylab("GENECONV global P-values") +
		scale_x + my_colors + my_theme
	p2 <- ggplot(data) + TM_rect + TM_text +
		geom_segment(aes(x = pos, y = ident, xend = pos2, yend = ident2, linetype = Prot, color = Group)) +
		ylab("Pairwise identity") +
		scale_x + my_theme + my_colors	
} else {
	p1 <- ggplot() + theme_void()
	p2 <- ggplot() + theme_void()
}

trees <- lapply(json$trees, `[[`, "newickString") %>%
	lapply(paste0, ";") %>%
	lapply(function(x) read.newick(text = x)) %>%
	lapply(as_tibble) %>%
	lapply(separate, col = "label", into = c("Gene_dash", "Domain"), sep = "__", fill = "left", remove = F) %>%
	lapply(left_join, y = metadata, by = "Gene_dash") %>%
	lapply(function(x) mutate(x, label = format_label(Gene, Domain, Completeness))) %>%
	lapply(`class<-`, c("tbl_tree", "data.frame")) %>%
	lapply(as.treedata) %>%
	lapply(function(tree) {
		ggtree(tree, layout = "ape") +
			geom_tiplab2(aes(color = Group), size = 3, hjust = -0.1) +
			geom_treescale(offset = 0.05, x = 1.8) +
			my_colors +
			theme(legend.position = "none")
	})

p123 <- annotate_figure(
	ggarrange(p1, p2, p3, ncol = 1, align = "v", labels = c("A", "B", "C"), hjust = 2),
	bottom = "Codon alignment coordinate"
)

p45 <- ggarrange(plotlist = trees, nrow = 1)

g <- ggarrange(p123, p45, ncol = 1, heights = c(3,2), labels = c("", "D"), hjust = 2) + theme(plot.margin = margin(0,0,0,20))
ggsave(output_file, g, width = 12, height = 12)
