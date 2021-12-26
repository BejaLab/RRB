
library(dplyr)
library(tidyr)

input  <- snakemake@input
output <- snakemake@output

log.files  <- input$log_files
tree.files <- input$tree_files

best.log.file  <- output$log_file
best.tree.file <- output$tree_file

best.index <- lapply(log.files, readLines) %>%
	lapply(as.data.frame) %>%
	lapply(setNames, "line") %>%
	setNames(log.files) %>%
	bind_rows(.id = "fname") %>%
	filter(grepl("BEST SCORE FOUND", line)) %>%
	extract(line, into = "score", regex = ": (.+)", convert = T) %>%
	arrange(-score) %>%
	head(n = 1) %>%
	pull(fname) %>%
	match(log.files)

ok.log  <- file.copy(log.files[best.index],  best.log.file)
ok.tree <- file.copy(tree.files[best.index], best.tree.file)
