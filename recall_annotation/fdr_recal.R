library(ggplot2)
library(plotly)
library(htmlwidgets)
library(rjson)

args = commandArgs(trailingOnly=TRUE)

num_gene <- 60664
num_exons <- 765571

parameters <- fromJSON(file="../parameters_local.json")

html_path <- parameters$HTML_PATH
genome_alignment_dir <- args[1]

genome_alignment <- tail(strsplit(genome_alignment_dir, "/", fixed=TRUE)[[1]], n=1)

RNAcode_overlap <- read.delim(paste0(genome_alignment_dir, "RNAcode_overlap_BT.tsv"), header=FALSE)
RNAcode_no_overlap <- read.delim(paste0(genome_alignment_dir, "RNAcode_no_overlap_BT.tsv"), header=FALSE)
gene_recall <- read.delim(paste0(genome_alignment_dir, "gene_recall.tsv"), header=FALSE)

colnames(RNAcode_no_overlap) <- c("p_value", "sum")
RNAcode_no_overlap$type <- "not annotated"

colnames(RNAcode_overlap) <- c("p_value", "sum")
RNAcode_overlap$type <- "annotated"
exon_dat <- rbind(RNAcode_overlap, RNAcode_no_overlap)
exon_dat$percent <- (exon_dat$sum / num_exons) * 100
exon_dat$log_p <- -log(exon_dat$p_value, base = 10)

colnames(gene_recall) <- c("p_value", "sum")
gene_recall$percent <- (gene_recall$sum / num_gene) * 100
gene_recall$log_p <- -log(gene_recall$p_value, base = 10)

g_exon <- ggplot(exon_dat, aes(y=percent, x=log_p, color=type)) +
  geom_line()
# g_exon
g_gene <- ggplot(gene_recall, aes(y=percent, x=log_p)) +
  geom_line()
# g_gene

saveWidget(as_widget(ggplotly(g_exon)), paste0(html_path, genome_alignment, "_exon.html"), selfcontained = TRUE)
saveWidget(as_widget(ggplotly(g_gene)), paste0(html_path, genome_alignment, "_gene.html"), selfcontained = TRUE)
