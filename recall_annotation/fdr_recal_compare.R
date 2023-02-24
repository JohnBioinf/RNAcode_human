library(ggplot2)
library(rjson)

library(Cairo)

num_gene <- 20090

parameters <- fromJSON(file="../parameters_local.json")

html_path <- parameters$HTML_PATH
data_dir_raw <- parameters$DATA_DIR_OLD
data_dir_concat <- parameters$DATA_DIR

fdr_both <- data.frame()
gene_recall_both <- data.frame()

for(type in c("raw", "concat")){
    if(type == "raw"){
        data_dir = data_dir_raw
    }else{
        data_dir = data_dir_concat
    }
    genome_alignment_dir <- paste0(data_dir, "/multiz100way") 
    RNAcode_overlap <- read.delim(paste0(genome_alignment_dir, "/RNAcode_overlap.tsv"), header=FALSE)
    RNAcode_no_overlap <- read.delim(paste0(genome_alignment_dir, "/RNAcode_no_overlap.tsv"), header=FALSE)
    gene_recall <- read.delim(paste0(genome_alignment_dir, "/gene_recall.tsv"), header=FALSE)
    colnames(RNAcode_no_overlap) <- c("p_value", "sum_FP")
    colnames(RNAcode_overlap) <- c("p_value", "sum_TP")
    fdr <- cbind(RNAcode_overlap, RNAcode_no_overlap["sum_FP"])
    fdr$fdr <- (fdr$sum_FP / (fdr$sum_FP + fdr$sum_TP)) * 100
    fdr$log_p <- -log(fdr$p_value, base = 10)
    fdr$type <- type
    colnames(gene_recall) <- c("p_value", "sum")
    gene_recall$percent <- (gene_recall$sum / num_gene) * 100
    gene_recall$log_p <- -log(gene_recall$p_value, base = 10)
    gene_recall$type <- type
    if(gene_recall[1,4] == Inf){
        gene_recall[1,4] <- gene_recall[2,4] + 1
    }
    fdr_both <- rbind(fdr_both, fdr)
    gene_recall_both <- rbind(gene_recall_both, gene_recall)
}

g_fdr <- ggplot(fdr_both, aes(y=fdr, x=log_p, color=type)) +
  geom_line() +
  scale_y_continuous(name="fdr", limits=c(0,100)) +
  scale_x_continuous(name="-log(p-val)", limits=c(0,17))

g_recall <- ggplot(gene_recall_both, aes(y=percent, x=log_p, color=type)) +
  geom_line() +
  scale_y_continuous(name="recall", limits=c(0,100)) +
  scale_x_continuous(name="-log(p-val)", limits=c(0,17))

svg(paste0(html_path, "/fdr_WGA_RNAcode_human.svg"))
g_fdr
dev.off()

svg(paste0(html_path, "/recall_WGA_RNAcode_human.svg"))
g_recall
dev.off()
