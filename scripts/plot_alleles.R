#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()



parser$add_argument("--gene_coords", help="input gene coords file", required=TRUE, nargs=1)
parser$add_argument("--allelic_fraction_matrix", help="allelic fraction matrix", required=TRUE, nargs=1)
parser$add_argument("--output_image_filename", help="output image file name", required=TRUE, nargs=1)
parser$add_argument("--cell_annots_file", help="cell type annotation file", required=TRUE, nargs=1)

args = parser$parse_args()

gene_coords_file = args$gene_coords
allelic_fraction_matrix_file = args$allelic_fraction_matrix
output_image_filename = args$output_image_filename
cell_annots_file = args$cell_annots_file

library(tidyverse)
library(cowplot)
library(reshape2)

## ###############################
## Expression plot
## ###############################

dotsize=0.3

## parsing annotation coordinates
gencode_gene_pos = read.table(gene_coords_file, header=F, row.names=NULL, stringsAsFactors=F)
colnames(gencode_gene_pos) = c('genename', 'seqnames', 'start', 'end')

gencode_gene_pos$chr = str_replace(string=gencode_gene_pos$seqnames, pattern="chr", replacement="")
gencode_gene_pos = gencode_gene_pos %>% filter(chr %in% 1:22)
gencode_gene_pos = gencode_gene_pos %>% mutate(chr = ordered(chr, levels=1:22))
gencode_gene_pos$pos = as.numeric( (gencode_gene_pos$start + gencode_gene_pos$end)/2)


## get chr bounds for plotting later.
chr_maxpos = gencode_gene_pos %>% group_by(chr) %>% summarize(maxpos = max(pos))
chr_maxpos$minpos = 1


## parse cell annotations
message("-getting cell annotations")
cell_annots = read.table(cell_annots_file, header=F, row.names=NULL, sep="\t", stringsAsFactors=F)
colnames(cell_annots) = c('cell', 'celltype')

## cluster cells according to expr pattern:
#gene_expr_matrix = data %>% select(genename, cell, exp.norm.smoothed) %>% spread(key=cell, value=exp.norm.smoothed)

#rownames(gene_expr_matrix) = gene_expr_matrix$genename
#gene_expr_matrix = gene_expr_matrix[,-1]
#h = hclust(dist(t(gene_expr_matrix)))
#ordered_cells = colnames(gene_expr_matrix)[h$order]

#data$cell = ordered(data$cell, levels=ordered_cells)

## set up base plot


## ####################################
## SNPs
## ####################################

message("-building snp plot")

message("parsing matrix: ", allelic_fraction_matrix_file);
allele_matrix = read.table(allelic_fraction_matrix_file, header=T, row.names=1, check.names=F)

message("melting matrix.")
r_melt = melt(as.matrix(allele_matrix))
colnames(r_melt) = c('chrpos', 'cell', 'AF')

## include chr and postion separately.
data = r_melt %>% separate(chrpos, "::", into=c('seqnames', 'pos'), remove=FALSE)

data$chr = str_replace(string=data$seqnames, pattern="chr", replacement="")
data = data %>% filter(chr %in% 1:22)

data = data %>% mutate(chr = ordered(chr, levels=1:22))
data$pos = as.numeric(data$pos)


## get chr bounds for plotting later.
chr_maxpos = data %>% group_by(chr) %>% summarize(maxpos = max(pos))
chr_maxpos$minpos = 1

## select snps where the fraction of cells containing an allelic site is between 0.1 and 0.9 of cells with read coverage.
#het_snps = data %>% group_by(chrpos) %>%  mutate(l=sum(alleleCov>0), n.bulk=sum(totCov>0), E=l/n.bulk) %>% filter(E>0.1 & E<0.9) %>% pull('chrpos')
#data = data %>% filter(chrpos %in% het_snps)


#cell_groupings = split(cell_annots, cell_annots$celltype)

#for (cell_grouping in cell_groupings) {

#cells_want = cell_groupings[cell_grouping][[1]]$cell

#    data = data[,colnames(data) %in% cells_want]
    
    ## filter snps, require at least 3 cells have coverage
    min.cells = 3
    snps_min_cells = data %>%  filter(AF > 0) %>% group_by(chrpos) %>% tally() %>% filter(n>=min.cells) %>% pull(chrpos)
    data = data %>% filter(chrpos %in% snps_min_cells)
    
    data = data %>% mutate(mBAF = pmax(AF, 1-AF))
    
    ## further restrict to het snps based on allele frequency data
    het_snps = data %>% filter(AF > 0.25 & AF < 0.75) %>% group_by(chrpos) %>% tally() %>% filter(n>=3) %>% pull(chrpos)
    
    data = data %>% filter(chrpos %in% het_snps)
    
    
    data = data %>% filter(AF > 0)  ## just here
    
    
    ## set up base plot
    ## define chr boundaries based on max coordinates for now.
    
    
    snps_plot = ggplot(data=data) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
        theme_bw() +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
              ) +

    geom_vline(data=chr_maxpos, aes(xintercept=minpos), color=NA) +
    geom_vline(data=chr_maxpos, aes(xintercept=maxpos), color=NA) +
    
    geom_point(aes(x=pos, y=cell, color=mBAF), alpha=0.6, size=dotsize) + scale_radius()
    
    
    message("-writing image")
    
    png("test.alleles.png")
    plot(snps_plot)
    dev.off()

    stop("debug stop")
    
#}



#pg = NULL
#if (! is.null(infercnv_plot)) {
#    pg = plot_grid(expr_plot, infercnv_plot, snps_plot, ncol=1, align='v')
#} else {
#    pg = plot_grid(expr_plot, snps_plot, ncol=1, align='v')
#}

#ggsave (output_image_filename, pg, width = 13.33, height = 7.5, units = 'in', dpi = 300)





