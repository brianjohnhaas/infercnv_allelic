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

gcinfo(TRUE)

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


malignant_cells_idx = grep("malignant", cell_annots$celltype, value=F)
normal_cells = cell_annots$cell[-malignant_cells_idx]
malignant_cells = cell_annots$cell[malignant_cells_idx]



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


## filter snps, require at least x number of cells have evidence of each allelic variant
min.cells = 3
cells_w_ref_allele = rowSums(allele_matrix != 0 & allele_matrix < 0.5)
cells_w_alt_allele = rowSums(allele_matrix != 0 & allele_matrix > 0.5)

allele_matrix = allele_matrix[(cells_w_ref_allele >= min.cells & cells_w_alt_allele >= min.cells), ]

num_snps_all = nrow(allele_matrix)
message("-num het snps used: ", num_snps_all)

write.table(allele_matrix, file=paste0(allelic_fraction_matrix_file, num_snps_all, ".af.matrix"), quote=F, sep="\t")


## reset the alt allele fraction to the cell-population minor allele frequency

message("-setting alt allele fraction to the cell-population minor allele fraction")
mAF_allele_matrix = apply(allele_matrix, 1, function(x) {
    nonzero_val_idx = which(x>0)
    nonzero_vals = x[nonzero_val_idx]
    frac_high = sum(nonzero_vals>0.5)/length(nonzero_vals)
    if ( frac_high > 0.5) {
        x[x==1] = 0.999
        x[nonzero_val_idx ] = 1 - x[nonzero_val_idx]
    }
    x
})

allele_matrix = t(mAF_allele_matrix)



message("melting matrix.")
datamelt = melt(as.matrix(allele_matrix))
colnames(datamelt) = c('chrpos', 'cell', 'AF')


## include chr and postion separately.
datamelt = datamelt %>% separate(chrpos, "::", into=c('seqnames', 'pos'), remove=FALSE)
gc()

datamelt$chr = str_replace(string=datamelt$seqnames, pattern="chr", replacement="")
gc()

datamelt = datamelt %>% filter(chr %in% 1:22)

datamelt = datamelt %>% mutate(chr = ordered(chr, levels=1:22))
datamelt$pos = as.numeric(datamelt$pos)


## get chr bounds for plotting later.
chr_maxpos = datamelt %>% group_by(chr) %>% summarize(maxpos = max(pos))
chr_maxpos$minpos = 1



#cell_groupings = split(cell_annots, cell_annots$celltype)

#for (cell_grouping in cell_groupings) {

#cells_want = cell_groupings[cell_grouping][[1]]$cell

#    datamelt = datamelt[,colnames(datamelt) %in% cells_want]


datamelt = datamelt %>% filter(AF > 0)  ## if AF==0, then means we have no coverage.


#    datamelt = datamelt %>% mutate(mBAF = pmax(AF, 1-AF))

## further restrict to het snps based on allele frequency datamelt
#het_snps = datamelt %>% filter(AF > 0.25 & AF < 0.75) %>% group_by(chrpos) %>% tally() %>% filter(n>=3) %>% pull(chrpos)

#datamelt = datamelt %>% filter(chrpos %in% het_snps)


normal_datamelt = datamelt %>% filter(cell %in% normal_cells)


normal_snps_plot = ggplot(data=normal_datamelt) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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

    geom_point(aes(x=pos, y=cell, color=AF), alpha=0.6, size=dotsize) + scale_radius() +
    scale_color_gradient(low="blue", high="red")


malignant_datamelt = datamelt %>% filter(cell %in% malignant_cells)

malignant_snps_plot = ggplot(data=malignant_datamelt) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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

    geom_point(aes(x=pos, y=cell, color=AF), alpha=0.6, size=dotsize) + scale_radius() +
    scale_color_gradient(low="blue", high="red")


message("-writing image")

num_normal_cells = length(normal_cells)
num_malignant_cells = length(malignant_cells)

pg = plot_grid(normal_snps_plot, malignant_snps_plot, ncol=1, align='v', rel_heights=c(num_normal_cells, num_malignant_cells))


ggsave (output_image_filename, pg, width = 13.33, height = 7.5, units = 'in', dpi = 300)





