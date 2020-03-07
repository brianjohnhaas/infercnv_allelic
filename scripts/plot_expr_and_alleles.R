#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()



parser$add_argument("--gene_coords", help="input gene coords file", required=TRUE, nargs=1)
parser$add_argument("--rnaseq_counts_matrix", help="rnaseq counts matrix", required=TRUE, nargs=1)
parser$add_argument("--allelic_fraction_matrix", help="allelic fraction matrix", required=TRUE, nargs=1)
parser$add_argument("--output_image_filename", help="output image file name", required=TRUE, nargs=1)
parser$add_argument("--cell_annots_file", help="cell type annotation file", required=TRUE, nargs=1)
parser$add_argument("--malignant_cell_types", help="malignant cell types, comma-delimited (rest considered normals)", required=TRUE, nargs=1)

parser$add_argument("--infercnv_obj", help="infercnv object for co-plotting", required=FALSE, default=NULL, nargs=1)


args = parser$parse_args()

gene_coords_file = args$gene_coords
rnaseq_counts_matrix_file = args$rnaseq_counts_matrix
allelic_fraction_matrix_file = args$allelic_fraction_matrix
output_image_filename = args$output_image_filename
cell_annots_file = args$cell_annots_file

infercnv_obj = args$infercnv_obj

malignant_cell_types = strsplit(args$malignant_cell_types, "\\s*,\\s*")[[1]]

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

cell_annots = read.table(cell_annots_file, header=F, row.names=NULL, sep="\t", stringsAsFactors=F)
colnames(cell_annots) = c('cell', 'celltype')

malignant_cells = cell_annots$cell[cell_annots$celltype %in% malignant_cell_types]
if (length(malignant_cells) == 0) {
    stop("Error, cannot extract malignant cells from ", cell_annots_file)
}



message("-building expression plot")

data = read.table(rnaseq_counts_matrix_file, header=T, sep="\t")
data = t( t(data)/colSums(data)) * 1e6 # cpm normalization
log2data = log2(data+1)
zscaled_data = t(scale(t(data), scale=T, center=T))
zscaled_data[is.nan(zscaled_data)] = 0

gexp.norm = zscaled_data

## subtract normals:

malignant_cells_idx = which(colnames(gexp.norm) %in% malignant_cells)

normal_reference_data = gexp.norm[,-malignant_cells_idx]

gexp.norm = gexp.norm[,malignant_cells_idx]


## subtract tumor from normal:
normal_gene_means = rowMeans(normal_reference_data)

for (i in 1:length(normal_gene_means) ) {
    gexp.norm[i,] = gexp.norm[i,]  - normal_gene_means[i]
}


gexp.melt = melt(gexp.norm)
colnames(gexp.melt) = c('genename', 'cell', 'exp.norm')

data = left_join(gexp.melt, gencode_gene_pos, key='genename')

data = data %>% filter(! is.na(chr))

## do smoothing by cell and chr.

data = data %>% mutate(cellchr = paste(cell, seqnames, sep=":"))
splitdata = split(data, data$cellchr)

smoother = function(df) {
    df = df %>% arrange(pos)
    df$exp.norm.smoothed = caTools::runmean(df$exp.norm,k=101, align="center")
    ## recenter cells

    return(df)
}

data = do.call(rbind, lapply(splitdata, smoother))

## mean center cells again after smoothing
splitdata = split(data, data$cell)
data = do.call(rbind, lapply(splitdata, function(x) {
    m = mean(x$exp.norm.smoothed)
    x$exp.norm.smoothed = x$exp.norm.smoothed - m
    return(x)
    } ) )


## cluster cells according to expr pattern:
gene_expr_matrix = data %>% select(genename, cell, exp.norm.smoothed) %>% spread(key=cell, value=exp.norm.smoothed)

rownames(gene_expr_matrix) = gene_expr_matrix$genename
gene_expr_matrix = gene_expr_matrix[,-1]
h = hclust(dist(t(gene_expr_matrix)))
ordered_cells = colnames(gene_expr_matrix)[h$order]

data$cell = ordered(data$cell, levels=ordered_cells)

## set up base plot
## define chr boundaries based on max coordinates for now.

q = quantile(data$exp.norm.smoothed, c(0.025, 0.975))

data$exp.norm.smoothed[data$exp.norm.smoothed < q[1] ] = q[1]
data$exp.norm.smoothed[data$exp.norm.smoothed > q[2] ] = q[2]

noise_range = 0.2
data$exp.norm.smoothed[data$exp.norm.smoothed > -noise_range & data$exp.norm.smoothed < noise_range] = 0


expr_plot = ggplot(data=data) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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

    geom_point(aes(x=pos, y=cell, color=exp.norm.smoothed), alpha=0.6, size=dotsize) +

    scale_colour_gradient2(low = "blue",
                           mid = "white",
                           high = "red",
                           midpoint = 0,
                           space = "Lab",
                           guide = "colourbar",
                           aesthetics = "colour")




########################
## InferCNV integration
########################

infercnv_plot = NULL

if (! is.null(infercnv_obj)) {
    library(infercnv)

    infercnv_obj = readRDS(infercnv_obj)
    expr.data = infercnv_obj@expr.data
    expr.data = expr.data[, unlist(infercnv_obj@observation_grouped_cell_indices)] # just the tumor cells

    h = hclust(dist(t(expr.data)))
    ordered_cells = colnames(expr.data)[h$order]
    expr.data = expr.data[,ordered_cells]

    expr.data = melt(expr.data)

    colnames(expr.data) = c('genename', 'cell', 'value')

    q = quantile(expr.data$value, c(0.025, 0.975))
    expr.data$value[expr.data$value < q[1] ] = q[1]
    expr.data$value[expr.data$value > q[2] ] = q[2]

    expr.data = left_join(expr.data, gencode_gene_pos, key='genename')

    infercnv_plot = ggplot(data=expr.data) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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

        geom_point(aes(x=pos, y=cell, color=value), alpha=0.6, size=dotsize) +

        scale_colour_gradient2(low = "blue",
                               mid = "white",
                               high = "red",
                               midpoint = 1,
                               space = "Lab",
                               guide = "colourbar",
                               aesthetics = "colour")

}






## ####################################
## SNPs
## ####################################

message("-building snp plot")


allele_matrix = read.table(allelic_fraction_matrix_file, header=T, row.names=NULL)

r_melt = melt(allele_matrix)
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


## filter snps, require at least 3 cells have coverage
min.cells = 3
snps_min_cells = data %>%  filter(AF > 0) %>% group_by(chrpos) %>% tally() %>% filter(n>=min.cells) %>% pull(chrpos)
data = data %>% filter(chrpos %in% snps_min_cells)

data = data %>% mutate(mBAF = pmax(AF, 1-AF))

## further restrict to het snps based on allele frequency data
het_snps = data %>% filter(AF > 0.25 & AF < 0.75) %>% group_by(chrpos) %>% tally() %>% filter(n>=3) %>% pull(chrpos)

data = data %>% filter(chrpos %in% het_snps)  ## note, want het_snps from normal!!!  but dont have it here.



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

pg = NULL
if (! is.null(infercnv_plot)) {
    pg = plot_grid(expr_plot, infercnv_plot, snps_plot, ncol=1, align='v')
} else {
    pg = plot_grid(expr_plot, snps_plot, ncol=1, align='v')
}

ggsave (output_image_filename, pg, width = 13.33, height = 7.5, units = 'in', dpi = 300)





