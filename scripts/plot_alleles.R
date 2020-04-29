#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()

parser$add_argument("--gene_coords", help="input gene coords file", required=TRUE, nargs=1)
parser$add_argument("--allelic_fraction_matrix", help="allelic fraction matrix", required=TRUE, nargs=1)
parser$add_argument("--output_prefix", help="output file prefix", required=TRUE, nargs=1)
parser$add_argument("--cell_annots_file", help="cell type annotation file", required=TRUE, nargs=1)
parser$add_argument("--cluster_cells", help="cluster cells", required=FALSE, default=FALSE, action='store_true')
parser$add_argument("--smooth", help="smooth by chr", required=FALSE, default=FALSE, action='store_true')
parser$add_argument("--smK", help="smooth by k window size", required=FALSE, type='integer', default=31, nargs=1)
parser$add_argument("--trend_smK", help="smooth k for trend lines", required=FALSE, type='integer', default=31, nargs=1)
parser$add_argument("--center", help="center by cell", required=FALSE, default=FALSE, action='store_true')
parser$add_argument("--sample_size", help="sample counts of normal and malignant cells", required=FALSE, default=-1)
parser$add_argument("--colorscheme", help="color scheme. options: BlueRed, BlueYellowRed  (default: BlueRed)", required=FALSE, default="BlueRed")
parser$add_argument("--incl_each_chr_sep", help="include each chromosome plot separately", required=FALSE, default=FALSE, action='store_true')


args = parser$parse_args()

gene_coords_file = args$gene_coords
allelic_fraction_matrix_file = args$allelic_fraction_matrix
output_prefix = args$output_prefix
cell_annots_file = args$cell_annots_file
cluster_flag = args$cluster_cells
smooth_flag = args$smooth
smK = args$smK
trend_smK = args$trend_smK
center_flag = args$center
sample_size = args$sample_size
colorscheme = args$colorscheme
each_chr_sep_flag = args$incl_each_chr_sep

CELL_POINT_ALPHA = 0.6


library(tidyverse)
library(cowplot)
library(reshape2)

## ###############################
## Expression plot
## ###############################

dotsize=0.3

#gcinfo(TRUE)

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

if (sample_size > 0) {
    message("-sampling ", sample_size, " cells from normal and malignant sets")
    if (length(normal_cells) < sample_size) {
        normal_cells = sample(x=normal_cells, size=sample_size, replace=FALSE)
    }
    if (length(malignant_cells) < sample_size) {
        malignant_cells = sample(x=malignant_cells, size=sample_size, replace=FALSE)
    }
}


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

if (sample_size > 0) {
    cells_want = union(normal_cells, malignant_cells)
    allele_matrix = allele_matrix[, colnames(allele_matrix) %in% cells_want]
}

num_snps_all = nrow(allele_matrix)
message("-num het snps used: ", num_snps_all)

if (sample_size > 0) {
    write.table(allele_matrix, file=paste0(allelic_fraction_matrix_file, num_snps_all, ".af.matrix"), quote=F, sep="\t")
}

## reset the alt allele fraction to the cell-population minor allele

message("-setting alt allele fraction to the tumor cell-population minor allele")
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


if (cluster_flag) {
    ## define cell ordering.
    message("clustering cells")
    h = hclust(dist(t(allele_matrix)))
    ordered_cells = colnames(allele_matrix)[h$order]
}

message("melting matrix.")
datamelt = melt(as.matrix(allele_matrix))
colnames(datamelt) = c('chrpos', 'cell', 'AF')

if (cluster_flag) {
    datamelt$cell = ordered(datamelt$cell, levels=ordered_cells)
}

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


datamelt = datamelt %>% filter(AF > 0)  ## if AF==0, then means we have no coverage.

datamelt = datamelt %>% mutate(cellchr = paste(cell, seqnames, sep=":"))


if (smooth_flag) {
    splitdata = split(datamelt, datamelt$cellchr)

    ## smooth AFs across chromosomes.
    message("-smoothing")
    smoother = function(df) {
        df = df %>% arrange(pos)
        df$AF = caTools::runmean(df$AF, k=smK, align="center")
        return(df)
    }

    datamelt = do.call(rbind, lapply(splitdata, smoother))
}

if (center_flag) {

    ## center by cell:
    splitdata = split(datamelt, datamelt$cell)

    ## smooth AFs across chromosomes.
    message("-centering by cell")
    centerbycell = function(df) {

        df$AF = df$AF - mean(df$AF)

        return(df)
    }

    datamelt = do.call(rbind, lapply(splitdata, centerbycell))
}


midpt = mean(datamelt$AF)



make_plots = function(dataToPlot, output_image_filename) {

    ## normal cell plot

    normal_snps_plot = NULL
    num_normal_cells = length(normal_cells)

    dataToPlot$sample_type = "tumor"



    if (num_normal_cells > 0) {

        dataToPlot$sample_type[ dataToPlot$cell %in% normal_cells ] = "normal"

        message("-making normal plot")

        normal_dataToPlot = dataToPlot %>% filter(cell %in% normal_cells)


        normal_snps_plot = ggplot(data=normal_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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

            geom_point(aes(x=pos, y=cell, color=AF), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()


        if (colorscheme == "BlueRed") {
            normal_snps_plot = normal_snps_plot + scale_color_gradient(low="blue", high="red")
        } else {
            normal_snps_plot = normal_snps_plot + scale_color_gradient2(low="blue", mid='yellow', high="red", midpoint=midpt)
        }
    }

    ## malignant cell plot

    message("-making malignant plot")

    num_malignant_cells = length(malignant_cells)

    malignant_dataToPlot = dataToPlot %>% filter(cell %in% malignant_cells)

    malignant_snps_plot = ggplot(data=malignant_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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

        geom_point(aes(x=pos, y=cell, color=AF), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()

    if (colorscheme == "BlueRed") {
        malignant_snps_plot = malignant_snps_plot + scale_color_gradient(low="blue", high="red")
    } else {
        malignant_snps_plot = malignant_snps_plot + scale_color_gradient2(low="blue", mid='yellow', high="red", midpoint=midpt)
    }


    message("-making allele freq plot")


    allele_freq_means = dataToPlot %>%
        group_by(chrpos,sample_type) %>%
        mutate(grp_pos_mean_AF = mean(AF)) %>% select(chrpos, chr, pos, sample_type, grp_pos_mean_AF) %>% unique()



    ## smooth mean AFs across chromosomes for each sample type.
    message("-smoothing sample means for trend lines")

    allele_freq_means = allele_freq_means %>% mutate(sampleTypeChr = paste(sample_type, chr, sep=":"))
    splitdata = split(allele_freq_means, allele_freq_means$sampleTypeChr)

    smoother = function(df) {
        df = df %>% arrange(pos)
        df$grp_pos_mean_AF_sm = caTools::runmean(df$grp_pos_mean_AF, k=trend_smK, align="center")
        return(df)
    }

    allele_freq_means = do.call(rbind, lapply(splitdata, smoother))

    allele_freq_plot = allele_freq_means %>%
        ggplot(aes(x=pos, y=grp_pos_mean_AF)) +
        facet_grid (~chr, scales = 'free_x', space = 'fixed') +
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
        geom_point(aes(color=sample_type), alpha=0.2, size=0.2)


    allele_freq_plot_w_trendlines = allele_freq_plot +
        geom_line(data=allele_freq_means,
                  aes(x=pos, y=grp_pos_mean_AF_sm, color=sample_type), size=0.5, alpha=1)

    message("-writing heatmap image")

    if (num_normal_cells > 0) {

        ratio_normal_cells = max(0.25, num_normal_cells/(num_normal_cells + num_malignant_cells))

        pg = plot_grid(normal_snps_plot, allele_freq_plot_w_trendlines, malignant_snps_plot, ncol=1, align='v', rel_heights=c(ratio_normal_cells, 0.25, 1-ratio_normal_cells))

        ggsave (output_image_filename, pg, width = 13.33, height = 7.5, units = 'in', dpi = 300)

    } else {
        ggsave (output_image_filename, malignant_snps_plot, width = 13.33, height = 7.5, units = 'in', dpi = 300)
    }

}


#output_image_filename = paste0(output_prefix, ".genome.png")
#make_plots(datamelt, output_image_filename)

if (each_chr_sep_flag) {
    chrs = datamelt %>% select(chr) %>% unique() %>% pull('chr')
    for (chr_select in chrs) {
        chrdata = datamelt %>% filter(chr==chr_select)
        chr_image_filename = paste0(output_prefix, ".chr", chr_select, ".png")
        message("-making image for chr: ", chr_select)
        make_plots(chrdata, chr_image_filename)
    }

}

message("done")

