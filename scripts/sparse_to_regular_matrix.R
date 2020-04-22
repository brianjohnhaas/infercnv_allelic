#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))


parser = ArgumentParser()
parser$add_argument("--malt_ratios", help="malt ratios sparse matrix", required=TRUE, nargs=1)
parser$add_argument("--cell_names", help="cell names file", required=TRUE, nargs=1)
parser$add_argument("--chrpos_names", help="chr pos names file", required=TRUE, nargs=1);
parser$add_argument("--output", help="output matrix file name", required=TRUE, nargs=1);

args = parser$parse_args()


malt_ratios_file = args$malt_ratios
cell_names_file = args$cell_names
chrpos_names_file = args$chrpos_names
output_filename = args$output

library(Matrix)
data = readMM(file=malt_ratios_file)

cell_annots = read.table(cell_names_file)[,1]

chrpos = read.table(chrpos_names_file)[,1]

colnames(data) = cell_annots
rownames(data) = chrpos



data = as.matrix(data)

write.table(data, file=output_filename, quote=F, sep="\t")

message("done")

