#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))


parser = ArgumentParser()
parser$add_argument("--alt_ratios", help="alt ratios sparse matrix", required=TRUE, nargs=1)
parser$add_argument("--cell_names", help="cell names file", required=TRUE, nargs=1)
parser$add_argument("--chrpos_names", help="chr pos names file", required=TRUE, nargs=1);
parser$add_argument("--output", help="output matrix file name", required=TRUE, nargs=1);

args = parser$parse_args()


alt_ratios_file = args$alt_ratios
cell_names_file = args$cell_names
chrpos_names_file = args$chrpos_names
output_filename = args$output

library(Matrix)
data = readMM(file=alt_ratios_file)

cell_annots = read.table(cell_names_file, header=F, stringsAsFactors=F, check.names=F)[,1]
cell_annots = gsub("-", "_", cell_annots)

chrpos = read.table(chrpos_names_file, header=F, stringsAsFactors=F, check.names=F)[,1]


data = as.matrix(data)

colnames(data) = cell_annots
rownames(data) = chrpos


write.table(data, file=output_filename, quote=F, sep="\t")

message("done")

