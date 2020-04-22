#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()

parser$add_argument("--allelic_fraction_matrix", help="allelic fraction matrix", required=TRUE, nargs=1)
parser$add_argument("--cell_annots_file", help="cell type annotation file. format tsv with cols: 'cellname', 'celltype', 'patient'", required=TRUE, nargs=1)

args = parser$parse_args()


allelic_fraction_matrix_file = args$allelic_fraction_matrix
cell_annots_file = args$cell_annots_file

message("-reading cell annots file: ", cell_annots_file)
cell_annots_df = read.table(cell_annots_file, header=F, sep="\t", stringsAsFactors=F, check.names=F, row.names=NULL)
colnames(cell_annots_df) = c('cellname', 'celltype', 'patient')

message("-reading allelic fraction matrix: ", allelic_fraction_matrix_file)
data_AF_matrix = read.table(allelic_fraction_matrix_file, header=T, row.names=1, sep="\t", check.names=F)

message("\n* beginning partitioning by patient.")
patient_groupings = split(cell_annots_df, cell_annots_df$patient)
for (patient in names(patient_groupings)) {

    message("-processing patient: ", patient)
    
    patient_grouping_df = patient_groupings[[patient]]
    cells = patient_grouping_df$cellname

    patient_cell_annots_file = paste0(patient, ".cell_annots.tsv")
    message("-writing cell annots file: " , patient_cell_annots_file)
    write.table(patient_grouping_df[,1:2], file=patient_cell_annots_file, col.names=F, row.names=F, sep="\t", quote=F)

    patient_allele_matrix_file = paste0(patient, ".AF.matrix")
    message("-writing allele matrix file: ", patient_allele_matrix_file)            
    patient_AF_matrix = data_AF_matrix[, colnames(data_AF_matrix) %in% cells]
    write.table(patient_AF_matrix, file=patient_allele_matrix_file, quote=F, sep="\t", row.names=T)
    
}


message("done")

