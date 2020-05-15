#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

###############################################################################
#
#  --allelic_counts <str>     file containing list of *.covered_sites.tsv
#                             format:
#                                cell_name (tab) /path/to/cell_covered_sites.tsv
#
#  --cell_types <str>         file containing cell types and patient identity
#                             format:
#                                cell_name (tab) cell_type (tab) patient
#                             note: <no header line>
#
###############################################################################


__EOUSAGE__

    ;


my $help_flag;
my $allelic_counts_files;
my $cell_types_file;

&GetOptions ( 'h' => \$help_flag,
              'allelic_counts=s' => \$allelic_counts_files,
              'cell_types=s' => \$cell_types_file,
    );

if ($help_flag) {
    die $usage;
}

unless ($allelic_counts_files && $cell_types_file) {
    die $usage;
}



main: {
    
    my %cell_to_allelic_file = &parse_cell_to_allelic_file($allelic_counts_files);
    
    my %patient_to_cell_types = &parse_cell_type_info($cell_types_file);
    
    foreach my $patient (keys %patient_to_cell_types) {
        
        print STDERR "-processing $patient\n";
        open(my $patient_allele_files_ofh, ">$patient.allelic_files") or die "Error, cannot write to file $patient.allelic_files";
        open(my $patient_cell_types_ofh, ">$patient.cell_types") or die "Error, cannot write to file $patient.cell_types";
        
        my $cell_types_href = $patient_to_cell_types{$patient};
        
        my $missing_counter = 0;
        my $found_counter = 0;
        foreach my $cell_name (keys %$cell_types_href) {
            my $cell_type = $cell_types_href->{$cell_name};
            print $patient_cell_types_ofh join("\t", $cell_name, $cell_type) . "\n";
            my $allelic_file = $cell_to_allelic_file{$cell_name};
            unless ($allelic_file) {
                ## try abridged cell name
                $missing_counter++;
                print STDERR "[missing: $missing_counter]   Error, no allelic file for cell: $cell_name of patient $patient\n";
                next;
            
            }
            print $patient_allele_files_ofh join("\t", $cell_name, $allelic_file) . "\n";
            $found_counter++;
        }
        
        close $patient_cell_types_ofh;
        close $patient_allele_files_ofh;
        
        if ($missing_counter / ($missing_counter+$found_counter) > 0.05) {
            die "Error, missing too many cells";
        }
        
        my $cmd = "$FindBin::Bin/util/make_allelic_count_matrices.pl --cell_types $patient.cell_types --allelic_counts_file_list_file $patient.allelic_files --output_prefix $patient ";
        &process_cmd($cmd);
        
        $cmd = "$FindBin::Bin/util/sparse_to_regular_matrix.R --alt_ratios $patient.alt_ratios.matrix --cell_names $patient.cells.names --chrpos_names $patient.chrpos.names --output $patient.AF.matrix";
        &process_cmd($cmd);
    }
   
    print STDERR "done.\n";
    
    exit(0);
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd $cmd died with ret $ret";
    }
    return;
}


####
sub parse_cell_type_info {
    my ($cell_types_file) = @_;

    my %cell_type_info;

    open(my $fh, $cell_types_file) or die "Error, cannot open file: $cell_types_file";
    while(<$fh>) {
        chomp;
        my ($cell_name, $cell_type, $patient) = split(/\t/);
        
        $cell_type_info{$patient}->{$cell_name} = $cell_type;
        
    }
    close $fh;

    return(%cell_type_info);
    
}

####
sub parse_cell_to_allelic_file {
    my ($allelic_counts_file_list_file) = @_;

    my %cell_name_to_file;

    open(my $fh, $allelic_counts_file_list_file) or die "Error, cannot open file: $allelic_counts_file_list_file";
    while(<$fh>) {
        chomp;
        my ($cell_name, $allelic_file) = split(/\t/);
        $cell_name_to_file{$cell_name} = $allelic_file;
    }
    
    close $fh;

    return(%cell_name_to_file);
}
