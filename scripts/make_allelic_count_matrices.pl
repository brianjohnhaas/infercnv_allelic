#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use List::Util qw(min max);

my $usage = "\n\n\tusage: $0 allelic_files.tsv.file.listing\n\n";

my $tsv_list_file = $ARGV[0] or die $usage;

my $MIN_TOT_COV = 10;


main: {

    my %data;
    my %cell_names;

    my $counter = 0;
    my $num_values = 0;

    open(my $fh, $tsv_list_file) or die $!;
    while(<$fh>) {
        chomp;
        my $tsv_file = $_;
        $counter++;
        print STDERR "[$counter] $tsv_file ";
    
        my $cell_name = basename($tsv_file);
        my @pts = split(/_qc./, $cell_name);
        $cell_name = $pts[0] or die "Error, no cell name from $tsv_file";
        
        $cell_names{$cell_name} = 1;
    
        &process_tsv_file($tsv_file, $cell_name, \%data, \$num_values);
    }
    close $fh;
    
    
    my @cells = sort keys %cell_names;
    
    my @chrpos = sort keys %data;

    my $num_cells = scalar(@cells);
    my $num_chrpos = scalar(@chrpos);

    print STDERR "-matrix represents $num_chrpos snps vs. $num_cells cells\n";

    {
        ## write the cell and chrpos names for use w/ the sparse matrices.
        
        {
            open (my $ofh, ">cells.names") or die $!;
            print $ofh join("\n", @cells) . "\n";
            close $ofh;
        }
        
        {
            open (my $ofh, ">chrpos.names") or die $!;
            print $ofh join("\n", @chrpos) . "\n";
            close $ofh;
        }
    }
    
    

    my $num_cell_cols = scalar(@cells);
    my $num_chrpos_rows = scalar(@chrpos);

    my $MM_header = "%%MatrixMarket matrix coordinate real general\n" 
        . join("\t", $num_chrpos_rows, $num_cell_cols, $num_values) . "\n";
    
    my $malt_allele_counts_matrix_file = "allele.malt.counts.matrix";
    open(my $ofh_malt, ">$malt_allele_counts_matrix_file") or die $!;
    print $ofh_malt $MM_header;
    
    my $tot_counts_matrix_file = "allele.tot.counts.matrix";
    open(my $ofh_tot, ">$tot_counts_matrix_file") or die $!;
    print $ofh_tot $MM_header;
    
    my $allele_malt_freq_matrix_file = "allele.malt_freqs.matrix";
    open(my $ofh_freq, ">$allele_malt_freq_matrix_file") or die $!;
    print $ofh_freq $MM_header;
    
    my $snp_counter = 0;
    foreach my $snp (@chrpos) {
        $snp_counter++;
        
        my $cell_counter = 0;
        foreach my $cell (@cells) {
            $cell_counter++;
            
            if (exists $data{$snp}->{$cell}) {
                my $alt_allele_count = $data{$snp}->{$cell}->{ALT};
                my $tot_count = $data{$snp}->{$cell}->{TOT};
                unless (defined $alt_allele_count) {
                    die "Error, undefined alt allele count for $snp, $cell ";
                }
                unless ($tot_count) {
                    die "Error, no tot count for $snp, $cell";
                }
                
                my $malt_allele_count = max($alt_allele_count, $tot_count - $alt_allele_count);
                my $malt_allele_freq = sprintf("%.3f", $malt_allele_count / $tot_count);
                
                ## output matrix elements

                print $ofh_malt join("\t", $snp_counter, $cell_counter, $malt_allele_count) . "\n";
                print $ofh_tot join("\t", $snp_counter, $cell_counter, $tot_count) . "\n";
                print $ofh_freq join("\t", $snp_counter, $cell_counter, $malt_allele_freq) . "\n";
                
            }
        }
    }
    
    close $ofh_malt;
    close $ofh_tot;
    close $ofh_freq;
    
    print STDERR "-done\n";

    exit(0);
}


####
sub process_tsv_file {
    my ($tsv_file, $cell_name, $data_href, $num_values_sref) = @_;
            
    print STDERR "-processing $cell_name\n";
    
    open(my $fh, $tsv_file) or die "Error, cannot open file: $tsv_file";
    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my ($CONTIG, $POSITION, $REF_COUNT, $ALT_COUNT, $REF_NUCLEOTIDE, $ALT_NUCLEOTIDE) = split(/\t/);
        
        my $chrpos = join("::", $CONTIG, $POSITION);
        my $tot_cov = $REF_COUNT + $ALT_COUNT;
        
        if ($tot_cov >= $MIN_TOT_COV) {

            $data_href->{$chrpos}->{$cell_name}->{ALT} = $ALT_COUNT;
            $data_href->{$chrpos}->{$cell_name}->{TOT} = $tot_cov;
            
            $$num_values_sref++;
        }
    }
    close $fh;
        
}


