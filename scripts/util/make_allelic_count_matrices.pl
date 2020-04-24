#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use List::Util qw(min max);

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $help_flag;

my $allelic_counts_file_list_file;
my $output_prefix;


my $min_variant_coverage = 3;  # initial scanning of counts files

my $min_cells_with_coverage = 3;  # filtering of variants
my $min_het_cells_per_snp = 1;

my $min_covered_variants_per_cell = 50;  # filtering of cells.

my $low_range_het_snp_ratio = 0.25;
my $high_range_het_snp_ratio = 0.75;

my $SMALL_NONZERO_VALUE = 1e-3;


my $usage = <<__EOUSAGE__;

##################################################################################
#
#  * required:
#
#  --allelic_counts_file_list_file <string>   file containing the list of allelic counts files, one file per cell
#
#  --output_prefix <string>                   prefix for output filenames
#
#  * options w/ defaults
# 
#  Capturing variants:
#  --min_variant_coverage <int>               min read coverage at a site to be captured from allelic counts file.
#                                             (default: $min_variant_coverage)
#
#  Bulk filtering of variants:
#
#  --min_cells_with_coverage <int>            min number of cells with detected read coverage
#                                             (default: $min_cells_with_coverage)
#
#  --min_het_cells_per_snp <int>              minimum cells that must have evidence as heterozygous
#                                             (default: $min_het_cells_per_snp)
#
#  Bulk filtering of cells:
#
#  --min_covered_variants_per_cell <int>      minimum number of covered variants per cell
#                                             (default: $min_covered_variants_per_cell)
#
#  Het site requirements:
#
#  --low_range_het_snp_ratio <float>          low range for het snp ratio
#                                             (default: $low_range_het_snp_ratio)
#  
#  --high_range_het_snp_ratio <float>         high range for het snp ratio
#                                             (default: $high_range_het_snp_ratio)
#
#####################################################################################



__EOUSAGE__

    ;



&GetOptions ( 'h' => \$help_flag,
              
              'allelic_counts_file_list_file=s' => \$allelic_counts_file_list_file,
              
              'min_variant_coverage=i' => \$min_variant_coverage,
              'min_cells_with_coverage=i' => \$min_cells_with_coverage,

              'min_het_cells_per_snp=i' => \$min_het_cells_per_snp,
              
              'min_covered_variants_per_cell=i' => \$min_covered_variants_per_cell,
              
              'low_range_het_snp_ratio=f' => \$low_range_het_snp_ratio,
              'high_range_het_snp_ratio=f' => \$high_range_het_snp_ratio,

              'output_prefix=s' => \$output_prefix,

    );


if ($help_flag) {
    die $usage;
}

unless ($allelic_counts_file_list_file) {
    die $usage; 
}

main: {

    my %data;
    my %cell_names;

    my $counter = 0;
    
    print STDERR "## -parsing files.\n";
    open(my $fh, $allelic_counts_file_list_file) or die "Error, cannot open file: $allelic_counts_file_list_file";
    while(<$fh>) {
        chomp;
        my $tsv_file = $_;
        $counter++;
        print STDERR "[$counter] $tsv_file ";
    
        my $cell_name = basename($tsv_file);
        my @pts = split(/_qc./, $cell_name);
        $cell_name = $pts[0] or die "Error, no cell name from $tsv_file";
        
        $cell_names{$cell_name} = 1;
    
        &process_tsv_file($tsv_file, $cell_name, \%data, $min_variant_coverage);
        
    }
    close $fh;


    print STDERR "## -filtering variants\n";
    &filter_variants(\%data, $min_cells_with_coverage, $min_het_cells_per_snp, $low_range_het_snp_ratio, $high_range_het_snp_ratio);

    print STDERR "## -filtering cells\n";

    my @good_cells;
    my $num_values = &filter_cells(\%data, $min_covered_variants_per_cell, \@good_cells);
        
    my @chrpos = sort keys %data;

    my $num_cells = scalar(@good_cells);
    my $num_chrpos = scalar(@chrpos);

    print STDERR "-matrix represents $num_chrpos snps vs. $num_cells cells\n";

    {
        ## write the cell and chrpos names for use w/ the sparse matrices.
        
        {
            open (my $ofh, ">$output_prefix.cells.names") or die $!;
            print $ofh join("\n", @good_cells) . "\n";
            close $ofh;
        }
        
        {
            open (my $ofh, ">$output_prefix.chrpos.names") or die $!;
            print $ofh join("\n", @chrpos) . "\n";
            close $ofh;
        }
    }
    
    

    my $num_cell_cols = scalar(@good_cells);
    my $num_chrpos_rows = scalar(@chrpos);

    my $MM_header = "%%MatrixMarket matrix coordinate real general\n" 
        . join("\t", $num_chrpos_rows, $num_cell_cols, $num_values) . "\n";
    
    my $alt_allele_counts_matrix_file = "$output_prefix.alt.counts.matrix";
    open(my $ofh_alt, ">$alt_allele_counts_matrix_file") or die $!;
    print $ofh_alt $MM_header;
    
    my $tot_counts_matrix_file = "$output_prefix.tot.counts.matrix";
    open(my $ofh_tot, ">$tot_counts_matrix_file") or die $!;
    print $ofh_tot $MM_header;
    
    my $allele_alt_freq_matrix_file = "$output_prefix.alt_ratios.matrix";
    open(my $ofh_freq, ">$allele_alt_freq_matrix_file") or die $!;
    print $ofh_freq $MM_header;
    
    my $snp_counter = 0;
    foreach my $snp (@chrpos) {
        $snp_counter++;
        
        my $cell_counter = 0;
        foreach my $cell (@good_cells) {
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
                unless ($alt_allele_count > 0) {
                    $alt_allele_count = $SMALL_NONZERO_VALUE; # make it small but nonzero to differentiate from true zeros when converting to a full matrix from sparse matrix
                }
                
                my $alt_allele_freq = sprintf("%.3f", max($alt_allele_count / $tot_count, $SMALL_NONZERO_VALUE));
                
                ## output matrix elements

                print $ofh_alt join("\t", $snp_counter, $cell_counter, $alt_allele_count) . "\n";
                print $ofh_tot join("\t", $snp_counter, $cell_counter, $tot_count) . "\n";
                print $ofh_freq join("\t", $snp_counter, $cell_counter, $alt_allele_freq) . "\n";
                
            }
        }
    }
    
    close $ofh_alt;
    close $ofh_tot;
    close $ofh_freq;
    
    print STDERR "-done\n";

    exit(0);
}


####
sub process_tsv_file {
    my ($tsv_file, $cell_name, $data_href, $min_variant_coverage) = @_;
            
    print STDERR "-processing $cell_name\n";
    
    open(my $fh, $tsv_file) or die "Error, cannot open file: $tsv_file";
    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my ($CONTIG, $POSITION, $REF_COUNT, $ALT_COUNT, $REF_NUCLEOTIDE, $ALT_NUCLEOTIDE) = split(/\t/);
        
        my $chrpos = join("::", $CONTIG, $POSITION);
        my $tot_cov = $REF_COUNT + $ALT_COUNT;
        
        if ($tot_cov >= $min_variant_coverage) {
            
            $data_href->{$chrpos}->{$cell_name}->{ALT} = $ALT_COUNT;
            $data_href->{$chrpos}->{$cell_name}->{TOT} = $tot_cov;
            
        }
    }
    close $fh;
        
}


####
sub filter_variants {
    my ($data_href, 
        $min_cells_with_coverage, $min_het_cells_per_snp, 
        $low_range_het_snp_ratio, $high_range_het_snp_ratio) = @_;

    
    my %snps_to_purge;

    foreach my $snp (keys %$data_href) {
        
        my @cells = keys %{$data_href->{$snp}};
        
        ## see if sufficient number of cells
        my $num_cells = scalar(@cells);
        if ($num_cells < $min_cells_with_coverage) {
            $snps_to_purge{$snp} = 1;
            next;
        }

        # count number of het sites
        my $num_het_sites = 0;
        
        my $num_ref_only = 0;
        my $num_alt_only = 0;
        
        foreach my $cell (@cells) {
            my $alt_allele_count = $data_href->{$snp}->{$cell}->{ALT};
            my $tot_count = $data_href->{$snp}->{$cell}->{TOT};

            my $ratio = $alt_allele_count / $tot_count;
            
            if ($ratio >= $low_range_het_snp_ratio && $ratio <= $high_range_het_snp_ratio) {
                $num_het_sites++;
            }
            elsif ($ratio < $low_range_het_snp_ratio) {
                $num_ref_only++;
            }
            elsif ($ratio > $high_range_het_snp_ratio) {
                $num_alt_only++;
            }
        }
        
        if ($num_ref_only > 0 && $num_alt_only > 0) { 
            $num_het_sites += min($num_ref_only, $num_alt_only); # err on the conservative side.
        }

        
        if ($num_het_sites < $min_het_cells_per_snp) {
            $snps_to_purge{$snp} = 1;
            next;
        }
        
    }
    
    print STDERR "--purging " . scalar(keys(%snps_to_purge)) . "\n";
    foreach my $snp (keys %snps_to_purge) {
        delete $data_href->{$snp};
    }

    print STDERR "--remaining variants: " . scalar(keys(%$data_href)) . "\n";
    
    return;
}

####
sub filter_cells {
    my ($data_href, $min_covered_variants_per_cell, $good_cells_aref) = @_;
    

    my %cell_to_variant_counter;
    
    foreach my $snp (keys %$data_href) {
        
        my @cells = keys %{$data_href->{$snp}};
        
        foreach my $cell (@cells) {
            $cell_to_variant_counter{$cell}++;
        }
    }

    
    my $num_data_points = 0;
    foreach my $cell (keys %cell_to_variant_counter) {
        my $num_covered_variants = $cell_to_variant_counter{$cell};

        if ($num_covered_variants >= $min_covered_variants_per_cell) {
            push(@$good_cells_aref, $cell);
            $num_data_points += $num_covered_variants;
        }
    }
    
    return($num_data_points);
    
}

