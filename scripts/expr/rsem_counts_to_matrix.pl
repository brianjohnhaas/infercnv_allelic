#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 file.list.file\n\n";

my $file_list_file = $ARGV[0] or die $usage;

main: {

    my %data;

    my %cells;

    my $counter = 0;
    open(my $fh, $file_list_file) or die $!;
    while(<$fh>) {
        chomp;
        my $file = $_;
        if ($file =~ /\/([^\/]+)_rsem.genes/) {
            my $cellname = $1;
            
            $counter++;
            print STDERR "-[$counter] processing $file\n";
            &add_to_expr_matrix($cellname, $file, \%data);
            $cells{$cellname}++;
        }
        else {
            die "Error, cannot decipher $file";
        }

        #if ($counter > 10) { last;  }
    }
    close $fh;
    

    ## output matrix:
    my @gene_ids = sort keys %data;
    my @cellnames = sort keys %cells;
    
    open(my $ofh, ">cells.rsem.gene_counts.matrix") or die $!;
    
    print $ofh "\t" . join("\t", @cellnames) . "\n";

    foreach my $gene_id (@gene_ids) {
        my @vals = ($gene_id);
        
        foreach my $cellname (@cellnames) {
            my $expr_val = $data{$gene_id}->{$cellname};
            if (! defined $expr_val) {
                die "Error, no expr value for $gene_id in cell $cellname ";
            }
            push (@vals, $expr_val);
        }

        print $ofh join("\t", @vals) . "\n";
    }

    close $ofh;

    print STDERR "-done\n";
    
    exit(0);
        

}

####
sub add_to_expr_matrix {
    my ($cellname, $file, $data_href) = @_;
    
    open(my $fh, $file) or die "Error, cannot open file: $file";
    my $header = <$fh>;
    unless($header =~ /^gene_id/) {
        die "Error, file $file has unrecognized formatting";
    }
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $gene_id = $x[0];
        my $count = $x[7];
        $data_href->{$gene_id}->{$cellname} = $count;
    }
    close $fh;

    return;
}


    
