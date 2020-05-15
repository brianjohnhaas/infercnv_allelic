#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

my $usage = "usage: $0 input.mtx\n\n";

my $input_matrix = $ARGV[0] or die $usage;

my %ensg_to_sym;
{
    open(my $fh, "$FindBin::Bin/ensg_gene_symbol_ids.tsv") or die $!;
    while(<$fh>) {
        chomp;
        my ($ensg_id, $symbol) = split(/\t/);
        $ensg_to_sym{$ensg_id} = $symbol;
    }
    close $fh;
}


my %seen;
open(my $fh, $input_matrix) or die $!;
my $header = <$fh>;
$header =~ s/-/_/g;
print $header;
while(<$fh>) {
    my @x = split(/\t/);
    my $ensg_id = $x[0];
    my $gene_symbol = $ensg_to_sym{$ensg_id} or die "Error, no gene symbol for $ensg_id";
    
    if ($seen{$gene_symbol}) {
        $seen{$gene_symbol}++;
        $gene_symbol = $gene_symbol . ".dupsym" . $seen{$gene_symbol};
        print STDERR "-duplicate: $gene_symbol\n";
    }
    else {
        $seen{$gene_symbol} = 1;
    }

    $x[0] = $gene_symbol;
    print join("\t", @x);
    

}

exit(0);

