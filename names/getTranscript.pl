#!/usr/bin/perl 

# get the single transcript from multi-fasta when given a list of transcript names
# usage: getTranscript.pl <names.tab> <multifasta> > file.fa
# outputs to STFOUT unless redirected

use strict;
use warnings;

my $names = $ARGV[0]; 
open(FH, '<', $names) or die $!;
chomp( my @names = <FH>);

my $fasta = $ARGV[1];
open(FA, '<', $fasta) or die $!;

sub find {
    my $match = $_[0];
    while (<FA>) {
        if (/(>$match)/s) {
            print $1;
        }
    }
}


foreach my $name (@names) {
    find($name);
}
