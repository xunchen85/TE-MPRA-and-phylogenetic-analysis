#!/usr/bin/perl

use strict;
use warnings;

my %seqs;
$/ = "\n>"; # read fasta by sequence, not by lines

while (<>) {
	s/>//g;
	my ($seq_id, @seq) = split (/\n/, $_);
	my $seq = uc(join "", @seq); # rebuild sequence as a single string
	my $len = length $seq;
	my $numA = $seq =~ tr/A//; # removing A's from sequence returns total counts
	my $numC = $seq =~ tr/C//;
	my $numG = $seq =~ tr/G//;
	my $numT = $seq =~ tr/T//;
	my $numN = $seq =~ tr/N//;
	my $num_ = $seq =~ tr/-//;
	my $GC_perC = ($numC + $numG)/($len - $numN);
	my $num_nt = $len-$num_-$numN;
	if ($num_nt >0) {
            print ">$seq_id\n$seq\n"
	}
	#print "$seq_id: Size=$len non_nt=$num_nt GC=$GC_perC A=$numA  C=$numC  G=$numG  T=$numT N=$numN gap=$num_\n";
	    }
