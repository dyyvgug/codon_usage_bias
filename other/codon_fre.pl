#!usr/bin/perl -w
use strict;

my $geneName; my %CDS;
open (hand1, "CDS_DNA.fa") or die $!
while (<hand1>) {
	chomp;
	if (/>/)  {
		$geneName = $_;
		next}
	$CDS{$geneName} = $_; }

close hand1;

my %codonFrq;

open (hand2, "top200.txt") or die $!;
while (<hand2>) {
	chomp;
	my @a = split(\t);
	if (exists $CDS{$a[0]} ) {
		my $seq = $CDS{$a[0]};
			for (my $i=0, $i < length($seq), $i += 3) {
				my $codon = substr($seq, $i, 3);
				$codonFrq{$codon} ++; }
				}
				}


open (hand3, ">codon_frq_top200.txt");
foreach my $i (keys %codonFrq)   {
	print hand3 "$i\t$codonFrq{$i}\n"; }
	
	close hand2;
	close hand3;
	
	
	
