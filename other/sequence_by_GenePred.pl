#!/usr/bin/perl -w
use List::Compare;
use strict;
# ===========================================================================================================================================
# This original script came from Yunkun Dang, it's aimed at get the CDS from genome FASTA file and annotation gtf file. Since the gtf file is 
# out of fashion, so Yingying Dong modified this script, we can get the CDS from GFF3 files.
# ===========================================================================================================================================

my $ge ; my %genome; my $species = "Apis_mellifera" ;  

open (hand2, "./$species/genome.fa") or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ge=$_;
    $ge=~ s/^>//;
    $ge=~ s/ .+$//; 
                             
    print "$ge\n";           
    next;}
    $genome{$ge} .= $_; }

close hand2;

my %cds ; my %utr ; my %cod ;
open (REF, "./$species/ref.txt") or die $!;
open (Seq, "./$species/CDS");  

while (<REF>)   {
	
   chomp;
   my @a = split /\t/;
   my @ins = split /,/, $a[8];
   my @ine = split /,/, $a[9];
   my $cds_count = 0;
   my $utr_count = 0;
   print "processing $a[0]\n";
 
   my @cds = ();
    if ($a[6] - $a[5] > 0)   {                               # cdsEnd-cdsStart,coding region.
	   
		for my $i ( 0..$a[7]-1 )  {                  # 0~exonCount-1
			
			for my $j ($ins[$i]..$ine[$i]-1)  {  # Add elements from exonEnds the first column to exonEnds the zero column,column++
				                             # until exonCount is 0.
				push @cds, $j; }  }          # Push: add an element from the end of the array.
		my $full = @cds;
		#print "all=".join(",",@cds)."\n";	     # Rough view the all Coding region.

		my @utrL = ();
                my @utrR = ();
		
		@utrL = $a[3]..$a[5]-1 if $a[3] - $a[5] != 0;
		@utrR = $a[6]..$a[4] if $a[4] - $a[6] != 0;                   

		my $lcl = List::Compare->new(\@cds, \@utrL); # List::Compare - Compare elements of two or more lists.
		@cds = $lcl->get_unique;		    
		@cds = sort {$a <=> $b} @cds;
                
		my $lcr = List::Compare->new(\@cds, \@utrR);
		@cds = $lcr->get_unique;
		@cds = sort {$a <=> $b} @cds;
		
		my $CDS = @cds;
		my $cds_intron = $a[6]-$a[5]-$CDS;

		my $seq ; 					     
		my %cds_exon ;
		my $position = $cds[0];                                                  # Position is the first element of codon region sequence.
		my $ex_number = 1;

		foreach my $i (@cds)  {
			if ($i <= $position + 1)  {                                      #If()loop function,count the number of exon in coding region.
				#print "$i,$position\t";
				$cds_exon{$ex_number} .= substr($genome{$a[1]}, $i, 1) ; #Substr($string,offset,length)function,intercept the chromosome name 
				$position ++;  }
			else {
				$ex_number ++;
				$cds_exon{$ex_number} .= substr($genome{$a[1]}, $i, 1) ;
				$position = $i;  }					 
			
			$seq .= substr($genome{$a[1]}, $i, 1) ;
			$seq =~ tr/atcg/ATCG/;
			}                                                              
			
		if ($a[2] eq '-')  {                                                     
			foreach my $id (keys %cds_exon)  {                               # Take out the keys from %cds_exon.
				$cds_exon{$id} = reverse $cds_exon{$id};                 # Get the exon in coding region of anti-chain.
				$cds_exon{$id} =~ tr/ATCG/TAGC/;  }

			$seq = reverse $seq;						
			$seq =~ tr/ATCG/TAGC/;  }

		print Seq ">$a[0]\n$seq\n";
	}  }
		
close REF; close Seq; 
