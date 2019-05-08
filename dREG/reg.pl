#!usr/bin/perl -w
use strict;

open (hand1, "/media/hp/disk2/DYY/dREG/0325out/0325peak_closer_sort.txt") or die $!;
open (hand2, ">/media/hp/disk2/DYY/dREG/0325out/dREG_peak_analysis/enhencer_re.txt");
while (<hand1>)  {
	chomp;
	if (/[CDS]/)     {

		
		print hand2 "";  }


}

close hand1; close hand2;
