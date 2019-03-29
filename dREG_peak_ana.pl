#！usr/bin/perl -w
use strict;
=head1 the purpose
提取出dREG包计算出的peak对应的fasta序列，该peak预示着转录起始区域。该组peak来源于人慢粒白血病细胞K562的GRO-seq，GSM号为GSM1480325。
=cut

mkdir "dREG_peak_analysis";
my $ge ; my %genome; my $species = "hg19" ;
open (hand2, "/media/hp/disk1/DYY/reference/genome/$species/genome.fa") or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;       # Remove"Enter" of per line.
    if (/^>/)       {        # If the ">" symbol is found in reference gennome file.
    $ge=$_;
    $ge=~ s/^>//;            # Remove the ">" symbol.
    $ge=~ s/ .+$//;          # Remove the space bar and the content behind it.Beacause the reference genome file from NCBI,The naming rule for each line of 
                             #  the fasta file is that there is other information besides the id information. In order to extract the information later, only    
                             #  the id information is extracted. 
    print "$ge\n";           
    next;}
    $genome{$ge} .= $_; }
close hand2;
#=========================================================================================================
# This part of bedgraph is for PRO-seq.
#=========================================================================================================

my %plus;
open (hand3, "../GSM1480325_K562_GROseq_plus.bedGraph") or die $!;
while (<hand3>)   {
	$_=~ s/\s+$//;       # Remove"Enter" of per line.
	my ($chr,$start,$end,$val) = split /\t/;
	next if $val == 0;
	for my $i ($start..$end-1)  {
		    	$plus{$chr."_".$i} = $val; }
         }
close hand3;

my %minus;
open (hand4, "../GSM1480325_K562_GROseq_minus.bedGraph") or die $!;
while (<hand4>)   {
	$_=~ s/\s+$//;       # Remove"Enter" of per line.
	my ($chr,$start,$end,$val) = split /\t/;
	next if $val == 0;
	for my $i ($start..$end-1)  {
		    	$minus{$chr."_".$i} = $val; }
         }
close hand4;

##=========================================================
# This part of bedgraph is for GRO-seq
##=========================================================
my %gro_plus;
open (hand5, "../GSM1480325_K562_GROseq_plus.bedGraph") or die $!;
while (<hand5>)   {
	$_=~ s/\s+$//;       # Remove"Enter" of per line.
	my ($chr,$start,$end,$val) = split /\t/;
	next if $val == 0;
	for my $i ($start..$end-1)  {
		    	$gro_plus{$chr."_".$i} = $val; }
         }
close hand5;

my %gro_minus;
open (hand6, "../GSM1480325_K562_GROseq_minus.bedGraph") or die $!;
while (<hand6>)   {
	$_=~ s/\s+$//;       # Remove"Enter" of per line.
	my ($chr,$start,$end,$val) = split /\t/;
	next if $val == 0;
	for my $i ($start..$end-1)  {
		    	$gro_minus{$chr."_".$i} = $val; }
         }
close hand6;

##=========================================================
# This part is peak of dREG
##=========================================================

open (dREG, "0325peak_closer_sort.txt") or die $!;
open (F, ">dREG_peak_analysis/peak_analysis.txt");
#print F "transcript_id\tdREG_score\tpromoter_signal\tpromoter_signal_density\tGRO_seq_gene_body_density\n";

# output file should contain gene|transcript name, dREG score, dREG region average signal, which use the reference file(bedgraph) to calcualate the density of signal 
while (<dREG>)   {
	   chomp;     
	next if $_ !~ /transcript_type=protein_coding/;
   	my @a = split /\t/;                                  
   	if ( $a[11] eq "transcript" )  {                      	  # ~ end positions.
		my $reg_start = $a[1];
		my $reg_end = $a[2];
		my $chr = $a[0];
		my $reg_score = $a[3];
		my $gene_id = $a[7];
		my $transcript_start = $a[5];
		my $transcript_end = $a[6];
		
		my $distance;
		if ($a[9] eq "+") {
			$distance = abs($reg_start - $transcript_end);}
		else {
			$distance = abs($reg_end - $transcript_start);}

		next if $distance > 500; # with this creteria we should known that the region is not as enhancer.
#====================================================================================
# this part is PRO-seq. 
#====================================================================================

		my $signal = 0;
		for my $i ($reg_start..$reg_end-1)   {
			$signal += $plus{$chr."_".$i} if exists $plus{$chr."_".$i};
			$signal += $minus{$chr."_".$i} if exists $minus{$chr."_".$i};  }

		my $density = $signal/($reg_end - $reg_start);


#=====================================================================================
# this part is for GRO-seq elongation region
#=====================================================================================
		my $gro_elo = 0;
		if ($a[9] eq "-")  {
			for my $j ($transcript_start..$transcript_end-120)   {
				$gro_elo += $minus{$chr."_".$j} if exists $minus{$chr."_".$j};  } }

		else {
			for my $j ($transcript_start+120..$transcript_end)   {
				$gro_elo += $plus{$chr."_".$j} if exists $plus{$chr."_".$j}; }}

		my $elo_density = $gro_elo/($transcript_end - $transcript_start);
			
        print F "$chr\t$gene_id\t$reg_score\t$signal\t$density\t$gro_elo\t$elo_density\n"; 
				
				}
}
close dREG; close F;


