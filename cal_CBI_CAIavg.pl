#!/usr/bin/perl
use strict;
use warnings;
use Number::Format qw/round/;
use Statistics::Basic qw(correlation);
use Statistics::Descriptive; 

my $stat = Statistics::Descriptive::Full->new();
# create a table by calculating the CBI with self defined random and optimized codons for mouse and human. 
# average CAI of each gene
# report the length of ORF and gene 


&cbi ('Apis_mellifera');
#&cbi ('mouse.GRCm38');

sub cbi                           {

   my ($sample) = @_;

   my %cai; my %score;
   open (hand1, "CBI_ref.txt") or die $!;   #(K\tAAA\t0.25\t1)
   while (<hand1>) {
       next if /^AA/;
       $_=~ s/\s+$// ;
       my @a = split /\t/;
       #print join(",",@a)."\n";
       $cai{$a[1]} = $a[2];
       $score{$a[1]} = $a[3];  
                   }
   close hand1;

   my %seq; my $trans; my $gene;
   open (hand2, "CDS_DNA.fa") or die $!;
   while (<hand2>)   {
       chomp;
       if (/^>/)   {
             $gene = $_; 
             $gene =~ s/^>//;
             next  }
        
       $seq{$gene} .= $_;  
                     }
    close hand2;

   my @cai=(); my @cbi = (); my %uni_seq; 

    open (hand3, ">$sample\_CBI_CAIavg.txt");
 
    print hand3 "transcription_id\tlength\tGC\tavg_CAI\tCBI\n";
    
    foreach my $id (sort keys %seq)   {
        my $si = length($seq{$id});
        print hand3 "$id\t$si";
        
        #count gc content
        my $c = () = $seq{$id} =~ /[cC]/g;
        my $g = () = $seq{$id} =~ /[gG]/g;
        my $gc = round(($c+$g)/$si,3);
        print hand3 "\t$gc";
   
        my $sum_opt = 0; my $sum_ran = 0; my @cai_gene; my $count = 0;
        for (my $i=0; $i<length($seq{$id}); $i+=3)   {
           my $cod = substr($seq{$id},$i,3);
           next if $cod =~ /TAG|TAA|TGA|ATG|TGG/;
           $count++;
           push @cai_gene, $cai{$cod} if exists $cai{$cod};
           $sum_opt++ if exists $score{$cod} and $score{$cod} == 3; 
           $sum_ran++ if exists $score{$cod} and $score{$cod} == 2; 
                                                      }
	$stat->add_data(@cai_gene);
        my $cai_avg = round($stat->geometric_mean(), 3);
	$stat->clear();
        my $cbi = round(($sum_opt-$sum_ran)/($count-$sum_ran),3);
        push @cai, $cai_avg;
        push @cbi, $cbi;
        print hand3 "\t$cai_avg\t$cbi\n";    }
     close hand3;
     my $cor = correlation(\@cai,\@cbi);
     print "\n\n$sample avg_CAI to CBI Pearson's r = $cor\n\n";

	
                                 }
