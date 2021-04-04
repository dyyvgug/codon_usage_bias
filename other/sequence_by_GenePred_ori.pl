#£¡usr/bin/perl -w
use List::Compare;
use strict;
#use Bio::Util::codonUsage qw(translate);

=head1 the purpose
This script is to extract the sequence from GenePred file, which can be created with gtf or gff3 with UCSC tools gtfToGenePred or gff3ToGenePred.
It is still unknown why some sequence give wrong coding region, but these sequences, or id will be stored in a translate_warning files, with original DNA sequence and translation. 
As a refrence, this warning file will also be taken to excluding them in the further analyses of PRF.
For codon frequency, these warning data are not included as well.  
=cut
#---------------------------------------------------------------------------------------------------------------------------
# Step 0.Assign numerical values to species,ge(gene),files1(FASTA file).
# Hash %genome arrays are used to store genetic information. 
# Set the species name as a variable to facilitate scripts to run in large batches.
#---------------------------------------------------------------------------------------------------------------------------
mkdir "ref";
my $ge ; my %genome; my $species = "Apis_mellifera" ;  
#---------------------------------------------------------------------------------------------------------------------------
#Step 1.Remove the first line of fasta files,and remove "enter".
#---------------------------------------------------------------------------------------------------------------------------
open (hand2, "/media/hp/disk1/DYY/reference/genome/$species/genome.fa") or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
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
#---------------------------------------------------------------------------------------------------------------------------
#Step 2.Declaring hash array to cds(Coding region ),utr(Untranslated Region),cod().
#Use GenePred reference file the annotate the peaks, either at 3' end or internal CDS.
#---------------------------------------------------------------------------------------------------------------------------	
	
my %cds ; my %utr ; my %cod ;
open (REF, "/media/hp/disk1/DYY/reference/annotation/$species/ref.txt") or die $!;
open (Seq, ">ref/CDS_DNA.fa");  
open (P, ">ref/CDS_pep.fa");
open (W, ">ref/translate_warning.txt");
open (Size, ">ref/CDS_intron_size.txt");
print Size "gene_id\ttranscript\tCDS\tCDS_intron\n";
open (cds_exon, ">ref/CDS_exons.txt");
#-----------------------------------------------------------------------------------------------------------------------------
# Step 3.According to the GenePred annotation file, count the number of accumulated coding areas.
#-----------------------------------------------------------------------------------------------------------------------------

while (<REF>)   {
	
   chomp;     #Remove the "enter" at the end of the string.
   my @a = split /\t/;                                       # GenePred files separated by "Tab".
   my @ins = split /,/, $a[8];                               # exonStarts(Exon start positions)separated by commas.
   my @ine = split /,/, $a[9];                               # exonEnds(Exon End positions)separated by commas.
   my $cds_count = 0;
   my $utr_count = 0;
   print "processing $a[0]\n";                               #View progress,$a[0] is geneName.
   # define CDS region, first get all exon.
   my @cds = ();
    if ($a[6] - $a[5] > 0)   {                               # cdsEnd-cdsStart,coding region.
	   
		for my $i ( 0..$a[7]-1 )  {                  # 0~exonCount-1
			
			for my $j ($ins[$i]..$ine[$i]-1)  {  # Add elements from exonEnds the first column to exonEnds the zero column,column++
				                             # until exonCount is 0.
				push @cds, $j; }  }          # Push: add an element from the end of the array.
		my $full = @cds;
		#print "all=".join(",",@cds)."\n";	     # Rough view the all Coding region.
#---------------------------------------------------------------------------------------------------------------------------------------------
# Step 4.According to the GenePred annotation file, calculate the position of the left untranslated region and the position of the right 
# untranslated region.
#---------------------------------------------------------------------------------------------------------------------------------------------	
		my @utrL = ();                               # Declaring array to utrL(the left of untranslated region).
                my @utrR = ();			             # Declaring array to utrR(the right of untranslated region).
		
		@utrL = $a[3]..$a[5]-1 if $a[3] - $a[5] != 0;# When txStart > cdsStart,
		                   			     # utrL = txStart(transcription start)~cdsStart(Coding region start).
		@utrR = $a[6]..$a[4] if $a[4] - $a[6] != 0;  # When txEnd > cdsEnd,                     
						             # utrR = cdsEnd(Coding region end)~txEnd(transcription end).
#----------------------------------------------------------------------------------------------------------------------------------------------
# Step 5.Compare elements of two lists (@cds and @utr)by the List::Compare module,get those items which appear only in the @cds list.
# That is calculating the number of coding region excludes the number of untranslated region.
# Then calculating the number of introns in the coding region.
#----------------------------------------------------------------------------------------------------------------------------------------------
		my $lcl = List::Compare->new(\@cds, \@utrL); # List::Compare - Compare elements of two or more lists.
		@cds = $lcl->get_unique;		     # Get those items which appear (at least once) only in the @cds list.
		@cds = sort {$a <=> $b} @cds;                # Sort by numerically.
		#print join(",",@cds)."\n";   
                
		my $lcr = List::Compare->new(\@cds, \@utrR);
		@cds = $lcr->get_unique;
		@cds = sort {$a <=> $b} @cds;
		#print join(",",@cds)."\n";                   # Check the coding region to remove the left and right untranslated region.
		
		my $CDS = @cds;                              # The array @cds is assigned to the variable $CDS to get the length of the array.
		my $cds_intron = $a[6]-$a[5]-$CDS;	     # The variable $cds_intron is cdsEnd - cdsStart - the count of coding region.
		print Size "$a[0]\t$full\t$CDS\t$cds_intron\n"; 
# geneName\t the_number_of_all_coding_region\t the_number_of_coding_region_excludes_the_number_of_untranslated_region\t the_number_of_introns_ # # # ## in_the_coding_region.
#-------------------------------------------------------------------------------------------------------------------------------------------------
# Step 6.Calculating the number of exons in the coding region.
#-------------------------------------------------------------------------------------------------------------------------------------------------
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
			#print "\n\n";
			$cds_exon{$ex_number} .= substr($genome{$a[1]}, $i, 1) ;
			$position = $i;  }					 #
		
		$seq .= substr($genome{$a[1]}, $i, 1) ;
		$seq =~ tr/atcg/ATCG/;                                           # FASTA file,some sequence are in lower case,
		}                                                                # but perl language distinguishes between uppercase and   
                                                                                 # lowercase letters.
		
	if ($a[2] eq '-')  {                                                     # eq "==",eq usually comparing strings.When - for strand.
		foreach my $id (keys %cds_exon)  {                               # Take out the keys from %cds_exon.
			$cds_exon{$id} = reverse $cds_exon{$id};                 # Get the exon in coding region of anti-chain.
			$cds_exon{$id} =~ tr/ATCG/TAGC/;  }

		$seq = reverse $seq;						 # Get the anti-chain.
		$seq =~ tr/ATCG/TAGC/;  }

	
	
	if ( $a[2] eq "-")  {                                                    # When - for DNA strand. 
		my $exon_order = 1;                                              
		foreach my $id (sort {$b <=> $a} keys %cds_exon)  {              # Arrange keys from small to large.
			print cds_exon ">$a[0]\_exon_$exon_order\n$cds_exon{$id}\n";
			$exon_order ++; } 
	}
	else { 								         # When + for DNA stamd.
		foreach my $id (sort {$a <=> $b} keys %cds_exon)  {              # Arrange keys from large to small. 
			print cds_exon ">$a[0]\_exon_$id\n$cds_exon{$id}\n"; }
	}
# Output cds_exon form: geneName\ Corresponding exon \ Corresponding conding region.
#----------------------------------------------------------------------------------------------------------------------------------------------
# Step 7.Calculate the amino acid chain that the current DNA strand will translate according to the central rule.
#----------------------------------------------------------------------------------------------------------------------------------------------
	my $protein; 
	for ( my $i = 0; $i < length($seq)-1; $i += 3 )    {
		my $aa = translate(substr($seq, $i, 3));                         # substr()funcution gets codons(adjacent three nucleotides).
		$protein .= $aa; }                                               # We all know that three nucleotides correspond to one amino acid.

	if ($protein !~ /X|^.+\*.+$/ and $protein =~ /^M.+\*$/)  {               # When three codons do not contain X or *, but start with M.
	
		print Seq ">$a[0]\n$seq\n";                                      # Output Seq form: geneName\sequence
		print P ">$a[0]\n$protein\n";                                    # Output P form:geneName\corresponding protein
		# count codons in proteins that are 
		# 1, started with ATG and stop with stop codons
		# 2, no internal * and no irrugular X, here X are not right codon, either cotain N or not 3 nt.
		for ( my $i = 0; $i < length($seq)-1; $i += 3 )    {
			my $cod = substr($seq, $i, 3);                           # Get codons.
			$cod{$cod} ++; }      }
	else {	
		print W ">$a[0]\n$protein\n$seq\n";}                             # Until sequence end,output geneName\corresponding protein\seq
	}  }
		
close REF; close Seq; close cds_exon;	

open (hand4, ">ref/codon_frequency.txt");

foreach my $c (keys %cod)   {                                                    # Take out the keys from %cod(codons hash).
	my $aa = translate($c);
	print hand4 "$c\t$aa\t$cod{$c}\n";}                                      

close hand4;




sub translate {                                                 # In the translation process, adjacent three nucleotides correspond to one amino acid.
 my ($cod) = shift;
 
 my (%codon2aa) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
 
 if (exists $codon2aa{$cod}) {					# Some codons cannot correspond to the standard codon table.
  	return $codon2aa{$cod};  } 
 else { 
	print "warning, $cod is not a codon\n";
	return "X";}
              }
