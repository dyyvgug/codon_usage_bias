#！usr/bin/perl -w

use strict;

# first step, load the TPM value of ribo-seq. Each gene has its own expression level

my $flank = 100;my $ge ; my %genome; my %val;my %cod;my %trna;

open (hand1,"/home/hp/Desktop/other_riboseq/C_elegans_Ensl_WBcel235/experiment2/aligned/SRR1804340_abund.out") or die $!;
while (<hand1>)   {
	chomp;
	my @a = split /\t/; 
	$val{$a[0]} = $a[8];
}
close hand1;


# open cDNA fasta file 
open (hand2, "/media/hp/disk1/DYY/reference/annotation/C_elegans_Ensl_WBcel235/ref/CDS_DNA.fa") or die $!;
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


# calculate , open the file contain genes that is highly transcribed and highly translated. This should come from your analyses
open (hand3,"/home/hp/Desktop/other_riboseq/C_elegans_Ensl_WBcel235/experiment2/aligned/ribo_num/SRR1804340_hE_hT_only_geneID.txt") or die $!;
# assume only 1 column
while (<hand3>)   {
	chomp;
	#if (exists genome{$_} and exists val{$_} )   {
		my $seq = $genome{$_};
		for ( my $i = 3 + $flank ; $i < length($seq)-1 - $flank; $i += 3 )    {
			my $cod = substr($seq, $i, 3); 
                          # Get codon frequency;
			$cod{$cod} ++; 
			$trna{$cod} += $val{$_}; }      #}
close hand3;

open (hand4, ">40_codon_frequency.txt");

foreach my $c (keys %cod)   {                                                    # Take out the keys from %cod(codons hash).
	my $aa = translate($c);
	print hand4 "$c\t$aa\t$cod{$c}\n";}                                      

close hand4;

open (hand4, ">40_tRNA_demand.txt");

foreach my $c (keys %trna)   {                                                    # Take out the keys from %cod(codons hash).
	my $aa = translate($c);
	print hand4 "$c\t$aa\t$trna{$c}\n";}                                      

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
}
