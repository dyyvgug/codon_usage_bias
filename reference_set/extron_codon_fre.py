#!/usr/bin/python
#coding:utf-8

# This program calculates the codon frequency in DNA sequence.
aa_codon = {
    'A':['GCT','GCC','GCA','GCG'],'C':['TGT','TGC'],
    'D':['GAT','GAC'],'E':['GAA','GAG'],'F':['TTT','TTC'],
    'G':['GGT','GGC','GGA','GGG'],'H':['CAT','CAC'],
    'K':['AAA','AAG'],'I':['ATT','ATC','ATA'],
    'L':['TTA','TTG','CTT','CTC','CTA','CTG'],'M':['ATG'],
    'N':['AAT','AAC'],'P':['CCT','CCC','CCA','CCG'],
    'Q':['CAA','CAG'],'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
    'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
    'Y':['TAT','TAC'],'T':['ACT','ACC','ACA','ACG'],
    'V':['GTT','GTC','GTA','GTG'],'W':['TGG'],
    'STOP':['TAG','TAA','TGA']
}
codon_count = {
    'GCT':0,'GCC':0,'GCA':0,'GCG':0,'CGT':0,'CGC':0,'CGA':0,'CGG':0,
    'ACT':0,'ACC':0,'ACA':0,'ACG':0,'GTT':0,'GTC':0,'GTA':0,'GTG':0,
    'GGT':0,'GGC':0,'GGA':0,'GGG':0,'TCT':0,'TCC':0,'TCA':0,'TCG':0,
    'CTT':0,'CTC':0,'CTA':0,'CTG':0,'CCT':0,'CCC':0,'CCA':0,'CCG':0,
    'ATT':0,'ATC':0,'ATA':0,'ATG':0,'GAT':0,'GAC':0,'GAA':0,'GAG':0,
    'CAT':0,'CAC':0,'CAA':0,'CAG':0,'AAT':0,'AAC':0,'AAA':0,'AAG':0,
    'AGT':0,'AGC':0,'AGA':0,'AGG':0,'TTT':0,'TTC':0,'TTA':0,'TTG':0,
    'TGT':0,'TGC':0,'TGG':0,'TAT':0,'TAC':0,'TAA':0,'TAG':0,'TGA':0
}

# Write the frequency of each codon to a file.
def calc_freq(codon_count,out_file):
    count_tot = {}
    for aa in aa_codon.keys():
        n = 0
        for codon in aa_codon[aa]:
            n = n + codon_count[codon]
        count_tot[aa] = float(n)
    for aa in aa_codon.keys():
        for codon in aa_codon[aa]:
            if count_tot[aa] != 0.0:
                freq = codon_count[codon] / count_tot[aa]
            else:
                freq = 0.0
            out_file.write('{:>4}\t{:>5}\t{:>4}\t{:>5.3}\n'.format(aa,codon,codon_count[codon],freq))

in_file = open('40_lE_lT_seq_209.fa','r')
out_file = open('40_lE_lT_codon_fre.txt','w')

# Reads the DNA sequence into a single string.
dna = ''
for line in in_file:
    if not line.startswith('>'):
        dna = dna + line.strip()

# Scans the sequence frame by frame,counts the number of occurrences of each codon,and stores it in codon_count dictionary.
# Then calls calc_freq()
out_file.write(' AA\tcodon\thits\tfrequency\n')
prot = ''
for i in range(0,len(dna),3):
    codon = dna[i:i+3]
    if codon in codon_count:
        codon_count[codon] = codon_count[codon] + 1
calc_freq(codon_count,out_file)

in_file.close()
out_file.close()
