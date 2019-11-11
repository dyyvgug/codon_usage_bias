#!/usr/bin/python
# coding:utf-8
# =========================================================================================================
# Yingying Dong.2019-11-11. mCAI pipeline v3. The weight from genes that product ribosomal protein.The genes
#   ID from formatted GFF3 annotation file that already have CDS header.
# =========================================================================================================
import os
import sys
import argparse
import rpy2.robjects as robjects
from scipy import stats

parser = argparse.ArgumentParser(description='mCAI pipeline v3.The annotation file is formatted GFF3.', prog='mCAI_ribo_gff', usage='%(prog)s [options]')
parser.add_argument('--spe', nargs='?', type=str, help='species name')
parser.add_argument('--spA', nargs='?', type=str, help='species name abbreviation')
parser.add_argument('--ann', nargs='*', type=str, default='ref.gff', help='input gtf format file')
parser.add_argument('--inp', nargs='*', type=str, default='CDS_DNA.fa',
                    help='expected FASTA file of calculations RSCU and weight')
args = parser.parse_args()

### Get ribosomal gene ID from the gff3 annotation file.At the same time,extract their DNA sequence.

os.chdir('/media/hp/disk1/DYY/reference/annotation/{}/'.format(args.spe))
gff_path = '/media/hp/disk1/DYY/reference/annotation/{}/'.format(args.spe)

gff = open('{}{}'.format(gff_path, args.ann), 'r')
os.system('cp ref.gff refbackup.gff')
os.system('sed -i \'s/ID.*Parent=//g\' ref.gff')
os.system('sed -i \'s/;Dbxref.*product=/\t/g\' ref.gff')
os.system('sed -i \'s/;protein_id.*//g\' ref.gff')
os.system('sed -i \'/ribosomal/!d\' ref.gff')

r_script1 = '''
gff = read.table('ref.gff',sep = '\t',header = F,quote = '')
gff = unique(gff)
write.table(gff$V9,file = 'ribo_geneID_u.txt', quote = FALSE,row.names = F, col.names = F)
'''
robjects.r(r_script1)

id_file = open('ribo_geneID_u.txt', 'r')

DNA_path = '/media/hp/disk1/DYY/reference/annotation/{}/ref/'.format(args.spe)
DNA = open('{}CDS_DNA.fa'.format(DNA_path), 'r')
ex_seq = open('{}ribo_seq.fa'.format(DNA_path), 'w')

id_table = []
for i in id_file:
    id_table.append(i.strip())

database = {}
for m in DNA:
    if m.startswith('>'):
        keys = m.lstrip('>').strip()
        database[keys] = []
    else:
        database[keys].append(m.strip())
for line in id_table:
    if line in database.keys():
        ex_seq.write('>' + str(line) + '\n' + str(''.join(database[line])) + '\n')
    else:
        print(line + ' Its sequence may not be coding sequence')

gff.close()
id_file.close()
ex_seq.close()
### Calculates codon frequency,RSCU value,weight in the gene FASTA sequence file.

aa_codon = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'],
    'K': ['AAA', 'AAG'], 'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'],
    'N': ['AAT', 'AAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'Y': ['TAT', 'TAC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'],
    'STOP': ['TAG', 'TAA', 'TGA']
}
codon_count = {
    'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0,
    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0,
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
    'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
    'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
    'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0,
    'TGT': 0, 'TGC': 0, 'TGG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 'TGA': 0
}
aa_num = {
    'F': 2, 'Y': 2, 'C': 2, 'H': 2, 'Q': 2, 'N': 2, 'K': 2, 'D': 2, 'E': 2,
    'P': 4, 'T': 4, 'V': 4, 'A': 4, 'G': 4,
    'L': 6, 'R': 6, 'S': 6,
    'W': 1, 'M': 1, 'STOP': 1,
    'I': 3,
}


# Write the frequency of each codon to a file.
def calc_freq(codon_count, out_file):
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
                RSCU = codon_count[codon] / count_tot[aa] * aa_num[aa]
            else:
                freq = 0.0
                RSCU = 0.0
            out_file.write('{}\t{}\t{}\t{:.3f}\t{:.3f}\n'.format(aa, codon, codon_count[codon], freq, RSCU))


in_file = open('{}ribo_seq.fa'.format(DNA_path), 'r')
out_file = open('ribo_codon_fre_RSCU.txt', 'w')

# Reads the DNA sequence into a single string.
dna = ''
for line in in_file:
    if not line.startswith('>'):
        dna = dna + line.strip()

# Scans the sequence frame by frame,counts the number of occurrences of each codon,and stores it in codon_count dictionary.
# Then calls calc_freq()
out_file.write(' AA\tcodon\thits\tfrequency\tRSCU\n')

for i in range(0, len(dna), 3):
    codon = dna[i:i + 3]
    if codon in codon_count:
        codon_count[codon] = codon_count[codon] + 1
calc_freq(codon_count, out_file)

in_file.close()
out_file.close()

r_script = '''
gene_fre = read.table("ribo_codon_fre_RSCU.txt",header = T,sep = '\t',quote = "")
df <- gene_fre
df$Weights <- ave(df$RSCU,df$AA,FUN=function(x) x/max(x)) # calculate weight
df = df[,-c(1,3,4,5)]
df = df[-c(30,61,62,63,64),]
write.table(df,file = "ribo_RSCU_weight.txt",sep = '\t',quote = F,row.names = F,col.names = F) # for calculate CAI
'''
robjects.r(r_script)

### Calculates CAI value.

fa_file = open('{}{}'.format(DNA_path, args.inp), 'r')
weight_file = open("ribo_RSCU_weight.txt", 'r')
CAI_file = open('{}{}_mCAI_ribo.txt'.format(DNA_path, args.spA), 'w')

weight_table = []
for line in weight_file:
    weight_table.append(line.strip().split('\t'))
codon_weight = {}
for i in weight_table:
    codon_weight[i[0]] = float(i[1])
# print (codon_weight)

CAI_file.write('gene_id\tmCAI_value\n')
# Reads the DNA sequence into a single string.
dna = ''
weight_list = []
for line in fa_file:
    if line.startswith('>') and dna == '':
        # The sequence title is now obtained.
        header = line.strip()
    elif not line.startswith('>'):
        dna = dna + line.strip()
    elif line.startswith('>') and dna != '':
        for j in range(0, len(dna), 3):
            codon = dna[j:j + 3]
            if codon in codon_weight:
                weight_list.append(codon_weight[codon])
        #print(header)
        # print(weight_list)
        CAI = stats.gmean(weight_list)
        CAI_file.write('{}\t{}\n'.format(header, CAI))
        header = line.strip()
        dna = ''
        weight_list = []

fa_file.close()
CAI_file.close()
weight_file.close()

os.chdir(DNA_path)
os.system('sed -i \'s/>//g\' {}_mCAI_ribo.txt'.format(args.spA))