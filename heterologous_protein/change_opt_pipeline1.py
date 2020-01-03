#!/usr/bin/python
# coding:utf-8
# ===========================================================================================================
# Yingying Dong.2019-12-10.Modified date: 2019-12-31.Change DNA sequence to improve expression pipeline v1.0.
#  Change the sequence codons to optimal,and avoid polyN and restriction sites.
# ===========================================================================================================
import os
import sys
import decimal
import argparse
parser = argparse.ArgumentParser(description='Change DNA sequence to improve expression.v1.0', prog='change_opt', usage='%(prog)s [options]')
parser.add_argument('--spe', nargs='?', type=str, help='species name')
parser.add_argument('--spA', nargs='?', type=str, help='species name abbreviation')
parser.add_argument('--inDNA', nargs='*', type=str, default=0, help='FASTA file for DNA sequences of genes that wish to increase expression')
parser.add_argument('--inPro', nargs='*', type=str, default=0, help='FASTA file for protein sequences of genes that wish to increase expression')
parser.add_argument('--out', nargs='*', type=str, default='optimized_seq.fa', help='file name of output optimized sequence')

args = parser.parse_args()

if (args.inDNA != 0):
    DNAfile = str(args.inDNA)
    DNAfile = DNAfile.lstrip('[\'')
    DNAfile = DNAfile.rstrip('\']')
    codon_table = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', \
                   'CGG': 'R', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', \
                   'GTG': 'V', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', \
                   'TCG': 'S', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', \
                   'CCG': 'P', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', \
                   'GAG': 'E', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', \
                   'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'TTT': 'F', 'TTC': 'F', 'TTA': 'L', \
                   'TTG': 'L', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', \
                   'TGA': 'STOP'}

    with open(DNAfile, 'r+') as f:
        f.seek(0, 2)
        f.write('>')
    f.close()

    dna_file = open(DNAfile, 'r')
    aa_file = open('aa.txt', 'w')

    dna = ''
    for line in dna_file:
        if line.startswith('>') and dna == '':
            # The sequence title is now obtained.
            header = line
            print(header)
        elif not line.startswith('>'):
            dna = dna + line.strip()
        elif line.startswith('>') and dna != '':
            prot = ''
            for i in range(0, len(dna), 3):
                codon = dna[i:i + 3]
                if codon in codon_table:
                    if codon_table[codon] == 'STOP':
                        prot = prot + '*'
                    else:
                        prot = prot + codon_table[codon]
                else:
                    # handle too short codons
                    prot = prot + '-'
            aa_file.write(header)
            j = 0
            while j < len(prot):
                print(prot[j:j + 48])
                aa_file.write(prot[j:j + 48] + '\n')
                j = j + 48
            # aa_file.write(header + prot + '\n')
            dna = ''
            header = line
    aa_file.seek(0, 2)
    aa_file.write('>')

    dna_file.close()
    aa_file.close()

if(args.inPro != 0):
    Profile = str(args.inPro)
    Profile = Profile.lstrip('[\'')
    Profile = Profile.rstrip('\']')
    print(Profile)
    with open(Profile, 'r+') as f2:
        f2.seek(0, 2)
        f2.write('>')
    f2.close()

    aa_seq = open(Profile,'r')
    AA_file = open('aa.txt','w')
    AA_seq = ''
    for line in aa_seq:
        if line.startswith('>') and AA_seq == '':
            # The sequence title is now obtained.
            header = line
            AA_file.write(header)
        elif not line.startswith('>'):
            AA_seq = AA_seq + line.strip()
        elif line.startswith('>') and AA_seq != '':
            AA_seq = AA_seq.upper()
            j = 0
            while j < len(AA_seq):
                print(AA_seq[j:j + 48])
                AA_file.write(AA_seq[j:j + 48] + '\n')
                j = j + 48
            AA_seq = ''
            header = line
    AA_file.seek(0, 2)
    AA_file.write('>')
    aa_seq.close()
    AA_file.close()

init_seq = open('aa.txt', 'r')
opt_seq_file = open('bac_Kp_opt.fa', 'w')
RSCU = open('ribo_codon_fre_RSCU.txt', 'r')
# opt_rare = open('opt_rara_sub_list','w')

aa_codon_fre = {}
RSCU_table = []
for i in RSCU:
    RSCU_table.append(i.strip().split('\t'))
for j in RSCU_table:
    aa_codon_fre[j[1]] = {j[0]: j[4]}
# print(aa_codon_fre)

optimalA = max(decimal.Decimal(x['A']) for x in aa_codon_fre.values() if
               'A' in x)  # First determine if it exists, otherwise a keyerror will occur.
keyA = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'A': str(optimalA)})]
opt_dic = {
    'A': keyA}  # This dictionary stores the optimal codon and corresponding the amino acid.In order to design the amino acid sequence set as the optimal codon DNA sequence.
optimale = max(decimal.Decimal(x['STOP']) for x in aa_codon_fre.values() if 'STOP' in x)
keye = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'STOP': str(optimale)})]

rareA = min(decimal.Decimal(x['A']) for x in aa_codon_fre.values() if 'A' in x)
raKeyA = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'A': str(rareA)})]
sub_dic = {keyA: raKeyA}  # This dictionary stores the optimal codon and rarest codon for an amino acid.


# print(opt_dic)

def findopt(aa):
    optimal = max(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
    optKey = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(optimal)})]
    opt_dic[aa] = optKey
    rare = min(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
    raKey = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(rare)})]
    sub_dic[optKey] = raKey
    return opt_dic, sub_dic


list_aa = ['C', 'D', 'E', 'F', 'G', 'H', 'K', 'I', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'Y', 'T', 'V', 'W']
for a in list_aa:
    findopt(a)
opt_dic['*'] = keye
# print(opt_dic)
# print(sub_dic)

prot = ''
for line in init_seq:
    if line.startswith('>') and prot == '':
        # The sequence title is now obtained.
        header = line
        print(header)
    elif not line.startswith('>'):
        prot = prot + line.strip()
    elif line.startswith('>') and prot != '':
        opt_dna = ''
        for i in range(0, len(prot), 1):
            j = prot[i]
            if j in opt_dic.keys():
                opt_dna = opt_dna + opt_dic[j]
            else:
                print(j + ' It is not amino acid')
        opt_seq_file.write(header)
        j = 0
        while j < len(opt_dna):
            opt_seq_file.write(opt_dna[j:j + 48] + '\n')
            j = j + 48
        prot = ''
        header = line
opt_seq_file.seek(0, 2)
opt_seq_file.write('>')

RSCU.close()
init_seq.close()
opt_seq_file.close()

# ===============================================================================================
# Remove sequence fragments with many repeated bases and restriction site , then replaced
#  with to rare codons.
# ===============================================================================================
full_opt = open('bac_Kp_opt.fa', 'r')
re_rep = open('bac_rm_rep.fa', 'w')


def findrep(base,opt_seq):
    pos = 0
    while pos >= -1:
        pos = opt_seq.find(base)
        if pos == -1:
            print('These sequences no longer have this restriction site or consecutive same nucleotides {}'.format(base))
            break
        else:
            if pos % 3 == 0:
                if opt_seq[pos:pos+3] in sub_dic.keys():
                    print(opt_seq[pos:pos + 3])
                    opt_seq = opt_seq.replace(opt_seq[pos:pos + 3], sub_dic[opt_seq[pos:pos+3]])
                    pos = opt_seq.find(base)
                    continue
            elif pos % 3 == 1:
                if opt_seq[pos-1:pos+2] in sub_dic.keys():
                    print(opt_seq[pos-1:pos+2])
                    opt_seq = opt_seq.replace(opt_seq[pos-1:pos+2], sub_dic[opt_seq[pos-1:pos+2]])
                    pos = opt_seq.find(base)
                    continue
            elif pos % 3 == 2:
                if opt_seq[pos+1:pos+4] in sub_dic.keys():
                    print(opt_seq[pos+1:pos+4])
                    opt_seq = opt_seq.replace(opt_seq[pos+1:pos+4], sub_dic[opt_seq[pos+1:pos+4]])
                    pos = opt_seq.find(base)
                    continue
            else:
                print('Something wrong')

    return opt_seq


opt_seq = ''
for line in full_opt:
    if line.startswith('>') and opt_seq == '':
        header = line
        print(header)
        re_rep.write(header)
    elif not line.startswith('>'):
        opt_seq = opt_seq + line.strip()
    elif line.startswith('>') and opt_seq != '':
        rep_posA = findrep('AAAAA',opt_seq)
        rep_posC = findrep('CCCCC',opt_seq)
        rep_posG = findrep('GGGGG',opt_seq)
        rep_posT = findrep('TTTTT',opt_seq)
        rep_EcoRI = findrep('GAATTC',opt_seq)
        rep_SalI = findrep('GTCGAC',opt_seq)
        rep_NdeI = findrep('GATATG',opt_seq)
        rep_XhoI = findrep('CTCGAG',opt_seq)

        j = 0
        while j < len(opt_seq):
            print(opt_seq[j:j + 48])
            re_rep.write(opt_seq[j:j + 48] + '\n')
            j = j + 48
        opt_seq = ''
        header = line
full_opt.close()
re_rep.close()
