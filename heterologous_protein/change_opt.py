#!/usr/bin/python
# coding:utf-8
# =========================================================================================================
# Yingying Dong.2019-12-10. Change the sequence codons to optimal.
# =========================================================================================================
import os
import sys
import decimal

init_seq = open('GFP_Kp_aa.txt', 'r')
opt_seq = open('GFP_Kp_opt.fa', 'w')
RSCU = open('ribo_codon_fre_RSCU.txt', 'r')

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
opt_dic = {'A': keyA}
optimale = max(decimal.Decimal(x['STOP']) for x in aa_codon_fre.values() if 'STOP' in x)
keye = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'STOP': str(optimale)})]


# print(opt_dic)

def findopt(aa):
    optimal = max(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
    opt_dic[aa] = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(optimal)})]
    return opt_dic


list_aa = ['C', 'D', 'E', 'F', 'G', 'H', 'K', 'I', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'Y', 'T', 'V', 'W']
for a in list_aa:
    findopt(a)
opt_dic['*'] = keye
# print(opt_dic)

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
        for i in range(0, len(prot),1):
            j = prot[i]
            if j in opt_dic.keys():
                opt_dna = opt_dna + opt_dic[j]
            else:
                print(j + ' It is not amino acid')
        opt_seq.write(header)
        j = 0
        while j < len(opt_dna):
            print(opt_dna[j:j + 48])
            opt_seq.write(opt_dna[j:j + 48] + '\n')
            j = j + 48
        prot = ''
        header = line

RSCU.close()
init_seq.close()
opt_seq.close()
