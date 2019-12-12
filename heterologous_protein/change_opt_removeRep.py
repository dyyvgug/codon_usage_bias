#!/usr/bin/python
# coding:utf-8
# =========================================================================================================
# Yingying Dong.2019-12-10.Modified date: 2019-12-12. Change the sequence codons to optimal.
# =========================================================================================================
import os
import sys
import decimal

init_seq = open('GFP_Kp_aa.txt', 'r')
opt_seq_file = open('GFP_Kp_opt.fa', 'w')
RSCU = open('ribo_codon_fre_RSCU.txt', 'r')
#opt_rare = open('opt_rara_sub_list','w')

aa_codon_fre = {}
RSCU_table = []
for i in RSCU:
    RSCU_table.append(i.strip().split('\t'))
for j in RSCU_table:
    aa_codon_fre[j[1]] = {j[0]: j[4]}
# print(aa_codon_fre)

optimalA = max(decimal.Decimal(x['A']) for x in aa_codon_fre.values() if 'A' in x)  # First determine if it exists, otherwise a keyerror will occur.
keyA = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'A': str(optimalA)})]
opt_dic = {'A': keyA}       # This dictionary stores the optimal codon and corresponding the amino acid.In order to design the amino acid sequence set as the optimal codon DNA sequence.
optimale = max(decimal.Decimal(x['STOP']) for x in aa_codon_fre.values() if 'STOP' in x)
keye = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'STOP': str(optimale)})]

rareA = min(decimal.Decimal(x['A']) for x in aa_codon_fre.values() if 'A' in x)
raKeyA = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'A': str(rareA)})]
sub_dic = {keyA: raKeyA}    # This dictionary stores the optimal codon and rarest codon for an amino acid.
# print(opt_dic)

def findopt(aa):
    optimal = max(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
    optKey = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(optimal)})]
    opt_dic[aa] = optKey
    rare = min(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
    raKey = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(rare)})]
    sub_dic[optKey] = raKey
    return opt_dic,sub_dic

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
        for i in range(0, len(prot),1):
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
opt_seq_file.seek(0,2)
opt_seq_file.write('>')

RSCU.close()
init_seq.close()
opt_seq_file.close()

#===============================================================================================
# Remove sequence fragments with many repeated bases and replaced with to rare codons.
#===============================================================================================
full_opt = open('GFP_Kp_opt_test.fa','r')
re_rep = open('GFP_Kp_opt_rm_rep.fa','w')

def findrep(base):
    remove = ''
    pos = opt_seq.find(base)
    if pos == -1:
        print('These sequences no longer have 5 or more consecutive same nucleotides {}'.format(base))
    elif pos % 3 == 0:
        if opt_seq[pos:pos+3] in sub_dic.keys():
            remove = opt_seq.replace(opt_seq[pos:pos+3],sub_dic[opt_seq[pos:pos + 3]])
    elif pos % 3 == 1:
        if opt_seq[pos-1:pos+2] in sub_dic.keys():
            remove = opt_seq.replace(opt_seq[pos-1:pos+2],sub_dic[opt_seq[pos-1:pos+2]])
    elif pos % 3 == 2:
        if opt_seq[pos+1:pos+4] in sub_dic.keys():
            remove = opt_seq.replace(opt_seq[pos+1:pos+4],sub_dic[opt_seq[pos+1:pos+4]])
    else:
        print('Something wrong')
    j = 0
    while j < len(remove):
        print(remove[j:j + 48])
        re_rep.write(remove[j:j + 48] + '\n')
        j = j + 48

opt_seq = ''
for line in full_opt:
    if line.startswith('>') and opt_seq == '':
        header = line
        print(header)
        re_rep.write(header)
    elif not line.startswith('>'):
        opt_seq = opt_seq + line.strip()
    elif line.startswith('>') and opt_seq != '':
        rep_posA = findrep('AAAAA')
        rep_posC = findrep('CCCCC')
        rep_posG = findrep('GGGGG')
        rep_posT = findrep('GGT')
        opt_seq = ''
        header = line
full_opt.close()
re_rep.close()
