#!/usr/bin/python
# coding:utf-8
# Yingying Dong.2019-10-25.Modified date:2019-10-28.
# This program calculates CAI,weights from reference gene set that high RNAseq TPM and high Ribo-seq TPM.
import sys
from scipy import stats

weight_path = '/home/hp/Desktop/other_riboseq/C_elegans_Ensl_WBcel235/experiment2/aligned/ribo_num/'
fa_file = open("CDS.fa", 'r')
weight_file = open("{}hE_hT_RSCU_weight.txt".format(weight_path), 'r')
CAI_file = open('Ce_mCAI.txt', 'w')

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
		#print(weight_list)
		print(header)
		CAI = stats.gmean(weight_list)
		CAI_file.write('{}\t{}\n'.format(header, CAI))
		header = line.strip()

fa_file.close()
CAI_file.close()
weight_file.close()
