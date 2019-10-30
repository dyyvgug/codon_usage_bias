#!/usr/bin/python
#coding:utf-8
#==========================================================================================================
# Yingying Dong.2019-10-25.Modified date:2019-10-29.
# This program calculates CAI value for given FASTA sequence,
#  weight from the reference gene set that high RNAseq TPM and high Ribo-seq TPM.
#==========================================================================================================
import sys,os
import argparse
from scipy import stats

parser = argparse.ArgumentParser(description='calculates CAI',prog = 'cal_CAI',usage = '%(prog)s [options]') 
parser.add_argument('--spe', nargs='?',type=str,help='species name' )
parser.add_argument('--spA', nargs='?',type=str,help='species name abbreviation and experiment number' )
parser.add_argument('--exp',nargs='?', type=int,help='experiment number')
parser.add_argument('--inp', nargs='?',type=str, help='expect DNA sequence of FASTA files that calculate modified CAI values') 
args = parser.parse_args()

os.chdir('/media/hp/disk1/DYY/reference/annotation/{}/ref'.format(args.spe))
weight_path = '/home/hp/Desktop/other_riboseq/{}/experiment{}/aligned/ribo_num/'.format(args.spe,args.exp)
fa_file = open(args.inp, 'r')
weight_file = open("{}hE_hT_RSCU_weight.txt".format(weight_path), 'r')
CAI_file = open('{}_mCAI.txt'.format(args.spA), 'w')

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
		print(header)
		#print(weight_list)
		CAI = stats.gmean(weight_list)
		CAI_file.write('{}\t{}\n'.format(header, CAI))
		header = line.strip()
		dna = ''
		weight_list = []

fa_file.close()
CAI_file.close()
weight_file.close()

os.system('sed -i \'s/>//g\' {}_mCAI.txt'.format(args.spA))
