#!/usr/bin/python
#coding:utf-8
#=========================================================================================================
# Yingying Dong.2019-10-30.Get gene sequence of high expression high translation from CDS sequence.
#=========================================================================================================
import os
import argparse

parser = argparse.ArgumentParser(description='find DNA sequence',prog='find seq',usage='%(prog)s [options]')
parser.add_argument('--spe', nargs='?',type=str,help='species name')
parser.add_argument('--exp',nargs='?', type=int,help='experiment number')
args = parser.parse_args()


os.chdir('/media/hp/disk1/DYY/reference/annotation/{}/ref'.format(args.spe))
DNA = open('CDS_DNA.fa','r')
hE_hT_path = '/home/hp/Desktop/other_riboseq/{}/experiment{}/aligned/ribo_num/'.format(args.spe,args.exp)
hE_hT = open('{}hE_hT_only_geneID.txt'.format(hE_hT_path),'r')
hE_hT_seq = open('{}hE_hT_seq.fa'.format(hE_hT_path),'w')

hE_hT_table = []
database = {}
for i in hE_hT.readlines():
    hE_hT_table.append(i.strip())
for j in DNA:
    if j.startswith('>'):
        keys = j.lstrip('>').strip()
        database[keys] = []
    else:
        database[keys].append(j.strip())
for line in hE_hT_table:
    hE_hT_seq.write('>' + str(line) + '\n' + str(''.join(database[line])) + '\n')

hE_hT.close()
DNA.close()
hE_hT_seq.close()