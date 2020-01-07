#!/usr/bin/python
#coding:utf-8
#=========================================================================================================
# Yingying Dong.2019-10-30.Get gene sequence of high expression high translation from CDS sequence.
#=========================================================================================================
import os
import argparse

parser = argparse.ArgumentParser(description='find DNA sequence',prog='find seq',usage='%(prog)s [options]')
parser.add_argument('--spe', nargs='?',type=str,help='species name')
args = parser.parse_args()


os.chdir('/media/hp/disk1/DYY/reference/annotation/{}/ref'.format(args.spe))
DNA = open('CDS_DNA.fa','r')
hE_path = '/media/hp/Katniss/DYY/aligned/Saccharomyces_cerevisiae/experiment3/'
hE = open('{}SRR5422018_low30_exp_only_name.txt'.format(hE_path),'r')
hE_seq = open('{}SRR5422018_low30_seq.fa'.format(hE_path),'w')

hE_table = []
database = {}
for i in hE.readlines():
    hE_table.append(i.strip())
for j in DNA:
    if j.startswith('>'):
        keys = j.lstrip('>').strip()
        database[keys] = []
    else:
        database[keys].append(j.strip())
for line in hE_table:
    if line in database.keys():
        hE_seq.write('>' + str(line) + '\n' + str(''.join(database[line])) + '\n')

hE.close()
DNA.close()
hE_seq.close()
