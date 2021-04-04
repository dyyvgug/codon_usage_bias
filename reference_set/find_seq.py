#!/usr/bin/python
#coding:utf-8
#=========================================================================================================
# Yingying Dong.2019-10-30.Modified date:2020-1-7.Get gene sequence of high expression high translation from CDS sequence.
#=========================================================================================================
import os
import argparse

parser = argparse.ArgumentParser(description='find DNA sequence',prog='find seq',usage='%(prog)s [options]')
parser.add_argument('--samp', nargs='?',type=str,help='sample name')
args = parser.parse_args()


DNA = open('CDS_DNA.fa','r')
hE = open('{}.txt'.format(args.samp), 'r')
hE_seq = open('{}_seq.fa'.format(args.samp),'w')

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