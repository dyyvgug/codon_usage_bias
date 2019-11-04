#!/usr/bin/python
#coding:utf-8
#=========================================================================================================
# Yingying Dong.2019-11-4.Get gene ID from the gene name.At the same time,extract their DNA sequence.
#=========================================================================================================
import os
import argparse

parser = argparse.ArgumentParser(description='gene name to gene ID',prog = 'nameToID',usage = '%(prog)s [options]')
parser.add_argument('--spe', nargs='?',type=str,help='species name' )
parser.add_argument('--spA', nargs='?',type=str,help='species name abbreviation' )
parser.add_argument('--inp', nargs='*',type=str,default='ref.gtf',help='input gtf format file')
args = parser.parse_args()

os.chdir('/media/hp/disk2/GO_KEGG/GO/{}/CC'.format(args.spA))
gtf_path = '/media/hp/disk1/DYY/reference/annotation/{}/'.format(args.spe)
name_file = open('CC_only_symbol.txt','r')
gtf = open('{}{}'.format(gtf_path,args.inp),'r')
id_file = open('CC_geneID','w')

DNA_path = '/media/hp/disk1/DYY/reference/annotation/{}/ref/'.format(args.spe)
DNA = open('{}CDS_DNA.fa'.format(DNA_path),'r')
ex_seq = open('{}CC_seq.fa'.format(DNA_path),'w')

name_table = []
name_id = {}
gtf_table = []
id_table = []
for i in name_file.readlines():
    name_table.append(i.strip())
for j in gtf:
    gtf_table.append(j.strip().split('\t'))
for k in gtf_table:
    #print(k[9])
    name_id[k[9]] = k[8]
for l in name_id.keys():
    if l in name_table:
        #print(name_id[l])
        id_table.append(name_id[l])
        id_file.write(name_id[l] + '\n')

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

name_file.close()
gtf.close()
id_file.close()
ex_seq.close()