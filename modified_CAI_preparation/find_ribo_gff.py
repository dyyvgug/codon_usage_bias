#!/usr/bin/python
# coding:utf-8
# =========================================================================================================
#  Yingying Dong.2019-11-8.Find ribosome genes from gff.
# =========================================================================================================
import os
import re
import sys
import argparse
import rpy2.robjects as robjects

parser = argparse.ArgumentParser(description='find ribosome genes from gff.', prog='find_ribo',
                                 usage='%(prog)s [options]')
parser.add_argument('--spe', nargs='?', type=str, help='species name')
parser.add_argument('--ann', nargs='*', type=str, default='ref.gff', help='input gtf format file')
args = parser.parse_args()

os.chdir('/media/hp/disk1/DYY/reference/annotation/{}/'.format(args.spe))
gff_path = '/media/hp/disk1/DYY/reference/annotation/{}/'.format(args.spe)
gff = open('{}{}'.format(gff_path, args.ann), 'r')
id_file = open('ribo_geneName.txt', 'w')


gff_table = []
id_table = []
keyword = re.compile(r'ribosomal', re.IGNORECASE)
for i in gff:
    gff_table.append(i.strip().split('\t'))

for k in gff_table:
    riboLine = keyword.search(k[9])
    if riboLine:
        print(k[8] +'\t' + k[9] +'\n')
        id_file.write(k[8] +'\t' + k[9] +'\n')

gff.close()
id_file.close()

r_script = '''
gff = read.table('ribo_geneName.txt',sep = '\t',header = F,quote = '')
gff = unique(gff)
write.table(gff$V1,file = 'ribo_geneName_u.txt', quote = FALSE,row.names = F, col.names = F)
'''
robjects.r(r_script)
