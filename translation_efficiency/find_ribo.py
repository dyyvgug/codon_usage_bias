#!/usr/bin/python
#coding:utf-8
# Yingying Dong.
import os
import re

os.chdir("/home/hp/Desktop/other_riboseq/mouse_mm10/experiment1/aligned/ribo_num/")
other_ribo = open("SRR8728404_other_ribo_gene.txt",'r')
hE_hT = open("SRR8728404_hiE_ht_gene.txt",'r')
other_ribo_filter = open("other_ribo_filter.txt",'w')
hE_hT_filter = open("hE_ht__ribo_filter.txt",'w')
hE_hT_no_ribo = open("hE_ht-ribo.txt",'w')
keyword = re.compile(r'rpl|rps',re.IGNORECASE)
def findRibo(ribokey):
	for line in other_ribo:
		riboLine = ribokey.search(line)
		if riboLine:
			other_ribo_filter.write(line)
	for line in hE_hT:
		riboLine = ribokey.search(line)
		if riboLine:
			hE_hT_filter.write(line)
		else:
			hE_hT_no_ribo.write(line)
#findRibo(keyword1)
findRibo(keyword)
other_ribo_filter.close()
hE_hT_filter.close()
hE_hT_no_ribo.close()
