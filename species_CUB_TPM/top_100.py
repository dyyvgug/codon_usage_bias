import os
#提取TPM最大的200对应的id，生成新表
TPM = open("TPM.txt","r").read()
a = TPM.replace('\t',',').replace('\n',',')
ls = a.split(',')
rna = ls[::2]
tpm = ls[1::2]
del rna[0]
del tpm[0]
tpm = [ float(x) for x in tpm ]
dic = dict(zip(rna,tpm))
items = list(dic.items())
#print(items)

items.sort(key=lambda x:x[1], reverse=True)
with open('top100.txt', 'w') as ff:
    for i in range(100):
        rna_, tpm_ = items[i]
        #print ("{0}\t{1}".format(rna_, tpm_))
        ff.write("\n{0}\t{1}".format(rna_, tpm_))
        #or ff.write('\n'.join("{0}\t{1}".format(rna_, tpm_)))
ff.close()
# 提取TPM的TOP200的ref信息，CDS区，正反链信息等。
fo = open("top100.txt")
txt = fo.read()
a1 = txt.replace('\t',',').replace('\n',',')
ls1 = a1.split(",")
del ls1[0]
print(ls1)
rna1 = ls1[::2]
tpm1 = ls1[1::2]
#dic1 = dict(zip(rna1,tpm1))
#print(dic)
fp = open('ref.txt','r')
with open('top100_out.txt', 'w') as fg:
    for line in fp:
        string = line.split('\t')
        if string[0] in rna1:
            print(string[0])
            fg.writelines(line)
            #ff.write('\t'.join([string[0], string[1],string[2],string[3]])+)
fo.close()
fp.close()
fg.close()
