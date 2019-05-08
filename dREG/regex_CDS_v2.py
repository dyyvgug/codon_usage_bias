#regex_CDS_v2.py
closer = open('/media/hp/disk2/DYY/dREG/0325out/0325peak_closer_sort.txt','r')
for line in closer.readlines():
    if 'HAVANA\tCDS' in line:
        with open('/media/hp/disk2/DYY/dREG/0325out/dREG_peak_analysis/CDS_enhancer.txt', 'a') as fi:
            fi.writelines(line)
    else:
        print('This line no CDS')
closer.close()
fi.close()