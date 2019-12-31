# ==============================================================================================================
# Yingying Dong.2019-12-30. Remove sequence fragments with many repeated bases and restriction site ,
#  then replaced with to rare codons.
# ==============================================================================================================
full_opt = open('bac_Kp_opt.fa', 'r')
re_rep = open('bac_rm_rep.fa', 'w')


def findrep(base,opt_seq):
    pos = 0
    while pos >= -1:
        pos = opt_seq.find(base)
        if pos == -1:
            print('These sequences no longer have this restriction site or consecutive same nucleotides {}'.format(base))
            break
        else:
            if pos % 3 == 0:
                if opt_seq[pos:pos+3] in sub_dic.keys():
                    print(opt_seq[pos:pos + 3])
                    opt_seq = opt_seq.replace(opt_seq[pos:pos + 3], sub_dic[opt_seq[pos:pos+3]])
                    pos = opt_seq.find(base)
                    continue
            elif pos % 3 == 1:
                if opt_seq[pos-1:pos+2] in sub_dic.keys():
                    print(opt_seq[pos-1:pos+2])
                    opt_seq = opt_seq.replace(opt_seq[pos-1:pos+2], sub_dic[opt_seq[pos-1:pos+2]])
                    pos = opt_seq.find(base)
                    continue
            elif pos % 3 == 2:
                if opt_seq[pos+1:pos+4] in sub_dic.keys():
                    print(opt_seq[pos+1:pos+4])
                    opt_seq = opt_seq.replace(opt_seq[pos+1:pos+4], sub_dic[opt_seq[pos+1:pos+4]])
                    pos = opt_seq.find(base)
                    continue
            else:
                print('Something wrong')

    return opt_seq


opt_seq = ''
for line in full_opt:
    if line.startswith('>') and opt_seq == '':
        header = line
        print(header)
        re_rep.write(header)
    elif not line.startswith('>'):
        opt_seq = opt_seq + line.strip()
    elif line.startswith('>') and opt_seq != '':
        rep_posA = findrep('AAAAA',opt_seq)
        rep_posC = findrep('CCCCC',opt_seq)
        rep_posG = findrep('GGGGG',opt_seq)
        rep_posT = findrep('TTTTT',opt_seq)
        rep_EcoRI = findrep('GAATTC',opt_seq)
        rep_SalI = findrep('GTCGAC',opt_seq)
        rep_NdeI = findrep('GATATG',opt_seq)
        rep_XhoI = findrep('CTCGAG',opt_seq)

        j = 0
        while j < len(opt_seq):
            print(opt_seq[j:j + 48])
            re_rep.write(opt_seq[j:j + 48] + '\n')
            j = j + 48
        opt_seq = ''
        header = line
full_opt.close()
re_rep.close()
