#tran_lower.py
txt = open("FUS.fa","r").read()
txt = txt.upper()

with open('FUS.fa', 'w') as ff:
    ff.writelines(txt)

ff.close()