#!/usr/bin/python

fafile = open("CDS_DNA.fa","r")
fafilter = open("filter.fa","w")

i=0
for line in fafile:
        #print i,i%4,line
        if i%2==0:
                seqID=line.strip("\n")
        elif i%2==1:
                sequence=line.strip("\n")                
                if len(sequence)>=300 :
                        fafilter.write(seqID+"\n"+sequence+"\n")
                        print (seqID)
        i+=1

fafile.close()
fafilter.close()

