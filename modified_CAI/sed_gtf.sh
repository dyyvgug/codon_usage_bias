#!/bin/bash
sed -i 's/ geneID ".*";//g' ref.gtf
sed -i 's/"; gene_id \".*\"; /\t/g' ref.gtf
sed -i 's/"; /\t/g' ref.gtf
sed -i 's/transcript_id "//g' ref.gtf
sed -i 's/gene_name "//g' ref.gtf
sed -i 's/";//g' ref.gtf

#-------------2-----------------
sed -i 's/transcript_id \".*\"; gene_id//g' ref.gtf
sed -i 's/"; gene_name "/\t/g' ref.gtf
sed -i 's/";//g' ref.gtf
sed -i 's/"//g' ref.gtf
sed -i 's/ //g' ref.gtf
