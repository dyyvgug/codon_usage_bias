#!/bin/bash
#--------------1----------------- Ec
sed -i 's/ geneID ".*";//g' ref.gtf
sed -i 's/"; gene_id \".*\"; /\t/g' ref.gtf
sed -i 's/"; /\t/g' ref.gtf
sed -i 's/transcript_id "//g' ref.gtf
sed -i 's/gene_name "//g' ref.gtf
sed -i 's/";//g' ref.gtf

#-------------2----------------- Sc Dm includeCe etc.  Most species.
sed -i 's/transcript_id \".*\"; gene_id//g' ref.gtf
sed -i 's/"; gene_name "/\t/g' ref.gtf
sed -i 's/";//g' ref.gtf
sed -i 's/"//g' ref.gtf
sed -i 's/ //g' ref.gtf

#-------------3----------------- Ce
sed -i '/transcript_id \"rna-NC.*/d' ref.gtf
sed -i 's/transcript_id \".*\"; gene_id//g' ref.gtf
sed -i 's/"; gene_name "/\t/g' ref.gtf
sed -i 's/";//g' ref.gtf
sed -i 's/"//g' ref.gtf
sed -i 's/ //g' ref.gtf

#-------------4------------------ Mm
sed -i 's/\"; gene_type \".*\"; gene_name \"/\t/g' ref.gtf
sed -i '/exon/d' ref.gtf
sed -i '/transcript/d' ref.gtf
sed -i 's/\"; level.*//g' ref.gtf
sed -i 's/gene_id "//g' ref.gtf
