step 1. download SRA files.

step 2. quality control. 

step 3. remove adaptor.

step 4. use bowtie2 to map.
  Firstly, create a index. For instance,
  
    $ bowtie2-build eschColi_K_12_MG1655-mature-tRNAs.fa tRNA
  

step 5. to bam.

step 6. to bed.

step 7. R script: tRNAseq_anticodon_condon.R
