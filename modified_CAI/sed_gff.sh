#!/bin/bash 
#get gene product

sed -i 's/ID.*locus_tag=//g' ref.gff
sed -i 's/;orig_transcript_id.*product=/\t/g' ref.gff
sed -i '/partial=true/d' ref.gff
sed -i 's/;.*//g' ref.gff
sed -i '/region/d' ref.gff
sed -i '/#/d' ref.gff
