#!/bin/bash

#download from NCBI
wget -o ./in_files/assembly_summary_bacteria.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
#cp assembly_summary.txt ./in_files/assembly_summary_bacteria.txt


#generate wget links for each genome
awk -F "\t" '$12=="Complete Genome"  && $8 ~ /Pseudomonas/ {print $20}' ./assembly_summary.txt > ./ftp_folder_bacteria.txt
awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder_bacteria.txt > faa_files_bacteria.sh


#download faa files from NCBI
echo downloading...
bash faa_files_bacteria.sh >download_log
#mv *.gz ./in_files
gzip -d *.gz
echo merging indexed faa...
#add genome index to each faa
awk '/>/ {$0=$0" "FILENAME}1' *protein.faa > indexed_merged_protein.faa
#rm -r GCA*.faa
#rm assembly*
echo Finished!!!
