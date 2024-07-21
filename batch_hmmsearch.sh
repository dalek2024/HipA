#!/bin/bash
while getopts q:o:d: flag
do
    case "${flag}" in
        q) hmmfile_dir=${OPTARG};; ### the directory to hmm files
        o) out_dir=${OPTARG};; ###the directory to save outputs
        d) database=${OPTARG};; ###database should be the indexed_merged_protein.faa
    esac
done

function linebreaks_removal {
        awk '!/^#/' $1 | tr -d "\n" | awk '{gsub(".faa", "\n", $0)}1'
        echo  $0
}

#index database
esl-sfetch --index $database

FILES=$hmmfile_dir/*.hmm
for f in $FILES
do
        file=$(basename $f)
        echo processing $file
        tbl_file=$out_dir"/tbl_"$file".out"
        std_output=$out_dir"/"$file"_e0.001.out"
        hmmsearch --cpu 40 -E 0.001 --textw 300 --tblout $tbl_file $f $database > $std_output
        linebreaks_removal $tbl_file > $out_dir"/tbl_"$file"_processed.out"
        grep -v '^#' $tbl_file | awk '{print $1}' > $out_dir"/"$file"_hit_ids.txt"
        esl-sfetch -f $database $out_dir"/"$file"_hit_ids.txt" > $out_dir"/"$file"_hits.fasta"
done

rm $out_dir/*_e0.001.out

# an exmaple to execute 'sh batch_hmmsearch.sh -q ./hmm_files -o ./hmm_out/hmm_raw -d ./in_files/indexed_merged_protein.faa'
