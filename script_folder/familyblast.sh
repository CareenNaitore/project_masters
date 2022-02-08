cd .././results.folders
cd  ./filtered_sequences
cd  ./miRNA_filtered
module load ncbi-blast/2.4.0+ 
sed 1d "miRNA_family2"|sort -u|awk '{print ">"$1"-"$3"\n"$2}' > miRNAfamily2.fa 

#created database using availables sequences 
makeblastdb -in miRNAfamily2.fa -out miRNAfamily2 -parse_seqids -dbtype nucl

for i in allmiRNA.fa ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in allmiRNA.fa ;
do
  blastn -query "$i" -db miRNAfamily2 -out all.miRNAblast2 -num_threads 16  -task blastn-short  -qcov_hsp_perc 80  -max_target_seqs 3  -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" # blast against nucleotide database
done

