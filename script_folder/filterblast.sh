cd .././results.folders
cd  ./filtered_sequences
cd  ./miRNA_filtered
#created database using availables sequences 
makeblastdb -in mydb.fsa -parse_seqids -dbtype nucl
#convert fastq to fasta 
for i in *_30 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *_30 ;
do
fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
done
#collapse to unique reads 
for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fa2 ;
do
collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
done
#blastn command 
for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
  blastn -query "$i" -db miRNA -out "$i".miRna.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" # blast against nucleotide database
done
#count function 
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
   cat "$i" | wc -l > $i.miRNAcount #count number for rRNA 
done
#selecting first column 
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cut -f1 "$i" > "$i".filtered
done
#filtering sequnces from fasta format 
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep -A1 -f "$a".fasta.miRna.blast.outfmt6.filtered "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".miRNA
done
#RIBOSOMAL RNA FILTERING 
#change folder
cd .././rRNA_filtered
for i in *_30 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *_30 ;
do
fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
done
# collapse the  fasta file to unique reads 
for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fa2 ;
do
collapse_reads_md.pl "$i" seq > "$i".fasta #collapses the reads to unique
done
#blast function
for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
 blastn -query "$i" -db rRNA -out "$i".rRNA.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" #blast command for rRNA
done
# count function
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do 
 cat "$i" | wc -l > $i.rRNAcount #count number for rRNA 
done 
# select first column

for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cut -f1 "$i" > "$i".rRNA.filtered
done
# filter sequence from fasta filr
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep -A1 -f "$a".fasta.blast.outfmt6.rRNA.filtered "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".others #filter sequences that are others 
done
#TRANSFER RNA FILTERING 
#changing folder 
cd .././tRNA_filtered
for i in *_30 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *_30 ;
do
fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
done
#collapse into unique reads 
for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fa2 ;
do
collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
done
#blast function
for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
 blastn -query "$i" -db tRNA -out "$i".tRNA.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" # count number for tRNA
done
#count function
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cut -f1 "$i" > "$i".tRNA.filtered
done
#select first column 
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
 cat  "$i"| wc -l  > $i.tRNAcount #count number for rRNA 
done 
#filter sequences 
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep -A1 -f "$a".fasta.blast.outfmt6.tRNA.filtered "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".tRNA
done
#OTHERTYPES OF small NON-CODING RNA IN RFAM
#cd  directory    
cd .././others_filtered
#converting fastq to fastafile
for i in *_30 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *_30 ;
do
fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
done
#collapsing to unique read
for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fa2 ;
do
collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
done
#blast function
for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
  blastn -query "$i" -db othersfiltered -out "$i".blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" #blast command to identify 
done
#count function 
for i in *.outfmt6; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
  grep -c "$i" > "$i".otherscount #count number for others
done 
#select first column
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6;
do
cut -f1 "$i" > "$i".others.filtered # identifies the file that the command 
done
#filtering from fasta file 
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep -A1 -f "$a".fasta.blast.outfmt6.others.filtered "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".others  
done
#for i in *.others ;
#do
#a=$(basename "$i" .others)
 #echo "$i" 
 #; done # identifies the file that the command is work on 
#for i in `cat *.others|sort -u`;
#do awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (seq != "$i") {print header, seq, qheader, qseq}}' "$a".fastq   > "$a".others.filt;  done
# awk '{ print "\""$0"\""}' read62.others |awk 'BEGIN{FS=OFS="\""} {gsub(/[[:space:]]/,"",$2)} 1' # adding double quotation and removing space                          

