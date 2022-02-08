#checking for quality of  fastq file using fastqc 
cd ..
cd ./results.folders # changing in a directory for results 
cd ./high_quality_sequences

mkdir fastqc_analysis
for i in /home/cnaitore/mirnaexpresionprofile/"May 2018 Careen miRNA sequencing results"/input.folders/Raw_sequence/*.fastq; do echo "$i" | uniq ; done # identifies files you have 
for i in /home/cnaitore/mirnaexpresionprofile/"May 2018 Careen miRNA sequencing results"/input.folders/Raw_sequence/*fastq; do
    fastqc --o fastqc_analysis --noextract "$i" #command run for fastqc file 
 done 

#Trimming the adapters 
cd ..
cd ..
cd ./input.folders # change directory 
cd ./Raw_sequence  
#no need to load the software on the cluster 
for i in *.fastq ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fastq ;
do  
a=$(basename "$i" .fastq) #identifies the index for for each library
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC  -a "$a".index -a CTCGTATGCCGTCTTCTGCTTG -e 0.25 -m 15 --match-read-wildcards  --untrimmed-output "$i"_FILE --too-short-output "$i"_tooshort -o "$i"_trim "$i"  #removes adapters 
 done 
#for i in *.fastq ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fastq ; do 
#trimmomatic SE -threads 36 \
#"$i" ./Trimmomatics/"$i".trim \
#ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
 #      LEADING:30 TRAILING:30 SLIDINGWINDOW:4:25 MINLEN:17 CROP:35
#done
#MOVE THE FILE IN THE OUTPUT FOLDERS 
for i in *_FILE ; do 
 mv "$i" /home/cnaitore/mirnaexpresionprofile/"May 2018 Careen miRNA sequencing results"/results.folders/high_quality_sequences # MOVE FILES  
 done 
#MOVE THE FILE IN THE$ OUTPUT FOLDERS 
for i in *_trim ; do
 mv "$i" /home/cnaitore/mirnaexpresionprofile/"May 2018 Careen miRNA sequencing results"/results.folders/high_quality_sequences # MOVES FILES 
 done 
#MOVE THE FILE IN THE OUTPUT FOLDERS 
for i in *_tooshort ; do
 mv "$i" /home/cnaitore/mirnaexpresionprofile/"May 2018 Careen miRNA sequencing results"/results.folders/high_quality_sequences # MOVES FILES 
 done 

# filtering poor quality reads
cd ..
cd ..
cd ./results.folders
cd ./high_quality_sequences 
module load fastx
for i in *_trim ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *_trim ; do 
fastq_quality_filter -q 20 -p 80 -i "$i" -o "$i"_30 -Q33 -v 

done 
cd ./Trimmomatics
#removes a overexpressed  sequence 
for i in *.trim ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.trim ;
do
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (seq != "TGCTTGGACTACATATGGTTGAGGGTTGTA") {print header, seq, qheader, qseq}}' "$i" > "$i".filt
 done 
for i in *.filt ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.filt ;
do
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (seq != "CTGCTTGGACTACATATGGTTGAGGGTTGTA") {print header, seq, qheader, qseq}}' "$i" > "$i".filt2
done
#trim sequence length
module load fastx
for i in *.filt2; do echo "$i" | uniq ; done
for i in *.filt2; do
 fastx_trimmer -i "$i" -o "$i"_trim10 -f 1 -l 35  -m 18 -Q 33 -v #trimmed the sequence to length 35
  done 

#check the length distribution of sequences
for i in *_trim10 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *_trim10 ;
do
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' "$i" > "$i".length  # to determine what length remains after trimming 

done

#checking for quality of  fastq file using fastqc 
for i in *.trim ; do echo "$i" | uniq ; done # identifies files you have 
for i in *.trim; do
  fastqc --o fastqc_analysis --noextract "$i" #command run for fastqc file 
 done
#to splits big jobs and use more processing speed
cat Adfem6.filtered.fastq |parallel --gnu --pipe --block 60000000 'grep -wEf outfile -A2 -B1 > parallel_test_{#}'cat Adfem6.filtered.fastq |parallel --gnu --pipe --block 60000000 'grep -wEf outfile -A2 -B1 > parallel_test_{#}'
