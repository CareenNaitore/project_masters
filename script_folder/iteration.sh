#WQWQtriming
i=Adfem2_1.fastq
trimmomatic SE -threads 36 \
"$i" ./Trimmomatic/"$i" \
ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
               LEADING:3 TRAILING:0 SLIDINGWINDOW:4:20 MINLEN:17 CROP:45                      
cd ./Trimmomatic
##Remove errant sequence
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (seq != "TGCTTGGACTACATATGGTTGAGGGTTGTA") {print header, seq, qheader, qseq}}' "$i" > "$i".filt
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (seq != "CTGCTTGGACTACATATGGTTGAGGGTTGTA") {print header, seq, qheader, qseq}}' Adfem2_1.fastq.filt > "$i".filt2
#quality check 
mkdir Qc                           
fastqc --o Qc  --noextract -t 36 "$i"            
fastqc --o Qc  --noextract -t 36 "$i".filt2
