#GENOME MAPPING(ANNOTATION OF THE SEQUENCE USING THE INSECT GENOME  )
# change directory 
cd ..
cd ./results.folders
cd ./small_noncodingRNA.true
#remove whitespaces 
sed -i 's/ /_/g' GCA_001077435.1_ASM107743v1_genomic.fna
#edit the name 
cat GCA_001077435.1_ASM107743v1_genomic.fna | cut -d_ -f 1,7 |sed 's/,$//'  > GCA_001077435.1_ASM107743v1_genomic.fasta
#build index 
bowtie-build GCA_001077435.1_ASM107743v1_genomic.fasta  Glossina-morsitans # create index files 
#Aligning to the genome command 
for i in *.txt ; do echo "$i" | uniq ; done
for i in *.txt ; do   
   mapper.pl "$i" -d -e -h -l 18 -r 8 -i -j -m -q -p Glossina-morsitans -s "$i"_gmm.fa -t "$i"_gmm.arf -v    # alignment for all genome sequence      
        
                      done  

#ANNOTATION OF SMALL NON-CODING RNA IN GLOSSINA PALLIDIPES GENOME
#change directory  
cd ./miRNA_anlysis
#software to the envirnment 
module load infernal
module load easel
#create index
cmpress Rfam.cm #indexes files 
#command for small RNA annotation
for i in *_ASM107743v1_genomic.fasta ; do echo "$i" | uniq ; done
for i in *_ASM107743v1_genomic.fasta ; do
      cmscan --cpu 8 --cut_ga --rfam --nohmmonly --tblout "$i".tblout --fmt 2  --clanin Rfam.clanin Rfam.cm  "$i" > "$i".cmscan #function to annotate the genome with non coding RNA
 grep -v "=" GCA_001077435.1_ASM107743v1_genomic.fasta.tblout > Glossina-morsitans.deoverlapped.tblout  # deoverlaps the sequence(uniquereads) 
      done
#ARRANGING THE GENOME FILE INTO A BED FILE 
awk '{ print $4, $10, $11, $2, $12}'  OFS='\t' Glossina-genome.deoverlapped.tblout >Gpa.bed 
awk '{ print $4, $10, $11, $2, $12}'  OFS='\t' Glossina-morsitans.deoverlapped.tblout >Gmm.bed # building bed files 
cat Gmm.bed > Gmm.bed.2
cat Gmm.bed.2| awk '{ a=$2; b=$3; if (a<b)  {print $1, $2, $3, $4, $5} else if (a>b) {print $1, $3, $2, $4, $5} }' OFS='\t' >Gmm.bed.3 # interchanges values in column2 which is greater than in column three
awk '{ print $1}' Gpa.bed.3 | cut -d_ -f 1,7 | sed 's/,$//' > Gpa.bed.4 #editting the scaffold name 
awk '{ getline v < "Gpa.bed.4"; split( v, a ); print a[1], $2, $3 ,$4 ,$5 ,$6} ' OFS='\t' Gpa.bed.3 > Gpa.bed.5 # editting the scaffold name
sort -k4,4 -k5,5 Gmm.bed.3|groupBy -g 1,4,5 -c 4,2,3,5 -o count,min,max,distinct |awk -v OFS='\t' '{print $1, $5, $6, $2, $3, $4}' > Gmm.sorted.bed.3 # sorting using bed tools 
awk '{print $1, $2, $3, $4,".", ".", $5, $6}' OFS='\t'  Gmm.sorted.bed.3 >  Gmm.sorted.bed.5
#ARRANGING THE aLIGNMENT FILES FROM MAPPING PROCESS
for i in *_gmm.arf; do echo "$i" | uniq ; done
for i in *_gmm.arf; do
awk '{ print $6, $8, $9, $1 ,$2 ,$11}' OFS='\t'  "$i" > $i.Bed.5 #building bed files 
     done
for i in *_gmm.arf.Bed; do echo "$i" | uniq ; done
for i in *_gmm.arf.Bed; do
awk '{ print $1}'  OFS='\t' "$i" | cut -d_ -f 1,7 | sed 's/,$//' > "$i".4   
 done
#BEDTOOLS ANNOTATION 
#conda environment 
for i in *_gmm.arf.Bed.5; do echo "$i" | uniq ; done
for i in *_gmm.arf.Bed.5; do
bedtools intersect -wa -wb -s -a "$i" -b Gmm.sorted.bed.5| sort | uniq  > "$i".noncodingRNA #annotates aligned sequences                    
       done
for i in *_gmm.arf.Bed.5; do echo "$i" | uniq ; done
for i in *_gmm.arf.Bed.5; do 
bedtools intersect -v  -a "$i" -b Gmm.sorted.bed.5| sort | uniq  > "$i".unaligned #annotates unaligned sequences 
       done
#sorts the aligned files aligned 
for i in *_gmm.arf.Bed.5.noncodingRNA; do echo "$i" | uniq ; donerf.Bed.5
for i in *_gmm.arf.Bed.5.noncodingRNA; do
sort -k4,4 -k5,5 "$i"|groupBy -g 1,4,6,10,5 -c 4,2,3,5,6 -o count,min,max,distinct,distinct |awk -v OFS='\t' '{print $1, $7, $8, $2, $4, $3 ,$5 ,$6, $9}' > "$i".3 #sorts aligned sequences 
 done 
for i in *_gmm.arf.Bed.5.unaligned; do echo "$i" | uniq ; done
for i in *_gmm.arf.Bed.5.unaligned; do
 sort -k4,4 -k5,5 "$i"|groupBy -g 1,4,6,5,10 -c 4,2,3,5,6 -o count,min,max,distinct,distinct |awk -v OFS='\t' '{print $1, $6, $5, $2, $4, $3 ,$7}' > "$i".3  #sorts unaligned sequences 
 done
#Results Generated
for i in *_gmm.arf.Bed.5.noncodingRNA.3; do echo "$i" | uniq ; done
for i in *_gmm.arf.Bed.5.noncodingRNA.3; do
   awk -v OFS='\t' '{print $1, $2, $3, $4, $4, $5 ,$6 ,$7 ,$8}' "$i" | awk '{gsub("[A-Z]*[0-9]*_[0-9]*_x","",$5)}1' |awk '{gsub("_[0-9]*_x[0-9]*","",$4)}1' OFS='\t' > "$i".sorted 
      done
for i in *_gmm.arf.Bed.5.unaligned.3 ; do echo "$i" | uniq ; done
for i in *_gmm.arf.Bed.5.unaligned.3 ; do
      awk -v OFS='\t' '{print $1, $2, $3, $4, $4, $5 ,$6 ,$7}' "$i" | awk '{gsub("[A-Z]*[0-9]*_[0-9]*_x","",$5)}1' |awk '{gsub("_[0-9]*_x[0-9]*","",$4)}1' OFS='\t' > "$i".sorted
done
