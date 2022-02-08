# MICRORNA TARGET PREDICTION 
# !/bin/sh
# identifying the target of the microRNA using the software miranda
cd ../results.folders/targetscan_analysis/
#awk '/^>/ { getline seq } seq != "Sequence unavailable" { print $0 "\n" seq}' Gpa_utr.fa > GPA_UTR.fa #filtering out unwanted sequences/unavailable  

### softwre loaded on the conda environment 
#cd ../results.folders/targetscan_analysis
#for i in GPA_UTR.fa ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in GPA_UTR.fa ;
#do  
#a=$(basename "$i"  _UTR.fa) #identifies the index for for each library
#miranda  miRNA.target2  "$i" -sc 150 -en 20 -strict  -keyval -out "$a".targetscan4 # command for target prediction 
 #done
## editing results  
for i in *.targetscan4 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.targetscan4 ;
do
a=$(basename "$i"  _UTR.fa) #identifies the index for for each library
grep -A 2 "Score for this Scan:" "$i" | sort | grep '>>' >"$a".result4 #extracting information 
  done
## identifying the target of the microRNA using the software RNAhybrid
# identifying the target of the microRNA using the software RNAhybrid 
module load rnahybrid
for i in GPA_UTR.fa ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in GPA_UTR.fa ;
do  
a=$(basename "$i"  _UTR.fa) #identifies the index for for each library
RNAhybrid -s 3utr_fly -e −20 -p 0.05 -q  miRNA.target2  -t "$i" -c >"$a".rnahybridfile2 # command for target prediction 
 done
