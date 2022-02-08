cd ../results.folders
cd ./miRNA_novelquality_true

module load mirdeep2
#extract_miRNAs.pl mature.fasta aga > mature_ref3.fasta
extract_miRNAs.pl mature.fasta afu,ssr,api,bdo,bib,tca,bmo,cqu,dan,der,dyr,dmo,dpe,dps,dqu,dse,dsi,dvi,dwi,dya,hme,imi,mse,ngi,nlo,pca,pxy,cel,aae,aga,hsa,mmu > mature_ref4.fasta # extract mature miRNA 
extract_miRNAs.pl hairpin.fasta afu,ssr,api,bdo,bib,tca,bmo,cqu,dan,der,dyr,dmo,dpe,dps,dqu,dse,dsi,dvi,dwi,dya,hme,imi,mse,ngi,nlo,pca,pxy,cel,aae,dme > hairpin_ref2.fasta # extract mature miRNA 
#extract_miRNAs.pl hairpin.fasta aga > hairpin3.fasta # extract hairpin miRNA 

#create txt files using the replicates for exaple txt file contain (a1,a2,a3) #

#Annotation with the GENOME 
module load mirdeep2
bowtie-build Glossina-pallidipes-IAEA_SCAFFOLDS_GpalI1.fa   Glossina-pallidipes 

for i in *.txt ; do echo "$i" | uniq ; done

for i in *.txt ; do   
    
mapper.pl "$i" -d -e -h -l 18 -r 8 -i -j -m -q -p Glossina-Pallidipes -s "$i"_gpa.fa -t "$i"_gpa.arf -v          
        
                    #done  


# IDENTIFICATION OF MIRNA USING MIRDEEP2 
for i in *.fa ; do echo "$i" |uniq ; done 
for i in *.fa ;
do 
C=$(basename "$i" .fa)
miRDeep2.pl "$i" GCA_000688715.1_Glossina_pallidipes-1.0.3_genomic.fasta  "$C".arf  novel.fasta  mature_ref4.fasta novel_hairpin.fasta -r "$i".2  > adfem2.log 
done

#QUANTIFICATION OF MIRNA IN SAMPLE USING MIRDEEP2 
for i in *.fa ; do echo "$i" |uniq ; done
for i in *.fa ;
do
quantifier.pl -p novel_hairpin.fasta -m novel.fasta -r "$i"  -y"$i" -k -j
 
 done 
