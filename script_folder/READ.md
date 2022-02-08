# MiRNA TRANSCRIPTOME ANALYSIS IN TSESTE  USING BIOINFORMATIC SOFTWARE
## introduction
The pipeline involves 

1.quality assements

2.genome annotation

3.filtering out other types of small non coding RNA 

4.miRNA analysis(miRNA identification and miRNA differential expression)

5.Target prediction and annotation

## Quality assesment
For *quality assesment*: the script **adapter_removal.sh** is suited for this step, it is a bash script i designed to help in improving the quality of my sequence, this scripts involves various procedure, 
               
                1. quality check (pre-triming quality check using fastqc software)
                
                2. triming the adapters(using cutadapt and fastx toolkitsoftware)
                
                3. removal of contaminants(using a command awk command)
               
                4. Quality check (post-trimming quality check using fastqc softwre)
                
                5. Determing the length of my sequences using (awk command)
### Length distribution
This involves calculation of the length distribution in each sequnce library. The scripts used are in bash scripting and r scripts format: **lengthdistribution.sh* and *lengthdistribution_R2.R**.After determining the length of the sequences using the awk command.The *lengthdistribution.sh* reshapes the data for analysis in r programming to generates distribution in length. 

## Genome annotation 
*Genome annotation* involves aligning to the reference genome. The script used in the analysis is  **noncodingRNA.sh**, also the script combine both perl and bash scripting involves:
1. Building indexes for the genome ( Most aligners require one to build indexes for better alignment using bowtie softwares) 

2. alignment of the sequences to the build indexes(using the bowtie software and mapper.pl command)

## Filtering out unwanted sequences(other types of small non coding RNA)
*Filtering* process involves removing other types of sRNA, The libraries contains sequences (length 18-35bps), currently we are interested in identifying miRNA.The procedure involves a database known as [Rfam](https://academic.oup.com/nar/article/46/D1/D335/4588106) and their availble software [infernal](http://eddylab.org/infernal/). the script used is  **noncodingRNA.sh** currentily has about six procedures:
1.aligning to the genome(mapper.pl script from mirdeep2)

2.Creating indexes for alignment in Rfam(infernal command)

3.annotating the genome with small noncoding genes(creating a GTF file)(uing the command for infernal and the database available)

4.generated bedfiles(using bedtools)

5.Generated overlaping features

6.uing linux command filtered in **filterbash.sh** out the unwanted sequence and then used the once (matched miRNA and unmatched sequences)

## MiRNA analysis
### miRNA identification and characterisation
miRNA identifcation involves annotating miRNA genes on the sequences of the libraries. It involves the software [miRdeep2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245920/) and the database used for the analysis is [miRBase](http://www.mirbase.org/)(v 22). The bash script developed includes **miRDeep2.sh**:
1. Developing index file [bowtie software](http://bowtie-bio.sourceforge.net/index.shtml) 
2. aligning sequence [mapper.pl](https://github.com/rajewsky-lab/mirdeep2)
3. miRNA identification analysis [miRdeep2.pl](https://github.com/rajewsky-lab/mirdeep2)
4. quantifies the number of reads [quantifier.pl](https://github.com/rajewsky-lab/mirdeep2)

#### miRNA differential expression
miRNA expression profiling is analyzed by two different software [Persus software](http://www.perseus-framework.org) and [R programming](https://cran.r-project.org/bin/windows/base/). The expression profiles were analyzed by first examining the correlation between the replicates, if there was a high correlation between the replicates stastical analysis was done using the [edgeR package](https://bioconductor.org/packages/release/bioc/html/edgeR.html) applying negative binomial test comparing between consecutive life stages and sexes.The expression profiles was determined by the false discovery rate(adjusted p value) and the log fold change between the replicates.

**codes for differential expression and plots** 

1.[Persus software](http://www.perseus-framework.org), a machine-learning based module  integrated with a collection of statistical tools  used for protemic down stream analysis, but can also be used in the study of other omics data analyis for covering normalization, pattern recogntion, time series analyis, cross omics comparison and multiple hypothesis testing   was used for analyis of variation within the replicates and also correlation with the miRNA relative abundance using pearson correlation coefficient module.
The software has a powerful visulization platform for the expression profile of miRNA, i was able to generate a heat map graph using the software  

2.**heatmap.ipynb**:this script is in python and is used to generate the heat map and the expression profile from data, it works with seaborn package in matlibplot for analysis of my sequence the scripts contains codes in jupyter lab and comments for each code hence it is easy to flow. 

3. **volcano_plots.ipynb**: the
script is in python and is used to generate volcanopots for differential expression visualization in consecutive life stages currently the scripts intergrates with python packages packages  including:
numpy
pandas 
matplotlab 


## Target prediction and functional annotation
miRNA target prediction is analyzed by two different software [miRanda software](http://www.mirtoolsgallery.org/miRToolsGallery/node/1055) and [RNAhybrid](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1538877/). The miRNA genes from the consecutive life stages and sexes and the 3 untranslated mRNA regions from [vectorbase](https://www.vectorbase.org) was used to determine the function of the miRNA in our species. We applied the **targetprediction.sh** script: it has four steps:
 1. Determining the target using miRanda

 2. Determining the target using RNAhybrid

 3. Extracting information from boith results of miRanda and RNA hybrid
 
 4 determining how many target genes has been identified in both software

 Functional annotation was conducted by BLAST2GO SOFTWARE (Visual basic software)