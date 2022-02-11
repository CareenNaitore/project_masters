# MiRNA TRANSCRIPTOME ANALYSIS In TSESTE (*Glossina pallidipes*) USING BIOINFORMATIC SOFTWARE.
This projects entails applying the different programming languages of my master's  project, involving identification of developmental miRNA genes in an insect vector (*G. Pallidipes* ) for african trypanosomiasis. One of the Neglected tropical diseases  (ntd) in Kenya.

## Introduction

MicroRNAs (miRNAs) are small non-coding RNA of 21-24 nucleotide long (Lucas & Raikhel, 2013)⁠. Different studies suggest that miRNA’s function is to regulate gene expression in almost all physiological and biological process in both eukaryotes and prokaryotes (Carthew & Sontheimer, 2009; Valencia-Sanchez et al., 2006). These genes are crucial for mRNA stability and posttranslation process within an organism. This led to the need of understanding how miRNA genes play a crucial role in the developmental cycle of our insect vector. Eventually, providing information that could be used as a new innovative way to control these vectors. Hence my study involved using both molecular and computational analysis to provide this information. In this plaform, i am providing the code and data used in my study: i applied languages such as `Bash language`,`python`, `R programming language`. This experiment was undertaken in ICIPE Kenya and my master's degree offered in Jomo Kenyatta University of Agriculture and Technology. 

Small RNA transcriptome data for tsetse (G.pallidipes) was generated using illumina hiseq 2500 from Korea.Data was analysed using softwares such as  [miRdeep2] (https://github.com/rajewsky-lab/mirdeep2)  for identification of novel and Known miRNA genes in our species, [R programming] (https://cran.r-project.org/bin/windows/base/) for statistical analysis (differential miRNA gene expression), [miRanda software](http://www.mirtoolsgallery.org/miRToolsGallery/node/1055) and [RNAhybrid](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1538877/) were applied to identify target prediction for our miRNA genes. The reads for the datasets are from developmental stages (larva,pupa,teneral adult(young adult recently hatched from pupa),non-teneral adult (mature reproductive adults).The datasets are in raw forms of fastq libraries. which have three replicates each. Data is contained the **Data folder** and the script for my analysis is contained in **script_folder**

### PIPELINE

The methodology used in this study to analysis miRNA includes:

    1.Quality assements
    2.Genome annotation
    3.Filtering out other types of small non coding RNA
    4.MiRNA analysis(miRNA identification and miRNA differential expression)
    5.Target prediction

Project Task

    1. Reproduce the pipeline by setting up the workspace in the system you will be using
    2. Perform downstream analysis for miRNA analysis
    
Through the project, you need to demonstrate collaborative research skills, informative visualization, and report writing.
