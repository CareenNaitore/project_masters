# curl https://www.ebi.ac.uk/ebisearch/ws/rest/rfam?query=entry_type:%22Family%22%20AND%20rna_type:%22IRES%22&format=json&size=100&start=0 # to get the different family 
  grep -Eo '\bRF[0-9]*\b' tRNA |sort|uniq > tRNA.index # identify the index 
  cat tRNA.index |wc -l #counting the number of id you have 
  cat tRNA.index >tRNA.txt #copy in text file 
  sed 's/\(.*\)/\1.fa.gz/'  tRNA.txt > tRNA_file #converting into a fasta file
  sed 's/^/ftp:\/\/ftp.ebi.ac.uk\/pub\/databases\/Rfam\/14.0\/fasta_files\//' tRNA_file > tRNA.file # adding the URL
   wget -i tRNA.file # get all the files in download 
   gunzip RF[0-9]*.fa.gz | cat RF[0-9]*.fa  > tRNA.fa # making An rNA database
   cat *.fa* | grep '^>' | sort | uniq -d  #sorts and removes duplicates  
   rm -r  RF[0-9]*.fa # removing all the RF* data


