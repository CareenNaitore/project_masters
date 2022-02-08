#joining the first length distribution with the second one 
for i in *1_1.fastq_trim_30.filt.filt2.length; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *1_1.fastq_trim_30.filt.filt2.length;
do  
a=$(basename "$i"  1_1.fastq_trim_30.filt.filt2.length) #identifies the index for for each library
awk 'NR == FNR {
  _[$1] = $2; next
  }
{
  print $0, _[$1]
  }' "$i" "$a"2_1.fastq_trim_30.filt.filt2.length |sort -u > "$a".sorted
done
#sorting the third length distribution file  
for i in *3_1.fastq_trim_30.filt.filt2.length; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *3_1.fastq_trim_30.filt.filt2.length;
do
a=$(basename "$i"  3_1.fastq_trim_30.filt.filt2.length) #identifies the index for for each library
 sort -u "$i" > "$a"sorted.3
done
 #joining the third file to the other groups 
for i in *.sorted; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.sorted;
do
a=$(basename "$i"  .sorted) #identifies the index for for each library
join "$i" "$a"sorted.3 >"$a".joined
done

#formating the file for r analysis
for i in *.joined; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.joined;
do
a=$(basename "$i"  .joined) #identifies the index for for each library
cat "$i" | awk '{print $1, $3, $2 ,$4}' OFS="\t" > "$a"_lengthdistribution.tsv
done
