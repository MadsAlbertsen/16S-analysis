###################################################################################################
#
#  Version 4.1
#
#  This workflow script generates OTU tables from raw V13 16S amplicon data.
#
#  It is currently only supported for internal use on Aalborg University,
#  but feel free to use any parts that you find usefull.
#
#  Mads Albertsen, 2014
#
###################################################################################################

clear
echo ""
echo "Running: 16S V13 workflow version 4.1"
date

echo ""
echo "Finding your samples and copying them to the current folder"
while read samples
do
a="_";
NAME=$samples$a;
find /space/sequences/ -name $NAME*R1* 2>/dev/null -exec gzip -cd {} \; | head -n 1  | sed 's/\@/>/' >> id.txt
find /space/sequences/ -name $NAME*R1* 2>/dev/null -exec gzip -cd {} \; | head -n 200000 >> forward.fastq
find /space/sequences/ -name $NAME*R2* 2>/dev/null -exec gzip -cd {} \; | head -n 200000 >> reverse.fastq
done < samples

paste -d "\t" id.txt samples > sampleid.txt
date

echo ""
echo "Removing bad quality reads"
trimmomatic-0.32.jar PE forward.fastq reverse.fastq forward_qs.fastq.gz s1.fastq.gz reverse_qs.fastq.gz s2.fastq.gz SLIDINGWINDOW:1:3 MINLEN:275 2> temp.log
date

echo ""
echo "Merging reads"
flash -m 25 -M 200 forward_qs.fastq.gz reverse_qs.fastq.gz -o temp > flash.log
perl /space/users/malber06/16S-analysis/scripts/trim.fastq.length.pl -i temp.extendedFrags.fastq -o merged_l.fastq -m 425 -x 525 > temp.log
date

echo ""
echo "Removing potential phiX contamination"
split -l 10000000 merged_l.fastq splitreads
for sreads in splitreads*
do
b=".temp1";
c=".temp2";
NAMEb=$sreads$b;
NAMEc=$sreads$c;
usearch7 -fastq_filter $sreads -fastaout $NAMEb -quiet
usearch7 -usearch_global $NAMEb -db /space/databases/phix/phix.fasta -strand both -id 0.97 -notmatched $NAMEc -quiet
done

cat *.temp2 > merged_ls.fasta
rm splitreads*
date

echo ""
echo "Dereplicating reads"
perl /space/users/malber06/16S-analysis/scripts/uparse.to.dereplicate.pl -i merged_ls.fasta -o uniques.fa -r reads.fa
date

echo ""
echo "Clustering into OTUs"
usearch7 -cluster_otus uniques.fa -otus otus.fa -quiet
python /space/users/malber06/python_scripts/fasta_number.py otus.fa OTU_ > otus_named.fa
date

echo ""
echo "Mapping reads to the OTUs"
split -l 10000000 reads.fa splitreads

for sreads in splitreads*
do
a=".uc";
NAME=$sreads$a;
usearch7 -usearch_global $sreads -db otus_named.fa -strand plus -id 0.97 -uc $NAME -quiet
done

cat *.uc > readmap
rm splitreads*
mv readmap readmap.uc
date

echo ""
echo "Classifying the OTUs"
parallel_assign_taxonomy_rdp.py -c 0.8 -i otus_named.fa -o taxonomy -r /space/databases/midas/MiDAS_S119_1.16.0.fasta -t /space/databases/midas/MiDAS_S119_1.16.0.tax -O 10 --rdp_max_memory 30000 --rdp_classifier_fp /space/users/mni/software/rdp_classifier_2.2/rdp_classifier-2.2.jar
date

echo ""
echo "Comparing to MiDAS OTUs"
usearch7 -usearch_global otus_named.fa -db /space/databases/midas/MiDAS_S119_1.16.0.otus.fa -strand both -id 0.95 -userout MiDAS_S119_1.16.0_hits.txt -userfields query+target+id -quiet
date

echo ""
echo "Making an OTU table"
perl /space/users/malber06/16S-analysis/scripts/uparse.to.otutable.pl -s sampleid.txt -u readmap.uc -o otutable.txt -t taxonomy/otus_named_tax_assignments.txt
date


echo ""
echo "Calculating stats"
echo "ReadID\tSampleID\tRaw\tTrimmed\tMerged\tPhiXScreened\tMapped" > stats.txt
sed 's/ /:/' sampleid.txt | cut -d ":" -f1,2,3,8,9,10,11 | sort -k1,1 | sed 's/>//' > temp.stats1.txt
sed -n 1~4p forward.fastq | sed 's/\@//' | sed 's/ /:/' | cut -d ":" -f1,2,3,8,9,10,11 | sort | uniq -c > raw.log
gzip -cd forward_qs.fastq | sed -n 1~4p | sed 's/\@//' | sed 's/ /:/' | cut -d ":" -f1,2,3,8,9,10,11 | sort | uniq -c > trim.log
sed -n 1~4p merged_l.fastq | sed 's/\@//' | sed 's/ /:/' | cut -d ":" -f1,2,3,8,9,10,11 | sort | uniq -c > merged.log
grep ">" merged_ls.fasta | sed 's/ /:/' |  sed 's/>//' | cut -d ":" -f1,2,3,8,9,10,11 | sort | uniq -c > phix.log
grep "^H" readmap.uc | cut -f9 | sed 's/ /:/' | cut -d ":" -f1,2,3,8,9,10,11 | sort | uniq -c > mapped.log

join -a1 -a2 -e "0" -1 1 -2 2 -o 0 1.2,2.1 temp.stats1.txt raw.log > temp.stats2.txt
join -a1 -a2 -e "0" -1 1 -2 2 -o 0 1.2,1.3,2.1 temp.stats2.txt trim.log > temp.stats3.txt
join -a1 -a2 -e "0" -1 1 -2 2 -o 0 1.2,1.3,1.4,2.1 temp.stats3.txt merged.log > temp.stats4.txt
join -a1 -a2 -e "0" -1 1 -2 2 -o 0 1.2,1.3,1.4,1.5,2.1 temp.stats4.txt phix.log > temp.stats5.txt
join -a1 -a2 -e "0" -1 1 -2 2 -o 0 1.2,1.3,1.4,1.5,1.6,2.1 temp.stats5.txt mapped.log >> stats.txt
date

echo ""
echo "Removing temp files"
rm otus.fa
mv otus_named.fa otus.fa
rm *.log
rm readmap.uc
rm sampleid.txt
rm reads.fa
rm uniques.fa
rm -r jobs/
rm -r taxonomy/
rm *fastq*
rm temp*
rm merged*
rm id.txt
date

echo ""
echo "Done. Enjoy."
date
echo ""
