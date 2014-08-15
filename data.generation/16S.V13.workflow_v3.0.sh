###################################################################################################
#  V3.0
#  This workflow script generates otu_tables from raw V13 16S amplicon data.
#
#  It is currently only supported for internal use on Aalborg University,
#  but feel free to use any parts if needed.
#
#  Mads Albertsen 2014
#
#  V3.0
#  Local install of RDP classifier requried. Remember to insert correct file path to RDP excutable
#  in the script.
#
###################################################################################################

clear
echo ""
echo "Running: 16S V13 workflow version 3.0"

echo ""
echo "Finding your samples and copying them to your current directory "
while read samples
do
a="_";
NAME=$samples$a;
find /space/sequences/ -name $NAME* -exec cp -t . {} \;
done < samples

echo ""
echo "Unpacking all data and keeping the first 50.000 reads"
head -q -n 200000 *R1* > forward.fastq
head -q -n 200000 *R2* > reverse.fastq

echo ""
echo "Removing bad quality reads"
java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar PE forward.fastq reverse.fastq p1.fastq.gz s1.fastq.gz p2.fastq.gz s2.fastq.gz SLIDINGWINDOW:1:3 MINLEN:275

echo ""
echo "Removing potential phiX contamination"
bowtie2 -x /space/databases/phix/phix -1 p1.fastq.gz -2 p2.fastq.gz -S screened.sam --un-conc screened.fastq --al-conc phix.fastq -p 60 --reorder

echo ""
echo "Merging reads"
flash -r 300 -f 475 -s 50 -m 25 -M 200 screened.1.fastq screened.2.fastq -o merged
perl /space/users/malber06/16S-analysis/scripts/trim.fastq.length.pl -i merged.extendedFrags.fastq -o merged_screened.fastq -m 425 -x 525

echo ""
echo "Dereplicating reads"
perl /space/users/malber06/16S-analysis/scripts/uparse.to.dereplicate.pl -i merged_screened.fastq -o uniques.fa -r reads.fa

echo ""
echo "Clustering into OTUs"
usearch7 -cluster_otus uniques.fa -otus otus.fa -uc clusters.uc
python ~/python_scripts/fasta_number.py otus.fa OTU_ > otus_named.fa

echo ""
echo "Mapping reads to the OTUs"
split -l 10000000 reads.fa splitreads

for sreads in splitreads*
do
a=".uc";
NAME=$sreads$a;
usearch7 -usearch_global $sreads -db otus_named.fa -strand plus -id 0.97 -uc $NAME
done

cat *.uc > readmap
rm splitreads*
mv readmap readmap.uc

echo ""
echo "Generating a sample id file"
head -n 1 *R1* | sed -n '1~3p' | cut -f2 -d " " | cut -f1 -d "_" | cat -n > sample.name.txt
head -n 1 *R1* | sed -n '2~3p' | sed 's/\@/>/' | cat -n > sample.header.txt
join sample.header.txt sample.name.txt | sed 's/ /,/' | cut -f2-3 -d , | rev | sed 's/ /,/' | rev | sed 's/,/\t/' > sampleid.txt

echo ""
echo "Making an OTU table"
perl /space/users/malber06/16S-analysis/scripts/uparse.to.otutable.pl -s sampleid.txt -u readmap.uc -o otutable.txt

echo ""
echo "Classifying the OTUs"
parallel_assign_taxonomy_rdp.py -c 0.8 -i otus_named.fa -o taxonomy -r /space/databases/green_gene/gg_13_5_otus/rep_set/97_otus.fasta -t /space/databases/green_gene/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt -O 10 --rdp_max_memory 30000 --rdp_classifier_fp /space/users/mni/software/rdp_classifier_2.2/rdp_classifier-2.2.jar
