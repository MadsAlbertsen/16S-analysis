###################################################################################################
#
#  This workflow script generates otu_tables from raw V13 16S amplicon data.
#
#  It is currently only supported for internal use on Aalborg University,
#  but feel free to use any parts if needed.
#
#  Mads Albertsen 2013
#
###################################################################################################
# PARAMETERS
# take the first n reads from each raw file
maxseqs=50000

clear
echo ""
echo "Running: 16S V13 workflow version 2.0"

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
nlines=$((maxseqs*4))
gunzip *.gz
head -q -n $nlines *R1* > r1.fastq
head -q -n $nlines *R2* > r2.fastq

echo ""
echo "Screening for phiX contamination"
bowtie2 --version | head -n 1
bowtie2 -x /space/databases/phix/phix -1 r1.fastq -2 r2.fastq -S screened.sam --un-conc screened.fastq --al-conc phix.fastq -p 40 --reorder

echo ""
echo "Merging read 1 and 2 using FLASH"
flash -r 300 -f 475 -s 50 -m 50 -M 200 screened.1.fastq screened.2.fastq -o merged

echo ""
echo "Removing reads outside the length criteria"
perl /space/users/malber06/16S-analysis/scripts/trim.fastq.length.pl -i merged.extendedFrags.fastq -o merged_screened.fastq -m 425 -x 525

echo ""
echo "Generating pseudo merged file"
sed -n '1~4s/^@/>/p;2~4p' merged_screened.fastq | sed 's/ 1:N:0//' > merged_screened_pseudo.fasta

echo ""
echo "Generating a sample id file"
head -n 1 *R1* | sed -n '1~3p' | cut -f2 -d " " | cut -f1 -d "_" | cat -n > sample.name.txt
head -n 1 *R1* | sed -n '2~3p' | sed 's/\@/>/' | cat -n > sample.header.txt
join sample.header.txt sample.name.txt | sed 's/ /,/' | cut -f2-3 -d , | rev | sed 's/ /,/' | rev | sed 's/,/\t/' > sampleid.txt

echo ""
echo "Formatting the reads for qiime"
perl /space/users/malber06/16S-analysis/scripts/merged.to.qiime.pl -i merged_screened_pseudo.fasta -s sampleid.txt -r -u 2 -m 10000

echo ""
echo "Qiime version and dependencies"
print_qiime_config.py

echo ""
echo "Qiime - de novo clustering, making OTUs"
pick_otus.py -i seqs.fna -m uclust -s 0.97

echo ""
echo "Qiime - pick representative set of sequences"
pick_rep_set.py -i uclust_picked_otus/seqs_otus.txt -f seqs.fna -o rep_set.fna -m most_abundant

echo ""
echo "Qiime - assign taxonomy - using MiDAS taxonomy"
parallel_assign_taxonomy_rdp.py -c 0.8 -i rep_set.fna -o rdp_assigned_taxonomy -r /space/databases/green_gene/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -t /space/databases/midas/midas_taxonomy.txt --rdp_max_memory 10000 -O 10

#GreenGenes: -r /space/databases/green_gene/gg_13_5_otus/rep_set/97_otus.fasta -t /space/databases/green_gene/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt
#MiDAS:      -r /space/databases/green_gene/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -t /space/databases/midas/midas_taxonomy.txt
#Silva:      -r /space/databases/silva/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta -t /space/databases/silva/Silva_111_post/taxonomy/97_Silva_111_taxa_map_RDP_6_levels.txt

echo ""
echo "Qiime - generate biom file"
make_otu_table.py -i uclust_picked_otus/seqs_otus.txt -t rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o otu_table.biom

echo ""
echo "Qiime - generate otu table from biom file"
convert_biom.py -i otu_table.biom -o otu_table.txt -b --header_key taxonomy

echo ""
echo "Cleaning your directory"

rm *.fastq
rm histograms.txt
rm *merged*
rm sampleid.txt
rm reordered.fa
rm reduced.stats.txt 
rm screened.sam
rm sample.*
rm map.txt
rm -r jobs/
rm seqs.fna

echo ""
echo "Completed - Enjoy"
echo ""
