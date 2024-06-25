!/bin/bash
  
############################
## Parameter setting starts

ncore=

srr=
# SRR id
# e.g., SRR9940470

fastq_dir=
# Fastq file directory
# e.g., path/to/SRR9940470

cutout_dir=
# Output directory for clean reads

align_dir=
# Output directory for alignment results

hisat3n_index=
# HISAT3N index directory

ref_dir=

fasta_file=
# e.g., Homo_sapiens.GRCh38.dna.primary_assembly.fa

## Parameter setting ends
############################

echo "Start analysis for DART-seq sequencing data ......"
echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

echo "\nFastp"
echo "\nfastp -i ${fastq_dir}/${srr}_1.fastq -I ${fastq_dir}/${srr}_2.fastq -o ${cutout_dir}/${srr}_clean_1.fastq -O ${cutout_dir}/${srr}_clean_2.fastq -q 20 -D -w ${ncore_fastp}\n"
fastp \
	-i ${fastq_dir}/${srr}_1.fastq -I ${fastq_dir}/${srr}_2.fastq \
	-o ${cutout_dir}/${srr}_clean_1.fastq -O ${cutout_dir}/${srr}_clean_2.fastq \
	-q 20 -D -w ${ncore}



echo "\nHISAT-3N alignment"
hisat-3n -x ${hisat3n_index} -q -1 ${cutout_dir}/${srr}_clean_1.fastq -2 ${cutout_dir}/${srr}_clean_2.fastq --base-change C,T --repeat | \
	samtools sort - -@ ${ncore} -o ${align_dir}/align.sort.sam -O sam



echo "\nGenerate HISAT-3N table\n"
hisat-3n-table \
        -p ${ncore} \
        --alignments ${align_dir}/align.sort.sam \
        --ref ${ref_dir}/${fasta_file} \
        --output-name ${align_dir}/output.tsv \
        --base-change C,T

echo "\nConvert Sam file to Bam file"
echo 'samtools view -@ ${ncore} -Shb ${align_dir}/align.sort.sam -o ${align_dir}/align.sort.bam\n'
samtools view -@ ${ncore} -Shb ${align_dir}/align.sort.sam -o ${align_dir}/align.sort.bam
rm ${align_dir}/align.sort.sam

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
