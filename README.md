# m6AConquer

This repository holds the scripts used to process the sequencing data in the m6AConquer database. The scripts are organized into two main directories:

1. **IndexBuilding**: Contains scripts used to build alignment indexes for Hisat-3N and STAR.
2. **ProcessingScript**: Contains the actual scripts used in various sequencing data processing pipelines.

## Table of Contents

- [IndexBuilding](#indexbuilding-scripts)
  - [Example: HISAT3N_A2G.txt](#hisat3n_a2gtxt)
- [ProcessingScripts](#processing-scripts)
  - [Example: m6ACE-seq.txt](#m6ace-seqtxt)
- [Usage](#usage)
- [Contributing](#contributing)

## IndexBuilding

### Example: HISAT3N_A2G.txt

This script is used to build the HISAT3N (A-to-G) index.

#### Parameters

- `reference_dir`: Path to your reference directory.
- `gft_file`: GTF file, e.g., `Homo_sapiens.GRCh38.110.gtf`.
- `genome_file`: Genome file, e.g., `Homo_sapiens.GRCh38.dna.primary_assembly.fa`. This file should be indexed using `samtools faidx genome.fa`.
- `hisat3n_tool_dir`: Path to the Hisat3n tool directory, e.g., `path/to/Hisat3n/hisat-3n`.
- `genome_index_name`: Name for the genome index, e.g., `GRCh38_tran`.

#### Example

```bash
#!/bin/bash

############################
## Parameter setting starts

reference_dir=
# Path to your reference directory

gft_file=
# e.g., Homo_sapiens.GRCh38.110.gtf

genome_file=
# e.g., Homo_sapiens.GRCh38.dna.primary_assembly.fa
# should be indexed
# samtools faidx genome.fa

hisat3n_tool_dir=
# e.g., path/to/Hisat3n/hisat-3n

genome_index_name=
# e.g., GRCh38_tran

## Parameter setting ends
############################

echo "\nStart building HISAT3N (A-to-G) index ......\n"
echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

echo "\nStart genome index building......"

echo "Generate splice site index and exon index"
$hisat3n_tool_dir/hisat2_extract_splice_sites.py $reference_dir/$gft_file > $reference_dir/genome.ss
$hisat3n_tool_dir/hisat2_extract_exons.py $reference_dir/$gft_file > $reference_dir/genome.exon
echo "Hisat3n index building for genome"
hisat-3n-build --base-change A,G --repeat-index --ss $reference_dir/genome.ss --exon $reference_dir/genome.exon $reference_dir/$genome_file $reference_dir/$genome_index_name

echo "Finish genome index building......\n"
```

## ProcessingScripts

### Example: m6ACE-seq.txt

This script is used to process m6ACE-seq data.

#### Parameters

- `ncore_fastp`: Number of cores for `fastp`.
- `ncore_star`: Number of cores for `STAR`.
- `srr`: SRR ID, e.g., `SRR8380776`.
- `fastq_dir`: Fastq file directory, e.g., `path/to/SRR8380776`.
- `cutout_dir`: Output directory for clean reads.
- `align_dir`: Output directory for alignment results.
- `star_index`: STAR index directory.
- `ref_dir`: Reference directory.
- `gtf_file`: GTF file, e.g., `gencode.v45.annotation.gtf`.

#### Example

```bash
!/bin/bash

############################
## Parameter setting starts

ncore_fastp=

ncore_star=

srr=
# SRR id
# e.g., SRR8380776

fastq_dir=
# Fastq file directory
# e.g., path/to/SRR8380776

cutout_dir=
# Output directory for clean reads

align_dir=
# Output directory for alignment results

star_index=
# STAR index directory

ref_dir=

gtf_file=
# e.g., gencode.v45.annotation.gtf

## Parameter setting ends
############################

echo "Start trimming for m6ACE-seq sequencing data ......"
echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

if [ ! -d "${cutout_dir}" ];then mkdir -p "${cutout_dir}" ;fi

echo "\nfastp\n"
echo "\nfastp -i ${fastq_dir}/${srr}_1.fastq -I ${fastq_dir}/${srr}_2.fastq -o ${cutout_dir}/${srr}_clean_1.fastq -O ${cutout_dir}/${srr}_clean_2.fastq -g -U --umi_loc=read1 --umi_len=8 -D -w ${ncore_fastp}\n"
fastp \
	-i ${fastq_dir}/${srr}_1.fastq -I ${fastq_dir}/${srr}_2.fastq \
	-o ${cutout_dir}/${srr}_clean_1.fastq -O ${cutout_dir}/${srr}_clean_2.fastq \
	-g -q 20 -U --umi_loc=read1 --umi_len=8 -D -w ${ncore_fastp}


echo "\nSTAR --runThreadN ${ncore_star} --genomeDir ${star_index} --readFilesIn ${cutout_dir}/${srr}_clean_1.fastq ${cutout_dir}/${srr}_clean_2.fastq --sjdbGTFfile ${ref_dir}/${gtf_file} --outFileNamePrefix ${align_dir}/\n"
STAR --runThreadN ${ncore_star} --genomeDir ${star_index} --readFilesIn ${cutout_dir}/${srr}_clean_1.fastq ${cutout_dir}/${srr}_clean_2.fastq --sjdbGTFfile ${ref_dir}/${gtf_file} --outFileNamePrefix ${align_dir}/

echo "\nsamtools view -F 3852 -b ${align_dir}/Aligned.out.sam -o ${align_dir}/reads.bam\n"
samtools view -F 3852 -b ${align_dir}/Aligned.out.sam -o ${align_dir}/reads.bam

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
```

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/m6AConquer-scripts.git
   cd m6AConquer-scripts
   ```

2. Navigate to the appropriate directory (IndexBuilding or ProcessingScript).

3. Fill in the required parameters in the script files with your data or reference paths.

4. Run the scripts using bash:

   ```bash
   bash scriptname.txt
   ```

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request for any improvements or bug fixes.
