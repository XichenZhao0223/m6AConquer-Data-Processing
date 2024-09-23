# m6AConquer: Data Processing


## Table of Contents
 1. [Overview][1]
 2. [Index Building][2]
    [Example: HISAT3N_A2G.txt][3]
 4. [Processing Scripts][4]
    [Example: m6ACE-seq.txt][5]
 6. [Usage][6]
 7. [Contributing][7]
 
## 1. Overview
This repository holds the scripts used to process the sequencing data in the m6AConquer database. The scripts are organized into two main directories:

- **Index Building**: Contains scripts used to build alignment indexes for HISAT-3N and STAR.
- **Processing Scripts**: Contains the actual scripts used in various sequencing data processing pipelines.

An overview table of the scipts included in the repository is:

|  **Type**   | **Technique/Software**  | **File** |
|  :----:  | :----:  | :----: | :----:|
| Index Building  | HISAT3N | [HISAT3N_A2G.txt][8] |
| Index Building  | HISAT3N | [HISAT3N_C2T.txt][9] |
| Index Building | STAR | [STAR.txt][10] |
| Processing Scripts | DART-seq |[Bulk-scDART-seq.txt][11]<br>[DART-seq.txt][12]|
| Processing Scripts | eTAM-seq | [0.index_building.txt][13]<br>[1.trim_dedup.txt][14]<br>[2.map_pileup.txt][15]
| Processing Scripts | GLORI | [0.index_building.txt][16]<br>[1.trim_dedup.txt][17]<br>[2.run_GLORI.txt][18]<br>[change_reference.py][19]
| Processing Scripts | m6ACE-seq | [m6ACE-seq.txt][20]
| Processing Scripts | m6A-REF-seq |[0.prepare_reference_file.txt][21]<br>[1.cutadapt_align.txt][22]<br>[2.count_motif.txt][23]
| Processing Scripts | m6A-SAC-seq | [1.trim_dedup_paired.txt][24]<br>[1.trim_dedup_single.txt][25]<br>[2.star_align_paired.txt][26]<br>[2.star_align_single.txt][27]<br>[3.bam_pileup.txt][28]
| Processing Scripts | MAZTER-seq | [MAZTER-seq.txt][29]|
| Processing Scripts | Oxford-nanopore (m6Anet) | [0.fast5_basecall.txt][30]<br>[1.m6Anet.txt][31]

## 2. Index Building

### Example: HISAT3N_A2G.txt

This script is used to build the HISAT3N (A-to-G) index.

#### Parameters

- `reference_dir`: Path to your reference directory.
- `gft_file`: GTF file, e.g., `Homo_sapiens.GRCh38.110.gtf`.
- `genome_file`: Genome file, e.g., `Homo_sapiens.GRCh38.dna.primary_assembly.fa`. This file should be indexed using `samtools faidx genome.fa`.
- `hisat3n_tool_dir`: Path to the Hisat3n tool directory, e.g., `path/to/Hisat3n/hisat-3n`.
- `genome_index_name`: Name for the genome index, e.g., `GRCh38_tran`.

#### Codes

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

## 3. Processing Scripts

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

#### Codes

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

## 4. Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/XichenZhao0223/m6AConquer-Data-Processing
   cd m6AConquer-Data-Processing
   ```

2. Navigate to the appropriate directory (IndexBuilding or ProcessingScript).

3. Fill in the required parameters in the script files with your data or reference paths.

4. Run the scripts using bash:

   ```bash
   bash scriptname.txt
   ```

## 5. Contributing

Contributions are welcome! Please fork the repository and submit a pull request for any improvements or bug fixes.


  [1]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#overview
  [2]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#index-building
  [3]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#example-hisat3n_a2gtxt
  [4]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#processing-scripts
  [5]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#example-m6ace-seqtxt
  [6]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#usage
  [7]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing?tab=readme-ov-file#contributing
  [8]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/IndexBuilding/HISAT3N_A2G.txt
  [9]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/IndexBuilding/HISAT3N_C2T.txt
  [10]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/IndexBuilding/STAR.txt
  [11]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/DART-seq/Bulk-scDART-seq.txt
  [12]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/DART-seq/DART-seq.txt
  [13]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/eTAM-seq/0.index_building.txt
  [14]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/eTAM-seq/1.trim_dedup.txt
  [15]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/eTAM-seq/2.map_pileup.txt
  [16]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/GLORI/0.index_building.txt
  [17]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/GLORI/1.trim_dedup.txt
  [18]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/GLORI/2.run_GLORI.txt
  [19]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/GLORI/change_reference.py
  [20]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6ACE-seq/m6ACE-seq.txt
  [21]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-REF-seq/0.prepare_reference_file.txt
  [22]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-REF-seq/1.cutadapt_align.txt
  [23]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-REF-seq/2.count_motif.txt
  [24]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-SAC-seq/1.trim_dedup_paired.txt
  [25]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-SAC-seq/1.trim_dedup_single.txt
  [26]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-SAC-seq/2.star_align_paired.txt
  [27]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-SAC-seq/2.star_align_single.txt
  [28]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/m6A-SAC-seq/3.bam_pileup.txt
  [29]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/MAZTER-seq/MAZTER-seq.txt
  [30]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/Oxford-nanopore%28m6Anet%29/0.fast5_basecall.txt
  [31]: https://github.com/XichenZhao0223/m6AConquer-Data-Processing/blob/main/ProcessingScripts/Oxford-nanopore%28m6Anet%29/1.m6Anet.txt
