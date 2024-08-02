# Genome annotation

## software requirements (version)
```
samtools v1.9
TrimGalore v0.6.4_dev
Trinity v2.14
```
## Prepare Reads
#### trim adapters with [trimgalore](https://github.com/FelixKrueger/TrimGalore)
```
trim_galore --fastqc --gzip -j 8 --paired --max_n 0 \
RNA_end1.fastq.gz \
RNA_end2.fastq.gz

trim_galore --fastqc --gzip -j 8 --paired --max_n 0 \
RNA_R1.fastq.gz \
RNA_R2.fastq.gz
```
## RNA-Seq _de novo_ Assembly - using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
```
Trinity --seqType fq --max_memory 100G \
--left RNA.trimmed_end1.fastq.gz \
--right RNA.trimmed_end2.fastq.gz \
--CPU 40 --output trinity_v2.14
```
#### get trinity statistics
```
util/TrinityStats.pl trinity_v2.14.Trinity.fasta > Trinity_v214_Stats.txt
```
#### Index Trinity.fasta file
```
samtools faidx trinity_v2.14.Trinity.fasta -o TrinityIndexed.fasta
```
#### check initial completeness of trinity assembly with [BUSCO](https://github.com/WenchaoLin/BUSCO-Mod)
```
busco -m transcriptome -c 15 -l eukaryota_odb10 --update-data -o busco_euk -i trinity_v2.14.Trinity.fasta
busco -m transcriptome -c 15 -l bacteria_odb10 --update-data -o busco_euk -i trinity_v2.14.Trinity.fasta
busco -m transcriptome -c 15 -l archaea_odb10 --update-data -o busco_euk -i trinity_v2.14.Trinity.fasta
```
