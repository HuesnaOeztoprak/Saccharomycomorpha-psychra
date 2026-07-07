# Genome annotation
## expected software version
	EDTA v1.9.8 	
	Bedtools Version: v2.26.0
	STAR v2.5.1a
	Trinity v2.1.1
	PASA 2.5.2
	stringtie v2.2.0
	TransDecoder 5.5.0
	EVidenceModeler 1.1.1
	BRAKER2 v2.1.6
	gffread v0.12.1
	emapper 2.1.6
	InterProScan version 5.61-93.0
	clusterProfiler v4.6.2

## software requirements (version)
```
samtools v1.9
TrimGalore v0.6.4_dev
Trinity v2.14
fastp 0.21.0
```
##  GenomeMask
### 1 Hardmask
	  EDTA.pl --genome polished.fasta --sensitive 1 --anno 1  --threads 50 --overwrite 1
### 2 Softmask
	  bedtools maskfasta -fi unpolished.fasta -fo polished.softmask.fasta -bed polished.fa.mod.EDTA.TEanno.gff3 -soft
   
## Trim raw RNA reads
#### trim adapters with [trimgalore](https://github.com/FelixKrueger/TrimGalore)
```
trim_galore --fastqc --gzip -j 8 --paired --max_n 0 \
RNA_end1.fastq.gz \
RNA_end2.fastq.gz

trim_galore --fastqc --gzip -j 8 --paired --max_n 0 \
RNA_R1.fastq.gz \
RNA_R2.fastq.gz
```
#### trim adapters with [FastP](https://github.com/OpenGene/fastp)
```
fastp -q 20 -u 5 -l 25 -y -n 0 --detect_adapter_for_pe --dont_overwrite --thread 16 \
--in1 RNA_end1.fastq.gz --in2 RNA_end2.fastq.gz \
--out1 RNA.trimmed_end1.fastq.gz --out2 RNA.trimmed_end2.fastq.gz \
-h RNA_FastP.html -j RNA_FastP.json
```

## [BRAKER](https://github.com/Gaius-Augustus/BRAKER)
### 
```
hisat2 -x final.softmasked.fasta -p 40 -1 RNA.trimmed_end1.fastq.gz -2 RNA.trimmed_end2.fastq.gz | samtools view -b -@ 30 | samtools sort -@ 30 -o final.bam
```
### use braker with GENEMARK, PROTHINT, AUGUSTUS ...
```
braker.pl --species Saccharomycomorpha_strain255 --genome final.softmasked.fasta --gff3 --UTR off --bam  final.bam --threads 40 --useexisting --PROTHINT_PATH /gmes_linux_64/ProtHint/bin/
```
### keep longest isoform
```
agat_sp_keep_longest_isoform.pl --gff ${PWD}/braker.gff3 -o ${PWD}/final.braker.dedup.gff3
```
### generate protein fasta
```
funannotate util gff2prot -g final.dedup.gff3 -f ../final.softmasked.fasta > final.dedup.protein.fasta
```
### eggNOG mapper version 2.1.4 
```
eggnog-mapper-2.1.4/emapper.py --cpu 30  --data_dir /RAID/Data/databases/eggnog-mapper-data \
 --score 0.01   --seed_ortholog_evalue 0.01 \
 -i final.dedup.protein.fasta -o EggNOG.001
```
### Interproscan - using [Galaxy](https://usegalaxy.eu/jobs/)

### Assign functional annotation to gene predictions
```
funannotate annotate -i funannotate_out --cpus 60 --eggnog out.emapper.annotations --iprscan InterProScan-5.54-87.0.xml --busco_db funannotate_db/protists --force
```
### Add InterproScan annotations to GFF3
```
agat_sp_manage_functional_annotation.pl \
  --gff braker.dedup.gff3 \
  --i InterProScan-5.59-91.0.xml \
  -o final.gff3
```
