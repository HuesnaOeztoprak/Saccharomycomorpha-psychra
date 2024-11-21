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
## [Funannotate](https://github.com/nextgenusfs/funannotate/tree/master)
### RNA-seq mediated training of Agustus/GeneMArk
```
funannotate train -i Spa255.polished.alt1.softmasked.fasta -o funannotate_out -l Rna.FastP_R1.fastq.gz -r Rna.FastP_R2.fastq.gz --cpus 40 --trinity /home/hoeztopr/Data/hoeztopr/Spa/Transcriptome/fastP/trinity_v2.14_${i}_fastP.Trinity.fasta --no_trimmomatic --no_normalize_reads --species "Saccharomycomorpha"
```
### GENEMARK
```
gmes_linux_64/gmes_petap.pl --ES --max_intron 3000 --soft_mask 2000 --cores 40 --sequence Spa255.polished.alt1.softmasked.fasta
```
### Gene prediction
```
/NVME/Software/funannotate-docker predict -i Spa255.polished.alt1.softmasked.fasta -o funannotate_out \
-s "Saccharomycomorpha" --cpus 40 --organism other --busco_db protists --optimize_augustus --weights glimmerhmm:0 snap:0 --genemark_gtf genemark.gtf
```
### eggNOG mapper - [web application](http://eggnog-mapper.embl.de/))
emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_pu2cxinq --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_pu2cxinq --override -m diamond --dmnd_ignore_warnings -i /emapper_web_jobs/emapper_jobs/user_data/MM_pu2cxinq/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel > /emapper_web_jobs/emapper_jobs/user_data/MM_pu2cxinq/emapper.out 2> /emapper_web_jobs/emapper_jobs/user_data/MM_pu2cxinq/emapper.err

### Interproscan - using [Galaxy](https://usegalaxy.eu/jobs/)

### Assign functional annotation to gene predictions
```
funannotate annotate -i funannotate_out --cpus 60 --eggnog out.emapper.annotations --iprscan InterProScan-5.54-87.0.xml --busco_db funannotate_db/protists --force
```
