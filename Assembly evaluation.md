
# Assembly evaluation

[KAT](https://github.com/TGAC/KAT) version 2.4.2
```sh
kat comp -o kat_comp_hap0 hifi_reads.fastq.gz assembly.fasta
```

[BUSCO](https://busco.ezlab.org/) version 5.0.0

```sh
busco -i assembly.fasta -m genome -c 20 -o busco_out_eukaryota_odb10 -l eukaryota_odb10
busco -i assembly.fasta -m genome -c 20 -o busco_out_bacteria_odb10 -l bacteria_odb10
busco -i assembly.fasta -m genome -c 20 -o busco_out_fungi_odb10 -l fungi_odb10

```
### reads vs assembly
[minimap2](https://github.com/lh3/minimap2) version 2.24r1122
[SAMtools](https://github.com/samtools/samtools) version 1.11
```sh
minimap2 -ax map-hifi assembly.fasta hifi_reads.fastq.gz | samtools view -b | samtools sort -o minimap2_hifi.bam
```
### contaminants in assembly
```sh
kat sect -t 16 -o cov_GC% assembly.fasta hifi_reads.fastq
```
#### not done
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) version 2.6.0
```sh
blastn -query assembly.fasta -db nt -outfmt "6 qseqid staxids bitscore std sscinames scomnames" \
	-max_hsps 1 -evalue 1e-25 -out blast.out
```

[BlobTools2](https://blobtoolkit.genomehubs.org/blobtools2/) version 2.3.3
```sh
blobtools add --fasta assembly.fasta --cov minimap2_hifi.bam --hits blast.out \
	--busco busco_out_eukaryota_odb10/run_eukaryota_odb10/full_table.tsv \
	--taxdump taxdump --create hap0_BLOBDIR
blobtools view hap0_BLOBDIR
```

<img src="./fig/final_scaffolds_BLOBDIR.blob.circle.svg" width=600>
