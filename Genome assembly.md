# Phased assembly pipeline

## *De novo* assembly

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.25
```sh
hifiasm -o assembly -l 0 hifi_reads.fastq.gz 
```
	##### The hifiasm assembly was manually curated using the graphical assembly
[PECAT](https://github.com/lemene/PECAT)) version 
```sh
pecat.pl config cfg
  #### include cfg script?
```

## [Purging](https://github.com/dfguan/purge_dups)
```sh
minimap2 -x map-hifi -t 30 ${i}.fasta hifi_reads.fastq.gz | gzip -c - > minimap2_${i}.paf.gz
purge_dups/bin/split_fa draft_assembly.fasta > ${i}.split
minimap2 -xasm5 -DP ${i}.split ${i}.split | gzip -c - > ${i}.split.self.paf.gz
purge_dups/bin/pbcstat minimap2_${i}.paf.gz

        #check hifiasm.purge_${i}.png (k-mer spectra) to set cutoffs
       
purge_dups/bin/calcuts -l 1 -m 300 -u 1000 PB.stat > cutoffs 2>calcults.log
purge_dups/scripts/hist_plot.py -c cutoffs PB.stat hifiasm.purge_${i}.png
        #purge haplotigs and overlaps
purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov ${i}.split.self.paf.gz> dups.bed 2>hifiasm.purge_${i}.log
        #get purged primary and haplotig seq from draft assembly
/purge_dups/bin/get_seqs dups.bed draft_assembly.fasta
```
```sh
mv purged.fa hifiasm.purged.alt1.fasta
mv hap.fa hifiasm.purged.alt2.fasta
```
## Scaffolding
[RagTag](https://github.com/malonge/RagTag)
	#### The polished pecat asembly was used to scaffold the curated hifiasm assembly
```sh
ragtag.py scaffold ref.fasta query.fasta
```

?[instaGRAAL](https://github.com/koszullab/instaGRAAL) version 0.1.6 no-opengl branch

```sh
instagraal --level 5 --cycles 100 hicstuff_out assembly.fasta instagraal_out
```

```sh
instagraal-polish -m polishing -f assembly.purged.fasta -j NNNNNNNNNN \
	-i instagraal_out/hicstuff_out/test_mcmc5/info_frags.txt \
	-o assembly.hic_scaffolds.fasta
```

?## Polishing 

[HyPo](https://github.com/kensung-lab/hypo) v1.0.3
[minimap2](https://github.com/lh3/minimap2) version 2.24r1122
[SAMtools](https://github.com/samtools/samtools) version 1.11
```sh
minimap2 --secondary=no --MD -ax map-hifi gap_filled.fasta hifi_reads.fastq.gz | samtools view -Sb - > mapped-ccs.bam
samtools sort -o mapped-ccs.sorted.bam mapped-ccs.bam
samtools index mapped-ccs.sorted.bam

hypo -d gap_filled.fasta -r hifi_reads.fastq.gz -s 200m -c 100 -b mapped-ccs.sorted.bam \
	-o polished.fasta
```
