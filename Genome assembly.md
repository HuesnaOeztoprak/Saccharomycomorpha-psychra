# Assembly pipeline

## *De novo* assembly

[hifiasm](https://github.com/chhylp123/hifiasm) version 0.25
```sh
hifiasm -o assembly -l 0 hifi_reads.fastq.gz 
```
	##### The hifiasm assembly was manually curated using the graphical assembly
[PECAT](https://github.com/lemene/PECAT)
```sh
pecat.pl config cfg
```

## Scaffolding
[RagTag](https://github.com/malonge/RagTag) version 2.1.0
```sh
ragtag.py scaffold hifiasm_assembly.fasta pecat_assembly.fasta
```
##### The polished pecat asembly was used to scaffold the curated hifiasm assembly
