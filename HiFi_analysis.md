
# HiFi reads evaluation

## software requirements (version)
```
kat
Genomescope - jellyfish
Smudgeplot
```
## Genome structure
#### [KAT](https://github.com/TGAC/KAT) version 2.4.2
```
kat hist -o kat_hist hifi_reads.fastq.gz
kat gcp -o kat_gcp hifi_reads.fastq.gz
```
<img>

#### [Smudgeplot](https://github.com/KamilSJaron/smudgeplot) version 0.2.5
[KMC](https://github.com/tbenavi1/KMC) version 
```
mkdir tmp_smudge
ls hifi_reads.fastq.gz > FILES
kmc -k27 -ci1 -cs10000 @FILES kmcdb tmp_smudge
kmc_tools transform kmcdb histogram kmcdb_k27.hist -cx10000

L=$(smudgeplot.py cutoff kmcdb_k27.hist L)
U=$(smudgeplot.py cutoff kmcdb_k27.hist U)
echo $L $U

kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
```
![spa255_smudgeplot_smudgeplot](https://github.com/user-attachments/assets/ad0e8696-6bba-4d47-9a7f-6ba06e548b86)

#### [GenomeScope](https://github.com/tbenavi1/genomescope2.0)
![Genomescope_Spa255](https://github.com/user-attachments/assets/028df1fc-77cf-41f7-b9f6-0efb97e83463)

