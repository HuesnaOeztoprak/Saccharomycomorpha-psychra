
# HiFi reads evaluation

## Genome structure
#### [KAT](https://github.com/TGAC/KAT) version 2.4.2
```
kat hist -o kat_hist hifi_reads.fastq.gz
kat gcp -o kat_gcp hifi_reads.fastq.gz
```
![hist_Spa255 hifi_reads_x1000](https://github.com/user-attachments/assets/ca8b2f4a-c8c0-41fa-b7cc-2c2264813c1c) ![gcp_Spa255 hifi_reads mx](https://github.com/user-attachments/assets/a062062f-dfdf-4d58-b2d5-26bca81dd226)



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
