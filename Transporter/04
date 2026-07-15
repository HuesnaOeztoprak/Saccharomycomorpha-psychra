module load seqkit   # if needed

seqkit split -p 12 \
  /home/hoeztopr/Scratch/hoeztopr/spa/assembly/final/Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.fasta \
  -O split_fasta

ls split_fasta/
#-rw-rw-r-- 1 hoeztopr hoeztopr 896K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_001.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 853K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_002.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 875K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_003.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 912K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_004.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 850K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_005.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 848K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_006.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 928K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_007.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 982K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_008.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 973K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_009.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 852K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_010.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 869K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_011.fasta
#-rw-rw-r-- 1 hoeztopr hoeztopr 853K Feb 19 10:22 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_012.fasta


####
#!/bin/bash

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export TORCH_NUM_THREADS=4
        #run deeptmhmm on individual files - requires all available cores, very intensive

bash run_deeptmhmm_2.sh split_fasta/Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.part_XX.fasta


#combine individual output files
grep -h "Number of predicted TMRs" results/*/TMRs.gff3 \
| sed 's/# //' \
| awk '{print $1"\t"$6}' \
> deeptmhmm_tm.txt


#| nTM | interpretation                                             |
#| --- | ---------------------------------------------------------- |
#| 0   | soluble protein                                            |
#| 1   | single-pass membrane protein                               |
#| ≥2  | multi-pass membrane protein (often transporter candidates) |


#genes with tm helix
awk '$2>0' deeptmhmm_tm.txt | wc -l
2335
#with multi-pass
awk '$2>=2' deeptmhmm_tm.txt | wc -l
1327
