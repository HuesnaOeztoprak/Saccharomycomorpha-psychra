module load seqkit   # if needed

seqkit split -p 12 \
  protein.fasta \
  -O split_fasta

ls split_fasta/
protein.part_001.fasta
protein.part_002.fasta
protein.part_003.fasta
protein.part_004.fasta
protein.part_005.fasta
protein.part_006.fasta
protein.part_007.fasta
protein.part_008.fasta
protein.part_009.fasta
protein.part_010.fasta
protein.part_011.fasta
protein.part_012.fasta


####
#!/bin/bash

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export TORCH_NUM_THREADS=4
        #run deeptmhmm on individual files - requires all available cores, very intensive

bash run_deeptmhmm_2.sh split_fasta/protein.part_XX.fasta


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
