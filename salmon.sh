salmon index -t transcript.fa -i salmon_index

salmon quant \
  -i salmon_index \
  -l A \
  -1 Spa255.rna.FastP_R1.fastq.gz \
  -2 Spa255.rna.FastP_R2.fastq.gz \
  -p 20 \
  -o salmon_output
