# Homology Search against Transporter Database

blastp \
  -db tcdb_blast \
  -query Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.fasta \
  -out tcdb_hits.tsv \
  -evalue 1e-20 \
  -qcov_hsp_perc 50 \
  -max_target_seqs 5 \
  -outfmt 6 -num_threads 30

cut -f1 tcdb_hits.tsv | sort | uniq | wc -l
        #981 

awk '
{
split($2,a,"|")
tc=a[4]
split(tc,b,".")
print $1"\t"b[1]"."b[2]
}
' tcdb_hits.tsv | sort -u > tcdb_families.txt
