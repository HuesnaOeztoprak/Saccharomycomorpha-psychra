#To add all patched up accessions to the combined dataframe.
grep -v '^#' pfam.tbl | awk '{print $1"\t"$2}' | grep -E "^ABC_ATPase|^ABC_C\b|^ABCC10_N" \
  | awk '{acc=$2; sub(/\.[0-9]+$/,"",acc); print acc}' | sort | uniq

#this will give the PF accession numbers which need to be added
#PF09818
#PF21709
#PF24358

echo "PF00439" >> pfam_accessions_only.txt   # ABC_ATPase
echo "PF12848" >> pfam_accessions_only.txt   # ABC_C
echo "PF09818" >> pfam_accessions_only.txt   #
echo "PF21709" >> pfam_accessions_only.txt   #
echo "PF24358" >> pfam_accessions_only.txt   #

