!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# Paths
# ---------------------------

GFF="final.gff3"
GENOME="final.fasta.fai"
PREFIX="final_"

# ---------------------------
# 1. Clean scaffold names (remove spaces)
# ---------------------------
sed -E 's/[[:space:]]+/\t/g' "$GFF" > ${PREFIX}_clean.gff3

# ---------------------------
# 2. Extract genes
# ---------------------------
echo "Extracting genes..."
awk '$3=="gene"' ${PREFIX}_clean.gff3 > ${PREFIX}_genes.gff3

# ---------------------------
# 3. Convert GFF → BED (6 columns)
# ---------------------------
echo "Converting to 6-column BED..."
awk 'BEGIN{FS="\t";OFS="\t"} {
    id=".";
    n=split($9,a,";");
    for(i=1;i<=n;i++){
        if(a[i] ~ /^ID=/){
            sub(/^ID=/,"",a[i]);
            id=a[i];
        }
    }
    print $1, $4-1, $5, id, ".", $7
}' ${PREFIX}_genes.gff3 > ${PREFIX}_genes.bed

# ---------------------------
# 4. Sort BED
# ---------------------------
echo "Sorting gene BED..."
sortBed -i ${PREFIX}_genes.bed > ${PREFIX}_genes.sorted.bed

# ---------------------------
# 5. Prepare genome file (chrom + length)
# ---------------------------
echo "Preparing genome file..."
cut -f1,2 "$GENOME" > ${PREFIX}.genome

# ---------------------------
# 6. Create 100 kb windows
# ---------------------------
echo "Generating 100 kb windows..."
bedtools makewindows -g ${PREFIX}.genome -w 100000 > ${PREFIX}_windows100kb.bed

# ---------------------------
# 7. Sort windows
# ---------------------------
sortBed -i ${PREFIX}_windows100kb.bed > ${PREFIX}_windows100kb.sorted.bed

# ---------------------------
# 8. Compute gene basepair fraction per window
# ---------------------------
echo "Calculating gene basepair fraction per window..."
bedmap --echo --bases-uniq-f \
  ${PREFIX}_windows100kb.sorted.bed \
  ${PREFIX}_genes.sorted.bed \
  | awk -F"|" 'BEGIN{OFS="\t"} {if($2=="") $2=0; print $1,$2}' \
  > ${PREFIX}_gene_bp_fraction_100kb.bed

echo "Done: gene basepair fraction per 100 kb window calculated."
