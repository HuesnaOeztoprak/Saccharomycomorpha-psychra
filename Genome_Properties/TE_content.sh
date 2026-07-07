### To use our FasTE output gff to visualize TE distribution along scaffolds

TECSV=final.fasta.out_tidy_noasterisk_repeatmasker.out_RM_TRIPS.csv
GENOME=final.fasta.fai

#name and length
cut -f1,2 "$GENOME" > genome.txt

#$3=qry_id, $8=merged_qrystart, $9=merged_qryend
#remove quotes
awk -F',' 'BEGIN{OFS="\t"}
NR>1{
    gsub(/"/,"",$3)
    gsub(/"/,"",$8)
    gsub(/"/,"",$9)
    chr=$3
    start=$8-1
    end=$9
    if(start < 0) start=0
    if(chr != "" && end > start) print chr, start, end
}' "$TECSV" \
> TE.unsorted.bed

#sort + merge overlapping TEs
sortBed -i TE.unsorted.bed > TE.sorted.bed
bedops --merge TE.sorted.bed > TE.merged.bed

#create 100 kb windows
bedtools makewindows -g genome.txt -w 100000 \
> windows100kb.bed
sortBed -i windows100kb.bed > windows100kb.sorted.bed

#calculate TE fraction per window
bedmap --echo --bases-uniq-f \
  windows100kb.sorted.bed \
  TE.merged.bed \
  | tr "|" "\t" \
  > windows100kb.TE_fraction.bed

#add SynVisio header
echo -e "Chr\tStart\tEnd\tValue" | cat - windows100kb.TE_fraction.bed \
> windows100kb.TE_fraction.txt
