### Instead of counting how many exons are in 100kb, we want to calculate the exon base pairs per window

GFF="final.gff3"
GENOME="final.fasta.fai"
PREFIX="final_"

#extract exons from gff
awk '$3=="exon"' $GFF > exons_only.gff3


#to use awk, the spaces need to be replaced by tabs (FS)
sed -E 's/ +/\t/g' \
  exons_only.gff3 \
  > exons_only.tabbed.gff3

#convert to bed (awk instead of gff2bed) because we only have the augustus ggf3
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
}' exons_only.tabbed.gff3 \
>exons.bed

#sort
sortBed -i exons.bed > Sexons.sorted.bed

## extract "genome" (Scaffold name + length), to be able to know length for windows
cut -f1,2 $GENOME > final.fasta.genome.fai

#make windows
bedtools makewindows -g final.fasta.genome.fai -w 100000 > final_windows100kb.bed
sortBed -i final_windows100kb.bed > final_windows100kb.sorted.bed

#exon bases per window
#bases-uniq-f: Reports the fraction of distinct bases from ref-file ‘s element elements in map-file." (https://bedops.readthedocs.io/en/latest/content/summary.html)
#ref: 100kb windows
#map: exon coordinates
#-> how many bases in every 10kb window (=fraction) is covered by exons?

bedmap --echo --bases-uniq-f \
   windows100kb.sorted.bed \
    exons.sorted.bed \
    | tr "|" "\t" \
    > exon_bp_100kb.bed

echo "Calculated exon basepairs per 100 kb window."

#results are displayed as relative values, meaning: number exon counts per window/ windowsize 
