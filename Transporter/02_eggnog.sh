#Evidence line via EggNOG annotation

awk -F"\t" '$12 != "-" && $12 != "" {print $1, $12}' out.emapper.annotations > eggnog_KO.txt
echo eggnog_KO.txt
wc -l eggnog_KO.txt

#KEGG_TC
awk -F"\t" '$18 != "-" && $18 != "" {print $1, $18}' out.emapper.annotations > eggnog_TC.txt
echo eggnog_TC.txt
wc -l eggnog_TC.txt #526

#PFAM of the TC proteins
awk -F"\t" '$18 != "-" && $18 != "" {print $1, $21}' \
out.emapper.annotations \
> eggnog_TC+Pfam.txt
echo eggnog_TC+Pfam.txt
wc -l eggnog_TC+Pfam.txt #526
