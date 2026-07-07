# Functional motifs in secreted enzymes

### Adjust fasta file format
```
sed -E 's/^>([^ ]+).*/>\1/' \
final.braker.dedup.fasta \
> fixed.fasta
```

## Identify secreted enzymes
### [Phobius](https://phobius.sbc.su.se/) - A combined transmembrane topology and signal peptide predictor
#### Count secreted and transmembrane domains
```sh
python count_SP+TM_domains.py
```
#### Merging Phobius outout with KO terms
##### extract secreted proteins (5th coloum)
```sh
awk -F'\t' '$5 == "True"' Phobius_summary_Spa255.tsv > Secreted_Proteins.tsv
KOs.txt only with ~5k hits I will use KOs_emapper with 5k KOs
python merge_KO_data.py
```
##### isolate all the KO terms of the secreted proteins
```sh
#cut -f8 Secreted_Proteins_KO.tsv | sort | uniq > Unique_KO_Secreted_Values.txt
```
##### isolate Proein_IDs of secreted proteins
```sh
awk -F'\t' '$5 == "True" {print $1}' Phobius_summary_Spa255.tsv > secreted_protein_IDs.txt
```
##### extract corresponding protein sequences
```sh
python extract_secreted_proteins.py
```
### [SecretomeP-2.0](https://services.healthtech.dtu.dk/services/SecretomeP-2.0/) - ab initio predictions of non-classical i.e. not signal peptide triggered protein secretion
##### As max file size is limited to 100 seq
```sh
seqkit split -f secreted_proteins.fasta -s 100 -O fasta_chunks/
```
#### Run on website
##### Merge single outputs into one file
```sh
awk 'FNR==1 && NR!=1 {next} /^-+/ {next} 1' secretomeP2_results/secreted_proteins.part00*.txt > secretomeP2_results/SecretomeP_AllResults.tsv
```
##### Filter the secreted proteins
```sh
awk 'NR==1 || $4 >= 0.6' SecretomeP_AllResults.tsv > SecretomeP_Filtered.tsv
```
##### Extract secreted protein ID
```sh
cut -f5 SecretomeP_Filtered.tsv | tail -n +2 > Secreted_Protein_IDs.txt
```
### [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) - Prediction of Signal Peptides and their cleavage sites in all domains of life
##### As max file size is limited to 1000 proteins
```sh
seqkit split -s 1000 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.fasta -O split_fastas
```
#### Run on website
##### Merge single outputs into one file
```sh
tail -q -n +3 signalP6.0/output*/prediction_results.txt >> SignalP_AllResults.txt
```
##### Filter the secreted proteins
```sh
awk '$3 == "SP(Sec/SPI)" || NR==1' SignalP_AllResults.txt > SignalP_Secreted.txt
awk '$3 == "SP"' SignalP_AllResults.txt > SignalP_Secreted.txt
```
#### Run [DeepLoc2.1]()

##### Merge single outputs into one file
```
#tab delim
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
awk 'BEGIN {FS=","; OFS="\t"} {$1=$1; print}' output/results_${i}.csv > output/${i}_results_DeepLoc2.1.tsv
done

head output/1_results_DeepLoc2.1.tsv -n 2
#Protein_ID      Localizations   Signals Membrane types  Cytoplasm       Nucleus Extracellular   Cell membrane   Mitochondrion   Plastid Endoplasmic reticulum   Lysosome/Vacuole        Golgi apparatus Peroxisome Peripheral       Transmembrane   Lipid anchor >
#g1.t1   Cytoplasm|Nucleus               Soluble 0.6338000297546387      0.6726999878883362      0.047600001096725464    0.07289999723434448     0.17030000686645508     0.010099999606609344    0.191899999976158140.09520000219345093      0.10859999805688858  >

#Protein_ID      Localizations   Signals Membrane types  Cytoplasm       Nucleus Extracellular   Cell membrane   Mitochondrion   Plastid Endoplasmic reticulum      Lysosome/Vacuole        Golgi apparatus Peroxisome      Peripheral      Transmembrane   Lipid >
#FUN_000001-T1   Cytoplasm|Nucleus               Soluble 0.4844000041484833      0.567300021648407       0.2387000024318695      0.15080000460147858     0.454800009727478  0.013299999758601189    0.22869999706745148     0.07880000025033951     0.042199999094>
####

# Extract the first two header lines from the first file
head -n 1 output/1_results_DeepLoc2.1.tsv > DeepLoc2.1_AllResults.txt

# Merge all files, skipping headers in the rest
tail -q -n +2 output/*_results_DeepLoc2.1.tsv >> DeepLoc2.1_AllResults.txt
####
echo Extracellular
grep 'Extracellular' DeepLoc2.1_AllResults.txt | wc -l
echo Signal peptide
grep 'Signal peptide' DeepLoc2.1_AllResults.txt | wc -l
```


#### Run [NetGPI](https://services.healthtech.dtu.dk/services/NetGPI-1.1/)


## Homologous analysis
        ##### 1.     Blasting against Fungal and Bacterial and Protistan DB
#/home/jbast/anaconda3/envs/EDTA/bin/blastp
#/RAID/Data/databases/FunSecKB/Fungal_secretome.fasta

#### Against Fungal Secretomes
#makeblastdb -in /RAID/Data/databases/FunSecKB/Fungal_secretome.fasta -dbtype prot -out FunSecKB_DB

#blastp -query secreted_proteins.fasta -db /RAID/Data/databases/FunSecKB/FunSecKB_DB -outfmt 6 -evalue 1e-5 -num_threads 40 > fungal_hits.txt
        #done

#### Against Bacterial Secretomes (LAB-Secretome)
#blastp -query your_secreted.fasta -db LAB_Secretome.fasta -outfmt 6 -evalue 1e-5 > bacterial_hits.txt

#### Against ProtSecKB (protists)
#blastp -query your_secreted.fasta -db ProtSecKB.fasta -outfmt 6 -evalue 1e-3 > protist_hits.txt

#####      Parameters: -evalue 1e-5 for stringent fungal/bacterial comparisons && -evalue 1e-3 for broader protistan searches (higher diversity)


        ### 2.     Orthofinder
        #fungi to do!
#orthofinder -f your_secreted.fasta FunSecKB.fasta LAB_Secretome.fasta ProtSecKB.fasta
#mkdir orthofinder_input
#cp /home/hoeztopr/Data/databases/FunSecKB/*.fas ./orthofinder_input
#cp secreted_proteins.fasta ./orthofinder_input
#cp /home/hoeztopr/Data/databases/SecretomeP2.0/secretome_species_split/*.fasta ./orthofinder_input

#/NVME/Software/OrthoFinder_source/orthofinder.py -f orthofinder_input -t 40 -a 40  # Adjust threads (-t) and parallel analyses (-a) as needed

#/home/hoeztopr/Scratch/hoeztopr/spa/Transcriptome/secretion/orthologs
docker run --rm -u $(id -u) -v "$PWD":"$(pwd)" -w "$(pwd)" davidemms/orthofinder orthofinder -f orthofinder_input -t 40 -a 40



