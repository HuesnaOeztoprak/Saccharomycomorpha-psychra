# Functional motifs in secreted enzymes
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

#### Run [NetGPI](https://services.healthtech.dtu.dk/services/NetGPI-1.1/)
