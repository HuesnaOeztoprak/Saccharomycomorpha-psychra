## Workflow to identify functional motifs in secreted enzymes

#### if needed: adjust fasta file format (no terminal *)
```
sed -E 's/^>([^ ]+).*/>\1/' final.braker.dedup.fasta > fixed.fasta
```

### Run [Phobius](https://phobius.sbc.su.se/) - A combined transmembrane topology and signal peptide predictor
#### Run Phobius on website and then count secreted and transmembrane domains from Phobius output:
```
python 01_count_SP+TM_domains.py
```

```
import pandas as pd

# Input & Output Files
phobius_file = "phobius_short_output_Spa255.txt"
output_file = "Phobius_summary_Spa255.tsv"

# Read Phobius short output
with open(phobius_file, "r") as f:
    lines = f.readlines()[1:]  # Skip header

# Data storage
results = []

# Parse each line
for line in lines:
    parts = line.strip().split()
    if len(parts) < 4:  # Skip malformed lines
        continue

    seq_id, tm, sp, prediction = parts[0], int(parts[1]), parts[2], parts[3]

    # Determine localization
    is_secreted = sp == "Y" and tm == 0
    is_membrane = tm > 0
    is_cytoplasmic = sp == "0" and tm == 0

    results.append([seq_id, tm, sp, prediction, is_secreted, is_membrane, is_cytoplasmic])

# Convert to DataFrame
df = pd.DataFrame(results, columns=["Protein_ID", "TM_count", "Signal_Peptide", "Prediction", 
                                    "Secreted", "Membrane_Bound", "Cytoplasmic"])

# Save as TSV
df.to_csv(output_file, sep="\t", index=False)
print(f"Summary saved to {output_file}")

# Print counts
print("\nSummary:")
print(f"Total proteins analyzed: {len(df)}")
print(f"Secreted proteins: {df['Secreted'].sum()}")
print(f"Membrane-bound proteins: {df['Membrane_Bound'].sum()}")
print(f"Cytoplasmic proteins: {df['Cytoplasmic'].sum()}")
```

### Run [SecretomeP-2.0](https://services.healthtech.dtu.dk/services/SecretomeP-2.0/) - ab initio predictions of non-classical i.e. not signal peptide triggered protein secretion
##### As max file size is limited to 100 seq
```sh
seqkit split -f secreted_proteins.fasta -s 100 -O fasta_chunks/
```
#### Run on website and then:
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

### Run [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) - Prediction of Signal Peptides and their cleavage sites in all domains of life
##### As max file size is limited to 1000 proteins
```sh
seqkit split -s 1000 Spa255.hifiasm.v25.ragtag.pecat.braker.dedup.fasta -O split_fastas
```
#### Run on website
##### Merge single outputs into one file
```sh
tail -q -n +3 signalP6.0/output*/prediction_results.txt >> SignalP_AllResults.txt
```
##### Extract the secreted proteins
```sh
awk '$3 == "SP(Sec/SPI)" || NR==1' SignalP_AllResults.txt > SignalP_Secreted.txt
awk '$3 == "SP"' SignalP_AllResults.txt > SignalP_Secreted.txt
```

### Run [DeepLoc2.1]()
##### As max file size is limited to 500 seq
```sh
seqkit split -f secreted_proteins.fasta -s 500 -O split_fastas_deeploc/
```
##### Merge single outputs into one file
```sh
#to ensure tab delim
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
awk 'BEGIN {FS=","; OFS="\t"} {$1=$1; print}' output/results_${i}.csv > output/${i}_results_DeepLoc2.1.tsv
done

## this is how it is supposed to look like:
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
```


### Run [NetGPI](https://services.healthtech.dtu.dk/services/NetGPI-1.1/)
                #Step 1: Filter Secretory Pathway Proteins from signalP output: use SignalP_secreted_proteins.fasta

                #Step 2: Split SignalP_Secretory.fasta into 5000-Sequence Batches
#Use seqkit to split into 5000-sequence files:  
        #seqkit split -s 5000 ../SignalP_secreted_proteins.fasta -O split_fastas_netgpi

#Now the folder split_fastas_netgpi/ contains:
        #SignalP_Secretory.part_001.fasta
        #SignalP_Secretory.part_002.fasta
#Each file has ≤ 5000 sequences.

                #Step 3: Submit to NetGPI
#Manual Steps:
#Go to: NetGPI Server
#Upload a batch file (e.g., SignalP_Secretory.part_001.fasta).
#Choose "Short Output" (recommended for high-throughput).
#Run the Prediction & Download the Results.
#Repeat for all batch files.

        #output: GPI anchors that keep them attached to the outer face of the plasma membrane
output_protein_type.txt: list of anchored and not anchored GPI
output_mature.fasta: proteins which are GPI-Anchored (82)
output.gff3: annotation of set

### Generate a combined dataframe and filter final set of secreted proteins.

```sh
python 01_get_dataframe.py
```

```sh
python 02_filterin_dataframe.py
```

```sh
bash 03_extract_protein_sequences.sh
```

```sh
python 04_extract_protein_sequences.py
```
