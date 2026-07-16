import pandas as pd

### Extract Protein IDs from FASTA ###
fasta_file = "protein.fasta"
protein_ids = []

with open(fasta_file, "r") as f:
    for line in f:
        if line.startswith(">"):  
            protein_id = line.split()[0][1:]  # Remove ">" and take the first part
            protein_ids.append(protein_id)

print(f" Extracted {len(protein_ids)} protein IDs from {fasta_file}")

### Load Prediction Data ###

# SignalP results
signalp_file = "SignalP_AllResults.csv"
signalp_df = pd.read_csv(signalp_file, sep="\t", comment="#", header=None, names=["Protein_ID", "Prediction", "OTHER", "SP_Sec", "CS_Position"])
signalp_df["SP_SignalP"] = signalp_df["Prediction"].apply(lambda x: 1 if "SP" in x else 0)

# Phobius results
phobius_file = "Phobius.tsv"
phobius_df = pd.read_csv(phobius_file, sep="\t")
phobius_df["SP_Phobius"] = phobius_df["Signal_Peptide"].apply(lambda x: 1 if x == "Y" else 0)

# DeepLoc results (keep all columns)
deeploc_file = "DeepLoc2.1_AllResults.txt"
deeploc_df = pd.read_csv(deeploc_file, sep="\t")

# NetGPI results
netgpi_file = "output_protein_type.txt"
netgpi_df = pd.read_csv(netgpi_file, sep="\t", comment="#", header=None, 
                        names=["Protein_ID", "Seq_length", "Pred_GPI_Anchored", "Omega_site_pos", "Likelihood", "Amino_acid"])

# Create a binary column indicating whether a protein is GPI-anchored
netgpi_df["GPI_Anchored_NetGPI"] = netgpi_df["Pred_GPI_Anchored"].apply(lambda x: 1 if "GPI-Anchored" in x else 0)

### Merge Data ###
merged_df = pd.DataFrame({"Protein_ID": protein_ids})  # Use FASTA IDs as the base

# Merge all data using left joins
merged_df = merged_df.merge(signalp_df[["Protein_ID", "SP_SignalP"]], on="Protein_ID", how="left")
merged_df = merged_df.merge(phobius_df[["Protein_ID", "SP_Phobius", "Membrane_Bound", "Secreted"]], on="Protein_ID", how="left")
merged_df = merged_df.merge(deeploc_df, on="Protein_ID", how="left")  # Include all DeepLoc columns
merged_df = merged_df.merge(netgpi_df[["Protein_ID", "GPI_Anchored_NetGPI"]], on="Protein_ID", how="left")

# Fill missing values with 0
merged_df.fillna(0, inplace=True)

# Rename Phobius columns
merged_df.rename(columns={"Membrane_Bound": "Membrane_Bound_Phobius", "Secreted": "Secreted_Phobius"}, inplace=True)

### Save Merged Data ###
merged_df.to_csv("merged_secreted_proteins.csv", sep="\t", index=False)
print(f"\n Merged data saved to merged_secreted_proteins.csv")
print("\n Sample Merged Data:")
print(merged_df.head())














