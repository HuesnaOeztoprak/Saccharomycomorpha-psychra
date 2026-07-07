import pandas as pd

# Load dataset
df = pd.read_csv("merged_secreted_proteins.csv", sep="\t")

### Step 1: Signal Peptide Conflict Resolution ###
def classify_sp_conflict(row):
    # Binary SP predictions
    phobius_sp = row["SP_Phobius"] == 1
    signalp_sp = row["SP_SignalP"] == 1
    # DeepLoc-derived features
    deeploc_sp = row["Signals"] == "Signal peptide"
    deeploc_soluble = row["Membrane types"] == "Soluble"
    deeploc_extracellular = row["Extracellular"] >= 0.5
    # count functions
    sp_votes = sum([phobius_sp, signalp_sp, deeploc_sp])

    if phobius_sp and signalp_sp and deeploc_sp:
        return 1  # High-confidence SP, all agree
    elif phobius_sp and signalp_sp:
        return 1  # High confidence SP (DeepLoc neutral or noisy)
    elif deeploc_sp and deeploc_soluble and deeploc_extracellular:
        return 2  # Potential SP (non-classical or weak signal, extracellular)
    elif sp_votes >= 2:
        return 2  # Potential SP (2 out of 3 tools agree)
    else:
        return 3  # Likely no SP

df["SP_Conflict_Resolution"] = df.apply(classify_sp_conflict, axis=1)

        #adding labels to the numbers above
SP_CONFIDENCE_LABELS = {
    1: "High_confidence_SP",
    2: "Potential_SP",
    3: "No_SP"
}

df["SP_Confidence_Label"] = df["SP_Conflict_Resolution"].map(SP_CONFIDENCE_LABELS)

#keep only SP-positive proteins
df_filtered = df[df["SP_Conflict_Resolution"] <= 2].copy()


### Step 2: Transmembrane Conflict Resolution ###
def classify_tm_conflict(row):
    phobius_tm = row["Membrane_Bound_Phobius"]
    deeploc_tm = row["Transmembrane"] >= 0.5

    if not phobius_tm and not deeploc_tm:
        return 1  # 100% agreement (Non-TM)
    elif not phobius_tm and deeploc_tm:
        return 2  # Flag for validation (DeepLoc suggests TM)
    elif phobius_tm and not deeploc_tm:
        return 3  # Check for signal anchors
    else:
        return 0  # Ambiguous case

df_filtered["TM_Conflict_Resolution"] = df_filtered.apply(classify_tm_conflict, axis=1)
df_filtered = df_filtered[df_filtered["TM_Conflict_Resolution"] != 0]


### Step 3: GPI Anchoring Conflict Resolution ###
def classify_gpi_conflict(row):
    netgpi_score = row["GPI_Anchored_NetGPI"]
    deeploc_membrane = row["Transmembrane"] > 0.7
    deeploc_type = row["Membrane types"] == "Transmembrane"

    if netgpi_score > 0.5 and (deeploc_membrane or deeploc_type):
        return "High_confidence_GPI"
    elif netgpi_score > 0.5:
        return "Potential_GPI"
    else:
        return "No_GPI"

df_filtered["GPI_Conflict"] = df_filtered.apply(classify_gpi_conflict, axis=1)

# Remove strong GPI-anchored proteins
df_filtered = df_filtered[df_filtered["GPI_Conflict"] != "High_confidence_GPI"]

### Step 4: Competing Compartment Exclusion ###
competing_filter = (
    (df_filtered["Mitochondrion"] < 0.5) &
    (df_filtered["Nucleus"] < 0.5) &
    (df_filtered["Cytoplasm"] < 0.7) #signal peptide & partial cytosolic probabilty thus raised to from 0.5 to 0.65 (0.7?)
)

df_filtered["Compartment_Conflict"] = ~competing_filter  # True if flagged

# Apply competing_filter 
df_final = df_filtered[competing_filter]

# Save updated dataset
df_final.to_csv("final_secreted_proteins.csv", sep="\t", index=False)

print("Filtering complete! Wohoo! File saved as 'final_secreted_proteins.csv'.")
print(df_final[["Protein_ID", "SP_Conflict_Resolution", "TM_Conflict_Resolution", "GPI_Conflict", "Compartme
