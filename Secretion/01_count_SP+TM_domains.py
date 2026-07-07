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

