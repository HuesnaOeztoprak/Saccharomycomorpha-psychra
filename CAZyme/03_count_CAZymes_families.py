import pandas as pd
import re

# Load the original file (tab-separated)
df = pd.read_csv("CAZymes_high_confidence.tsv", sep="\t")

# Keep only relevant columns
df = df[["Gene ID", "HMMER"]]

# Remove rows where HMMER is missing or 'N'
df = df[df["HMMER"] != "N"]

# Extract CAZy family name (remove coordinates like (28-344))
df["CAZy_Family"] = df["HMMER"].str.extract(r"^([A-Z]+\d+)")

# Count occurrences of each family
family_counts = df["CAZy_Family"].value_counts().reset_index()
family_counts.columns = ["CAZy_Family", "Count"]

# Sort by count (descending)
family_counts = family_counts.sort_values(by="Count", ascending=False)

# Save to file (same format as your previous run)
family_counts.to_csv(
    "CAZy_Family_Counts.txt",
    sep=" ",
    index=False,
    header=False
)

print("✅ New CAZy_Family_Counts.txt generated successfully!")
print(family_counts.head())
