import pandas as pd

# Load data
df = pd.read_csv("../overview.txt", sep="\t")

# Step 1: High-confidence (≥2 tools)
high_conf = df[df["#ofTools"] >= 2]

# Step 2: Include HMMER-only hits
hmmer_only = df[(df["HMMER"] != "N") & (df["#ofTools"] == 1)]

# Combine both for a broader confidence set
combined_conf = pd.concat([high_conf, hmmer_only]).drop_duplicates()

# Save results
high_conf.to_csv("CAZymes_high_confidence.tsv", sep="\t", index=False)
combined_conf.to_csv("CAZymes_combined_confidence.tsv", sep="\t", index=False)

print(f"✅ High-confidence: {len(high_conf)} genes")
print(f"✅ Combined (includes HMMER-only): {len(combined_conf)} genes")

################
        #this is the identificationa and counting of most abundant families
from collections import Counter
import re

# Helper to extract CAZy family names from tool columns
def extract_families(entry):
    if entry == "N" or pd.isna(entry): return []
    # Split by '+' and extract CAZy terms like GH23 or CE9
    return re.findall(r"(CBM\d+|GH\d+(_\d+)?|GT\d+|CE\d+|PL\d+|AA\d+(_\d+)?|GT2_Glycos_transf_2)", entry)

# Apply to relevant tool columns
all_families = []

for col in ["HMMER", "Hotpep", "DIAMOND"]:
    for entry in combined_conf[col]:
        all_families.extend(extract_families(str(entry)))

# Count occurrences
family_counts = Counter(all_families)
pd.DataFrame(family_counts.items(), columns=["Family", "Count"]).sort_values(by="Count", ascending=False).to_csv("CAZyme_family_abundance.tsv", sep="\t", index=False)
