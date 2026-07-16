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

print(f" High-confidence: {len(high_conf)} genes")
print(f" Combined (includes HMMER-only): {len(combined_conf)} genes")
