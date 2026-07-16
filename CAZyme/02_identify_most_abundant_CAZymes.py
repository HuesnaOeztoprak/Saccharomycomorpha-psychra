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
