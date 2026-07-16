#✳️
!/usr/bin/env python3
"""
telomere detection with gap-tolerant merging and
scaffold classification into Complete / Single_telomere / No_telomere.
"""

import pandas as pd

# ── Input ────────────────────────────────────────────────────────────────────
TIDK   = "Spa255.hifiasm.v25.ragtag.pecat_sorted_search_AACCCT_w500_telomeric_repeat_windows.tsv"
KARYO  = "scaffold_lengths.txt"
WINDOW = 500          # bp per tidk window
MAX_REPEATS = 83      # ≈ 500 / 6  (AACCCT = 6 bp)
FRAC_THR   = 0.20     # ≥ 20 % of window must be telomeric
GAP_TOL    = 5_000    # merge telomeric windows separated by ≤ 5 kb
END_DIST   = 20_000   # telomere must be within 20 kb of scaffold end

# ── Load data ────────────────────────────────────────────────────────────────
df = pd.read_csv(TIDK, sep="\t")
df.columns = ["id", "window", "fwd", "rev", "motif"]
df["start"] = df["window"] - WINDOW
df["end"]   = df["window"]
df["score"] = df["fwd"].fillna(0) + df["rev"].fillna(0)
df["frac"]  = df["score"] / MAX_REPEATS

karyo = pd.read_csv(KARYO, sep="\t", header=None, names=["id", "length"])
lengths = dict(zip(karyo["id"], karyo["length"]))

# ── Filter by fraction threshold ─────────────────────────────────────────────
hits = df[df["frac"] >= FRAC_THR].copy()
hits = hits.sort_values(["id", "start"])

# ── Gap-tolerant merging ─────────────────────────────────────────────────────
merged = []
for chrom, grp in hits.groupby("id"):
    s = e = None
    for _, r in grp.iterrows():
        rs, re = int(r["start"]), int(r["end"])
        if s is None:
            s, e = rs, re
        elif rs <= e + GAP_TOL:      # ← KEY FIX: tolerate small gaps
            e = max(e, re)
        else:
            merged.append([chrom, s, e])
            s, e = rs, re
    if s is not None:
        merged.append([chrom, s, e])

tel = pd.DataFrame(merged, columns=["scaffold", "start", "end"])

# ── Classify each interval as LEFT / RIGHT end ──────────────────────────────
def classify_end(row):
    scaf_len = lengths.get(row["scaffold"], None)
    if scaf_len is None:
        return "unknown"
    if row["start"] < END_DIST:
        return "left"
    if row["end"] > scaf_len - END_DIST:
        return "right"
    return "internal"          # interstitial — flag for manual check

tel["end_type"] = tel.apply(classify_end, axis=1)

# ── Per-scaffold summary ─────────────────────────────────────────────────────
summary = (
    tel.groupby("scaffold")["end_type"]
       .apply(lambda x: set(x))
       .reset_index()
)
summary.columns = ["scaffold", "end_types"]

def classify_scaffold(end_types):
    has_left  = "left"  in end_types
    has_right = "right" in end_types
    if has_left and has_right:
        return "Complete"
    if has_left or has_right:
        return "Single_telomere"
    return "Internal_only"   # interstitial telomere — unusual, worth flagging

summary["class"] = summary["end_types"].apply(classify_scaffold)

# Add scaffolds with NO telomere at all
all_scaffolds = pd.DataFrame({"scaffold": list(lengths.keys())})
summary = all_scaffolds.merge(summary, on="scaffold", how="left")
summary["class"] = summary["class"].fillna("No_telomere")

# ── Chromosome-scale filter (length ≥ 1 Mb AND has ≥ 1 end telomere) ────────
summary["length"] = summary["scaffold"].map(lengths)
summary["chrom_scale"] = (
    (summary["length"] >= 1_000_000) &
    (summary["class"].isin(["Complete", "Single_telomere"]))
)

# ── Write outputs ─────────────────────────────────────────────────────────────
tel[tel["end_type"].isin(["left","right"])].to_csv(
    "telomeres_highlight.txt", sep="\t",
    columns=["scaffold","start","end"], header=False, index=False
)

summary.to_csv("scaffold_classification.tsv", sep="\t", index=False)

complete = summary[summary["class"] == "Complete"]["scaffold"]
complete.to_csv("complete_chromosomes.txt", index=False, header=False)

chrom_scale = summary[summary["chrom_scale"]]["scaffold"]
chrom_scale.to_csv("chrom_scale_scaffolds.txt", index=False, header=False)

# ── Print summary ─────────────────────────────────────────────────────────────
print(summary["class"].value_counts().to_string())
print(f"\nChromosome-scale (≥1 Mb + telomere): {summary['chrom_scale'].sum()}")
print(f"  of which Complete (T2T):            {(summary['chrom_scale'] & (summary['class']=='Complete')).sum()}")
print("\nOutput files:")
print("  telomeres_highlight.txt")
print("  scaffold_classification.tsv")
print("  complete_chromosomes.txt")
print("  chrom_scale_scaffolds.txt")
