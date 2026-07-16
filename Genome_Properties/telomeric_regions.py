#✳️
#!/usr/bin/env python3
import pandas as pd

f = "../tidk/AACCCT_w500_telomeric_repeat_windows.tsv"
df = pd.read_csv(f, sep="\t")

df.columns = ["id", "window", "forward_repeat_number", "reverse_repeat_number", "telomeric_repeat"]
df["start"] = df["window"] - 500
df["end"] = df["window"]
df["score"] = df["forward_repeat_number"].fillna(0) + df["reverse_repeat_number"].fillna(0)

thr = 20
hits = df[df["score"] >= thr].copy()
hits = hits.sort_values(["id", "start"])

merged = []
for chrom, grp in hits.groupby("id"):
    s = None
    e = None
    for _, r in grp.iterrows():
        if s is None:
            s, e = int(r["start"]), int(r["end"])
        elif int(r["start"]) <= e:
            e = max(e, int(r["end"]))
        else:
            merged.append([chrom, s, e])
            s, e = int(r["start"]), int(r["end"])
    if s is not None:
        merged.append([chrom, s, e])

out = pd.DataFrame(merged, columns=["chr","start","end"])
out.to_csv("telomeres_highlight.txt", sep="\t", header=False, index=False)

