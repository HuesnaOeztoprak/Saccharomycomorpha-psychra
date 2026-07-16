#✳️
#!/usr/bin/env python3
import pandas as pd

cov = pd.read_csv("coverage.bedgraph", sep="\t", header=None, names=["chr","start","end","cov"])
cov["win"] = (cov["start"] // 100000).astype(int)

rows = []
for (chrom, win), grp in cov.groupby(["chr","win"], sort=False):
    wstart = int(win * 100000)
    wend = int((win + 1) * 100000)
    total = ((grp["end"] - grp["start"]) * grp["cov"]).sum()
    bases = (grp["end"] - grp["start"]).sum()
    mean_cov = total / bases if bases else 0
    rows.append([chrom, wstart, wend, mean_cov])

out = pd.DataFrame(rows, columns=["chr","start","end","value"])
out.to_csv("coverage_100kb.txt", sep="\t", header=False, index=False, float_format="%.6f")
