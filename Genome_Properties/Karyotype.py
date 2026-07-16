#✳️
#!/usr/bin/env python3
import pandas as pd

fai = pd.read_csv("scaffold_lengths.txt", sep="\t", header=None, names=["chr","len"])
with open("karyotype.txt", "w") as out:
    for i, row in fai.iterrows():
        out.write(f"chr - {row['chr']} {row['chr']} 0 {int(row['len'])} grey\n")
