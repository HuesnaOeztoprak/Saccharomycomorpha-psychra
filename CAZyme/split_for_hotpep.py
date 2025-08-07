import os
from pathlib import Path
from Bio import SeqIO


input_fasta = "../funannotate_Spa255_out_old/annotate_results/Saccharomycomorpha.proteins.fa"  # ← Change this to your FASTA path
output_dir = "Hotpep/Input_Spa255"
num_chunks = 12  # Adjust based on how many orfsX.txt files you want

records = list(SeqIO.parse(input_fasta, "fasta"))
chunk_size = len(records) // num_chunks + 1

os.makedirs(output_dir, exist_ok=True)

for i in range(num_chunks):
    chunk = records[i*chunk_size:(i+1)*chunk_size]
    output_file = os.path.join(output_dir, f"orfs{i+1}.txt")
    with open(output_file, "w") as out_f:
        SeqIO.write(chunk, out_f, "fasta")

print(f"✅ Split into {num_chunks} files in {output_dir}")
