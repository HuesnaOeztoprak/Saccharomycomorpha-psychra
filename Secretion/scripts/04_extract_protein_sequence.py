from Bio import SeqIO

# Input files
fasta_file = "protein.fasta"
id_file = "secreted_Protein_IDs.txt"
output_fasta = "secreted_proteins.fasta"

# Read secreted protein IDs
with open(id_file) as f:
    secreted_ids = set(line.strip() for line in f)

# Extract matching sequences
count = 0
with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_id = record.id.split()[0]  # Extract only the first part of header
        if fasta_id in secreted_ids:
            SeqIO.write(record, out_f, "fasta")
            count += 1

print(f" Extracted {count} secreted proteins to {output_fasta}")
