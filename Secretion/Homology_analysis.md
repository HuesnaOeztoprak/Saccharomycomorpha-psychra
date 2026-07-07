## Homologous analysis - Blasting against Fungal and Bacterial DB

#### Against Fungal Secretomes
#makeblastdb -in /RAID/Data/databases/FunSecKB/Fungal_secretome.fasta -dbtype prot -out FunSecKB_DB
#blastp -query secreted_proteins.fasta -db /RAID/Data/databases/FunSecKB/FunSecKB_DB -outfmt 6 -evalue 1e-5 -num_threads 40 > fungal_hits.txt

#### Against Bacterial Secretomes (LAB-Secretome)
#blastp -query your_secreted.fasta -db LAB_Secretome.fasta -outfmt 6 -evalue 1e-5 > bacterial_hits.txt

#####      Parameters: -evalue 1e-5 for stringent fungal/bacterial comparisons && -evalue 1e-3 for broader protistan searches (higher diversity)
