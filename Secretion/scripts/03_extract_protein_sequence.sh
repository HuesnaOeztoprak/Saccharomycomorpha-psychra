#only the first protID coloum
awk '{print $1}' secreted_proteins.csv > secreted_Protein_IDs.txt
        #count 986
wc -l secreted_Protein_IDs.txt
