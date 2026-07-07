#only the first protID coloum
awk '{print $1}' final_secreted_proteins.csv > final_secreted_Protein_IDs.txt
        #count 986
wc -l final_secreted_Protein_IDs.txt
