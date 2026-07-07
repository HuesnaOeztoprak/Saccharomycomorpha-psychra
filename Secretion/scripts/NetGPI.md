      #How to NetGPI
                #Step 1: Filter Secretory Pathway Proteins from signalP output: use SignalP_secreted_proteins.fasta

                #Step 2: Split SignalP_Secretory.fasta into 5000-Sequence Batches
#Use seqkit to split into 5000-sequence files:  
        #seqkit split -s 5000 ../SignalP_secreted_proteins.fasta -O split_fastas_netgpi

#Now the folder split_fastas_netgpi/ contains:
        #SignalP_Secretory.part_001.fasta
        #SignalP_Secretory.part_002.fasta
#Each file has ≤ 5000 sequences.

                #Step 3: Submit to NetGPI
#Manual Steps:
#Go to: NetGPI Server
#Upload a batch file (e.g., SignalP_Secretory.part_001.fasta).
#Choose "Short Output" (recommended for high-throughput).
#Run the Prediction & Download the Results.
#Repeat for all batch files.

        #output: GPI anchors that keep them attached to the outer face of the plasma membrane
output_protein_type.txt: list of anchored and not anchored GPI
output_mature.fasta: proteins which are GPI-Anchored (82)
output.gff3: annotation of set
