#Usage: tidk explore [OPTIONS] <FASTA>
tidk explore --minimum 5 --maximum 12 assembly.fasta > explore.tsv
                        #likely telomeric region :AACCCT 
                        
#Usage: tidk search [OPTIONS] --string <STRING> --output <OUTPUT> --dir <DIR> <FASTA>
tidk search --string AACCCT --window 500 --dir . --output search_AACCCT_w500 --extension tsv assembly.fasta

        #plots the csv output as in svg
tidk plot --tsv search_AACCCT_w500_telomeric_repeat_windows.tsv --output tidk-plot_AACCCT
