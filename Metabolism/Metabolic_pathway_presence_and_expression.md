### KO terms of KofamKOALA and salmons transcript quantification enable pathway presence and expression visualization in ipath3

#### To visualize KEGG annotated genes and KO terms KofamKOALA webtool was utilized (https://www.genome.jp/tools/kofamkoala/; last accessed 31.03.2026)


#### Salmon is a tool to quantificate transcripts from RNA data in transcripts per million (tpm)

Script and input data as shown in:
```
Assembly_and_Annotation/salmon.sh
```


quant.sf and KofamKOALA output were then combined to produce a table of all genes found with their respective KO term and 
tpm value (expressed is equal/ greater than 1)



### In R: 

```
#Since KofamKOALA only allows input sequences up to 10000 entries, first our data needs to be split and afterwards manually combined
#library(seqinr)

#setwd()

#prots <- read.fasta("assembly.fasta")
#length(prots)


#split_at <- 8500

#write.fasta(
#  sequences = prots[1:split_at], 
#  names = names(prots)[1:split_at], 
#  file.out = "part1.fasta"
#)

#write.fasta(
#  sequences = prots[(split_at+1):length(prots)], 
#  names = names(prots)[(split_at+1):length(prots)], 
#  file.out = "part2.fasta"
#)

#cat("Part1:", length(prots[1:split_at]), "proteins\n")
#cat("Part2:", length(prots[(split_at+1):length(prots)]), "proteins\n")

#then combine manually after KofamKOALA run

#Script for ipath3

setwd()

library(dplyr)
library(readr)

##read salmon data and KofamKoala output to combine
quant <- read_tsv("quant.sf", col_types = cols())
kofam_raw <- read_lines("assemblyKofamKOALA_result_ko_combined.txt")

# KofamKOALA parsing
kofam_list <- list()
for(i in 1:length(kofam_raw)) {
  line <- kofam_raw[i]
  parts <- strsplit(line, "\\s+")[[1]]
  name <- parts[1]
  ko <- if(length(parts) > 1) parts[2] else NA_character_
  kofam_list[[i]] <- data.frame(Name = name, KO = ko, stringsAsFactors = FALSE)
}
#create list of gene with KO term
kofam <- bind_rows(kofam_list) %>% 
  mutate(KO = ifelse(is.na(KO) | KO == "", NA_character_, KO))

## all KOs present in our KofamKoala output recieve thin lines
all_kos <- kofam %>% 
  filter(!is.na(KO)) %>% 
  distinct(KO) %>% 
  pull(KO)

## only KOs that also have a tpm equal/ above 1 tpm recieve a thick line 
expr_kos <- quant %>% 
  filter(TPM >= 1) %>% 
  select(Name, TPM) %>% 
  inner_join(kofam %>% filter(!is.na(KO)), by = "Name") %>% 
  group_by(KO) %>% 
  summarise(TPM_sum = sum(TPM), .groups = "drop") %>%
  pull(KO)

## combining thick and thin line information in one data set
ipath_lines <- c()

# first all thin lines
for(ko in all_kos) {
  ipath_lines <- c(ipath_lines, paste0(ko, " W2"))  
}

# after all thick lines
for(ko in expr_kos) {
  ipath_lines <- c(ipath_lines, paste0(ko, " W10")) 
}

#save
writeLines(ipath_lines, "ipath3_KEGGmapper.txt")

