## Putative transporter protein identification workflow

### Homology searches against Transporter Classification Database (TCDB; Saier et al., 2021)
```sh
01_blast.sh
```
### EggNOG transporter classification and Pfam domain annotation
```sh
02_eggnog.sh
```
### DeepTMHMM (Hallgren et al., 2022) a deep learning-based topology predictor
```sh
04_DeepTMHmmer.sh
```
### Combining all results and validation of propoer integration of each result
```sh
05_combine.sh
```
```sh
05.2_adding_Pfams.sh
```
### Filtering all result-dataframe 
#### Proteins containing at least one predicted transmembrane helix and at least one transporter-related annotation were retained as candidates. 
#### three confidence tiers: high-confidence transporters were defined as proteins containing ≥3 evidence lines, medium confidence (n = 2), indicating agreement between two independent sources and low confidence (n =1 with ≥1 predicted TM helix), indicating a single transport-specific annotation supported by membrane topology. 
```sh
06_filter_dataframe.sh
```
