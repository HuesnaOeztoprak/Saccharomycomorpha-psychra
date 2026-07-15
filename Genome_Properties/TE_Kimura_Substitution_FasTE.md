### Overview of FasTE workflow:

Bell et al. 2021 (https://github.com/ellenbell/FasTE) "FasTE" workflow was utilized to create a *de novo* TE library. To further disclose used parameters and versions, steps were performed with the following:


(for not further disclosed steps, please refer to the FasTE workflow)

#### Step 1: 
The TE library generation using EDTA (v. 2.2.2; Ou et al., 2019) with parameters
```
“--sensitive 1 --anno 1” 
```
for whole-genome TE annotation


#### Step 2: 
Reclassification of TEs at superfamily level was enabled using DeepTE (python = v. 3.6;
tensorflow-gpu = v. 1.14.0; biopython; keras = v. 2.2.4; numpy = v. 1.16.0; Yan et al., 2020)


#### Step 3: 
Annotation against custom library using RepeatMasker (v. 4.2.2; http://www.repeatmasker.org/;
Tarailo‐Graovac & Chen, 2009) with parameters 
```
“-a -s -no_is”
```

#### Prior to downstream analysis,the output files were processed and filtered with the script “RM_Trips” in R (version 4.4.3; R Core Team, 2025) following the FasTE workflow


#### Visualization of Kimura Substitution Levels: 
Workflow as described in Öztoprak et. al 2025 (https://www.science.org/doi/10.1126/sciadv.adn0817) was performed.

