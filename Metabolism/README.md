### Overview of the workflow to investigate the metabolism of *S.* *Psychra*; the following scripts were utilized


#### To assess the metabolic pathway completeness by using EggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper) and KEGGaNOG script (https://github.com/iliapopov17/KEGGaNOG):
This results in a heatmap as overview (Figure S7)
```
Metabolic_pathway_completeness.md
```

and

#### To access metabolic pathway presence and expression via combining Salmon expression data (https://github.com/COMBINE-lab/salmon; TPM) and KofamKOALA (https://www.genome.jp/tools/kofamkoala/):
This results in an overview map (via ipath3; https://pathways.embl.de): Figure S8

```
Metabolic_pathway_presence_and_expression.md
```

This provides the file:
```
ipath3_Spa255_KEGGmapper_KofamKOALA.txt
```
This can be implemented into https://pathways.embl.de/ipath3.cgi?map=metabolic and then “select whole modules”. This provides an interactive map to click on KEGG modules and further investigate the metabolic map and its properties.



