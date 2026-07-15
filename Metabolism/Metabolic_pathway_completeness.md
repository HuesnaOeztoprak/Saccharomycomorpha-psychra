### To visualize pathway completeness, first EggNOG-mapper (v. 2.1.4) was used and then visualized with KEGGaNOG (v. 1.1.19)

##### EggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper) 
```
Functional annotation was done using EggNOG for details see Assembly_and_Annotation/RNA-seq de novo assembly.md
```

##### KEGGaNOG script (https://github.com/iliapopov17/KEGGaNOG)

```
KEGGaNOG -i EggNOG.001 -o KEGGaNOG_grouped -g
```


This provides a heatmap with the completness of metabolic pathways as output
