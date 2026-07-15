## CAZyme annotation

### [dbcan2](https://github.com/linnabrown/run_dbcan/blob/master/README.md) - annotation tool for automated CAZyme annotation
#### Run dbcan2 with Hmmer and DIAMOND
```sh
python run_dbcan.py \
    /home/hoeztopr/Scratch/hoeztopr/spa/Transcriptome/funannotate_Spa255_out_old/annotate_results/Saccharomycomorpha.proteins.fa protein \
    --out_dir /home/hoeztopr/Scratch/hoeztopr/spa/Transcriptome/Carbohydrate_active_enzymes/dbcan2_Spa255 \
    --tools 'hmmer' 'diamond' --dia_cpu 12 --hmm_cpu 12 \
    --cgc_sig_genes all \
    --db_dir /NVME/Software/dbCAN2/run_dbcan/db/
```
#### Run Hotpep seperately and combine with dbcan2
##### split files for hotpep as input file should be 1000 proteins
```sh
python split_for_hotpep.py
```
##### run hotpep with splitted files
```sh
python train_many_organisms_many_families.py Input_Spa255 4 5 5
```
###### Input_Spa255: your input folder 4: number of threads (should match orfsX.txt count) 5: minimum hits 5: minimum frequency

##### rename output.txt
```sh
cp Results/output.txt  dbcan2_Spa255/Hotpep.out
```
##### rerun run_dbcan.py with hotpep
##### run overview_table.py based on diamond_output.txt, hmmer_output.txt and hotpep_output.txt

## Filter dbcan data table 
```sh
python 01_filter_dbscan_overview_table.py
```
#### Identify high confidence CAZymes
```
python 02_identify_high_confidence_CAZymes.py
```
#### Count CAZyme families
```
python 03_count_CAZymes_families.py
```

