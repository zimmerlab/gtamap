For each mutation rate and scenatio there needs to be one big tsv file that 
contains the results of the array job where each mapper was called gene-wise 
(except for the last genome_wide_simulation_sim)



$i \in \left\lbrace 0,1,2,3,4,5 \right\lbrace$


**DEFAULT:**  
./pipeline_all_rna_default/pipeline_all_{i}/benchmarking_dataset_{i}.all_results.tsv

**SENS:**  
./pipeline_all_rna_fair_final/benchmarking_dataset_{i}.all_results.tsv

**SENS GENOME WIDE**:  
./genome_wide_simulation_sim/dataset_{i}.all_results.tsv


Thes tsv files can be created by:

1. Run grid job where each mappers has 17,575 jobs and after its execution
   gets evaluated for the current gene with compare_to_ground_truth.py 
   
   One mapper call in `jobs.txt` looks like this (for one mutation rate and setting):   
   `/usr/bin/time -v -o /mnt/raidbio/biocluster/projekte/BA_weyrich/timing/benchmarking_dataset_0/hisat2_ENSG00000170948.time /mnt/raidbio/biocluster/praktikum/genprakt-ws24/gobi_gta/other_mappers/HISAT/hisat2_local/hisat2-2.2.0/hisat2     -p 8 -x /mnt/studtemp/weyrich/index_hisat/ENSG00000170948/ENSG00000170948 -1 /mnt/raidbio/biocluster/projekte/BA_weyrich/data/simulated_reads/benchmarking_dataset_0/fw.fastq   -2 /mnt/raidbio/biocluster/projekte/BA_weyrich/data/simulated_reads/benchmarking_dataset_0/rw.fastq -S /mnt/raidbio/biocluster/projekte/BA_weyrich/results_mixed/benchmarking_dataset_0/hisat2_ENSG00000170948.sam --no-unal --very-sensitive --no-mixed ; rm /mnt/raidbio/biocluster/projekte/BA_weyrich/results_mixed/benchmarking_dataset_0/hisat2_ENSG00000170948.sam`

2. All results of the array job for all mappers get stored in `/results_mixed/benchmarking_dataset_0/MAPPER_GENE.sam`
3. Call `concat_all.sh` on the directory containing the results for the current mutation rate and scenario
4. As a result you should get one of the big tsv files used in `./pipeline_all_rna_default/pipeline_all_{i}/benchmarking_dataset_{i}.all_results.tsv`

Calling `compute_all.py`: 
```sh
compute_all.py pipeline_all_rna_fair_final/benchmarking_dataset_0.all_res pipeline_all_rna_fair_final/benchmarking_dataset_0  rna
```

on one big tsv also comutes each metric separately for further inspection
