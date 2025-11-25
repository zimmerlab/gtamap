
Uses three tsvs as input:
1. index_rna.tsv:  
   file    gene    user_time       system_time     cpu_percent     wall_time       wall_time_sec   max_rss_MB     exit_status     tool 
   
   Metrics of all index commands timed by the log of the array job
2. index_bowtie (contains same metric but only for bowtie)
3. index_bwa (contains same metric but only for bwa)

```sh
python3 analyze_index.py index_rna.tsv index_bowtie.tsv index_bwa.tsv  index_analysis
```


Call `logs_to_tsv_fast.sh` on a log dir of the corresponding grid job to create the tsv files

