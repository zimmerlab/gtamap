
Requires a combined log file of 17,575 calls per filter type with following columns


```sh
head ~/projects/bachelor_arbeit/analysis/filter_genelength/combined_logs_filter.tsv
timestamp       readsProcessed  readsAfterFiltering     readsMapped     confidentMappings       mappingLocations        percentConfidentOfMapped  percentFiltered percentMappedTotal      percentMappedFilter     meanLocationsPerRead    allocKB heapAllocKB     stackAllocKB    totalFileSize     bytesProcessed  percentFileProcessed    gene_id mutation_rate   type
2025-10-05 12:23:04     2632900 8446    2481    95      28220   3.83    0.32    0.09    29.37   11.37   18827   18827   1280    16806115601680611560      100.00  binned_ENSG00000000003  benchmarking_dataset_0_log      Binned
....
```
Also requires a "./combined_times.tsv" file:

```sh
file    tool    gene    user_time       system_time     cpu_percent     wall_time       wall_time_sec   max_rss_MB   type
gtamap_global_ENSG00000156172.time      gtamap  ENSG00000156172 11.92   0.85    80      0:15.84 15.84   129.73  Global
gtamap_binned_ENSG00000156172.time      gtamap  ENSG00000156172 12.93   0.55    91      0:14.72 14.72   178.99  Binned
```

where each mapper call of gtamap was timed with `time` across all 17,575 genes (x2 since both filters were tested)
