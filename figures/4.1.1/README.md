
The genome wide comp uses these three sam files:

```py
sam_files = [
    "/mnt/studtemp/weyrich/results_genome_wide/HISAT2_subset_all.sam",
    "/mnt/studtemp/weyrich/results_genome_wide/STAR_subset_all.sam",
    "/mnt/studtemp/weyrich/results_genome_wide/Minimap2_subset_all.sam",
]
```

They are the result of mapping the first 500,000 reads of SRR30429425
to the human reference genome with default params.

