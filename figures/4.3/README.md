Called `featureCounts` on genome wide SAM results
of 4.1.1/4.2:

```sh
featureCounts -p -T 8 -a /mnt/raidbio/biocluster/projekte/BA_weyrich/data/genome/gencode.v48.annotation.gtf -o counts_HISAT2.txt HISAT2_subset_all.sam
featureCounts -p -T 8 -a /mnt/raidbio/biocluster/projekte/BA_weyrich/data/genome/gencode.v48.annotation.gtf -o counts_STAR.txt STAR_subset_all.sam
featureCounts -p -T 8 -a /mnt/raidbio/biocluster/projekte/BA_weyrich/data/genome/gencode.v48.annotation.gtf -o counts_Minimap2.txt Minimap2_subset_all.sam
```


