Two directories are required:
Input fastq for mappers: 500,000 first read pairs of SRR30429425

1. /mnt/studtemp/weyrich/results_gene_wise_gtamap:
    - contains 17,575 SAM files of gene-wise runs of gtamap in format `gtamap_GENEID.sam`
2. /mnt/studtemp/weyrich/real_rna_benchmark_genes:
    - contains 17,575 SAM files x3 of gene-wise runs of external rna mappers when called in 
      default config (MAPPER_GENEID.sam)
