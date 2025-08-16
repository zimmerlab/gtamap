#### extract regions for gene, exon and cds from gtf
```bash
grep -vP "^#" Homo_sapiens.GRCh38.114.gtf | awk '($3=="gene" || $3=="exon" || $3=="CDS") && ($1 < 23 || $1=="X" || $1=="Y" || $1=="MT") {chr=($1 < 23 || $1 =="X" || $1=="Y" ? "chr"$1 : ($1=="MT" ? "chrM" : "chr"$1)); if (match($0, /(ENSG[0-9]+)/, a) && match($0, /gene_biotype "([^"]+)"/, b)) print chr","$4","$5","$3","a[1]","b[1]}' > ./genes.csv
```

#### extract gene id, gene name and biotype from gtf
```bash
grep -vP "^#" Homo_sapiens.GRCh38.114.gtf | awk '$3=="gene" {if (match($0, /gene_id "([^"]+)"/, id) && match($0, /gene_name "([^"]+)"/, name) && match($0, /gene_biotype "([^"]+)"/, b)) print id[1]","name[1]","b[1]}' > ./gene-names.csv
```
