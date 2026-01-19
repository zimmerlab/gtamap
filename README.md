# GTAMap
GTAMap is a locus-centric read mapping tool that inverts the traditional mapping paradigm. Rather than asking "where does this read belong?", GTAMap asks "what evidence exists for this locus?" This approach maximizes sensitivity for hypothesis-driven analysis of specific genomic regions, particularly genes with paralogs, pseudogenes, or repetitive elements.

## Key Features

- Targeted mapping to user-defined genomic regions
- Weighted evidence model based on k-mer uniqueness
- Support for both DNA-seq and RNA-seq data
- Improved recall for challenging genes (PMS2, CBS, SMN1/SMN2, etc.)
- Computational efficiency for locus-specific analysis

## Usage
### `index-pre`
Extract gene sequences from a genome
```
Usage:
  gtamap index-pre [flags]

Flags:
  -d, --downstream int           Number of bases to add downstream of the gene end position.
  -f, --fasta string             Fasta file (required)
      --fasta-file-name string   Output FASTA file name (within output directory) (only for --single-file) (default: genes.fa)
  -i, --fasta-index string       Fasta index file (default: [--fasta].fai)
  -l, --gene-ids strings         Gene IDs to extract (comma-separated). If not provided, all genes are extracted.
  -g, --gtf string               Genome annotation file (.gtf) (required)
  -h, --help                     help for index-pre
  -o, --output string            Output directory
  -s, --single-file              Write all gene sequences to a single fasta file
  -u, --upstream int             Number of bases to add upstream of the gene start position.

Global Flags:
      --config string     Path to config YAML
      --loglevel string   Log output level (ERROR, INFO, DEBUG) (default "INFO")
```

### `index`
Build the gtamap index (.gtai)
```
Usage:
  gtamap index [flags]

Flags:
  -f, --fasta string              Fasta file (required)
  -h, --help                      help for index
  -i, --index-file-name string    Output gtamap index file name (within output directory) (default: index.gtai)
  -k, --kmer-occurrences string   Genome-wide kmer occurrences (.tsv) to compute region-specific kmer frequencies
  -o, --output string             Output directory (required)
  -m, --regionmask string         Regionmask file (.bed) containing specific mismatch constraints per region
  -u, --use-fasta-file-name       Use the name of the fasta file (without extension) as index file name (within output directory)

Global Flags:
      --config string     Path to config YAML
      --loglevel string   Log output level (ERROR, INFO, DEBUG) (default "INFO")
```

### `map`
Run mapping
```
Usage:
  gtamap map [flags]

Flags:
  -h, --help                   help for map
  -i, --index string           Index file (*.gtai) (required)
  -o, --output string          Output directory (required)
  -p, --progress string        Progress file name (within output dir) or full path (default: mapping_progress.tsv)
  -r, --read-origin string     Specify read origin: 'dna' or 'rna' (required)
  -1, --reads-r1 string        FASTQ file containing the R1 reads (required)
  -2, --reads-r2 string        FASTQ file containing the R2 reads (if paired-end)
  -f, --sam-file-name string   Output SAM file name (within output directory) (default: aligned.sam)
  -e, --stagetime string       Progress stage time file name (within output dir) or full path (default: mapping_stage.tsv)
  -s, --stats string           Mapping statistics file name (within output dir) or full path (default: mapping_stats.tsv)
  -t, --threads int            Number of threads to use (default: all available)

Global Flags:
      --config string     Path to config YAML
      --loglevel string   Log output level (ERROR, INFO, DEBUG) (default "INFO")
```






