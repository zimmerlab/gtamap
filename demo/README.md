# Demo

Here we simulated all transcripts of TP53 (1000 per transcript) and additional 200 random transcripts from random genes (100 reads each)
using the `readSimulator.jar` (1% mutation rate).

```sh
java -jar readSimulator.jar
     -length 150
     -frlength 350
     -SD 80
     -mutationrate 1.0
     -seqerrrate 0.0
     -gtf /path/to/gencode.v48.annotation.gtf
     -fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa
     -fidx /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
     -readcounts ./data/readcounts.simulation
     -od ./data
```


## Usage

```sh
make clean
make all
```
