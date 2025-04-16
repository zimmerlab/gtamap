# GTA_Map


### Second Pass
During the initial pass, some reads may not be mappable with sufficient confidence. 
These reads are stored and re-analyzed in a second pass.

Reasons to be not mappable in the first pass:
- Maps to multiple locations within the same gene with sufficient confidence
- Bounds of kmers / read are not uniquely mappable (mismatches, indels, etc)

## Deletions and Introns

Deletions and introns are detected in a similar fashion as they both represent a
continous read sequence and a discontinuous reference sequence (gap in reference).

A gap is categorized as a deletion if the number of skipped bases is less than the 
configured number of intron bases.

### Left Normalization

There exists a decision problem when there are bases that can be assigned to both sides of the gap.

[!left-norm-1](resources/readme/left-norm.1.png)

In this example, the sequenced individual has a 9bp deletion.

The read sequence is `ATACAGT` which leaves the `C` base to be assigned to both sides of the gap.

The gap could either be `CGGAAAGGT` or `GGAAAGGTC`.

[!left-norm-2](resources/readme/left-norm.2.png)

Left normalization uses the leftmost beginning for the gap.