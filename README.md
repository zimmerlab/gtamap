# GTA_Map


### Second Pass
During the initial pass, some reads may not be mappable with sufficient confidence. 
These reads are stored and re-analyzed in a second pass.

Reasons to be not mappable in the first pass:
- Maps to multiple locations within the same gene with sufficient confidence
- Bounds of kmers / read are not uniquely mappable (mismatches, indels, etc)