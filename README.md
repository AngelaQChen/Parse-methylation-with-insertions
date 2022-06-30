# Parse methylation with insertions
## Parse CpG methylation from BAM files.

Existing tools CpGtools (PacBio) and modbam2bed (Nanopore) does not include CpG methylation within insertion sites in their output files. This script aims to fill that gap by parsing and summarizing all CpG methylation sites found on each read in the BAM file within the MM and ML tags. Insertion sites will be represented by insertion_start_site+num_bases_from_start_of_insertion (eg. 545253+24 where 545253 is the reference position of the insertion and the methylation position is 24 bases from the start of the insertion). Methylation on different haplotypes are separated into different files.

Requires: Python3 and pysam

```
usage: get_methylation.py [-h] [-b FILE] [-r STR] [-o STR]

Parse sorted, phased and indexed BAM file for methylation modification in specified region.

optional arguments:
  -h, --help            show this help message and exit
  -b FILE, --bam FILE   indexed and sorted BAM file with methylation data from primrose
  -r STR, --region STR  CHROM:START-END format for region
  -o STR, --out STR     prefix of the output files [output]
```

Output columns:
1. chrom
2. start
3. ref_pos
4. read_coverage
5. num_mod_calls
6. num_nonzero_mod_calls
7. num_high_qual_mod
8. avg_qual
9. avg_qual_percent
10. percent_high_qual_calls
11. mod_qual


Pleas note that the output is not sorted by position. This script should also be applicable for CpG methylation data obtained from Nanopore but it has not been tested.
