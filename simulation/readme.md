## Simulation
The scripts to generated simulated data are in this folder.
To simulate paired-end data, run `simu.hla.sh` like:
```
Simulation HLA capture reads for HLA Typing
USAGE: <PATH-TO>/simu.hla.sh -n <sample_size> -r <read_length> -1 <depth of haplotype1> -2 <depth of haplotype2>

-n        : sample size of simulation [required]

-r        : read length [required]

-1        : depth of haplotype1 [required]

-2        : depth of haplotype2 [required]
```
To simulate different sequencing protocols, use `sim_diff_platforms.py`.
