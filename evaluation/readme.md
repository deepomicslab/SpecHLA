## Evaluation
The scripts to evaluate the typing results are in this folder.

| Script | Description |
| --- | --- |
|cal_seq_accuracy.py | Evaluate the reconstructed HLA sequences |
|cal_resource.py | Calculate the computational resources  |
|assess_read_assign.py| Evaluate the read-binning accuracy in simulated data |
|eva_type_accuracy.py| Calculate G-group typing accuracy in real samples |
|visualize_trio/*| Generate the data format for visulization of trio-consistency |
|compare.1k.pl |Compare typing results in 1000G samples|
|get_truth_HSVC2_g.py| Infer G group HLA types from HGSVC2 samples as ground truth|
|get_HLA_alleles_from_assembly.py| Extract HLA alleles from HGSVC2 samples as ground truth|
|download_hgsvc_assemly.py| Download phased assemblies of HGSVC2 samples as ground truth |
|convert.allele.pl | map the HLA types to the most recent IMGT version|