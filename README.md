# Drug Interactions ANAlyzer (DIANA)

2017 Joaquim Aguirre-Plans

Structural Bioinformatics Laboratory

Universitat Pompeu Fabra



Drug Interactions ANAlyzer (DIANA) is a Python software which permits the prediction of drug combination using drug target and protein-protein interaction information

## Getting Started

### Summary

* **diana_v3.5.local.py**: Latest script version of DIANA to execute locally 
* **diana_v3.5.cluster.py**: Latest script version of DIANA to execute in a Slurm cluster
* **data/**: Folder where the data from the analysis is stored
* **results/**: Folder where the results of the analysis are stored
* **diana_analysis_results/**: Folder containing scripts to analyze the obtained results

### Execution

Execute DIANA locally:

```
python diana_v3.5.local.py 
  -d1 <drug name 1>
  -d2 <drug name 2>
  -sif <network>
  -tis <type of network codes>
  -ex <thresholds list>
  -db <database> # It may not be necessary if we have the Pickle files
  -prof # If introduced, it generates the profiles
  -comp # If introduced, it performs the comparison between profiles of the two drugs
```

Example:

```
python diana_v3.5.local.py 
  -d1 DCC0009
  -d2 DCC0001
  -sif  sif/human_eAFF_geneid_2017.sif
  -tis geneid
  -ex top_thresholds.list
  -db BIANA_JAN_2017
  -prof
```

