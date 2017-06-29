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

Before executing the program, download the following files:

```
wget http://geneontology.org/ontology/go-basic.obo
wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
```

Change the paths in line 1009-1010 of diana_v3.5.local.py.

#### Create the profiles of a drug:

```
python diana_v3.5.local.py 
  -d1 <drug name>
  -s1 <drug targets>
  -sif <network>
  -tis <type of network codes>
  -ex <thresholds list> # Cut-offs that will define the profiles 
  -db <database>
  -prof # If introduced, it generates the profiles
  -comp # If introduced, it performs the comparison between profiles of the two drugs
```
python diana_v3.5.local.py -d1 oltipraz -s1 seeds/oltipraz.seeds -sif sif/human_eAFF_geneid_2017.sif -tis geneid -ex small_analysis_thresholds.list -db BIANA_JAN_2017 -prof

Example:

```
python diana_v3.5.local.py 
  -d1 oltipraz
  -s1 seeds/oltipraz.seeds
  -sif sif/human_eAFF_geneid_2017.sif
  -tis geneid
  -ex top_thresholds.list
  -db BIANA_JAN_2017
  -prof # If introduced, it creates the profiles of the drug
```

#### Compare the profiles of a pair of drugs:

```
python diana_v3.5.local.py 
  -d1 <drug name 1>
  -d2 <drug name 2>
  -sif <network>
  -tis <type of network codes>
  -ex <thresholds list> # Cut-offs that will define the profiles 
  -db <database>
  -comp # If introduced, it performs the comparison between profiles of the two drugs
```

Example:

```
python diana_v3.5.local.py 
  -d1 oltipraz
  -d2 drug2
  -sif sif/human_eAFF_geneid_2017.sif
  -tis geneid
  -ex top_thresholds.list
  -db BIANA_JAN_2017
  -comp
```

