# Drug Interactions ANAlyzer (DIANA)

2017 Joaquim Aguirre-Plans

Structural Bioinformatics Laboratory

Universitat Pompeu Fabra



Drug Interactions ANAlyzer (DIANA) is a Python (v3) software which permits the prediction of drug combination using drug target and protein-protein interaction information

## Getting Started

### Summary

* **scripts**: Folder where the scripts to run the program are stored 
* **scripts/generate_profiles.py**: Script to generate the profiles of a drug 
* **scripts/compare_profiles.py**: Script to compare the profiles of two drugs 
* **scripts/download_go_files.py**: Script to download Gene Ontology files necessary for the analysis 
* **diana**: Folder where the source code is stored
* **workspace**: Folder where the results of the analysis are stored
* **workspace/profiles**: Folder containing the profiles generated
* **workspace/comparisons**: Folder containing the comparisons of profiles
* **config.ini**: Configuration file, containing the default parameters to run DIANA

### Prepare the configuration file

Before starting running DIANA, we need to modify the parameters in the configuration file as we want.

The configuration file is in the main directory, as **config.ini**

Here we can modify the parameters of BIANA database if we have it.

If we do not have BIANA, we can set use_cluster as True and use the downloaded data.

We also have to introduce the path to Guild executable. In case that you do not have it, you can download it and install it from https://github.com/emreg00/guild

### Execution

#### Download the GO files:

First, make sure that you have downloaded the Gene Ontology files.
To do it, run the following script:

```
python path/to/diana/scripts/download_go_files.py
```



#### Generate the profiles of a drug:

Now, create the profiles of the drugs that you want to compare.
To create the profile of one drug, run the following script:

```
python path/to/diana/scripts/generate_profiles.py 
  -j <job id> # Give an specific job identifier. If not, it will be created automatically
  -d <drug name>
  -t <drug targets (in Entrez geneID)>
  -pt <type IDs of the proteins (default is geneid)>
  -sif <network (in SIF format and geneID)>
  -th <thresholds list> # Cut-offs that will define the profiles generated
  -ws <workspace> # Directory where the results will be saved
```

Example of generate the profile having network and the targets:

```
python path/to/diana/scripts/generate_profiles.py -j metformin_with_targets -d metformin -t path/to/diana/workspace/targets/metformin.targets -sif path/to/diana/workspace/sif/human_eAFF_geneid_2017.sif
```

Example of generate the profile having network but not having targets.
The targets will be searched in BIANA.

```
python path/to/diana/scripts/generate_profiles.py -j metformin -d metformin -sif path/to/diana/workspace/sif/human_eAFF_geneid_2017.sif
```



#### Compare the profiles of a pair of drugs:

With this script we can compare the profiles of a pair of drugs:

```
python path/to/diana/scripts/compare_profiles.py 
  -j1 <job 1>
  -j2 <job 2>
  -sif <network (in SIF format and geneID)>
  -th <thresholds list> # Cut-offs that will define the profiles generated
  -ws <workspace> # Directory where the results will be saved
```

In the following example, we compare the profiles created for metformin with the profiles created for haloperidol:

```
python path/to/diana/scripts/compare_profiles.py -j1 metformin -j2 haloperidol -sif path/to/diana/workspace/sif/human_eAFF_geneid_2017.sif
```


#### Analyze the results:

There is a pipeline on how to analyze results in README_analysis_one_target.


