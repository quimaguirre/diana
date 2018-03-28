# Drug Interactions ANAlyzer (DIANA)

2017 Joaquim Aguirre-Plans

Structural Bioinformatics Laboratory

Universitat Pompeu Fabra



Drug Interactions ANAlyzer (DIANA) is a Python software which permits the prediction of drug combination using drug target and protein-protein interaction information

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
  -d <drug name>
  -t <drug targets (in Entrez geneID)>
  -pt <type IDs of the proteins (default is geneid)>
  -sif <network (in SIF format and geneID)>
  -th <thresholds list> # Cut-offs that will define the profiles generated
  -rad <radius> # The radius to generate the network of expansion if no network is introduced
  -tax <taxid> # The Taxonomy ID restriction to generate the network of expansion if no network is introduced
  -ws <workspace> # Directory where the results will be saved
```

Example of generate the profile having network and the targets:

```
python path/to/diana/scripts/generate_profiles.py -d 'metformin' -t path/to/diana/workspace/targets/metformin.targets -sif path/to/diana/workspace/sif/human_eAFF_geneid_2017.sif
```

Example of generate the profile having network but not having targets.
The targets will be searched in BIANA.

```
python path/to/diana/scripts/generate_profiles.py -d 'metformin' -sif path/to/diana/workspace/sif/human_eAFF_geneid_2017.sif
```

Example of generate the profile not having neither the network nor the targets.
The targets will be searched in BIANA.
The program will create a network of expansion using the targets as seeds and the options defined by the user as restrictions.
In this case, we restrict the network to only human interactions found at least by Yeast to Hybrid methods, and not further than 3 levels from the targets.

```
python path/to/diana/scripts/generate_profiles.py -d 'metformin' -rad 3 -tax 9606 -res Y2H
```



#### Compare the profiles of a pair of drugs:

With this script we can compare the profiles of a pair of drugs:

```
python path/to/diana/scripts/compare_profiles.py 
  -d1 <drug name 1>
  -d2 <drug name 2>
  -t1 <targets from drug 1 (in Entrez geneID)>
  -t2 <targets from drug 2 (in Entrez geneID)>
  -pt <type IDs of the proteins (default is geneid)>
  -sif <network (in SIF format and geneID)>
  -th <thresholds list> # Cut-offs that will define the profiles generated
  -ws <workspace> # Directory where the results will be saved
```

In the following example, we compare the profiles created for metformin using the targets in 'metformin.targets', with the profiles created for haloperidol using the targets from BIANA:

```
python path/to/diana/scripts/compare_profiles.py -d1 'metformin' -t1 path/to/diana/workspace/targets/metformin.targets -d2 'haloperidol' -sif path/to/diana/workspace/sif/human_eAFF_geneid_2017.sif
```



#### Analyze the results:

There is a pipeline on how to analyze results in README_analysis_one_target.


