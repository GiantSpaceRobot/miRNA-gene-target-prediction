# Predict miRNA gene targets

## Quickstart
* Run TargetScan, miRDB, and miRWalk on your favourite miRNA
* Run the following:
```
git clone https://github.com/GiantSpaceRobot/miRNA-gene-target-prediction.git
cd miRNA-gene-target-prediction
Rscript miRNA-gene-target-prediction-analysis.R TargetScan7.2.txt miRDB.tsv miRWalk_miRNA_Targets.csv My-Output yes 
```

### Parameters
* Argument 1: Targetscan txt output
* Argument 2: miRDB text output (formatted as .tsv)
* Argument 3: miRWalk csv output
* Argument 4: Output prefix
* Argument 5: Use all predicted genes for GO? {yes/no}

Note: Selecting "no" for argument 5 will instead use the genes predicted by at least two of the methods.

## Common miRNA gene targets
Another script is included in this repo for determining the common gene targets for three separate miRNAs.
In order to run this script, create a list of miRNA gene targets for each miRNA (one gene name per line).
Run the script as follows:
```
Rscript Common-Gene-Targets.R miRNA-1_all-predicted-gene-targets.tsv miRNA-2_all-predicted-gene-targets.tsv miRNA-3_all-predicted-gene-targets.tsv Common-gene-targets
```

### Parameters
* Argument 1: miRNA #1 gene names
* Argument 2: miRNA #2 gene names
* Argument 3: miRNA #3 gene names
* Argument 4: Output prefix

## Contributors
* Paul Donovan

## License
This project is licensed under the MIT License.

