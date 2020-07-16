# Predict miRNA gene targets

## Quickstart
* Run TargetScan, miRDB, and miRWalk on your favourite miRNA
* Run the following:
```
git clone https://github.com/GiantSpaceRobot/miRNA-gene-target-prediction.git
cd miRNA-gene-target-prediction
Rscript miRNA-gene-target-prediction-analysis.R TargetScan7.2.txt miRDB.tsv miRWalk_miRNA_Targets.csv My-Output yes 
```

## Parameters
* Argument 1: Targetscan txt output
* Argument 2: miRDB text output (formatted as .tsv)
* Argument 3: miRWalk csv output
* Argument 4: Output prefix
* Argument 5: Use all predicted genes for GO? {yes/no}

Note: Selecting "no" for argument 5 will instead use the genes predicted by at least two of the methods.

## Contributors
* Paul Donovan

## License
This project is licensed under the MIT License.

