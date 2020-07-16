# Predict miRNA gene targets

## Quickstart
* Run the following:
```
git clone https://github.com/GiantSpaceRobot/miRNA-gene-target-prediction.git
cd miRNA-gene-target-prediction
Rscript miRNA-gene-target-prediction-analysis.R miR-335-5p/TargetScan7.2__miR-335-5p.predicted_targets.txt miR-335-5p/miRDB.tsv miR-335-5p/miRWalk_miRNA_Targets.csv miR-335-5p/Results/miR-335-5p yes 
```

## Parameters
* Argument 1: Targetscan txt output
* Argument 2: miRDB text output (formatted as .tsv)
* Argument 3: miRWalk csv output
* Argument 4: Output prefix
* Argument 5: Use all gene for GO? yes/no

Note: Selecting no will use all genes predicted by at least two of the methods

## Contributors
* Paul Donovan

## License
This project is licensed under the MIT License.

