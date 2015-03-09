# Unnamed Biomarker Project

This is the developement repository for a biomarker selection method that currently has no name. It's code name is "modified lefse". A final "report" will be generated using R markdown - eventually...

## Contents

| Item                         | Type      | Description                                                                                                                       |
|------------------------------|-----------|-----------------------------------------------------------------------------------------------------------------------------------|
| Data                         | dir       | Contains the biome file and other data. Small enough to go onto github                                                            |
| BiomarkerSelectionMethods.R  | R script  | Main script.  Comparing to Schubert paper + selecting biomarkers with "modified lefse" and compares to schubert                   |
| SupportingBiomarkerMethods.R | R script  | The methodology guts. Contains the biomarker selection glmnet wrapper + also useful methods to reduce clutter in the main script. |
| UnnameBiomarkerProject.Rproj | R project | R project file 
| Reports | dir | contains abstracts and small summary reports of the findings|
|


## Status
Still haven't finished some analysis steps. Remaining to be done are the following:

1. Run binary response
2. Run continous response (using OTU 19 as proxy for amount of C.diff)
3. Get the AUC, compared to base model and combined model.

