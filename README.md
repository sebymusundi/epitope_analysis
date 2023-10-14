## Epitope analysis 
This script provides information on how to analyze predicted b-cell linear epitopes from Plasmodium species.

This repository has been primarily analyzed using bash and Rstudio. Most of the bash tools have been run on the conda environment. 

Rstudio and more specifically, the ```tidyverse``` package had been primarily used to filter information from important web-browser tools.



### Bepipred extractor in RStudio 
Bepipred 2.0 is a tool used to predict the presence of b-cell linear epitopes. 

Currently, Bepipred 2.0 is hosted at DTU University found here https://services.healthtech.dtu.dk/services/BepiPred-2.0/

Bepipred accepts a maximum of 50 protein sequences and  300,000 amino acids per submission with the length ranging from 10-300000.  The output for multiple protein sequences 
can be downloaded in the form of a CSV file. However, the output contains raw b-cell epitope sequences and takes time to process especially on epitopes that meet the defined 
threshold of 0.5. Herein I tried to extract the exposed epitopes that met the threshold score of 0.5. The final output file is in the form of a table containing 
the protein_ID, start and end position, epitope sequence, and length.  The first trial is run with 50 protein sequences. 

