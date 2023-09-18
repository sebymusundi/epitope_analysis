## Bepipred extractor in RStudio 
Bepipred 2.0 is a tool used to predict the presence of b-cell linear epitopes. 

Currently Bepipred 2.0 is hosted at DTU university found here https://services.healthtech.dtu.dk/services/BepiPred-2.0/

Bepipred accepts a maximum of 50 protein sequences and  300,000 amino acids per submission with the length ranging from 10-300000.  The output for multiple protein sequences 
can be downloaded in form of a csv file. However, the output contains raw b-cell epitope sequences and takes time to process especially on epitopes that meet the defined 
theshold of 0.5. Herein I tried to extract the exposed epitopes that met the theshold score of 0.5. The final output file is in form of a table containing 
the protein_ID, start and end position, epitope sequence and length.  The first trial is run with 50 protein sequences. 
