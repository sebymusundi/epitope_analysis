## Prediction of conserved cross-species B-cell linear epitopes

This repository contains information on the steps followed to identify conserved cross species b-cell linear epitopes associated with human malaria

### Requirements

Scripts included in this analysis were run on

- Bash
  - Conda environment instaled.
- R
- Nextflow
  - Conda and Docker used to run different programs
  - Whole genome sequencing data from field isolates used to analyze the conservation of individual predicted antigens containing epitopes of interest.
- Prediction tools including Deeploc, WoLF PSORT, Phobius, Bepipred 2.0 and Deeploc

### Steps

#### Downloading Plasmodium data and clustering redundant proteins

Plasmodium proteins from five species associated with human malaria were downloaded and CD-HIT was used to remove redundant proteins with a threshold of 90%.

#### Subcellular prediction

WOLF PSORT and Deeploc were used to predict the subcellular location of the proteins using P. falciparum as the reference organism. Intracellular proteins were also eliminated using their protein names.

#### Surfaceome prediction

The presence of a signal peptide and transmembrane helice predicted by using Phobius. Proteins with one signal peptide or one transmembrane helix were selected

#### Merozoite protein selection

We retrieved mass spectrometry data on proteins located on the late schizont-merozoite stage and cross selected against proteins located on the surface of the cell

#### B-cell linear epitope prediction

We utilized Bepipred 2.0 and Epidope to predict epitopes. Bepipred 2.0 is a tool used to predict the presence of b-cell linear epitopes. We then made an R script

Currently, Bepipred 2.0 is hosted at DTU University found here https://services.healthtech.dtu.dk/services/BepiPred-2.0/

Bepipred accepts a maximum of 50 protein sequences and 300,000 amino acids per submission with the length ranging from 10-300000. The output for multiple protein sequences can be downloaded in the form of a CSV file. However, the output contains raw b-cell epitope sequences and takes time to process especially on epitopes that meet the defined threshold of 0.5. Herein I tried to extract the exposed epitopes that met the threshold score of 0.5. The final output file is in the form of a table containing
the protein_ID, start and end position, epitope sequence, and length. The first trial is run with 50 protein sequences.
