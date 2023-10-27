#! /bin/bash 

VAR=$(cat final_cross_species_epitopes.tsv)



for i in ${VAR}
do 
grep ${i}  non_falciparum_vs_epitopes.tsv 
done 
