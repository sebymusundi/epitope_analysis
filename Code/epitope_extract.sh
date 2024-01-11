#! /bin/bash 

VAR=$(cat  final_predicted_epitopes.tsv)



for i in ${VAR}
do 
grep ${i}  non_falciparum_vs_epitopes.tsv 
done 
