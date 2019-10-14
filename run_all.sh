#!/bin/bash

classes=(0-10 10-30 30-50 50-90)
data_loc=$1
evt_tag=$2
theory_tag=$3
cent=$4
nev=$5
hard_type=$6

if [[ $cent == "true" ]]
then
   for class in "${classes[@]}"
   do
      mkdir -p output/${theory_tag,}/nev${nev}/${evt_tag}/${class}
      mkdir -p log/${theory_tag,}/nev${nev}/${evt_tag}/${class}
      
      ./run.sh ${data_loc}/${class} ${theory_tag,}/nev${nev}/${evt_tag}/${class}/${evt_tag}-${class} ${nev} ${hard_type} ${theory_tag^}
   done
else
   mkdir -p output/${theory_tag,}/nev${nev}/${evt_tag}
   mkdir -p log/${theory_tag,}/nev${nev}/${evt_tag}
   
   ./run.sh ${data_loc} ${theory_tag,}/nev${nev}/${evt_tag}/${evt_tag} ${nev} ${hard_type} ${theory_tag^}
fi
