#!/bin/bash

source ~/miniconda3/bin/activate ICL
source progress.sh



n=0
j=0

tasks_in_total=$(cat ../to_be_processed.txt |wc -l)

numfiles=1
while [[ $numfiles -le tasks_in_total ]]; do

        
    numpr=$(ps -u ankitsingh -o comm= | grep "^python$" | wc -l
    )

    
    if [ $numpr -lt 4 ]; then

    read snapdat
    
    echo $snapdat
    snap=$(basename "${snapdat%.*}")
    # Replace this with your processing code or command
    time python ~/hr5_metallicity/Code/All_data_gen.py $snap &
    numfiles=$((numfiles+1))
    show_progress $numfiles $tasks_in_total
    echo " "
    fi

done < ../to_be_processed.txt

