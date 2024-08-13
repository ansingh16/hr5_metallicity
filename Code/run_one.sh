#!/bin/bash

source ~/miniconda3/bin/activate ICL
source progress.sh




source ~/miniconda3/bin/activate ICL
source progress.sh

n=0
j=0

tasks_in_total=$(wc -l < ../to_be_processed.txt)

# for files in to_be_processed.txt run python script
for snapdat in $(cat ../to_be_processed.txt); do
    echo $snapdat
    snap=$(basename "${snapdat%.*}")
    # Replace this with your processing code or command
    python ~/hr5_metallicity/Code/All_data_gen.py $snap &
    PID=$!
    wait $PID
    j=$((j+1))
    show_progress $j $tasks_in_total
    echo " "
done
