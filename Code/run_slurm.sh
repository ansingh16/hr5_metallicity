#!/bin/bash

source /home/ankitsingh/miniconda3/bin/activate ICL

tasks_in_total=$(wc -l < /home/ankitsingh/hr5_metalicity/to_be_processed.txt)
numfiles=1
min_jobs_running=10

# Function to count the number of running Slurm jobs for the user
count_running_jobs() {
    squeue -u $(whoami) -h -o %j | grep "metal_" | wc -l
}

while [[ $numfiles -le tasks_in_total ]]; do

    # Count the number of currently running jobs
    numpr=$(count_running_jobs)
    
    if [ $numpr -lt $min_jobs_running ]; then
        # Read a line from the input file
        read -r snapdat
        
        # Extract the basename without extension
        snap=$(basename "${snapdat%.*}")
        echo "Submitting job for: $snap"

        # Create a temporary Slurm script
        temp_script=$(mktemp)
        cat <<EOT > $temp_script
#!/bin/bash
#SBATCH --job-name=metal_${snap}
#SBATCH --output=slurm_${snap}-%j.out
#SBATCH --error=slurm_${snap}-%j.err
#SBATCH --ntasks=1

source /home/ankitsingh/miniconda3/bin/activate ICL
python /home/ankitsingh/hr5_metalicity/Code/All_data_gen.py $snap
EOT

        # Submit the job and increment the counter
        sbatch $temp_script
        numfiles=$((numfiles + 1))

        # Clean up the temporary script
        rm $temp_script

    
    fi

    # Sleep for a short time before checking again
    sleep 5
done < /home/ankitsingh/hr5_metalicity/to_be_processed.txt
