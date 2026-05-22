#!/bin/bash

# This script contains common functions and trap settings for job handling.
# It is intended to be sourced by other main scripts.

# A variable to hold the log file path, must be defined in the main script
# before calling any functions from this utility.
# Example: logfile=${LOG_PATH}/${TIME}_${SAMPLE}.log

# Function to be executed at the start of a job.
start_job() {
    exec > "$logfile" 2>&1

    # Check if the script is running in a Slurm environment.
    if [ -n "$SLURM_JOB_ID" ]; then
        echo "*-----------------------------*"
        echo "$(date '+%Y-%m-%d %H:%M:%S') Job started"
        echo "SLURM JOB ID = $SLURM_JOB_ID"
        requested_mem_gb=$(( $SLURM_MEM_PER_NODE / 1024 ))
        echo "Requested memory: ${requested_mem_gb} GB"
        echo -e "*-----------------------------*\n"
    else
        echo "*-----------------------------*"
        echo "$(date '+%Y-%m-%d %H:%M:%S') Job started (Not in Slurm environment)"
        echo -e "*-----------------------------*\n"
    fi
    
}

# Function to handle job cancellation.
cancel_handler() {
    echo "=================="
    echo "Job was cancelled."
    echo "=================="
    exit 130 
}

# Trap SIGINT and SIGTERM signals to execute the cancel_handler function on job cancellation.
trap cancel_handler SIGINT SIGTERM

# Function to be executed at the end of the script.
function finish {
    EXIT_CODE=$?

    # Calculate total runtime.
    duration=$SECONDS
    hours=$((duration / 3600))
    minutes=$(( (duration % 3600) / 60 ))
    seconds=$(( duration % 60 ))

    # Determine job status based on the exit code (success, failure, or cancellation).
    if [ $EXIT_CODE -eq 0 ]; then
        echo -e "\n*------------------------------*"
        echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished"
        echo "*------------------------------*"
        echo -e "\nTotal runtime: ${hours} hours, ${minutes} minutes, ${seconds} seconds\n"
        echo "==========================="
        echo "Job completed successfully."
        echo "==========================="
    elif [ $EXIT_CODE -eq 130 ]; then
        echo -e "\nTotal runtime: ${hours} hours, ${minutes} minutes, ${seconds} seconds\n"
        echo "=================="
        echo "Job was cancelled."
        echo "=================="
    else
        echo -e "\nTotal runtime: ${hours} hours, ${minutes} minutes, ${seconds} seconds\n"
        echo "================================="
        echo "Job failed with exit code: $EXIT_CODE"
        echo "================================="
    fi
}

# Set a trap to execute the finish function upon script exit.
trap finish EXIT

