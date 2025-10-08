#!/bin/bash

# Default for optional value
threads=8

# Check if the script is being run with bash
if [ -z "$BASH_VERSION" ]; then
  echo "❌ This script must be run with bash, not sh."
  echo "run the script with bash run_RCTD.sh or ./run_RCTD.sh"
  exit 1
fi

while getopts i:o:r:t: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
        r) reference=${OPTARG};;
        t) threads=${OPTARG};;
    esac
done

# Check required arguments
if [[ -z "$input" || -z "$output" || -z "$reference" ]]; then
  echo "missing required input arguments. Run Rscript RCTD.R -h for a list of required input arguments"
  echo "See usage below:"
  echo ""
  Rscript RCTD.R -h
  exit 1
fi

# Set thread-related env vars
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run Rscript
Rscript RCTD.R -i "$input" -o "$output" -r "$reference" -t "$threads"

