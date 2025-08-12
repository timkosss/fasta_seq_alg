#!/bin/bash

#SBATCH -J make_big_csv
#SBATCH -o /wistar/auslander/Timothy/fasta_seq_alg/OUTFIles/output.txt
#SBATCH -e /wistar/auslander/Timothy/fasta_seq_alg/OUTFIles/output.err
#SBATCH --mem=32G
#SBATCH --ntasks=1

echo "Running On:"
srun hostname
srun uname -a

echo "Starting Python Script..."
python -u /wistar/auslander/Timothy/fasta_seq_alg/combine_csv.py 

#################################
#           end   script        #
#################################
