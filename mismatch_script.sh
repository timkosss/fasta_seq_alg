#!/bin/bash

#SBATCH -J mismatch_process
#SBATCH -o /wistar/auslander/Timothy/fasta_seq_alg/not_exact_matches/OUTFIles/output_%A_%a.txt
#SBATCH -e /wistar/auslander/Timothy/fasta_seq_alg/not_exact_matches/OUTFIles/output_%A_%a.err
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --array=1-25

echo "Running On:"
srun hostname
srun uname -a

LINES_PER_JOB=1000000
START_LINE=$(( (SLURM_ARRAY_TASK_ID - 1) * LINES_PER_JOB + 1 ))
END_LINE=$(( START_LINE + LINES_PER_JOB - 1 ))

sed -n "${START_LINE}, ${END_LINE}p" /wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt | python /wistar/auslander/Timothy/fasta_seq_alg/not_exact_matches/mismatch_seq_9mer.py


#################################
#           end   script        #
#################################