#!/bin/bash

#SBATCH -J mismatch_process
#SBATCH -o /wistar/auslander/Timothy/fasta_seq_alg/OUTFIles/output_%A_%a.txt
#SBATCH -e /wistar/auslander/Timothy/fasta_seq_alg/OUTFIles/output_%A_%a.err
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --array=1-500

echo "Running On:"
srun hostname
srun uname -a

LINES_PER_JOB=1000000
START_LINE=$(( (SLURM_ARRAY_TASK_ID - 1) * LINES_PER_JOB + 1 ))
END_LINE=$(( START_LINE + LINES_PER_JOB - 1 ))

TOTAL_LINES=$(wc -l < "/wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt")

# If next chunk is staring out of bounds, exit
if [ "$START_LINE" -gt "$TOTAL_LINES" ]; then
    echo "Job $SLURM_ARRAY_TASK_ID is out of bounds, Exiting."
    scancel ${SLURM_ARRAY_JOB_ID}
    exit 0
fi

# If end lin is out of bounds, make it the last available line
if [ "$END_LINE" -gt "$TOTAL_LINES" ]; then
    END_LINE=$TOTAL_LINES
fi

sed -n "${START_LINE}, ${END_LINE}p" /wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt | python /wistar/auslander/Timothy/fasta_seq_alg/ex_seq_better.py


#################################
#           end   script        #
#################################