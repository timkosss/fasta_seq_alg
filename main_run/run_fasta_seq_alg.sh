#!/bin/bash

#SBATCH -J fasta_seq_alg
#SBATCH -o /wistar/auslander/Timothy/microproteins_timka/outfiles/fasta_seq_%A.txt
#SBATCH -e /wistar/auslander/Timothy/microproteins_timka/outfiles/fasta_seq_%A.err
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00
#SBATCH --array=1-1

DB2=""
K=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --db2)
      DB2="$2"
      shift 2
      ;;
    --k)
      K="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "DB2 = $DB2"
echo "K   = $K"

python /wistar/auslander/Timothy/fasta_seq_alg/main_run/ex_seq_better.py "$DB2" "$K"
python /wistar/auslander/Timothy/fasta_seq_alg/main_run/mismatch_seq_9mer.py "$DB2" "$K"

