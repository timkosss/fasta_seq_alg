#!/bin/bash

#SBATCH -J create_dicts
#SBATCH -o /wistar/auslander/Timothy/microproteins_timka/outfiles/dicts_output_%A.txt
#SBATCH -e /wistar/auslander/Timothy/microproteins_timka/outfiles/dicts_output_%A.err
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00


DB1=""
K=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --db1)
      DB1="$2"
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

echo "DB1 = $DB1"
echo "K   = $K"

echo "Creating Dictionaries"
python /wistar/auslander/Timothy/fasta_seq_alg/main_run/create_dicts.py "$DB1" "$K"