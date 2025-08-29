import pandas as pd
import os
import sys
import pickle
from tqdm import tqdm # type: ignore

#================Configuration=====================#
os.chdir('/wistar/auslander/Timothy/fasta_seq_alg')

kmer_length = int(sys.argv[2])
database_2_filepath = str(sys.argv[1])
# fasta_a_filepath = sys.argv[1]
# fasta_b_filepath = sys.argv[2]
# fasta_b_filepath = "refseq_small.txt"
# fasta_b_filepath = "/wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt"

output_csv = (f"aa_matches.csv")

#================FASTA Reader=====================#
def fastaReader(filepath):
    
    # Generator that yields (sequence_id, sequence) from a FASTA file.
    
    print(f"Reading FASTA: {filepath}")
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                id = line
            elif id:
                yield(id, line)
                id = None
    print("Finished reading FASTA")

#================Kmer Dictionary=====================#
def build_kmer_dict(fasta_path, kmer_length):
    
    # Builds a dictionary mapping each unique kmer in FASTA A to the IDs of the sequences that contain them.
    # Ex. {ABCDEFGHIJ: [Protein1, Protein2]}
    
    kmer_dict = {}
    print("Building kmer dictionary...")
    for seq_id, seq in fastaReader(fasta_path):
        for i in range(len(seq) - kmer_length + 1):
            kmer = seq[i:i+kmer_length]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = []
            kmer_dict[kmer].append(seq_id)
    print(f"Kmer dictionary built with {len(kmer_dict)} entries.")
    return kmer_dict

#================Match Finder=====================#
def process_sequence(seqB_id, seq_b, kmer_length, kmer_dict):
    
    # Process a single sequence from FASTA B, checking for exact kmer matches in kmer_dict.
    # Returns a list of (seqA_id, seqB_id, matching_kmer)
    
    matches = []
    for i in range(len(seq_b) - kmer_length + 1):
        kmer = seq_b[i:i+kmer_length]
        if kmer in kmer_dict:
            for seqA_id in kmer_dict[kmer]:
                matches.append((seqA_id, seqB_id, kmer))
    return matches

def find_matches(kmer_length, kmer_dict):
    
    # Iterates through all the sequences in FASTA B and collects all matching kmers from FASTA A.
    
    try:
        matches = []
        print("Searching for matches...")
        for i, (seqB_id, seq_b) in enumerate(fastaReader(database_2_filepath), start= 1):
            matches.extend(process_sequence(seqB_id, seq_b, kmer_length, kmer_dict))
            print(f'{i} sequences from DataBase 2 completed.')
        return matches
    except KeyboardInterrupt:
        print("Process interrupted. Returning collected matches.")
        return matches

#================Main Execution=====================#
if __name__ == "__main__":
    with open('firstDict.pkl', 'rb') as f:
        kmer_dict = pickle.load(f)

    matches = find_matches(kmer_length, kmer_dict)

    df_matches = pd.DataFrame(matches, columns=['seqA_id', 'seqB_id', 'matching_kmer'])
    df_matches = df_matches.drop_duplicates()
    print(f"Found {len(df_matches)} unique matches.")

    df_matches.to_csv(output_csv, index=False)
    print(f"Matches saved to {output_csv}")
