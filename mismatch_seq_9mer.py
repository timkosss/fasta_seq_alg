import pandas as pd
import os
import pickle
import time
from tqdm import tqdm # type: ignore

#================Configuration=====================#
os.chdir('/wistar/auslander/Timothy/fasta_seq_alg')
kmer_length = 10

fasta_a_filepath = "/wistar/auslander/Bryant/projects/pdx_protein_seq/uniprotkb_AND_model_organism_9606_2024_05_06.fasta"
fasta_b_filepath = "refseq_small.txt"
# fasta_b_filepath = "/wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt"

first_dictionary_path = 'firstDict.pkl' #Dictionary that maps all unique kmers from FASTA A to the protein names that contain the kmer 
#Ex. {ABCDEFGH: [Protein 1, Protein 2]}

expanded_dictionary_path = 'superDict.pkl' #A list containing K amount of dictionaries; Each individul dictionary takes a kmer from first dictionary and deletes the Kth character of the sequence and maps to original kmer. 
#Ex. {BCDEFGH: [ABCDEFGH, BBCDEFGH]}

output_csv = "aa_matches_mismatch.csv"

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

#================Match Finder=====================#
def process_sequence(seqB_id, seq_b, kmer_length, expanded_dicts, first_dict):
    
    #Process a single sequence from Database B, checking for matches in Expanded Dict.
    #Returns a list of (seqA_id, matched_kmer_from_A, seqB_id, matched_kmer_from_B)
    
    matches = []
    for i in range(len(seq_b) - kmer_length + 1):
        matched_kmer_from_B = seq_b[i:i+kmer_length] #Gets substring of length K from the sequence
        for dict_number in range(kmer_length):
            kmer_with_deletion = matched_kmer_from_B[:dict_number] + matched_kmer_from_B[dict_number+1:] #Deletes AA at index dict_number from substring
            if kmer_with_deletion in expanded_dicts[dict_number]: #Compares it to dictionary of corresponding number
                for matched_kmer_from_A in expanded_dicts[dict_number][kmer_with_deletion]: #Look at every mapped kmer
                    for seqA_id in first_dict[matched_kmer_from_A]: #And then every mapped Sequence name to that kmer
                        matches.append((seqA_id, matched_kmer_from_A, seqB_id, matched_kmer_from_B))    
    return matches

def find_matches(fasta_b_filepath, kmer_length, expanded_dicts, first_dict):
    
    # Iterate through all the Sequences in FASTA B and collect all matches from that sequence.
    
    try:
        matches = []
        print("Searching for matches...")
        for seqB_id, seq_b in tqdm(fastaReader(fasta_b_filepath), desc="Scanning DB2", unit="seq"):
            matches.extend(process_sequence(seqB_id, seq_b, kmer_length, expanded_dicts, first_dict))
        return matches
    except KeyboardInterrupt:
        print(f"Process interrupted. Returning collected matches.")
        return matches

#================Main Execution=====================#
if __name__ == "__main__":
    print(f"Loading pickled dictionaries...")

    with open(first_dictionary_path, 'rb') as f:
        first_dict = pickle.load(f)    
    print("Loaded First Dictionary")

    with open(expanded_dictionary_path, 'rb') as f:
        expanded_dicts = pickle.load(f)
    print("Loaded Expanded Dictionaries")

    start_time = time.time()
    matches = find_matches(fasta_b_filepath, kmer_length, expanded_dicts, first_dict)
    elapsed_time = time.time() - start_time
    print(f"Matching completed in {elapsed_time:.2f} seconds.")


    df_matches = pd.DataFrame(matches, columns=['seqA_id', 'matched_kmer_from_A', 'seqB_id', 'matched_kmer_from_B'])
    df_matches = df_matches.drop_duplicates()
    print(f"Found {len(df_matches)} unique matches.")

    df_matches.to_csv(output_csv, index=False)
    print(f"Matches saved to {output_csv}")
