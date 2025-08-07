import pandas as pd
import os
import math
os.chdir('/wistar/auslander/Timothy/fasta_seq_alg')

#do it line by line so that with big databases you can get results constantly
        
def fastaReader(file):
    print("called fasta to dict")
    with open(file) as f:
        print("openned file")
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                id = line
            elif id:
                yield(id, line)
                id = None
        print("FINISHED")



def dict_aa(fasta, k):
    aa_dict = {}
    print("Dictionary In Progress")
    for seq_id, seq in fastaReader(fasta):
        for i in range(len(seq) - k + 1):
            aa = seq[i:i+k]
            if aa not in aa_dict:
                aa_dict[aa] = []
            aa_dict[aa].append(seq_id)
    print(f"Dictionary of {len(aa_dict)} Created")
    return aa_dict

def process_sequence(seq_id_b, seq_b, k, aa_index):
    matches = []
    for i in range(len(seq_b) - k + 1):
        aa = seq_b[i:i+k]
        if aa in aa_index:
            for seq_id_a in aa_index[aa]:
                
                matches.append((seq_id_a, seq_id_b, aa))    
    return matches

def find_matches(fasta_b_path, k, aa_index):
    try:
        matches = []
        print("Searching for matches")
        i = 0
        for seq_id_b, seq_b in fastaReader(fasta_b_path):
            matches.extend(process_sequence(seq_id_b, seq_b, k, aa_index))
            i+=1
            print(f"ONE DB2 DOWN!!!! WE HAVE GONE THROUGH {i} of DATABASE 2's SEQUENCES")
        return matches
    except KeyboardInterrupt:
        return matches
    


k = 10
fasta_a = "/wistar/auslander/Bryant/projects/pdx_protein_seq/uniprotkb_AND_model_organism_9606_2024_05_06.fasta"
#fasta_b = "/wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt"
fasta_b = "refseq_small.txt"

aa_index = dict_aa(fasta_a, k)

matches = find_matches(fasta_b, k, aa_index)

df_matches = pd.DataFrame(matches, columns=['seqA_id', 'seqB_id', 'matching_aa'])
df_matches = df_matches.drop_duplicates()
print(f"there are {len(df_matches)} entries")

output_csv = "aa_matches.csv"
df_matches.to_csv(output_csv, index=False)
print("Saved and FInished")