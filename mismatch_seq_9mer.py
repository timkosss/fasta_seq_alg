import pandas as pd
import os
import pickle
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



# def dict_aa(fasta, k):
#     aa_dict = {}
#     print("Dictionary In Progress")
#     for seq_id, seq in fastaReader(fasta):
#         for i in range(len(seq) - k + 1):
#             aa = seq[i:i+k]
#             if aa not in aa_dict:
#                 aa_dict[aa] = []
#             aa_dict[aa].append(seq_id)
#     print(f"Dictionary of {len(aa_dict)} Created")
#     return aa_dict

def dict_aa_9(index, k):
    pos_dicts = []
    for i in range(k): pos_dicts.append({})
    length = len(index)
    percent = 0
    for aa in index:
        percent+=1
        for i in range(k):
            mismatch = aa[:i] + aa[i+1:]
            if mismatch not in pos_dicts[i]:
                pos_dicts[i][mismatch] = []
            pos_dicts[i][mismatch].append(aa)
        print(f"we are {percent/length}% of the way done creating the mega dictionary")
    return pos_dicts

def process_sequence(seq_id_b, seq_b, k, super_dict, aa_index):
    matches = []
    for i in range(len(seq_b) - k + 1):
        aa = seq_b[i:i+k]
        for j in range(k):
            aa_9 = aa[:j] + aa[j+1:]
            if aa_9 in super_dict[j]:
                for aa_a in super_dict[j][aa_9]:
                    for seq_id_a in aa_index[aa_a]:
                        matches.append((seq_id_a,aa_a, seq_id_b, aa))    
    return matches

def find_matches(fasta_b_path, k, super_dict, aa_index):
    try:
        matches = []
        print("Searching for matches")
        i = 0
        for seq_id_b, seq_b in fastaReader(fasta_b_path):
            matches.extend(process_sequence(seq_id_b, seq_b, k, super_dict, aa_index))
            i+=1
            print(f"ONE DB2 DOWN!!!! WE HAVE GONE THROUGH {i} of DATABASE 2's SEQUENCES")
        return matches
    except KeyboardInterrupt:
        return matches
    


k = 10
fasta_a = "/wistar/auslander/Bryant/projects/pdx_protein_seq/uniprotkb_AND_model_organism_9606_2024_05_06.fasta"
fasta_b = "/wistar/auslander/prot/try0/refseq_db/bac_ref_seq.txt"
#fasta_b = "refseq_small.txt"

print("Reading FIRST dict ...")
with open('firstDict.pkl', 'rb') as f:
    aa_index = pickle.load(f)    
print("read in index")

# dict9 = dict_aa_9(aa_index, k)

# print("starting superDict Dump")
# with open('superDict.pkl', 'wb') as f:
#     pickle.dump(dict9, f)
# print("Dict Dumped ... Staring to read the dict...")

with open('superDict.pkl', 'rb') as f:
    super_dict = pickle.load(f)
print("read in super dictionary")

matches = find_matches(fasta_b, k, super_dict, aa_index)

df_matches = pd.DataFrame(matches, columns=['seqA_id', 'aa_a', 'seqB_id', 'aa_b'])
df_matches = df_matches.drop_duplicates()
print(f"there are {len(df_matches)} entries")

output_csv = "aa_matches_mismatch.csv"
df_matches.to_csv(output_csv, index=False)
print("Saved and FInished")