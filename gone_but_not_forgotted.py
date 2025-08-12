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

k=10
fasta = "/path/to/fasta/file/1"
aa_index = dict_aa(fasta, k)
dict9 = dict_aa_9(aa_index, k)
