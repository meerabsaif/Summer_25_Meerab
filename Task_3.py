import sys

def read_fasta(filename):
    seq_dict = {}
    key_id = None
    value_seq = ""
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()

                if not line:          
                    continue

                if line.startswith('>'):
                    if key_id is not None:
                        if valid_seq(value_seq):
                            seq_dict[key_id] = value_seq

                    key_id = line[1:]
                    value_seq = ""
                else:
                    value_seq += line.upper()

            if key_id is not None and valid_seq(value_seq):
                seq_dict[key_id] = value_seq
                
    except FileNotFoundError:
        print(f"File not found")
        sys.exit(1)
    
    return seq_dict

def valid_seq(seq):
    valid_nuc = {'A', 'T', 'C', 'G'}

    seq = seq.upper()
    for x in seq: 
         if x not in valid_nuc:
             return False 
    return True 


def filter_sequences(seq_dict, min_length):
    filtered = {}
    for yes, seq in seq_dict.items():
        if len(seq) >= min_length:
            filtered[yes] = seq
    return filtered


def write_fasta(seq_dict, output_file):
    try:
        with open(output_file, 'w') as wr:
            for id, seq in seq_dict.items():
                wr.write(f">{id}\n{seq}\n")
    except Exception as error:
        print(f"Error in writing: {error}")
        sys.exit(1)


def summary_stats(seqs):
    lengths = [len(seq) for seq in seqs.values()]
    avg_len = sum(lengths) / len(lengths) 
    max_len = max(lengths) 
    min_len = min(lengths) 
    return avg_len, max_len, min_len
  

if __name__ == "__main__":
    if len(sys.argv) != 2:
            print("To use this script type: python script.py <fasta_file>")
            sys.exit(1)
    
    file = sys.argv[1]
    try:
        min_length = int(input("Enter minimum seq length: "))
    except ValueError:
        print("Minimum length must be an integer.")
        sys.exit(1)

    seq_dict = read_fasta(file)
    total_sequences = len(seq_dict)
    filtered_dict = filter_sequences(seq_dict, min_length)
    filtered_count = len(filtered_dict)
    avg_len, max_len, min_len = summary_stats(filtered_dict)
    
    output_file = "filtered_sequences.fasta"
    write_fasta(filtered_dict, output_file)

    print(f"Total no. of sequences read: {total_sequences}")
    print(f"No. of Sequences that passed the length filter: {filtered_count}")
    print(f"Average length: {avg_len:.2f} bp")
    print(f"Longest sequence: {max_len} bp")
    print(f"Shortest sequence: {min_len} bp")
    print(f"Filtered sequences written to '{output_file}'")