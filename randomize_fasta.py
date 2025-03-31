import random
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta(file_path):
    try:
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt") as handle:
                return list(SeqIO.parse(handle, "fasta"))
        else:
            return list(SeqIO.parse(file_path, "fasta"))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []

def randomize_sequences(sequences, N, kmer_length=150):
    if not sequences:
        return []
    
    randomized_sequences = []
    
    for _ in range(N):
        seq = random.choice(sequences)
        sequence_length = len(seq.seq)
        
        if sequence_length < kmer_length:
            print(f"Sequence {seq.id} is shorter than kmer length of {kmer_length}. Skipping.")
            continue
        
        random_start = random.randint(0, sequence_length - kmer_length)
        kmer = seq.seq[random_start:random_start + kmer_length]
        randomized_sequences.append(SeqRecord(Seq(str(kmer)), id=f"{seq.id}_randomized", description="randomized kmer"))
    
    return randomized_sequences

def write_fasta(sequences, output_path):
    try:
        SeqIO.write(sequences, output_path, "fasta")
        print(f"Successfully wrote {len(sequences)} sequences to {output_path}")
    except Exception as e:
        print(f"Error writing FASTA file: {e}")

def write_paired_end_fasta(sequences, output_path_1, output_path_2, insert_size=300):
    forward_reads = []
    reverse_reads = []
    
    for seq in sequences:
        sequence_length = len(seq.seq)
        if sequence_length < insert_size:
            print(f"Sequence {seq.id} is shorter than insert size of {insert_size}. Skipping.")
            continue
        
        forward_read = seq.seq[:insert_size//2]
        reverse_read = seq.seq[-(insert_size//2):]
        
        forward_reads.append(SeqRecord(forward_read, id=f"{seq.id}_forward", description="paired-end forward read"))
        reverse_reads.append(SeqRecord(reverse_read, id=f"{seq.id}_reverse", description="paired-end reverse read"))
    
    write_fasta(forward_reads, output_path_1)
    write_fasta(reverse_reads, output_path_2)

def main(input_fasta, output_fasta, N, paired_end=False, insert_size=300):
    sequences = read_fasta(input_fasta)
    randomized_sequences = randomize_sequences(sequences, N)
    
    if paired_end:
        output_fasta_1 = output_fasta.replace(".fasta", "_1.fasta")
        output_fasta_2 = output_fasta.replace(".fasta", "_2.fasta")
        write_paired_end_fasta(randomized_sequences, output_fasta_1, output_fasta_2, insert_size)
    else:
        write_fasta(randomized_sequences, output_fasta)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Randomize sequences in a FASTA file by k-mer length.")
    parser.add_argument("input_fasta", help="Input FASTA (or gzipped FASTA) file")
    parser.add_argument("output_fasta", help="Output FASTA file")
    parser.add_argument("N", type=int, help="Number of sequences to randomize")
    parser.add_argument("--paired_end", action="store_true", help="Output paired-end FASTA files")
    parser.add_argument("--insert_size", type=int, default=300, help="Insert size for paired-end reads")
    args = parser.parse_args()
    
    main(args.input_fasta, args.output_fasta, args.N, args.paired_end, args.insert_size)