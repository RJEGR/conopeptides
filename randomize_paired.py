import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta(file_path):
    """Reads a FASTA file and returns a list of sequences."""
    try:
        return list(SeqIO.parse(file_path, "fasta"))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []

def generate_random_sequence(reference_sequences, insert_size):
    """Generates a random sequence of specified length from the reference."""
    seq = random.choice(reference_sequences)
    sequence_length = len(seq.seq)
    
    if sequence_length < insert_size:
        raise ValueError(f"Sequence {seq.id} is shorter than the insert size of {insert_size}.")
    
    start = random.randint(0, sequence_length - insert_size)
    return seq.seq[start:start + insert_size]

def create_paired_end_reads(sequence, insert_size):
    """Creates forward and reverse paired-end reads."""
    half_size = insert_size // 2
    
    if len(sequence) < insert_size:
        raise ValueError("Sequence is shorter than the insert size.")
    
    forward_read = sequence[:half_size]
    reverse_read = sequence[-half_size:].reverse_complement()
    
    return forward_read, reverse_read

def write_fastq(forward_reads, reverse_reads, output_path_1, output_path_2, phred_score=30):
    """Writes paired-end reads to FASTQ files with a specified PHRED quality score."""
    quality_string = chr(phred_score + 33) * len(forward_reads[0])
    
    with open(output_path_1, "w") as f1, open(output_path_2, "w") as f2:
        for i, (forward, reverse) in enumerate(zip(forward_reads, reverse_reads)):
            f1.write(f"@read{i}/1\n{forward}\n+\n{quality_string}\n")
            f2.write(f"@read{i}/2\n{reverse}\n+\n{quality_string}\n")

def main(input_fasta, output_prefix, N, insert_size, phred_score=30):
    """Main function to generate paired-end reads in FASTQ format."""
    sequences = read_fasta(input_fasta)
    forward_reads = []
    reverse_reads = []
    
    for _ in range(N):
        try:
            random_sequence = generate_random_sequence(sequences, insert_size)
            forward, reverse = create_paired_end_reads(random_sequence, insert_size)
            forward_reads.append(str(forward))
            reverse_reads.append(str(reverse))
        except ValueError as e:
            print(f"Skipping sequence: {e}")
    
    output_path_1 = f"{output_prefix}_1.fastq"
    output_path_2 = f"{output_prefix}_2.fastq"
    write_fastq(forward_reads, reverse_reads, output_path_1, output_path_2, phred_score)
    print(f"Paired-end FASTQ files written to {output_path_1} and {output_path_2}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Randomize sequences in a FASTA file and output paired-end FASTQ files.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_prefix", help="Output prefix for FASTQ files")
    parser.add_argument("N", type=int, help="Number of sequences to randomize")
    parser.add_argument("insert_size", type=int, help="Insert size for paired-end reads")
    parser.add_argument("--phred_score", type=int, default=30, help="PHRED quality score for FASTQ files (default: 30)")
    args = parser.parse_args()
    
    main(args.input_fasta, args.output_prefix, args.N, args.insert_size, args.phred_score)