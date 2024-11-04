from Bio import SeqIO
from Bio import AlignIO
import subprocess


def load_sequences(fasta_file):
    """ Load sequences from a list of FASTA files."""
    sequences = []
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record)
    if len(sequences) == 0:
        raise ValueError(f"No sequences found in {fasta_file}.")
    return sequences

def save_combined_fasta(sequences, output_file="combined_sequences.fasta"):
    """ Save multiple sequences into a single FASTA file for MSA."""
    SeqIO.write(sequences, output_file, "fasta")
    return output_file

def perform_msa(combined_fasta, tool="clustalw"):
    """ Perform MSA using ClustalW or Muscle."""
    if tool == "clustalw":
        cmd = ["conda", "run", "clustalw", "-align", f"-infile={combined_fasta}"]
        subprocess.run(cmd)
    elif tool == "muscle":
        cmd = ["conda", "run", "muscle", "-align", f"{combined_fasta}", "-output", "combined_sequences.aln"] # make sure to have it installed
        subprocess.run(cmd)
    else:
        raise ValueError("Unsupported tool. Use 'clustalw' or 'muscle'.")


def read_alignment(format="clustal"):
    """Read and print the alignment from the alignment file."""
    alignment = AlignIO.read("combined_sequences.aln", format)
    print(alignment)

if __name__ == "__main__":
    # List of FASTA files on local pc
    fasta_files = ["",
                   "",
                   ""] # add path of files here

    # Load and combine sequences into a single file
    sequences = load_sequences(fasta_files)
    combined_fasta = save_combined_fasta(sequences)

    # Perform MSA using selected tool
    tool = "clustalw" # change to "muscle"
    perform_msa(combined_fasta, tool)
    # Read and print the alignment
    read_alignment(format="clustal" if tool == "clustalw" else "fasta")

