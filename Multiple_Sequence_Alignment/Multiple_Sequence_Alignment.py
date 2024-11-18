from Bio import SeqIO
from Bio import AlignIO
import subprocess

class MSA():
    def __init__(self, fasta_file):
        """ Load sequences from a list of FASTA files."""
        self.sequences = []
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                self.sequences.append(record)
        if len(self.sequences) == 0:
            raise ValueError(f"No sequences found in {self}.")


    def save_combined_fasta(self, output_file="combined_sequences.fasta"):
        """ Save multiple sequences into a single FASTA file for MSA."""
        SeqIO.write(self.sequences, output_file, "fasta")
        return output_file


    def perform_msa(self, tool):
        """ Perform MSA using ClustalW or Muscle."""
        if tool == "clustalw": # make sure to have it installed
            cmd = ["conda", "run", "clustalw", "-align", f"-infile={combined_fasta}"]
            subprocess.run(cmd)
        elif tool == "muscle": # make sure to have it installed
            cmd = ["conda", "run", "muscle", "-align", f"{combined_fasta}", "-output", "combined_sequences.aln"]
            subprocess.run(cmd)
        else:
            raise ValueError("Unsupported tool. Use 'clustalw' or 'muscle'.")


    def read_alignment(self, format="clustal"):
        """Read and print the alignment from the alignment file."""
        alignment = AlignIO.read("combined_sequences.aln", format)
        print(alignment)


if __name__ == "__main__":
    # List of FASTA files on local pc
    fasta_files = [] # add path of files here

    # Load and combine sequences into a single file
    sequences = MSA(fasta_files)
    combined_fasta = sequences.save_combined_fasta()

    # Perform MSA using selected tool
    alignment = MSA(combined_fasta).perform_msa(tool="clustalw")  # or change to 'muscle'
    # Read and print the alignment
    tool = "clustalw"
    MSA(alignment).read_alignment(format="clustal" if tool == "clustalw" else "fasta")

