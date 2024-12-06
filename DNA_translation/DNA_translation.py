from Bio.Seq import Seq
from Bio import SeqIO

class DNA_sequence():
    def __init__(self, sequence):
        """ Load sequences from a list of FASTA files."""
        for record in SeqIO.parse(fasta_file, "fasta"):
            self.sequence = str(record.seq)
        if len(self.sequence) == 0:
            raise ValueError(f"No sequences found in {fasta_file}.")

    def translate(self) -> str:
        """Translate the nucleotide sequence in amino acids sequence using codon table"""
        if len(self.sequence) % 3 == 0 :
            protein_sequence = Seq.translate(Seq(self.sequence), "Standard")
        elif len(self.sequence) - 1 % 3 == 0:
            protein_sequence = Seq.translate(Seq(self.sequence)[:-1], "Standard")
        else:
            protein_sequence = Seq.translate(Seq(self.sequence)[:-2], "Standard")
        return protein_sequence


fasta_file = "/Users/giuse/pythonProject/Mycodes/samples/sequence1.txt" # add file path
sequence = DNA_sequence(fasta_file)
protein_sequence = sequence.translate()
print(f"The sequence is {protein_sequence}")



