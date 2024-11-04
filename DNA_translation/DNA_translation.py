from Bio.Seq import Seq
from Bio import SeqIO


def load_sequence(fasta_file):
    """ Load sequences from a list of FASTA files."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
    if len(sequence) == 0:
        raise ValueError(f"No sequences found in {fasta_file}.")
    return sequence

def translate(sequence):
    # translate sequence in aminoacids
    if len(sequence) % 3 == 0 :
        protein_sequence = Seq.translate(sequence, "Standard")
    elif (len(sequence) - 1) % 3 == 0:
        protein_sequence = Seq.translate(sequence[:-1], "Standard")
    else:
        protein_sequence = Seq.translate(sequence[:-2], "Standard")
    return protein_sequence


fasta_file = "" # add file path
sequence = load_sequence(fasta_file)
protein_sequence = translate(sequence)
print(f"The sequence is {protein_sequence}")



