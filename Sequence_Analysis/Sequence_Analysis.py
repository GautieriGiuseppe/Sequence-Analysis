import re
from typing import *
#from Bio import SeqIO

class ShortSeq:
    def __init__(self, startPos:int, endPos:int, seq:str) -> None:
        self.startPos, self.endPos, self.seq = startPos, endPos, seq

    def fromPattern(pattern:Pattern[str]|str, seq:str) -> list[Self]:
        return [ShortSeq(m.start, m.end, m.group()) for m in re.finditer(pattern, seq)]

    def __repr__(self) -> str:
        return f"  Position: {self.startPos}-{self.endPos}, Sequence: {self.seq}"

class Sequence:
    # Example patterns: TATA box, CAAT box, etc.
    STD_PROMOTER_PATTERNS = tuple(map(re.compile, (r'TATA[AT]A[AT]', r'CAAT', r'CGTACG')))

    def fromFasta(filePath:str) -> list[Self]:
        return list(map(Sequence.fromFastaRecord, SeqIO.parse(filePath, "fasta")))

    def fromFastaRecord(record) -> Self:
        return Sequence(str(record.seq), record.id, record.description)

    def __init__(self, seq:str, id = 0, descr = "Not provided.", promoterPatterns :tuple[Pattern[str], ...] = STD_PROMOTER_PATTERNS) -> None:
        self.id, self.sequence, self.description = id, seq, descr

        self.gc_content                = self.calculate_gc_content()
        self.promoters                 = self.findAllShortSeqsFromPatterns(promoterPatterns)
        self.nucleotide_counts         = self.calculate_nucleotide_counts_DNA()
        self.start_sites, self.dists   = self.find_transcription_start_site(self.promoters)

    def calculate_gc_content(self) -> float:
        """Calculate the GC content of a DNA sequence."""
        g = self.sequence.count('G')
        c = self.sequence.count('C')
        return (g + c) / len(self.sequence) * 100

    def calculate_nucleotide_counts_DNA(self) -> dict[str, int]:
        """Calculate the count of each nucleotide in a DNA sequence."""
        return {
            'A': self.sequence.count('A'),
            'T': self.sequence.count('T'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C')
        }

    def calculate_nucleotide_counts_RNA(self) -> dict[str, int]:
        """Calculate the count of each nucleotide in an RNA sequence."""
        return {
            'A': self.sequence.count('A'),
            'U': self.sequence.count('U'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C')
        }

    def findAllShortSeqsFromPatterns(self, patterns :tuple[Pattern[str]|str, ...]) -> list[ShortSeq]:
        """Identify promoter sequences in a DNA sequence based on given patterns."""
        return [ seq for pattern in patterns for seq in ShortSeq.fromPattern(pattern, self.sequence) ]

    def find_transcription_start_site(self, promoters:list[ShortSeq]) -> tuple[list[ShortSeq], list[int]]:
        """Identify the transcription start site using promoter distance and sequence."""
        startSites       = self.findAllShortSeqsFromPatterns([r"[CT][CT]A[ATGC][AT][CT][CT]"])
        distsToPromoters = [
            abs(transcrSeq.startPos - promoter.endPos)
            for transcrSeq in startSites
            for promoter   in promoters ]
        
        return startSites, distsToPromoters

    def __repr__(self) -> str:
        """Displays the results"""
        seqRepr = f"ID: {
            self.id
        }\nDescription: {
            self.description
        }\nLength: {
            len(self.sequence)
        }\nGC Content: {
            self.gc_content:.2f
        }%\n\nNucleotide Counts:"

        for nucleotide, count in self.nucleotide_counts.items(): seqRepr += f"  {nucleotide}: {count},"

        seqRepr += f"\nPromoter Sequences{':' if len(self.promoters) else ' not found.'}\n"
        for promoter in self.promoters: seqRepr += f"{promoter}\n"
        
        seqRepr += f"\nTranscription Start Sequences{':' if len(self.start_sites) else ' not found.'}\n"
        for transcrSeq in self.start_sites:
            d = self.dists[0] #problem here!!
            seqRepr += f"{transcrSeq}\nDistance: {d}, {
                'Abnormal' if d < 25 or d > 30 else 'Physiological'} promoter position"

        return seqRepr

if __name__ == "__main__":
    print(Sequence("GATCCACTACG"))