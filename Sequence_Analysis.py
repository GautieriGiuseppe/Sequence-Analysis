import re
from Bio import SeqIO


def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    g = sequence.count('G')
    c = sequence.count('C')
    return (g + c) / len(sequence) * 100


def calculate_nucleotide_counts_DNA(sequence):
    """Calculate the count of each nucleotide in a DNA sequence."""
    return {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }


def calculate_nucleotide_counts_RNA(sequence):
    """Calculate the count of each nucleotide in a DNA sequence."""
    return {
        'A': sequence.count('A'),
        'U': sequence.count('U'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }


def find_promoter_sequences(sequence, promoter_patterns):
    """Identify promoter sequences in a DNA sequence based on given patterns."""
    promoters = []
    for pattern in promoter_patterns:
        for match in re.finditer(pattern, sequence):
            promoters.append((match.start(), match.end(), match.group()))
    return promoters


def find_transcription_start_site(sequence, promoters):
    """Identify the transcription start site using promoter distance and sequence."""
    start_site_pattern = '[CT][CT]A[ATGC][AT][CT][CT]'
    start_site = []
    distance = 0
    for match in re.finditer(start_site_pattern, sequence):
        start_site.append((match.start(), match.end(), match.group()))
    for start_init, end_init, init in start_site:
        for start_prom, end_prom, promoter in promoters:
            distance = abs(start_init - end_prom)
    return start_site, distance


def analyze_genome(fasta_file, promoter_patterns):
    """Analyze a genome from a FASTA file, including GC content, nucleotide counts, and promoter sequences."""
    results = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)
        nucleotide_counts = calculate_nucleotide_counts_DNA(sequence)
        promoters = find_promoter_sequences(sequence, promoter_patterns)
        start_site = find_transcription_start_site(sequence, promoters)[0]
        distance = find_transcription_start_site(sequence, promoters)[1]
        results.append({
            'id': record.id,
            'description': record.description,
            'length': length,
            'gc_content': gc_content,
            'nucleotide_counts': nucleotide_counts,
            'promoters': promoters,
            'transcription_start_site': start_site,
            'distance_promoter': distance,
        })
    return results

fasta_file = "" # add path to fasta file here
promoter_patterns = [r'TATA[AT]A[AT]', r'CAAT', r'CGTACG']  # Example patterns: TATA box, CAAT box, etc.
analysis_results = analyze_genome(fasta_file, promoter_patterns)

# Output results
for result in analysis_results:
    print(f"ID: {result['id']}")
    print(f"Description: {result['description']}")
    print(f"Length: {result['length']}")
    print(f"GC Content: {result['gc_content']:.2f}%")
    print("Nucleotide Counts:")
    for nucleotide, count in result['nucleotide_counts'].items():
        print(f"  {nucleotide}: {count}")
    print("Promoter Sequences:")
    for start, end, promoter in result['promoters']:
        print(f"  Position: {start}-{end}, Sequence: {promoter}")
    print("Transcription start site: ")
    if result['transcription_start_site']:
        for start, end, start_site in result['transcription_start_site']:
            print(f"  Position: {start}-{end}, Sequence: {start_site}")
            if result['distance_promoter'] < 25 or result['distance_promoter'] > 30:
                print(f"  Distance: {result['distance_promoter']}, Abnormal promoter position")
            else:
                print(f"  Distance: {result['distance_promoter']}, Physiological positioning")
    else:
        print()
        print("  Transcription start site not found in this sequence")
    print()
    print(f"Protein sequence: {result['protein_sequence']}")
    print()