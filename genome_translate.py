from Bio import SeqIO, pairwise2
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
import sys
import glob
import os

SEQUENCE_DIR = "sequences"
GENOMES_DIR = f"{SEQUENCE_DIR}/genomes"
ORFS_DIR = f"{SEQUENCE_DIR}/orfs"
SPIKES_DIR = f"{SEQUENCE_DIR}/spikes"

REF_SPIKE_FILE = f"{SPIKES_DIR}/Wuhan-A.fasta"
START_SEQ = "MFVFLVLLP"

def levenshtein_distance(s1, s2):
    """Calculate the Levenshtein distance between two strings."""
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def find_approximate_start(sequence, start_seq, max_mutations=2):
    """Find the best approximate match for start_seq in sequence, allowing up to max_mutations."""
    best_pos = 0
    best_distance = float('inf')
    
    # Search in a sliding window
    for i in range(len(sequence) - len(start_seq) + max_mutations + 1):
        # Try different window sizes to account for insertions/deletions
        for window_size in range(len(start_seq) - max_mutations, len(start_seq) + max_mutations + 1):
            if i + window_size > len(sequence):
                continue
                
            window = sequence[i:i + window_size]
            distance = levenshtein_distance(start_seq, window)
            
            if distance <= max_mutations and distance < best_distance:
                best_distance = distance
                best_pos = i
                
        # If we found an exact or very close match, we can return early
        if best_distance == 0:
            break
    
    return best_pos if best_distance <= max_mutations else 0


def find_orfs_and_translate(seq, min_len=300):
    proteins = []
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate(to_stop=False)
            aa_seq = str(trans)
            current_protein = ""
            start = None
            for i, aa in enumerate(aa_seq):
                if aa == "*":
                    if len(current_protein) >= min_len:
                        start_pos = frame + i*3 if strand == 1 else len(seq) - frame - i*3
                        proteins.append(SeqRecord(Seq(current_protein),
                                                  id=f"orf_{strand}_{frame}_{start_pos}",
                                                  description=""))
                    current_protein = ""
                else:
                    current_protein += aa
    return proteins



def process_file(input_file):
    assert("Genome" in input_file), f"Input file must be a genome file: {input_file}"
    record = SeqIO.read(f"{GENOMES_DIR}/{input_file}", "fasta")
    orfs = find_orfs_and_translate(record.seq)
    orfs_file = f"{ORFS_DIR}/{input_file.replace('Genome', 'ORFs')}"
    output_file = f"{SPIKES_DIR}/{input_file.replace('-Genome', '')}"

    SeqIO.write(orfs, orfs_file, "fasta")
    print(f"Wrote {len(orfs)} protein sequences to {orfs_file}")

    ref_spike = SeqIO.read(REF_SPIKE_FILE, "fasta").seq

    alignments = []

    for record in SeqIO.parse(orfs_file, "fasta"):
        print(record)
        alignment = pairwise2.align.localms(ref_spike, record.seq, 2, -1, -0.5, -0.1, one_alignment_only=True)
        score = alignment[0].score
        aligned_seq = alignment[0].seqB
        print(f"Match with {record.id}, score: {score}")
        alignments.append((record.id, score, aligned_seq))

    alignments = sorted(alignments, key=lambda x: x[1], reverse=True)
    with open(output_file, "w") as out:
        start_pos = find_approximate_start(alignments[0][2], START_SEQ)
        found_sequence = alignments[0][2][start_pos:start_pos+len(START_SEQ)]
        distance = levenshtein_distance(START_SEQ, found_sequence)
        
        print(f"Looking for start sequence: {START_SEQ}")
        print(f"Found at position {start_pos}: {found_sequence}")
        print(f"Edit distance: {distance}")
        
        out.write(f">{input_file.replace('-Genome.fasta', '')}\n{alignments[0][2][start_pos:].replace("-", "")}\n")


def main(pattern):
    """Main function that processes files matching the given pattern."""
    # Handle glob patterns
    if '*' in pattern or '?' in pattern or '[' in pattern:
        # It's a glob pattern
        search_pattern = f"{GENOMES_DIR}/{pattern}"
        if not pattern.endswith('.fasta'):
            search_pattern += "*.fasta"
        
        matching_files = glob.glob(search_pattern)
        
        if not matching_files:
            print(f"No files found matching pattern: {pattern}")
            print(f"Searched in: {search_pattern}")
            return
        
        # Extract just the filenames from full paths
        input_files = [os.path.basename(f) for f in matching_files]
        
        print(f"Found {len(input_files)} files matching pattern '{pattern}':")
        for f in input_files:
            print(f"  - {f}")
        # Process each file
        for input_file in input_files:
            try:
                print(f"\n=== Processing {input_file} ===")
                process_file(input_file)
                print(f"✓ Completed processing {input_file}")
            except Exception as e:
                print(f"✗ Error processing {input_file}: {e}")
                continue
    else:
        # It's a single file
        if not pattern.endswith('.fasta'):
            pattern += '.fasta'
        process_file(pattern)


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) != 1:
        print("Usage: python genome_translate.py <input_file_or_pattern>")
        print("\nExamples:")
        print("  python genome_translate.py Beta-B.1.351-Sample-1-Genome.fasta  # Single file")
        print("  python genome_translate.py Beta*                               # All files starting with 'Beta'")
        print("  python genome_translate.py *-Genome                            # All genome files")
        sys.exit(1)
    input_pattern = args[0]
    print(input_pattern)
    main(input_pattern)
