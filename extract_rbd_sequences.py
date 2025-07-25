#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob

def extract_rbd_sequences(input_dir="sequences/spikes", output_dir="sequences/rbds"):
    """
    Extract RBD sequences (residues 319-541) from all FASTA files in the input directory.
    
    Args:
        input_dir (str): Directory containing input FASTA files
        output_dir (str): Directory to save output files (default: same as input)
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all .fasta files in the input directory
    fasta_files = glob.glob(os.path.join(input_dir, "*.fasta"))
    
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}")
        return
    
    print(f"Found {len(fasta_files)} FASTA files to process...")
    
    # RBD region: residues 319-541 (1-based) = positions 318-540 (0-based)
    rbd_start = 318  # 0-based index for residue 319
    rbd_end = 541    # 0-based end index (exclusive) for residue 541
    
    processed_count = 0
    
    for fasta_file in fasta_files:
        try:
            # Parse the input file
            records = list(SeqIO.parse(fasta_file, "fasta"))
            
            if not records:
                print(f"Warning: No sequences found in {fasta_file}")
                continue
                
            # Process each sequence in the file
            rbd_records = []
            
            for record in records:
                sequence = str(record.seq)
                
                # Check if sequence is long enough for RBD extraction
                if len(sequence) < rbd_end:
                    print(f"Warning: Sequence in {fasta_file} is too short ({len(sequence)} residues). Need at least {rbd_end} residues.")
                    continue
                
                # Extract RBD region
                rbd_sequence = sequence[rbd_start:rbd_end]
                
                # Create new record with RBD sequence
                rbd_record = SeqRecord(
                    Seq(rbd_sequence),
                    id=record.id + "-rbd",
                    description=record.description + f" RBD region (residues 319-541)"
                )
                
                rbd_records.append(rbd_record)
            
            if rbd_records:
                # Generate output filename
                base_name = os.path.basename(fasta_file)
                name_without_ext = os.path.splitext(base_name)[0]
                output_filename = f"{name_without_ext}-rbd.fasta"
                output_path = os.path.join(output_dir, output_filename)
                
                # Write RBD sequences to output file
                SeqIO.write(rbd_records, output_path, "fasta")
                
                print(f"Processed: {base_name} -> {output_filename} ({len(rbd_records)} sequence(s))")
                processed_count += 1
            
        except Exception as e:
            print(f"Error processing {fasta_file}: {str(e)}")
            continue
    
    print(f"\nCompleted! Processed {processed_count} files successfully.")
    print(f"RBD sequences extracted from residues 319-541 ({rbd_end - rbd_start} residues)")

def main():
    """Main function to run the RBD extraction."""
    
    # You can modify these paths as needed
    input_directory = "sequences/spikes"
    output_directory = "sequences/rbds"  # Save in the same directory
    
    # Check if input directory exists
    if not os.path.exists(input_directory):
        print(f"Error: Input directory '{input_directory}' does not exist!")
        return
    
    print("Extracting RBD sequences (residues 319-541) from spike protein sequences...")
    print(f"Input directory: {input_directory}")
    print(f"Output directory: {output_directory}")
    print("-" * 60)
    
    extract_rbd_sequences(input_directory, output_directory)

if __name__ == "__main__":
    main() 