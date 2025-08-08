#!/usr/bin/env python3

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
import sys
import os

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates"""
    return np.sqrt(np.sum((coord1 - coord2)**2))

def extract_ca_coordinates(structure, chain_id):
    """Extract CA coordinates from a structure for a specific chain"""
    coordinates = {}
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    res_id = residue.id[1]  # Get residue number
                    res_name = residue.get_resname()
                    
                    # Get CA atom
                    if 'CA' in residue:
                        ca_coord = residue['CA'].get_coord()
                        coordinates[res_id] = {
                            'residue_name': res_name,
                            'coordinates': ca_coord
                        }
    
    return coordinates

def parse_tmalign_alignment(alignment_file):
    """Parse TMalign alignment file to get residue pair mappings"""
    alignment_pairs = []
    
    if not os.path.exists(alignment_file):
        print(f"Warning: Alignment file {alignment_file} not found")
        return alignment_pairs
    
    with open(alignment_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('select') and ':A,' in line and ':B' in line:
                # Parse lines like: select  328:A,  11:B
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        # Extract residue numbers
                        a_part = parts[1].split(':')[0]
                        b_part = parts[2].split(':')[0]
                        res_a = int(a_part)
                        res_b = int(b_part.rstrip(','))
                        alignment_pairs.append((res_a, res_b))
                    except (ValueError, IndexError):
                        continue
    
    return alignment_pairs

def get_pairwise_distances(aligned_pdb_file, output_csv=None, alignment_file=None):
    """
    Extract pairwise distances between aligned residues from TMalign output
    
    Args:
        aligned_pdb_file (str): Path to the aligned PDB file from TMalign
        output_csv (str): Optional path to save results as CSV
        alignment_file (str): Optional path to TMalign alignment file
        
    Returns:
        pandas.DataFrame: DataFrame with residue pairs and distances
    """
    
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('aligned', aligned_pdb_file)
    
    # Extract coordinates for both chains
    chain_a_coords = extract_ca_coordinates(structure, 'A')  # First structure
    chain_b_coords = extract_ca_coordinates(structure, 'B')  # Second structure
    
    print(f"Chain A has {len(chain_a_coords)} residues")
    print(f"Chain B has {len(chain_b_coords)} residues")
    
    # Parse alignment file if provided
    aligned_pairs = []
    if alignment_file and os.path.exists(alignment_file):
        aligned_pairs = parse_tmalign_alignment(alignment_file)
        print(f"Found {len(aligned_pairs)} aligned pairs from TMalign alignment file")
    
    # If no alignment file or no pairs found, fall back to residue number matching
    if not aligned_pairs:
        print("Using residue number matching...")
        # Debug: Show residue number ranges
        if chain_a_coords:
            a_min, a_max = min(chain_a_coords.keys()), max(chain_a_coords.keys())
            print(f"Chain A residue range: {a_min} to {a_max}")
        
        if chain_b_coords:
            b_min, b_max = min(chain_b_coords.keys()), max(chain_b_coords.keys())
            print(f"Chain B residue range: {b_min} to {b_max}")
        
        # Find common residue numbers
        common_residues = set(chain_a_coords.keys()) & set(chain_b_coords.keys())
        aligned_pairs = [(res, res) for res in sorted(common_residues)]
    
    print(f"Using {len(aligned_pairs)} aligned residue pairs")
    
    # Calculate pairwise distances
    distance_data = []
    
    for res_a_num, res_b_num in aligned_pairs:
        # Check if both residues exist in the structures
        if res_a_num not in chain_a_coords or res_b_num not in chain_b_coords:
            continue
            
        res_a = chain_a_coords[res_a_num]
        res_b = chain_b_coords[res_b_num]
        
        # Calculate distance between CA atoms
        distance = calculate_distance(res_a['coordinates'], res_b['coordinates'])
        
        distance_data.append({
            'chain_a_residue_num': res_a_num,
            'chain_b_residue_num': res_b_num,
            'chain_a_residue': res_a['residue_name'],
            'chain_b_residue': res_b['residue_name'],
            'ca_distance': distance,
            'chain_a_coords': res_a['coordinates'],
            'chain_b_coords': res_b['coordinates']
        })
    
    # Create DataFrame
    df = pd.DataFrame(distance_data)
    
    # Add some statistics
    print(f"\nDistance Statistics:")
    print(f"Mean distance: {df['ca_distance'].mean():.3f} Å")
    print(f"Std distance:  {df['ca_distance'].std():.3f} Å")
    print(f"Min distance:  {df['ca_distance'].min():.3f} Å")
    print(f"Max distance:  {df['ca_distance'].max():.3f} Å")
    
    # Save to CSV if specified
    if output_csv:
        # Prepare simplified DataFrame for CSV (without coordinate arrays)
        csv_df = df[['chain_a_residue_num', 'chain_b_residue_num', 'chain_a_residue', 'chain_b_residue', 'ca_distance']].copy()
        csv_df.to_csv(output_csv, index=False)
        print(f"\nSaved results to {output_csv}")
    
    return df

def get_all_vs_all_distances(coords_dict):
    """Calculate all pairwise distances within a single structure"""
    residues = sorted(coords_dict.keys())
    distances = []
    
    for i, res1 in enumerate(residues):
        for j, res2 in enumerate(residues):
            if i <= j:  # Only calculate upper triangle + diagonal
                dist = calculate_distance(coords_dict[res1]['coordinates'], 
                                        coords_dict[res2]['coordinates'])
                distances.append({
                    'residue_1': res1,
                    'residue_2': res2,
                    'distance': dist,
                    'res1_name': coords_dict[res1]['residue_name'],
                    'res2_name': coords_dict[res2]['residue_name']
                })
    
    return pd.DataFrame(distances)

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_pairwise_distances.py <aligned_pdb_file> [output_csv]")
        print("Example: python extract_pairwise_distances.py aligned_output.pdb distances.csv")
        return
    
    aligned_pdb_file = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(aligned_pdb_file):
        print(f"Error: File {aligned_pdb_file} not found!")
        return
    
    print(f"Processing aligned structures from: {aligned_pdb_file}")
    print("=" * 60)
    
    # Try to find the alignment file (remove _all_atm_lig suffix)
    alignment_file = None
    if aligned_pdb_file.endswith('_all_atm_lig'):
        alignment_file = aligned_pdb_file[:-12]  # Remove '_all_atm_lig' (12 chars)
        if os.path.exists(alignment_file):
            print(f"Using alignment file: {alignment_file}")
        else:
            alignment_file = None
    
    # Extract pairwise distances between aligned residues
    df = get_pairwise_distances(aligned_pdb_file, output_csv, alignment_file)
    
    # Show some example distances
    print(f"\nSample pairwise distances:")
    print(df[['chain_a_residue_num', 'chain_b_residue_num', 'chain_a_residue', 'chain_b_residue', 'ca_distance']].head(10))
    
    # Identify residues with large structural differences
    threshold = 2.0  # Angstroms
    large_differences = df[df['ca_distance'] > threshold]
    
    if len(large_differences) > 0:
        print(f"\nResidues with CA distance > {threshold} Å:")
        print(large_differences[['chain_a_residue_num', 'chain_b_residue_num', 'chain_a_residue', 'chain_b_residue', 'ca_distance']])
    else:
        print(f"\nAll aligned residues have CA distance ≤ {threshold} Å")

if __name__ == "__main__":
    main() 