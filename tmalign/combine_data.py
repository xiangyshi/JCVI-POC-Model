RESIDUALS = [344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 392, 393, 394, 395, 396, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 467, 468, 469, 470, 471, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512]

import pandas as pd
import os
import glob
import numpy as np

def interpolate_distance(target_residue, residue_distances):
    """
    Interpolate distance for a missing residue using nearest neighbors
    
    Args:
        target_residue (int): The residue number to interpolate
        residue_distances (dict): Dictionary mapping residue numbers to distances
        
    Returns:
        float: Interpolated distance value
    """
    available_residues = sorted(residue_distances.keys())
    
    # If target residue exists, return its value
    if target_residue in residue_distances:
        return residue_distances[target_residue]
    
    # Find nearest lower and higher residues
    lower_residue = None
    higher_residue = None
    
    for res in available_residues:
        if res < target_residue:
            lower_residue = res
        elif res > target_residue and higher_residue is None:
            higher_residue = res
            break
    
    # Handle edge cases
    if lower_residue is None and higher_residue is None:
        # No available residues - this shouldn't happen
        return np.nan
    elif lower_residue is None:
        # Only higher residue available - use its value
        return residue_distances[higher_residue]
    elif higher_residue is None:
        # Only lower residue available - use its value
        return residue_distances[lower_residue]
    else:
        # Interpolate between lower and higher
        lower_dist = residue_distances[lower_residue]
        higher_dist = residue_distances[higher_residue]
        
        # Linear interpolation
        fraction = (target_residue - lower_residue) / (higher_residue - lower_residue)
        interpolated = lower_dist + (higher_dist - lower_dist) * fraction
        
        return interpolated

def process_csv_file(csv_path):
    """
    Process a single CSV file and extract/interpolate distances for target residues
    
    Args:
        csv_path (str): Path to the CSV file
        
    Returns:
        dict: Dictionary mapping residue numbers to distances
    """
    print(f"Processing: {os.path.basename(csv_path)}")
    
    # Read CSV
    df = pd.read_csv(csv_path)
    
    # Create dictionary mapping chain_a_residue_num to ca_distance
    residue_distances = {}
    for _, row in df.iterrows():
        residue_num = int(row['chain_a_residue_num'])
        distance = float(row['ca_distance'])
        residue_distances[residue_num] = distance
    
    print(f"  Found {len(residue_distances)} residues in CSV")
    print(f"  Residue range: {min(residue_distances.keys())} - {max(residue_distances.keys())}")
    
    # Extract/interpolate distances for target residues
    result = {}
    missing_count = 0
    
    for target_res in RESIDUALS:
        distance = interpolate_distance(target_res, residue_distances)
        result[target_res] = distance
        
        if target_res not in residue_distances:
            missing_count += 1
    
    print(f"  Interpolated {missing_count} missing residues out of {len(RESIDUALS)} total")
    
    return result

def main():
    """Main function to process all CSV files and combine data"""
    
    # Find all CSV files in tm-data directory
    csv_pattern = "tm-data/*.csv"
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"No CSV files found in tm-data directory")
        return
    
    print(f"Found {len(csv_files)} CSV files to process")
    print("=" * 60)
    
    # Process each CSV file
    all_data = {}
    
    for csv_file in sorted(csv_files):
        # Extract sequence name from filename
        basename = os.path.basename(csv_file)
        sequence_name = basename.replace('Wuhan-A-rbd_vs_', '').replace('-rbd_distances.csv', '')
        
        # Process the CSV
        distances = process_csv_file(csv_file)
        all_data[sequence_name] = distances
        
        print()
    
    # Create combined DataFrame
    print("Creating combined DataFrame...")
    
    # Convert to DataFrame format
    rows = []
    for sequence, distances in all_data.items():
        for residue, distance in distances.items():
            rows.append({
                'sequence': sequence,
                'residue': residue,
                'ca_distance': distance
            })
    
    combined_df = pd.DataFrame(rows)
    
    # Create final CSV with residues as columns and sequences as rows
    print("Creating final CSV with residues as columns...")
    
    final_data = []
    for sequence, distances in all_data.items():
        row = {'sequence': sequence}
        # Add each residue as a column
        for residue in RESIDUALS:
            row[f'residue_{residue}'] = distances[residue]
        final_data.append(row)
    
    final_df = pd.DataFrame(final_data)
    
    # Sort by sequence name for consistency
    final_df = final_df.sort_values('sequence').reset_index(drop=True)
    
    # Save the final result
    final_output = "residue_distances_by_sequence.csv"
    final_df.to_csv(final_output, index=False)
    print(f"Saved final result to: {final_output}")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"Total sequences processed: {len(all_data)}")
    print(f"Target residues: {len(RESIDUALS)}")
    print(f"Output CSV shape: {final_df.shape}")
    
    print(f"\nComparisons processed:")
    for sequence in sorted(all_data.keys()):
        print(f"  - {sequence}")
    
    print(f"\nResidues included as columns: {len(RESIDUALS)}")
    print(f"Residue range: {min(RESIDUALS)} - {max(RESIDUALS)}")
    
    # Show sample of final table
    print(f"\nSample of final table (first 5 sequences, first 10 residue columns):")
    sample_cols = ['sequence'] + [f'residue_{r}' for r in RESIDUALS[:10]]
    sample_df = final_df[sample_cols].head(5)
    
    # Format for better display
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    print(sample_df.round(3))

if __name__ == "__main__":
    main()

