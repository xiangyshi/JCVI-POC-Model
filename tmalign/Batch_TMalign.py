#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import glob
from pathlib import Path

def run_tmalign(pdb_ref, pdb_target, output_prefix):
    """Run TMalign between reference and target PDB"""
    cmd = ["./TMalign", pdb_ref, pdb_target, "-o", output_prefix]
    
    # Clean up any existing temp files first
    temp_files = glob.glob(f"{output_prefix}*")
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except:
            pass
    
    print(f"  Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"  TMalign completed successfully")
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        print(f"  TMalign failed for {pdb_target}: {e.stderr}")
        return False, None
    except Exception as e:
        print(f"  Error running TMalign: {e}")
        return False, None

def extract_distances(aligned_file, output_csv):
    """Extract pairwise distances using the extract_pairwise_distances.py script"""
    cmd = ["python", "extract_pairwise_distances.py", aligned_file, output_csv]
    print(f"  Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Distance extraction failed: {e.stderr}")
        return False
    except Exception as e:
        print(f"Error extracting distances: {e}")
        return False

def get_pdb_files(folder_path):
    """Get all PDB files from the specified folder"""
    pdb_extensions = ['*.pdb', '*.PDB']
    pdb_files = []
    
    for ext in pdb_extensions:
        pdb_files.extend(glob.glob(os.path.join(folder_path, ext)))
    
    return sorted(pdb_files)

def create_directories():
    """Create necessary directories"""
    os.makedirs("temp", exist_ok=True)
    os.makedirs("tm-data", exist_ok=True)

def get_base_name(pdb_path):
    """Get base name from PDB file path for naming output files"""
    return Path(pdb_path).stem

def main():
    parser = argparse.ArgumentParser(description='Batch TMalign analysis')
    parser.add_argument('pdb_ref', help='Reference PDB file')
    parser.add_argument('folder', help='Folder containing PDB files to compare')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.pdb_ref):
        print(f"Error: Reference PDB file '{args.pdb_ref}' not found!")
        sys.exit(1)
    
    if not os.path.isdir(args.folder):
        print(f"Error: Folder '{args.folder}' not found!")
        sys.exit(1)
    
    # Check if TMalign executable exists
    if not os.path.exists("./TMalign"):
        print("Error: TMalign executable not found in current directory!")
        sys.exit(1)
    
    # Check if extract_pairwise_distances.py exists
    if not os.path.exists("extract_pairwise_distances.py"):
        print("Error: extract_pairwise_distances.py not found in current directory!")
        sys.exit(1)
    
    # Create output directories
    create_directories()
    
    # Get all PDB files from folder
    pdb_files = get_pdb_files(args.folder)
    
    if not pdb_files:
        print(f"No PDB files found in folder '{args.folder}'")
        sys.exit(1)
    
    print(f"Reference PDB: {args.pdb_ref}")
    print(f"Target folder: {args.folder}")
    print(f"Found {len(pdb_files)} PDB files to process")
    print("-" * 60)
    
    successful = 0
    failed = 0
    
    for pdb_file in pdb_files:
        if args.verbose:
            print(f"\nProcessing: {pdb_file}")
        
        # Get base name for output files
        base_name = get_base_name(pdb_file)
        ref_name = get_base_name(args.pdb_ref)
        
        # Define output paths
        temp_prefix = "temp/temp"
        aligned_file = "temp/temp_all_atm_lig"
        output_csv = f"tm-data/{ref_name}_vs_{base_name}_distances.csv"
        
        # Skip if reference and target are the same file
        if os.path.abspath(args.pdb_ref) == os.path.abspath(pdb_file):
            if args.verbose:
                print(f"  Skipping self-comparison")
            continue
        
        # Run TMalign (always run fresh for each comparison)
        print(f"  Step 1: Running TMalign for {base_name}...")
        
        success, tmalign_output = run_tmalign(args.pdb_ref, pdb_file, temp_prefix)
        
        if not success:
            print(f"✗ TMalign failed for {base_name}")
            failed += 1
            continue
        
        # Check if aligned file was created
        if not os.path.exists(aligned_file):
            print(f"  Error: Expected TMalign output file '{aligned_file}' not found")
            print(f"  Available files: {glob.glob(f'{temp_prefix}*')}")
            failed += 1
            continue
        
        print(f"  Step 2: TMalign output file created: {aligned_file}")
        
        # Extract distances
        print(f"  Step 3: Extracting distances...")
        
        if extract_distances(aligned_file, output_csv):
            print(f"✓ {base_name} -> {output_csv}")
            successful += 1
        else:
            print(f"✗ Failed to extract distances for {base_name}")
            failed += 1
        
        # Clean up temporary files
        temp_files = glob.glob(f"{temp_prefix}*")
        for temp_file in temp_files:
            try:
                os.remove(temp_file)
            except:
                pass
    
    print("\n" + "=" * 60)
    print(f"BATCH PROCESSING COMPLETE")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total processed: {successful + failed}")
    print(f"Results saved in: tm-data/")
    print("=" * 60)

if __name__ == "__main__":
    main() 