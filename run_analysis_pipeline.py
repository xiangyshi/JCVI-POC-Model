#!/usr/bin/env python3
"""
Analysis Pipeline Driver Script

This script runs the complete structural analysis pipeline:
1. Execute tmalign/Batch_TMalign.py to generate pairwise distance CSV files
2. Execute tmalign/combine_data.py to combine and process the results

Usage:
    python run_analysis_pipeline.py <reference_pdb> <target_directory>
    
Example:
    python run_analysis_pipeline.py rbd-structures/Wuhan-A-rbd.pdb rbd-structures
"""

import subprocess
import sys
import os
import time

def run_command(command, description, working_dir=None):
    """
    Run a command and handle errors
    
    Args:
        command (list): Command to run as list of arguments
        description (str): Description of what the command does
        working_dir (str): Directory to run the command in
        
    Returns:
        bool: True if successful, False if failed
    """
    print(f"\n{'='*60}")
    print(f"STEP: {description}")
    print(f"{'='*60}")
    print(f"Running: {' '.join(command)}")
    if working_dir:
        print(f"Working directory: {working_dir}")
    
    try:
        start_time = time.time()
        result = subprocess.run(command, check=True, capture_output=True, text=True, cwd=working_dir)
        end_time = time.time()
        
        print(f"✓ SUCCESS - Completed in {end_time - start_time:.1f} seconds")
                
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ FAILED - Command returned non-zero exit status {e.returncode}")
        if e.stdout:
            print("STDOUT:")
            print(e.stdout)
        if e.stderr:
            print("STDERR:")
            print(e.stderr)
        return False
    except Exception as e:
        print(f"✗ ERROR - {str(e)}")
        return False

def check_prerequisites():
    """Check if required files and directories exist"""
    print("Checking prerequisites...")
    
    required_files = [
        "tmalign/Batch_TMalign.py",
        "tmalign/combine_data.py",
        "tmalign/extract_pairwise_distances.py",
        "tmalign/TMalign"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"✗ Missing required files: {missing_files}")
        return False
    
    # Check if TMalign is executable
    if not os.access("tmalign/TMalign", os.X_OK):
        print("✗ tmalign/TMalign is not executable")
        return False
    
    print("✓ All prerequisites satisfied")
    return True

def main():
    """Main pipeline execution"""
    
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python run_analysis_pipeline.py <reference_pdb> <target_directory>")
        print("Example: python run_analysis_pipeline.py rbd-structures/Wuhan-A-rbd.pdb rbd-structures")
        sys.exit(1)
    
    reference_pdb = sys.argv[1]
    target_directory = sys.argv[2]
    
    # Validate inputs
    if not os.path.exists(reference_pdb):
        print(f"✗ Reference PDB file not found: {reference_pdb}")
        sys.exit(1)
        
    if not os.path.isdir(target_directory):
        print(f"✗ Target directory not found: {target_directory}")
        sys.exit(1)
    
    print("STRUCTURAL ANALYSIS PIPELINE")
    print("="*60)
    print(f"Reference PDB: {reference_pdb}")
    print(f"Target directory: {target_directory}")
    print(f"Working directory: {os.getcwd()}")
    
    # Check prerequisites
    if not check_prerequisites():
        print("\n✗ PIPELINE FAILED - Prerequisites not satisfied")
        sys.exit(1)
    
    pipeline_start = time.time()
    
    # Step 1: Run Batch_TMalign.py (from tmalign directory)
    step1_success = run_command(
        ["python", "Batch_TMalign.py", f"../{reference_pdb}", f"../{target_directory}"],
        "Generate pairwise distance CSV files using TMalign",
        working_dir="tmalign"
    )
    
    if not step1_success:
        print("\n✗ PIPELINE FAILED - Step 1 (Batch_TMalign.py) failed")
        sys.exit(1)
    
    # Step 2: Run combine_data.py (from tmalign directory)
    step2_success = run_command(
        ["python", "combine_data.py"],
        "Combine and process distance data into final result table",
        working_dir="tmalign"
    )
    
    if not step2_success:
        print("\n✗ PIPELINE FAILED - Step 2 (combine_data.py) failed")
        sys.exit(1)
    
    
    # Pipeline completed successfully
    pipeline_end = time.time()
    total_time = pipeline_end - pipeline_start
    
    print(f"\n{'='*60}")
    print("🎉 PIPELINE COMPLETED SUCCESSFULLY! 🎉")
    print(f"{'='*60}")
    print(f"Total execution time: {total_time:.1f} seconds")
    
    # Show final results
    output_files = [
        "residue_distances_by_sequence.csv"
    ]
    
    print(f"\nOutput files generated in root directory:")
    for output_file in output_files:
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"  ✓ {output_file} ({file_size:,} bytes)")
        else:
            print(f"  ✗ {output_file} (not found)")
    
    print(f"\nMain result file: residue_distances_by_sequence.csv")
    print(f"This file contains:")
    print(f"  - Each row: one sequence comparison")
    print(f"  - Each column: one important residue distance")
    print(f"  - Missing residues interpolated using nearest neighbors")

if __name__ == "__main__":
    main() 