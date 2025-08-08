from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
import numpy as np
import math

def normalize(v):
    arr = np.asarray(v)
    norm = np.linalg.norm(arr)
    if norm < 1e-10:
        raise ValueError("Cannot normalize zero vector")
    return arr / norm

def get_residue_vector(residue, atom_name):
    return residue[atom_name].get_vector()

def build_local_frame(residue):
    """
    Build an orthonormal frame using CA as origin.
    X: CA→CB (or CA→N for glycine)
    Y: perpendicular to N–CA–C plane
    Z: completes right-handed frame
    """
    try:
        ca = get_residue_vector(residue, "CA")
        n = get_residue_vector(residue, "N")
        c = get_residue_vector(residue, "C")

        ca_arr = ca.get_array()
        n_vec = normalize(n.get_array() - ca_arr)
        c_vec = normalize(c.get_array() - ca_arr)

        resname = residue.get_resname().upper()
        if resname == "GLY":
            x_vec = n_vec
        else:
            cb = get_residue_vector(residue, "CB")
            x_vec = normalize(cb.get_array() - ca_arr)

        # Initial Y from N–CA–C plane
        cross_product = np.cross(n_vec, c_vec)
        cross_norm = np.linalg.norm(cross_product)

        if cross_norm < 1e-6:
            # Fallback: use X and C
            alt = np.cross(c_vec, x_vec)
            y_axis = normalize(alt) if np.linalg.norm(alt) > 1e-6 else np.array([0.0, 1.0, 0.0])
        else:
            y_axis = normalize(cross_product)

        z_axis = normalize(np.cross(x_vec, y_axis))
        y_axis = normalize(np.cross(z_axis, x_vec))  # re-orthogonalize

        return np.stack([x_vec, y_axis, z_axis], axis=0), ca_arr
    except Exception as e:
        raise RuntimeError(f"Failed to build frame for residue {residue.get_id()} ({residue.get_resname()}): {e}")

def compute_spherical_angles(vec_in_local_frame):
    x, y, z = vec_in_local_frame
    r = np.linalg.norm([x, y, z])
    if r < 1e-10:
        return 0.0, 0.0
    z_over_r = max(min(z / r, 1.0), -1.0)
    theta = math.acos(z_over_r) * 180 / math.pi  # elevation
    phi = math.atan2(y, x) * 180 / math.pi       # azimuth
    return theta, phi

def compute_chi1(residue):
    """Compute χ₁ torsion angle using standard definition N–CA–CB–CG (or variant)"""
    chi1_atoms = {
        "ARG": ["N", "CA", "CB", "CG"],
        "ASN": ["N", "CA", "CB", "CG"],
        "ASP": ["N", "CA", "CB", "CG"],
        "CYS": ["N", "CA", "CB", "SG"],
        "GLN": ["N", "CA", "CB", "CG"],
        "GLU": ["N", "CA", "CB", "CG"],
        "HIS": ["N", "CA", "CB", "CG"],
        "ILE": ["N", "CA", "CB", "CG1"],
        "LEU": ["N", "CA", "CB", "CG"],
        "LYS": ["N", "CA", "CB", "CG"],
        "MET": ["N", "CA", "CB", "CG"],
        "PHE": ["N", "CA", "CB", "CG"],
        "PRO": ["N", "CA", "CB", "CG"],
        "SER": ["N", "CA", "CB", "OG"],
        "THR": ["N", "CA", "CB", "OG1"],
        "TRP": ["N", "CA", "CB", "CG"],
        "TYR": ["N", "CA", "CB", "CG"],
        "VAL": ["N", "CA", "CB", "CG1"],
    }

    resname = residue.get_resname().upper()
    if resname not in chi1_atoms:
        return None
    try:
        a1, a2, a3, a4 = [get_residue_vector(residue, atom) for atom in chi1_atoms[resname]]
        angle_rad = calc_dihedral(a1, a2, a3, a4)
        return np.degrees(angle_rad)
    except KeyError:
        return None

def compute_sidechain_angles_from_residue(residue, terminal_atom="CB"):
    resname = residue.get_resname().upper()
    if terminal_atom not in residue:
        raise ValueError(f"{terminal_atom} not found in residue {resname} {residue.get_id()}")
    
    if resname == "GLY":
        return 0.0, 0.0, None  # No side chain

    frame, ca_coord = build_local_frame(residue)
    v = residue[terminal_atom].get_vector().get_array() - ca_coord
    local_coords = frame @ v
    theta, phi = compute_spherical_angles(local_coords)
    chi1 = compute_chi1(residue)
    return theta, phi, chi1

def compute_sidechain_angles(pdb_file, chain_id, res_id, terminal_atom="CB"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found")
    chain = model[chain_id]

    if res_id not in chain:
        available = [r[1] for r in chain.child_dict.keys()]
        raise ValueError(f"Residue {res_id} not found. Available: {available}")

    residue = chain[res_id]
    return compute_sidechain_angles_from_residue(residue, terminal_atom)

# ----------------------------
# Batch analysis of all PDB files
# ----------------------------

import os
import glob
import pandas as pd

terminal_atoms = {
    'ARG': 'CZ', 'LYS': 'NZ', 'ASP': 'CG', 'GLU': 'CD',
    'ASN': 'CG', 'GLN': 'CD', 'HIS': 'CG', 'PHE': 'CZ',
    'TYR': 'OH', 'TRP': 'CH2', 'SER': 'OG', 'THR': 'OG1',
    'CYS': 'SG', 'MET': 'SD', 'LEU': 'CD1', 'ILE': 'CD1',
    'VAL': 'CG1', 'PRO': 'CG', 'GLY': 'CA'
}

# Key residues to analyze (will be adjusted based on what's available in each file)
# These are important RBD residues, but we'll analyze whatever is available
key_residues = [417, 446, 449, 487, 489, 493, 500, 501, 502, 505]
key_residues = [res - 318 for res in key_residues]

# Get all PDB files in the pdbs directory
pdb_files = glob.glob("pdbs/*.pdb")
pdb_files.sort()  # Sort for consistent output

print(f"Found {len(pdb_files)} PDB files to analyze:")
for file in pdb_files:
    print(f"  - {os.path.basename(file)}")
print()

# Store all results for CSV output
all_results = {}

# Process each PDB file
for pdb_file in pdb_files:
    filename = os.path.basename(pdb_file)
    # Remove "-rbd.pdb" suffix for cleaner filenames
    clean_filename = filename.replace("-rbd.pdb", "").replace(".pdb", "")
    
    print(f"=== Analyzing {filename} ===")
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        
        # Check if chain A exists
        if "A" not in structure[0]:
            print(f"  Chain A not found in {filename}")
            print()
            continue
            
        chain = structure[0]["A"]
        
        # Get the range of residues in this structure
        residue_ids = list(chain.child_dict.keys())
        min_res = min(r[1] for r in residue_ids)
        max_res = max(r[1] for r in residue_ids)
        
        print(f"  Residue range: {min_res} - {max_res}")
        
        # Try to find key residues, but if not available, analyze all residues
        adjusted_residues = []
        for res in key_residues:
            adjusted_res = res - min_res + 1  # Adjust to 1-based indexing
            if adjusted_res in chain:
                adjusted_residues.append(adjusted_res)
        
        # If no key residues found, analyze all residues in the structure
        if not adjusted_residues:
            print(f"  Key residues not found, analyzing all residues in {filename}")
            adjusted_residues = list(chain.child_dict.keys())
            # Convert from (chain, res_id, insertion) format to just res_id
            adjusted_residues = [r[1] for r in adjusted_residues]
            # Limit to first 20 residues to avoid overwhelming output
            adjusted_residues = adjusted_residues[:20]
            print(f"  Analyzing first {len(adjusted_residues)} residues")
        
        # Initialize results for this file
        file_results = {}
        
        # Analyze each key residue
        results = []
        for res_id in adjusted_residues:
            try:
                residue = chain[res_id]
                resname = residue.get_resname().upper()
                terminal_atom = terminal_atoms.get(resname, 'CB')
                
                theta, phi, chi1 = compute_sidechain_angles_from_residue(residue, terminal_atom)
                chi1_str = f"{chi1:.2f}°" if chi1 is not None else "N/A"
                
                # Store results
                original_res = res_id + min_res - 1  # Convert back to original numbering
                results.append({
                    'residue': original_res,
                    'resname': resname,
                    'theta': theta,
                    'phi': phi,
                    'chi1': chi1
                })
                
                # Store in file_results for CSV
                file_results[f"{original_res + 318}_Theta"] = theta
                file_results[f"{original_res + 318}_Phi"] = phi
                file_results[f"{original_res + 318}_Chi1"] = chi1 if chi1 is not None else "N/A"
                
                print(f"  Residue {original_res} ({resname}): Theta={theta:.2f}°, Phi={phi:.2f}°, Chi1={chi1_str}")
                
            except Exception as e:
                original_res = res_id + min_res - 1
                print(f"  Residue {original_res}: Error - {e}")
        
        # Store file results
        all_results[clean_filename] = file_results
        
        # Summary statistics
        if results:
            valid_chi1 = [r['chi1'] for r in results if r['chi1'] is not None]
            if valid_chi1:
                avg_chi1 = sum(valid_chi1) / len(valid_chi1)
                print(f"  Average Chi1: {avg_chi1:.2f}°")
            print(f"  Successfully analyzed {len(results)} residues")
        
    except Exception as e:
        print(f"  Error processing {filename}: {e}")
    
    print()

# Create DataFrame and save to CSV
print("Creating CSV file...")
df = pd.DataFrame.from_dict(all_results, orient='index')
df.to_csv('sidechain_angles.csv')
print(f"Saved results to sidechain_angles.csv")
print(f"DataFrame shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
