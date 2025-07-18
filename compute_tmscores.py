#!/usr/bin/env python3
import os
import re
import subprocess
import csv
import sys

def parse_tmscore_output(output: str):
    """
    Parse TM-score stdout and return a dict with metadata.
    """
    data = {}
    # Structure1 line: extract length1
    m1 = re.search(r"Structure1:.*Length=\s*(\d+)", output)
    if m1:
        data['length1'] = int(m1.group(1))
    # Structure2 line: extract length2
    m2 = re.search(r"Structure2:.*Length=\s*(\d+)", output)
    if m2:
        data['length2'] = int(m2.group(1))
    # RMSD line
    m3 = re.search(r"RMSD of  the common residues=\s*([\d\.]+)", output)
    if m3:
        data['rmsd'] = float(m3.group(1))
    # TM-score line
    m4 = re.search(r"TM-score\s*=\s*([\d\.]+)", output)
    if m4:
        data['tm_score'] = float(m4.group(1))
    # MaxSub-score
    m5 = re.search(r"MaxSub-score=\s*([\d\.]+)", output)
    if m5:
        data['maxsub'] = float(m5.group(1))
    # GDT-TS-score
    m6 = re.search(r"GDT-TS-score=\s*([\d\.]+)", output)
    if m6:
        data['gdt_ts'] = float(m6.group(1))
    return data

def compute_tmscores(target_path: str, folder: str = "rbd-structures"):
    target = os.path.basename(target_path)
    results = []
    
    # Gather all pdb files in folder
    files = [f for f in os.listdir(folder) if f.lower().endswith((".pdb", ".ent"))]
    
    for subject in files:
        if subject == target:
            continue
        subject_path = os.path.join(folder, subject)
        # run TM-score
        cmd = ["./TMscore", target_path, subject_path]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            out = proc.stdout
            meta = parse_tmscore_output(out)
            row = {
                "target": target,
                "subject": subject,
                **meta
            }
            results.append(row)
        except subprocess.CalledProcessError as e:
            print(f"Error running TM-score on {subject}: {e}", file=sys.stderr)
    
    # write to CSV
    csv_name = f"{os.path.splitext(target)[0]}_tmscores.csv"
    with open(csv_name, "w", newline="") as csvfile:
        fieldnames = ["target", "subject", "length1", "length2", "rmsd", "tm_score", "maxsub", "gdt_ts"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"target": target, "subject": target, "length1": 1273, "length2": 1273, "rmsd": 0, "tm_score": 1, "maxsub": 1, "gdt_ts": 1})
        for row in results:
            writer.writerow(row)
    
    print(f"Results saved to {csv_name}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python compute_tmscores.py path/to/target.pdb", file=sys.stderr)
        sys.exit(1)
    target_file = sys.argv[1]
    compute_tmscores(target_file)
