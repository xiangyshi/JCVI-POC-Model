import subprocess
import sys

def main(variant):
    cat_cmd = f"cat sequences/spikes/Wuhan-A.fasta sequences/spikes/{variant}-* > align_in.fasta"
    align_cmd = f"mafft --clustalout align_in.fasta > align_out.fasta"
    alignf_cmd = f"mafft --auto --anysymbol align_in.fasta > alignf_out.fasta"
    aliview_cmd = f"aliview align_out.fasta"

    subprocess.run(cat_cmd, shell=True)
    subprocess.run(align_cmd, shell=True)
    subprocess.run(alignf_cmd, shell=True)
    subprocess.run(aliview_cmd, shell=True)

if __name__ == "__main__":
    variant = sys.argv[1]
    main(variant)