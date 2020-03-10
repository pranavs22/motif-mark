# motif-mark
Thie repository contains Python code to visualize motifs

About the script-
• First converts multiline fasta to single line fasta for ease of precessing

• Is compatible with Python3

• Uses argparse

• Accepts input as FASTA file and motifs file

• Outputs single figure

• Can handle multiple sequences and multiple motifs

• Can handle ambiguous motifs 

• Outputs svg image file

• Key/labeling


## Input-

1. A list of files containing sequences
2. A list of files containing motifs

## Output
It outputs a single vector-based image displaying all motifs present in the sequences.

## Example Run Command

python motif-3.py -f Figure_1.fasta -m Fig_1_motifs.txt -o C:/Bi624/motif-mark/motifs
