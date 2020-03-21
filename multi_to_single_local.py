#!usr/bin/env python3

#def get_args():
#    parser = argparse.ArgumentParser(description="Inputs to the script")
#    parser.add_argument("-f", "--FILE", help="FASTA FILE", required=True, type = str)
#    parser.add_argument("-m", "--OUTFILE", help="MOTIFs FILE", required=True, type = str)
#    return parser.parse_args()
#
#args = get_args()
#Fasta=args.FILE
#motifs=args.OUTFILE
import os

def parse_fasta(Fasta,outfile):
    single_dict={}
    with open(Fasta) as fa:
        for line in fa:
            line=line.strip()
            if line.startswith(">"):
                header=line
                single_dict[header]=""
            else:
                single_dict[header]+=line

    with open(outfile,"w") as o:
        for header,seq in single_dict.items():
            o.write(header)
            o.write("\n")
            o.write(seq)
            o.write("\n")
    return "Finished converting!"