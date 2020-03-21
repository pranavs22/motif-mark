#! usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:46:13 2020

@author: Pranav
"""
#Imports

import cairo
import argparse
import re
import numpy as np
import itertools
import matplotlib as mpl

def get_args():
    parser = argparse.ArgumentParser(description="Inputs to the script")
    parser.add_argument("-f", "--FILE", help="FASTA FILE", required=True, type = str)
    parser.add_argument("-m", "--MOTIF_FILE", help="MOTIFs FILE", required=True, type = str)
    parser.add_argument("-o", "--OUTFILE", help="OUTPUT IMAGE FILE", required=True, type = str)

    return parser.parse_args()

args = get_args()
Fasta=args.FILE
motifs=args.MOTIF_FILE
out=args.OUTFILE

#outfile='C:/Bi624/motif-mark/single_fasta.txt'
#parse_fasta(Fasta,outfile)
global context 
global x


surface=cairo.SVGSurface(out+".svg", 1000, 800)
context = cairo.Context(surface)
rgb_colors={}

for name , hex in mpl.colors.cnames.items():
    rgb_colors[name]=mpl.colors.to_rgb(hex)

single_dict={}
def parse_fasta(Fasta):
    with open(Fasta) as fa:
            for line in fa:
                line=line.strip()
                if line.startswith(">"):
                    header=line
                    single_dict[header]=""
                else:
                    single_dict[header]+=line
    single_fasta=out+"_single_fasta.fa"
    with open(single_fasta,"w") as sfa:
        for header,seq in single_dict.items():
            sfa.write(header)
            sfa.write("\n")
            sfa.write(seq)
            sfa.write("\n")
    return single_fasta
def parse_data(Fasta,motifs):
    print("Creating Image...")
    
###Variabes and declaraations
    single_fasta=parse_fasta(Fasta)
    
    motif_dict={'A':'[Aa]','T':'[TtUu]','G':'[Gg]','C':'[Cc]','U':'[TtUu]',
                'Y':'[CcTtUu]','y':'[CcTtUu]','W':'[AaTtUu]','S':'[GgCc]', 
                'M':'[AaCc]','K':'[GgTtUu]','R':'[AaGg]','Y':'[CcTtUu]',
                'B':'[CcGgTtUu]','D':'[AaGgTtUu]','H':'[AaCcTtUu]','V':'[AaCcGg]',
                'N':'[AaCcGgTt]'}
    
    m_motif=sequence=''
    motif_list=[]
    fasta={}
    motif_back=seq_pos=motif_color={}    
    x=10

    context.move_to(100,10+x)
    context.show_text("Motif-Mark")
    context.move_to(750,10+x)
    context.show_text("Legend ")
    context.move_to(750,30+x)
    context.show_text("Color      Gene")
    context.rectangle(740,x,160,150)

    context.stroke()            
##########Process file
    with open (single_fasta) as f,open(motifs) as m:
        for line in f:
            if line.strip().startswith(">"):
                header=line.strip()
                fasta[header]=''

            else:
                sequence=line.strip()
                fasta[header]=sequence
        motif=m.readlines()
        c=np.random.choice(list(rgb_colors.keys()),len(motif))    
        c_index=0
        
        for mot in motif:
            mot=mot.strip().upper()

            motif_color[mot.strip()]=rgb_colors[c[c_index]]
            c_index+=1
        x=20
        c_index=0
        
        for mot in motif:
            x+=20
            
            mot=mot.strip().upper()
                
            context.set_source_rgb(motif_color[mot][0],motif_color[mot][1],motif_color[mot][2])
            context.rectangle(750,x+15,40,20)
            context.fill()
            context.stroke()
            context.stroke()
            
            context.set_source_rgb(0,0,0)                  #Text color
            context.move_to(800,25+x)
            context.show_text(mot)
            context.stroke()   
            

            for letter in mot: 
                m_motif += motif_dict[letter]
            motif_list.append(m_motif)
            motif_back[m_motif]=mot
            m_motif=''
        c_index+=1

    exons=re.compile('[ATCG]+')
    seq=fasta.values()
    
    comb=list(itertools.product(seq,motif_list))    
    x=20
    for header,seq in fasta.items():

        e=exons.search(seq)
        start=e.span()[0]
        stop=e.span()[1]
        x+=70
#### Legend
        context.set_source_rgb(0,0,0)                  #Line color 

        context.move_to(50,10+x)
        context.show_text("Gene: "+ re.sub('>','',header)[:4]+"                          Location: "+re.sub('>','',header)[4:])
        context.stroke()            

        context.move_to(40,25+x)    
        seq_pos[seq]=25+x       #1
        context.line_to(len(seq)+40,25+x)
        context.set_line_width(1)
        context.set_source_rgb(0,0,0)                  #Line color 
        context.stroke()                                

        context.set_source_rgb(0,0,0)                  #Exon colour
        context.rectangle(start+40,15+x,stop-start,20)
        context.fill_preserve()
        context.stroke()    




    x=60
    for info in comb:
        x+=8
        
        make_image(motif_back,info,seq_pos,x,motif_color,fasta)
    context.stroke()  
############################### Draw EXON
#    print(motif)

    
    surface.finish()   
#    print(seq_pos)
    
    return "Done" 
def make_image(motif_back,info,seq_pos,x,motif_color,fasta):
    motif_count={}
    motif=info[1] 
    seq=info[0] 
    sp=re.finditer(motif,seq.upper())
    

    start=[i.start() for i in sp]
    motif_count[motif_back[motif]]=0

    context.set_source_rgb(motif_color[motif_back[motif]][0],motif_color[motif_back[motif]][1],motif_color[motif_back[motif]][2])   
    
    for i in start:
#        print(motif_back[motif],seq)
        motif_count[motif_back[motif]]+=1
        context.rectangle(i+40,seq_pos[seq]-10,4,20)
        context.fill()
        context.stroke()
    
    for k,v in fasta.items():
        if v==seq:
            print(k,motif_count)
                

    return motif_count
print(parse_data(Fasta,motifs))
