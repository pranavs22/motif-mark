# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:46:13 2020

@author: Pranav
"""
#Imports

import cairo
import argparse
import re
#from multi_to_single_local import *
# Parse the Fasta File
def get_args():
    parser = argparse.ArgumentParser(description="Inputs to the script")
    parser.add_argument("-f", "--FILE", help="FASTA FILE", required=True, type = str)
    parser.add_argument("-m", "--MOTIF_FILE", help="MOTIFs FILE", required=True, type = str)
    return parser.parse_args()

args = get_args()
Fasta=args.FILE
motifs=args.MOTIF_FILE
#outfile='C:/Bi624/motif-mark/single_fasta.txt'
#parse_fasta(Fasta,outfile)

global context 
global x
surface=cairo.SVGSurface("example.svg", 1000, 400)
context = cairo.Context(surface)
def parse_data(Fasta,motifs):

    motif_dict={'A':'[Aa]','T':'[TtUu]','G':'[Gg]','C':'[Cc]','U':'[TtUu]',
                'Y':'[CcTtUu]','y':'[CcTtUu]','W':'[AaTtUu]','S':'[GgCc]', 
                'M':'[AaCc]','K':'[GgTtUu]','R':'[AaGg]','Y':'[CcTtUu]',
                'B':'[CcGgTtUu]','D':'[AaGgTtUu]','H':'[AaCcTtUu]','V':'[AaCcGg]',
                'N':'[AaCcGgTt]'}
    
    m_motif=sequence=''
    motif_list=[]
    fasta={}
    x=10
    with open (Fasta) as f,open(motifs) as m:
        for line in f:

            if line.startswith(">"):
                header=line.strip()
            else:
                sequence=line.strip()
            fasta[header]=sequence
        motif=m.readlines()
        
        for mot in motif:
            mot=mot.strip().upper()
            for letter in mot: 
                m_motif += motif_dict[letter]
            motif_list.append(m_motif)
        m_motif=''
    for key,seq in fasta.items():
        for m in motif_list:
            x+=20
            print(make_image(seq,m,x))
    surface.finish()   
    return 
def make_image(seq,motif,x):

    num_motifs=0        
    exons=re.compile('[ATCG]+')
#    print(motif)
    e=exons.search(seq)
        
    start=e.span()[0]
    stop=e.span()[1]
    exon=seq[start:stop]
    intron_start=[]
    exon_start=[]
    intron1_motif=re.finditer(motif,seq)
    exon_motif=re.finditer(motif,exon)

    
    for i in intron1_motif:
        num_motifs+=1
        intron_start.append(i.start())
    
    for i in exon_motif:
        num_motifs+=1
        exon_start.append(i.start())
    if len(intron_start) != 0 or len(exon_start) != 0:
        context.set_source_rgb(0,0,0)                     #Line Colour
        context.move_to(10,25+x)                             #1
        context.line_to(len(seq)+x,25+x)
        context.set_line_width(1)
        context.set_source_rgb(0,0,0)                  #Text Colour
        context.move_to(700,25+x)
        context.show_text("".join(["Number of motifs"+str(num_motifs)]))
        
        context.stroke()

        context.stroke()            
        context.rectangle(start,15+x,stop-start,20)
        context.fill_preserve()
        context.stroke()    
 
    for i in intron_start:
        if motif=='[CcTtUu][Gg][Cc][CcTtUu]':
            context.set_source_rgb(0,0,0)                  #black colour
        if motif=='[Gg][Cc][Aa][TtUu][Gg]':
            context.set_source_rgb(255,160,122)                  #red colour
            context.stroke()    
            print('yes',motif)

        if motif=='[Cc][Aa][TtUu][Aa][Gg]':
            context.set_source_rgb(1,0,0)
            context.stroke()    
            print('yes',motif)

        if motif=='[CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu]':                  #Exon colour
            context.set_source_rgb(205,92,92)
            context.stroke()    
            print('yes',motif)

        context.rectangle(i,15+x,4,20)
        context.stroke()
        
    for i in exon_start:
        if str(motif)=='[CcTtUu][Gg][Cc][CcTtUu]':
            context.set_source_rgb(0,0,0)
            context.stroke()    
        if str(motif)=='[Gg][Cc][Aa][TtUu][Gg]':
            print('yes',motif)
            context.stroke()    
            context.set_source_rgb(255,160,122)                  #red colour
        if str(motif)=='[Cc][Aa][TtUu][Aa][Gg]':
            context.stroke()    
            context.set_source_rgb(1,0,0)
            context.stroke()    
        if str(motif)=='[CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu][CcTtUu]':                  #Exon colour
            context.set_source_rgb(205,92,92)

        context.rectangle(i+start,15+x,4,20)
        context.set_source_rgb(205,92,92)
        context.set_line_width(0.05)
        context.stroke()    
    
        context.set_source_rgb(0,0,0)                  #Exon colour
        context.rectangle(start,15+x,stop-start,20)
        context.fill_preserve()
        context.stroke()    


    context.fill_preserve()
    context.set_source_rgb(0,0,0)                  #Exon-motif colour
    context.stroke()    

    
    return 

    

print(parse_data(Fasta,motifs))
    