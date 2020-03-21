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
import matplotlib.pyplot as plt
import matplotlib as mpl
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
#yy=gcaug=catag=ygcy=0
surface=cairo.SVGSurface("example.svg", 1000, 800)
context = cairo.Context(surface)
rgb_colors={}

for name , hex in mpl.colors.cnames.items():
    rgb_colors[name]=mpl.colors.to_rgb(hex)

c=np.random.choice(list(rgb_colors.keys()),4)    

def parse_data(Fasta,motifs):

    
    motif_dict={'A':'[Aa]','T':'[TtUu]','G':'[Gg]','C':'[Cc]','U':'[TtUu]',
                'Y':'[CcTtUu]','y':'[CcTtUu]','W':'[AaTtUu]','S':'[GgCc]', 
                'M':'[AaCc]','K':'[GgTtUu]','R':'[AaGg]','Y':'[CcTtUu]',
                'B':'[CcGgTtUu]','D':'[AaGgTtUu]','H':'[AaCcTtUu]','V':'[AaCcGg]',
                'N':'[AaCcGgTt]'}
    
    m_motif=sequence=''
    motif_list=[]
    fasta={}
    motif_back=seq_pos={}    
    x=10

    context.move_to(100,10+x)
    context.show_text("Motif-Mark")
#    context.set_font_size(0.5)
    context.move_to(750,10+x)
    context.show_text("Legend ")
    context.move_to(750,30+x)
    context.show_text("Color      Gene")
    context.rectangle(740,x,160,150)

    context.stroke()            

    with open (Fasta) as f,open(motifs) as m:
        for line in f:
            if line.strip().startswith(">"):
                header=line.strip()
                fasta[header]=''

            else:
                sequence=line.strip()
                fasta[header]=sequence
        motif=m.readlines()
        
        print(len(motif))
        

        x=20
        c_index=0
        motif_color={}
        for mot in motif:
            x+=20
            c_index+=1
            mot=mot.strip()
            if mot=='ygcy':
                
                context.set_source_rgb(rgb_colors[c[c_index]][0],rgb_colors[c[c_index]][1],rgb_colors[c[c_index]][2])
                context.rectangle(750,x+15,40,20)
                context.fill()
                context.stroke()
            if mot=='GCAUG':
                context.set_source_rgb(1,0,0)
                context.rectangle(750,x+15,40,20)
                context.fill()
                context.stroke()
            if mot=='catag':
                context.set_source_rgb(0,0,1)
                context.rectangle(750,x+15,40,20)
                context.fill()
                context.stroke()
            if mot=='YYYYYYYYYY':
                
                context.set_source_rgb(0,1,0)
                context.rectangle(750,x+15,40,20)
                context.fill()

#                context.move_to(760,25+x)
                context.show_text(mot)
                context.stroke()
            
            context.set_source_rgb(0,0,0)                  #Text color
            context.move_to(800,25+x)
            context.show_text(mot)
            context.stroke()            
                                    
            mot=mot.upper()

            for letter in mot: 
                m_motif += motif_dict[letter]
            motif_list.append(m_motif)
            motif_back[m_motif]=mot
            m_motif=''
#    print(motif_list)
#    motif_seq={}
    exons=re.compile('[ATCG]+')
    seq=fasta.values()
    comb=list(itertools.product(seq,motif_list))    
    x=20
    for header,seq in fasta.items():

        e=exons.search(seq)
        start=e.span()[0]
        stop=e.span()[1]
        print(start)
    


        x+=70

        context.move_to(50,10+x)
        context.show_text("Gene: "+ re.sub('>','',header)[:4]+"                          Location: "+re.sub('>','',header)[4:])
        context.stroke()            



        context.move_to(1,25+x)    
        seq_pos[seq]=25+x       #1
        context.line_to(len(seq),25+x)
        context.set_line_width(1)
        context.set_source_rgb(0,0,0)                  #Line color 
        context.stroke()                                

        context.set_source_rgb(0,0,0)                  #Exon colour
        context.rectangle(start,15+x,stop-start,20)
        context.fill_preserve()
        context.stroke()    




    x=20
    for info in comb:
        x+=8
    
        make_image(motif_back,info,seq_pos,x)
#        print(key,"c",seq)
#        print(make_image(seq,motif_list,x))
        
#    print(motif_seq)
#    make_image(motif_seq,x)
    context.stroke()  
############################### Draw EXON
#    print(motif)

    
    surface.finish()   
#    print(seq_pos)
    print("DONE")
    return 
def make_image(motif_back,info,seq_pos,x):
    ygcy=yy=gcaug=catag=0
    motif=info[1] 
    seq=info[0] 
    sp=re.finditer(motif,seq.upper())
    start=[i.start() for i in sp]
    if motif_back[motif]=='YGCY':
#        print(motif_back)
        ygcy+=1
        context.set_source_rgb(0,255,255)
        
        for i in start:
            context.rectangle(i,seq_pos[seq]-10,4,20)
            context.fill()
            context.stroke()
    if motif_back[motif]=='GCAUG':
#        print(motif_back)
        context.set_source_rgb(1,0,0)
        gcaug+=1
        for i in start:
#            print(i)
            context.rectangle(i,seq_pos[seq]-10,4,20)
            context.fill()
            context.stroke()
    if motif_back[motif]=='YYYYYYYYYY':
#        print(start,seq)
        context.set_source_rgb(0,1,0)
        yy+=1
        for i in start:
#            print(i)
            context.rectangle(i,seq_pos[seq]-10,4,20)
            context.fill()
            context.stroke()
    if motif_back[motif]=='CATAG':
#        print(motif_back[motif])
        context.set_source_rgb(0,0,1)
#        print(start)
        catag+=1
        for i in start:
#            print(i)
            print(start)
            context.rectangle(i,seq_pos[seq]-10,4,20)
            context.fill()
            
        context.stroke()
            
            
            
    exons=re.compile('[ATCG]+')
#    print(motif)
    e=exons.search(seq)
#        
#    start=e.span()[0]
#    stop=e.span()[1]
#    exon=seq[start:stop]
#    intron_start=[]
#    exon_start=[]
##    if len(intron_start) != 0 or len(exon_start) != 0:
#################Draw Lines
#    x+=30
#    if len(start)!=None:
#        print('a')
#    if str(motif)=='[CcTtUu][Gg][Cc][CcTtUu]':
#        context.set_source_rgb(1,0,1)
#    if str(motif)=='[CcTtUu][Gg][Cc][CcTtUu][Gg][Cc][Aa][TtUu][Gg]':
#        print('yes',motif)
#        context.set_source_rgb(255,160,122)                  #red colour
#
#        
##        context.set_source_rgb(1,0,1)                  #Rectangle  color
#        
#        context.rectangle(start,15+x,stop-start,20)
#        context.fill_preserve()
#        context.stroke()    
##    for i in intron_start:
# 
#        
#    for i in exon_start:
#
#        context.rectangle(i+start,15+x,4,20)
#        context.set_source_rgb(205,92,92)
#        context.set_line_width(0.05)
#        context.stroke()    
#        
#        context.set_source_rgb(0,0,0)                  #Exon colour
#        context.rectangle(start,15+x,stop-start,20)
#        context.fill_preserve()
#        context.stroke()    
#
#
#        context.fill_preserve()
#        context.set_source_rgb(0,0,0)                  #Exon-motif colour
#        context.stroke()    
#
#    
    return motif_back

    
    

print(parse_data(Fasta,motifs))
    