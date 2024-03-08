#!/usr/bin/env python
'''My sweet tool documentation yo.'''

# This allows type checking to "see" stuff before they've been defined.
# Always make this the very first import of your file, top of the file.
from __future__ import annotations

import cairo
import math
import argparse
import re
#Global variables 
# to do - documentation???????
MARGIN_X = 100 + 50

COLOR_LIST=[
    (0,1, 1),   # no one knows?
    (1, 0, 0),  # Red
    (0, 1, 0),  # Green
    (0, 0, 1),  # Blue
    (1, 1, 0),  # Yellow
     
]

def get_args():
    parser = argparse.ArgumentParser(description="open four programs simultaneously")
    parser.add_argument("-f", "--fastafilename", type=str, help="firstfilename", required=True)
    parser.add_argument("-m", "--motiffilename", type=str, help="secondfilename", required=True)
    return parser.parse_args()

args=get_args()
fasta_file= args.fastafilename
motif=args.motiffilename

def oneline_fasta(filename: str) -> dict[str, str]:
    '''Reads a fasta file and returns a dictionary with header as keys and sequences as values'''
    # dict representation of a fasta file;
    # keys: headers, values: sequences
    fasta = {}

    with open(fasta_file, "r") as fh:
        header = None
        sequence = ""
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    fasta[header] = sequence
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        # Add the last sequence
        if header is not None:
            fasta[header] = sequence

    return fasta

# create a dictionary for motifs as key-motifs in file, values - identified motifs
#key- motifs in file
#value - regex expression

def motif_to_regex(motif: str) -> str:
    #ex- motif = ycgy --> "[CT]CG[CT]"
    motif_regex = ''
    for base in motif:
        if base.upper() == "A":
            motif_regex += "A"
        elif base.upper() == "C":
            motif_regex += "C"
        elif base.upper() == "G":
            motif_regex += "G"
        elif base.upper() == "T":
            motif_regex += "T"
        elif base.upper() == "U":
            motif_regex += "U"
        elif base.upper() == "R":
            motif_regex += "[AG]"
        elif base.upper() == "Y":
            motif_regex += '[CT]'
        elif base.upper() == "S":
            motif_regex += "[GC]"
        elif base.upper() == "W":
            motif_regex += "[AT]"
        elif base.upper() == "K":
            motif_regex += "[GT]"
        elif base.upper() == "M":
            motif_regex += "[AC]"
        elif base.upper() == "B":
            motif_regex += "[CGT]"
        elif base.upper() == "D":
            motif_regex += "[AGT]"
        elif base.upper() == "H":
            motif_regex += "[ACT]"
        elif base.upper() == "V":
            motif_regex += "[ACG]"
        elif base.upper() == "N":
            motif_regex += "[CTGAU]"
    return motif_regex

#print(motif_to_regex("ycgy"))

#key - a motif from the file : ygcy
#value - regex representation of motif: [CT]gc[CT]
motifs: dict[str, str] = {}

# key - a motif from the file : ygcy
#value - assigned color :  (1, 0, 0)
MOTIF_COLOR_MAP : dict[str, tuple[float, float, float]] = {}

# with open ("Fig_1_motifs.txt", "r") as file:
#     i=0
#     for motif in file:
#         my_motif= motif.strip().upper()
#         #setting motif- ygcy set as key of motifs dictionary
#         motifs[my_motif]= motif_to_regex(my_motif)
#         # set color to motifs
#         motif_color_map[my_motif]=COLOR_LIST[i]
#         i+=1

with open ("Fig_1_motifs.txt", "r") as file:
    for i, motif in enumerate(file):
        my_motif= motif.strip().upper()
        #setting motif- ygcy set as key of motifs dictionary
        motifs[my_motif]= motif_to_regex(my_motif)
        # set color to motifs
        MOTIF_COLOR_MAP[my_motif]=COLOR_LIST[i]
print(MOTIF_COLOR_MAP)

# def build_motifs(sequence: str) -> Motif:
#     motif_d = Motif(my_motif, match.start(), match.stop() )
#     for sequence in fasta_dict.items(): 
#         matches = re.finditer( motif_to_regex(sequence), sequence, re.IGNORECASE)
        
#         for match in matches:
#             return matches, match.start
        
def build_motifs(sequence: str, motifs: dict, gene_number: int) -> Motif:
    found_motifs = []
    for motif, regex_pattern in motifs.items():
        matches = re.finditer(regex_pattern, sequence, re.IGNORECASE)
        for match in matches:
            #motif_instance = Motif(motif, match.start(), len(sequence))
            motif_instance = Motif(motif, match.start(), match.end(), gene_number)
            #motif_instance = Motif(motif, match.span())
            found_motifs.append(motif_instance)
    return found_motifs
    #return motif_instance


class Motif:
    def __init__(self, the_type: str, the_start:int, the_end: int, the_gene_number:int):
        # Initialize instance attributes
        # type example - ycgy, yyyyyyyyy, (types of motifs)
        self.type = the_type
        self.start = the_start
        self.end = the_end
        self.gene_number = the_gene_number

    def __repr__(self) -> str:
        return f"Motif({self.type}, {self.start}, {self.end})"

    def motif_draw(self, context):
        # Set color and line width for drawing gene
        #color = self.get_color()
        print(f'debug 345 {self.start}')
        context.set_source_rgb(0, 0, 0)  # Black color
        context.set_line_width(30)
        context.move_to(self.start, self.gene_number * 100 + 50)
        context.line_to(self.end, self.gene_number * 100 + 50)
        # context.move_to(270, 28050)
        # context.line_to(280, 28050)
        # context.move_to(270, )
        # context.line_to(280, )
        context.stroke()

    def get_color(self) -> tuple:
            # Define a dictionary mapping motif types to RGB color values
            color_map = {
                "ygcy": (1, 0, 0),  # Red
                "GCAUG": (0, 1, 0),  # Green
                "catag": (0, 0, 1),  # Blue
                "YYYYYYYYYY": (1, 1, 0),  # Yellow
            }
            # If the motif type is in the color map, return its color; otherwise, return black
            return color_map.get(self.type, (0, 0, 0)) 

"""
dictionary with a key - motifs in file ex- 'YgcY' 
with value - regex expression of that motif EX- '[CT]gc[CT]

"""
#3 classes -class gene - store gene name, 
    #class seq- stors gene length, 
    #groups for saving information
def build_gene(header: str, sequence: str, gene_number: int) -> Gene:
    gene = Gene(header, gene_number, len(sequence))
    return gene

class Gene:
    def __init__(self, the_name: str, the_number: int, the_length: int):
        # Initialize instance attributes
        self.name = the_name 
        self.number = the_number
        self.length = the_length
    
    def __repr__(self) -> str:
        return f"Gene({self.name}, {self.number}, {self.length})"
    
    #def draw(self, the_context: what_type??):
    def gene_draw(self, context):
        # Set color and line width for drawing gene
        context.set_source_rgb(0, 0, 0)  # Black color
        context.set_line_width(2)
        context.move_to(0, self.number * 100 + 50)
        context.line_to(self.length, self.number * 100 + 50)
        context.stroke()
        # Draw gene label
        context.move_to(0, self.number * 100 + 40)
        context.show_text(self.name)

def build_exons(sequence: str, gene_number: int) -> list[Exon]:
    exons = [] # tuple that stores the start and end of exons
    exon_start = None
    # we need to loop over each character in the string to check for exons(capital bases)
    for i, base in enumerate(sequence):
        if base.isupper():
            if exon_start is None:
                exon_start = i # sets the start pos as the counter index i
        elif exon_start is not None: # at the end of the exon , append to list
            exon_stop = i - 1
            exon = Exon(exon_start, exon_stop, gene_number)
            exons.append(exon)
            #exons.append((exon_start / len(sequence), i / len(sequence)))  # Normalize positions in refrance to the length of the sequence 
            exon_start = None # set back to 0
    # if exon_start is not None: # continue to loop
        #exons.append((exon_start / len(sequence), len(sequence) / len(sequence)))  # Normalize positions
    return exons

class Exon:
    def __init__(self, the_start: int, the_stop: int, the_gene_number: int):
        # Initialize instance attributes
        self.start = the_start
        self.stop = the_stop
        self.gene_number = the_gene_number
    def __repr__(self) -> str:
        return f"Exon({self.start}, {self.stop})"

    #def draw_exon_positions(context, exon_positions, gene_number, sequence_length, gene_name):
    def draw_exon(self, context):
        # Draw vertical line for the length of the sequence
        context.move_to(self.start, self.gene_number * 100 + 50)
        context.line_to(self.stop, self.gene_number * 100 + 50)
        context.set_source_rgb(0, 0.5, 0)  # Black color
        context.set_line_width(10)
        context.stroke()

########
# Main #
########

# Define image width
WIDTH = 800
fasta_dict= oneline_fasta(fasta_file)
# Create ImageSurface for PNG
with cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, 100 * len(fasta_dict)) as surface:
    context = cairo.Context(surface)
    context.set_source_rgb(1, 1, 1)  # White background
    context.paint()

    # Iterate through the fasta_dict and draw exon positions for each     
    for gene_number, (gene_name, sequence) in enumerate(fasta_dict.items()):
            #draw motiffs
            motif_list = build_motifs(sequence, motifs, gene_number)
            for motif_instance in motif_list:
                print(motif_instance)
                motif_instance.motif_draw(context)
                #motif_instance.get_color(context)
            #draw build gene
            gene = build_gene(gene_name, sequence, gene_number)
            print(gene)
            gene.gene_draw(context)
            #draw build exon
            exons= build_exons(sequence, gene_number)
            for exon in exons:
                print(exon)
                exon.draw_exon(context)


    # Save the drawing to a PNG file
    surface.write_to_png('exon_positions.png')
