#!/usr/bin/env python

import cairo
import math
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="open four programs simultaneously")
    parser.add_argument("-f", "--fastafilename", type=str, help="firstfilename", required=True)
    parser.add_argument("-m", "--motiffilename", type=str, help="secondfilename", required=True)
    return parser.parse_args()

args=get_args()
fasta_file= args.fastafilename
motif=args.motiffilename

fasta = {}

def oneline_fasta(filename: str) -> dict[str, str]:
    '''Reads a fasta file and returns a dictionary with header as keys and sequences as values'''

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

motifs={}
# 
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

with open ("Fig_1_motifs.txt", "r") as file:
    for motif in file:
        my_motif= motif.strip().upper()
        #setting motif- ygcy set as key of motifs dictionary
        motifs[my_motif]= motif_to_regex(my_motif)

#with open("test.fa", "r") as fh2:
    #seq- aaaactcgaaa
for sequence in fasta.items(): 
    matches = re.finditer( motif_to_regex("catag"), sequence, re.IGNORECASE)
    for match in matches:
        print("-----")
        print(match)
        #print(type(match))
        #print(dir(match))
        # print(match.start())
        # print(match.end())
        print(match.span())

"""
dictionary with a key - motifs in file ex- 'YgcY' 
with value - regex expression of that motif EX- '[CT]gc[CT]

"""

#3 classes -class gene - store gene name, 
    #class seq- stors gene length, 
    #groups for saving information

class Gene:
    def __init__(self, the_name: str, the_number: int, the_length: int):
        # Initialize instance attributes
        self.name = the_name 
        self.number = the_number
        self.length = the_length
    
    def __repr__(self) -> str:
        return f"Gene({self.name}, {self.number}, {self.length})"
    
    #def draw(self, the_context: what_type??):
    def draw(self, context):
        # Set color and line width for drawing gene
        context.set_source_rgb(0, 0, 0)  # Black color
        context.set_line_width(5)

        # Calculate start and end positions for drawing the gene line
        start_x = self.number * 100
        end_x = start_x + self.length

        # Draw gene line
        context.move_to(start_x, 50)
        context.line_to(end_x, 50)
        context.stroke()

        # Draw gene label
        context.move_to(start_x, 30)
        context.show_text(self.name)

# Create instances of Gene
g1 = Gene("INSR", 1, 24)
g2 = Gene("MNSB", 2, 21)

# Create an SVG surface and draw genes on it
with cairo.SVGSurface("genes.svg", 500, 100) as surface:
    context = cairo.Context(surface)
    context.set_source_rgb(1, 1, 1)  # White background
    context.paint()

    # Draw genes
    g1.draw(context)
    g2.draw(context)

    surface.write_to_png('gene3.png')
    
def identify_exons(sequence: str) -> list:
    exon_positions = [] # tuple that stores the start and end of exons
    exon_start = None
    # we need to loop over each character in the string to check for exons(capital bases)
    for i, base in enumerate(sequence):
        if base.isupper():
            if exon_start is None:
                exon_start = i # sets the start pos as the counter index i
        elif exon_start is not None: # at the end of the exon , append to list
            exon_positions.append((exon_start / len(sequence), i / len(sequence)))  # Normalize positions in refrance to the length of the sequence 
            exon_start = None # set back to 0
    if exon_start is not None: # continue to loop
        exon_positions.append((exon_start / len(sequence), len(sequence) / len(sequence)))  # Normalize positions
    return exon_positions

class Exon:
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.exons = self.identify_exons(sequence)

    def draw_exon_positions(context, exon_positions, gene_number, sequence_length, gene_name):
        # Draw vertical line for the length of the sequence
        context.move_to(0, gene_number * 100 + 50)
        context.line_to(sequence_length, gene_number * 100 + 50)
        context.set_source_rgb(0, 0, 0)  # Black color
        context.set_line_width(2)
        context.stroke()

        # Set color and line width for drawing exon positions
        context.set_source_rgb(1, 0, 0)  # Red color
        context.set_line_width(6)  # Increased line width for better visualization

        # Calculate y position for drawing exon positions
        y_position = gene_number * 100 + 50  # Adjust as needed

        # Draw vertical lines for exon positions
        # parses through the tuple created 
        for start, end in exon_positions:
            # to find the exact position of the exon in refrance to the total length of the seq
            start_x = start * sequence_length
            end_x = end * sequence_length
            # draw the exons
            context.move_to(start_x, y_position)  
            context.line_to(end_x, y_position)  
            context.stroke()

        # Annotate the figure with gene names
        context.set_source_rgb(0, 0, 0)  # Black color
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(15)
        context.move_to(10, gene_number * 100 + 25)  # Adjust position for gene name
        context.show_text(gene_name)

    # Example usage:
    fasta_dict = {
        "Gene1": "atgCAGTcGACGAaaaaTCaaaTCG",
        "Gene2": "atCGaaaGATCGaaaaaaTCGATCG",
        "Gene3": "atCGaaaaGACGTACGTaaaaaTCG",
    }

    # Define image width
    WIDTH = 800

    # Create ImageSurface for PNG
    with cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, 100 * len(fasta_dict)) as surface:
        context = cairo.Context(surface)
        context.set_source_rgb(1, 1, 1)  # White background
        context.paint()

        # Iterate through the fasta_dict and draw exon positions for each sequence
        for i, (gene_name, sequence) in enumerate(fasta_dict.items()):
             #gene_name = header.split()[0]
             exon_positions = identify_exons(sequence)
             sequence_length = len(sequence) * 10  # Scale factor for length of sequence
             draw_exon_positions(context, exon_positions, i, sequence_length, gene_name)

        # Save the drawing to a PNG file
        surface.write_to_png('exon_positions.png')

    print("Exon positions saved as exon_positions.png")
