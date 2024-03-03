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

def oneline_fasta(filename: str) -> dict:
    '''Reads a fasta file and returns a dictionary with header as keys and sequences as values'''
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

#3 classes -class gene - store gene name, class seq- stors gene length, class motif- identifies motif positions
# class for drawing  
class Gene:
    def __init__(self, the_name, the_header):
        # Initialize instance attributes
        self.name = the_name
        self.header = the_header
        # You can initialize other attributes here

class Sequence:
    def __init__(self, the_name, the_sequence):
        # Initialize instance attributes
        self.name = the_name
        self.sequence = the_sequence
        # You can initialize other attributes here

class MOTIF:
    def __init__(self, motif_file):
        # Initialize instance attributes
        self.motif_file = motif_file
        # You can initialize other attributes here

class Drawing:
    def __init__(self, width= 2000.6, height = 200):
        self.width = width
        self.height = height
        self.surface = cairo.PDFSurface("motif_mark.png", width, height)
        self.context = cairo.Context(self.surface)
        
    def draw_line(self, x1, y1, x2, y2, line_width=1):
        self.context.set_line_width(line_width)
        self.context.move_to(x1, y1)
        self.context.line_to(x2, y2)
        self.context.stroke()
        
    def draw_rectangle(self, x, y, width, height):
        self.context.rectangle(x, y, width, height)
        self.context.fill()
