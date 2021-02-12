#!/usr/bin/env python
# Eric Tran
# Motif Mark Assignment
# Bi 625 Winter 2021

import re
import itertools
import cairo
import argparse

def argument_parser():
    '''Returns parser object'''
    # Create parser object.
    parser = argparse.ArgumentParser()

    # Create flags for motif file and fasta file.
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Absolute file path \
        for fasta file')
    parser.add_argument('-m', '--motif', type=str, required=True, help='Absolute file path \
        for motif file)')
    return parser.parse_args()

def get_motif_coordinates(motif, sequence):
    '''Returns a list of tuples specifying the motif coordinates'''
    coordinate_tuples = []

    for match in re.finditer(motif, sequence, re.IGNORECASE):
        current_motif_coordinate = (match.start(), match.end() - 1)
        coordinate_tuples.append(current_motif_coordinate)
    return coordinate_tuples

def get_exon_coordinates(sequence):
    '''Returns a dictionary of the exon coordinates'''
    capital_char_indices = []
    exon_coordinates = {}

    for count, char in enumerate(sequence):
        if char.isupper():
            capital_char_indices.append(count)
    exon_coordinates["exon"] = (capital_char_indices[0], capital_char_indices[-1])

    return exon_coordinates

def get_ambigious_coordinates(motif_list, sequence):
    '''Returns a dictionary of ambigious reads and their coordinates'''

    # A list of all ambigious motifs.
    ambiguous_motif_list = []
    pattern = re.compile('^[ATGCUatgcu]+$')

    for motif in motif_list:
        if re.search(pattern, motif):
            continue
        else:
            ambiguous_motif_list.append(motif)
    
    # Now generate a dictionary of ambigious motif combinations and their coordinates
    # Key: Ambigious Motif
    # Value: Dictionary (key: motif, value: list of coordinate tuples)
    vIUPAC = {
    "A":["A"            ],
    "C":[    "C"        ],
    "G":[        "G"    ],
    "T":[            "T"],
    "U":[            "U"],
    "W":["A",        "T"],
    "S":[    "C","G"    ],
    "M":["A","C"        ],
    "K":[        "G","T"],
    "R":["A",    "G",   ],
    "Y":[    "C",    "T"],
    "B":[    "C","G","T"],
    "D":["A",    "G","T"],
    "H":["A","C",    "T"],
    "V":["A","C","G",   ],
    "N":["A","C","G","T"],
    "Z":[               ],
    }

    ambiguous_motif_coordinates = {}

    # The following piece of code was referenced from stackoverflow and modified.
    for ambiguous_motif in ambiguous_motif_list:

        groups = itertools.groupby(ambiguous_motif.upper(), lambda char:char not in vIUPAC)
        splits = []
        for b,group in groups:
            if b:
                splits.extend([[g] for g in group])
            else:
                for nuc in group:
                    splits.append(vIUPAC[nuc])
        combinations = [''.join(p) for p in itertools.product(*splits)]

        temporary_motif_coordinates = {}

        for motif in combinations:
            temporary_motif_coordinates[motif] = get_motif_coordinates(motif, sequence)
        
        all_coordinates = []

        for coordinate_list in temporary_motif_coordinates.values():
            all_coordinates += coordinate_list
        
        ambiguous_motif_coordinates[ambiguous_motif] = all_coordinates

    return ambiguous_motif_coordinates

def get_all_coordinates(motif_list, sequence):
    '''Returns a dictionary of motifs as keys and a list of their corresponding 
    coordinates as tuples'''

    # Key: motif
    # Value: List of tuples. tuple(beginning index, end index)
    motif_coordinates = {}
    exon_coordinates = get_exon_coordinates(sequence)
    ambigious_coordinates = get_ambigious_coordinates(motif_list, sequence)

    for motif in motif_list:
        motif_coordinates[motif] = get_motif_coordinates(motif, sequence)
    
    # Merge motif and exon dictionaries into one
    all_coordinates = {**motif_coordinates, **exon_coordinates}

    # Add ambiguous coordinates.
    for ambig in ambigious_coordinates:
        all_coordinates[ambig] = ambigious_coordinates[ambig]

    return all_coordinates

def create_drawing(coordinates_dictionary, longest_intron, motif_list, sequence_dictionary):
    '''This function creates an image of the gene, and its exons, introns, and motifs'''

    # Calculate the longest motif
    longest_motif = 0
    for motif in motif_list:
        length = len(motif)
        if length > longest_motif:
            longest_motif = length
    
    # Calculate Width and Height of the surface
    longest_motif_font = 12 * longest_motif
    left_border = longest_motif_font + 37
    gene_count = len(coordinates_dictionary)
    height = gene_count * 150
    width = left_border + 30 + longest_intron

    with cairo.SVGSurface("motif_visualization.svg", width, height) as surface:
        cxt = cairo.Context(surface)

        # Generate a motif color pallete.
        color_palette = {}
        color_palette['ygcy'] = [1, 0.2, 0.2, 0.6] # Peach
        color_palette['GCAUG'] = [4, 0, 4, 0.5] # Purple
        color_palette['catag'] = [0.0, 1, 0.0, 1] # Green
        color_palette['YYYYYYYYYY'] = [0.0, 0.0, 1, 1] # Blue

        motif_index = 0
        # Create a legend for the motifs and a box representing their colors.
        for name in motif_list:
            if name == 'GCATG':
                continue
            if name == 'cauag':
                continue

            motif_height = motif_index * 35
            cxt.set_source_rgba(color_palette[name][0], color_palette[name][1], color_palette[name][2], color_palette[name][3])
            cxt.set_line_width(3)
            cxt.rectangle(10, motif_height, 15, 15)
            cxt.fill()

            cxt.set_source_rgb(0.1, 0.1, 0.1)
            cxt.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            cxt.set_font_size(15)
            cxt.move_to(28, motif_height + 10)
            cxt.show_text(name)
            motif_index += 1

        gene_index = 0
        # Draw an image for each gene from the FASTA input file.
        for key in sequence_dictionary:

            # Draw intron line (black line)
            sequence_length = len(sequence_dictionary[key])
            gene_height = 100 + (gene_index * 150)
            cxt.set_source_rgb(0.1, 0.1, 0.1)
            cxt.set_line_width(3)
            cxt.move_to(left_border, gene_height)
            cxt.line_to(left_border + sequence_length, gene_height)
            cxt.stroke()

            # Mark the exon location
            exon_coordinates_tuple = coordinates_dictionary[key]["exon"]
            exon_start = exon_coordinates_tuple[0]
            exon_end = exon_coordinates_tuple[1]
            cxt.set_source_rgb(0.1, 0.1, 0.1)
            cxt.set_line_width(3)
            exon_width = exon_end - exon_start
            cxt.rectangle(left_border + exon_end, gene_height - 20, exon_width, 40)
            cxt.stroke()

            # Create a title for each gene using its header line.
            header_height = 45 + (gene_index * 150)
            cxt.set_source_rgb(0.1, 0.1, 0.1)
            cxt.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            cxt.set_font_size(25)
            cxt.move_to(left_border, header_height)
            cxt.show_text(key)

            # Mark the motifs.
            motif_coordinates = coordinates_dictionary[key]
            for motif in motif_coordinates:
                if motif == "exon":
                    continue
                if motif == 'GCATG':
                    motif = 'GCAUG'
                if motif == 'cauag':
                    motif = 'catag'
                for motif_tuple in motif_coordinates[motif]:
                    cxt.set_source_rgba(color_palette[motif][0], color_palette[motif][1], color_palette[motif][2], color_palette[motif][3])
                    cxt.set_line_width(3)
                    cxt.move_to(left_border + motif_tuple[0], gene_height - 20)
                    cxt.line_to(left_border + motif_tuple[0], gene_height  + 20)
                    cxt.stroke()
            gene_index += 1

def main():

    # Retrieve parser object for fasta and motif file paths.
    args = argument_parser()
    fasta_file_path = args.fasta
    motif_file_path = args.motif

    motif_list = []
    motif_file = open(motif_file_path, 'r')

    # Read motifs into list.
    for motif in motif_file:
        motif = motif.strip()
        motif_list.append(motif)

    # Run through the motif list again. If there is a U base, create a motif
    # of the DNA sequence.
    ut_motifs = []
    for motif in motif_list:
        if 'u' in motif:
            ut_motifs.append(motif.replace('u', 't'))
        elif 't' in motif:
            ut_motifs.append(motif.replace('t', 'u'))
        elif 'U' in motif:
            ut_motifs.append(motif.replace('U', 'T'))
        elif 'T' in motif:
            ut_motifs.append(motif.replace('T', 'U'))
    
    # Update the original motif list
    motif_list += ut_motifs

    # Read in Fasta file information.
    fasta_file = open(fasta_file_path, 'r')

    # Dictionary to store header and sequences
    sequence_dictionary = {}

    # Variable to track the current header line.
    current_header = "none"

    for line in fasta_file:
        line = line.strip()
        
        if line[0] == '>':
            current_header = line
            sequence_dictionary[current_header] = ""
        else:
            sequence_dictionary[current_header] = sequence_dictionary[current_header] + line

    # Use nested dictionary to store keys (header) and values dictionary (motif: list
    # of tuple lists with coordinates)
    # Also calculate the exon coordinates (all capital letters in fasta)
    coordinates_dictionary = {}

    for header in sequence_dictionary:
        sequence = sequence_dictionary[header]
        coordinates_dictionary[header] = get_all_coordinates(motif_list, sequence)
    
    for key in coordinates_dictionary:
        print(key)

    gene_name = ""
    longest_intron = 0
    
    # Call the create_drawing function.
    for gene in sequence_dictionary:
        length = len(sequence_dictionary[gene])
        if length > longest_intron:
            longest_intron = length
            gene_name = gene
    create_drawing(coordinates_dictionary, longest_intron, motif_list, sequence_dictionary)
main()