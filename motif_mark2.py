#!/usr/env/bin python3

import os
import re 
import math
import cairo 
import random
import argparse
import numpy as np
from Bio import Seq 
from itertools import product
from collections import defaultdict

def get_arguments():
	parser = argparse.ArgumentParser(description="Motif Mark (MM) draws genes and motifs (to scale) given a FASTA and motif file. MM outputs a series of SVG images with introns, exons, and motifs clearly labeled. MM 1.0 takes in unlimited number of genes, but can only mark up to 16 unique motifs with the built-in color palette.")
	parser.add_argument("-f", "--file", type=str, help="Select a FASTA file with introns as lower-case letters and exons as upper-case letters.")
	parser.add_argument("-m", "--motif_file", type=str, help="Select a file with motifs separated by newlines. Ambiguous nucleotides are allowed.")
	return parser.parse_args()

args = get_arguments()
input = args.file
motif_file = args.motif_file

def CreateOutputFolder():
	if 'output' not in os.listdir():
		os.mkdir('./output') # creates output directory at the working directory
		print('output folder created in current working directory')
	else:
		print('output folder detected')

CreateOutputFolder()

# Calculate all combinations of ambiguous DNA sequences with biopython.
def extend_ambiguous_dna(seq):
	"""return list of all possible sequences given an ambiguous DNA input"""
	d = Seq.IUPAC.IUPACData.ambiguous_dna_values
	return(list(map("".join, product(*map(d.get, seq))))) # product() information... https://docs.python.org/2/library/itertools.html

# Initialize motifs
MotifDict = {} 
# Key: provided motif e.g. 'YGCY' | Value: ['TGCC', 'TGCT', 'CGCT', 'CGCC']
NumberOfMotifs = 0

def initialize_motifs(motif_file):
	with open(motif_file, "r") as motif_file:
		AmbiguousNucleotides = ('R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', 'Z', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n', 'z')
		global NumberOfMotifs # accesses global variable NumberOfMotifs
		for line in motif_file:
			line = str(line.strip().upper()) 
			NumberOfMotifs += 1
			# If motif contains any ambiguous bases, calculate all combinations and add to dict as list in values. elif motif is explicit, just add to dict. 
			if any(x in line for x in AmbiguousNucleotides):
				combos = extend_ambiguous_dna(line)
				MotifDict[str(line)] = extend_ambiguous_dna(line)
			else:
				MotifDict[line] = [line]
	return(MotifDict)

initialize_motifs(motif_file)

# Create FASTA dictionary
Gene = {}
# Key: gene header | Value: sequence of a gene in one line 
sequence = ''
LongestGene = int() # PyCairo canvas width
NumberOfGenes = 0 # PyCairo canvas height 

with open(input, "r") as fh:
	for line in fh:
		line = line.strip()
		# print(line)
		if line.startswith(">"):
			header = line[1:]
			Gene[header] = []
			NumberOfGenes += 1
			sequence = '' # resets sequence for next gene. 
		else:
			sequence = sequence + line
			if (len(sequence) > LongestGene): # iteratively compare and select the longest gene 
				LongestGene = len(sequence)
				Gene[header] = sequence
			else:
				Gene[header] = sequence

# Mark introns / exons / motifs
MotifStartSites = defaultdict(list) 
# Key: motif e.g. 'YGCY' | Value: [start sites] e.g. [40, 100, 193, 500]
def MarkGeneFeatures(seq):
	IntronSpan = [m.span() for m in re.finditer("[a-z]+", seq)]
	ExonSpan = [m.span() for m in re.finditer("[A-Z]+", seq)]
	for motif in MotifDict: # loop over motif 
		for submotif in MotifDict[motif]: # loop over sub-motifs within one motif (e.g. find 'tgcc' sub-motif within the 'ygcy' motif)
			Starts = [m.start() for m in re.finditer(submotif, seq, flags = re.IGNORECASE)] # return all start sites as list. one start site per motif
			MotifStartSites[motif].extend(Starts)
			# print('key', key, 'intron span', IntronSpan, 'exon span', ExonSpan, 'motif', motif, 'submotif', submotif, 'Starts', Starts)
	return(IntronSpan, ExonSpan, MotifStartSites)

# Prepare PyCairo instructions
GeneFeatures = defaultdict(list)
# Key: gene headers | Values: [reverse complementarity (0 = false, 1 = true), length of gene as int, [span of introns as pairs of tuples (start, end) using ints], [span of exons e.g. (start, end)], [motif start sites using ints]]
for key, seq in Gene.items():
	# Split workflow if gene is reverse complemented or not
	if 'reverse complement' in key:
		reverse = 1
		gene_len = len(seq)
		# key = key[:-22] # remove '(reverse complement)'
		GeneFeatures[key].append(reverse); GeneFeatures[key].append(gene_len)
		IntronSpan, ExonSpan, MotifStartSites = MarkGeneFeatures(seq)
		GeneFeatures[key].extend([IntronSpan]); GeneFeatures[key].extend([ExonSpan]); GeneFeatures[key].extend([MotifStartSites])
		MotifStartSites = defaultdict(list) # resets MotifStartSites so two genes don't get the same motif
	else:
		reverse = 0
		gene_len = len(seq)
		GeneFeatures[key].append(reverse); GeneFeatures[key].append(gene_len)
		IntronSpan, ExonSpan, MotifStartSites = MarkGeneFeatures(seq)
		GeneFeatures[key].extend([IntronSpan]); GeneFeatures[key].extend([ExonSpan]); GeneFeatures[key].extend([MotifStartSites])
		MotifStartSites = defaultdict(list) # resets MotifStartSites so two genes don't get the same motif

# Define a Color Palette 
MotifColors = {}
# Key: motif e.g. 'YGCY' | Value: (normalized rgb values) e.g. (0.3, 0.6, 0.1)
# Decided not to randomly generate colors (RGC) per motif because two RGC's may be too close in shade, so user must re-run the script again which is bad.
cols = [[0.7529411764705882, 0.2235294117647059, 0.16862745098039217], [0.9058823529411765, 0.2980392156862745, 0.23529411764705882], [0.6078431372549019, 0.34901960784313724, 0.7137254901960784], [0.5568627450980392, 0.22745098039215686, 0.6784313725490196], [0.1607843137254902, 0.5019607843137255, 0.7254901960784313], [0.20392156862745098, 0.596078431372549, 0.8588235294117647], [0.10196078431372549, 0.7372549019607844, 0.611764705882353], [0.08627450980392157, 0.6274509803921569, 0.5215686274509804], [0.15294117647058825, 0.6823529411764706, 0.3764705882352941], [0.1803921568627451, 0.8, 0.44313725490196076], [0.9450980392156862, 0.7686274509803922, 0.058823529411764705], [0.9529411764705882, 0.611764705882353, 0.07058823529411765], [0.9019607843137255, 0.49411764705882355, 0.13333333333333333], [0.8274509803921568, 0.32941176470588235, 0.0], [0.4980392156862745, 0.5490196078431373, 0.5529411764705883], [0.20392156862745098, 0.28627450980392155, 0.3686274509803922]]
# cols come from sixth row of shades from color chart at https://htmlcolorcodes.com/

if NumberOfMotifs > len(cols):
	print('ERROR: Number of Motifs exceeds built-in color palette (16)')
	exit()

temp = random.sample(cols, NumberOfMotifs) # randomly sample without replacement the number of motifs

# assign a motif a color set
i = 0
for header in MotifDict:
	MotifColors[header] = temp[i]
	i += 1

# Draw Introns + Exons + Motifs + Legends
width, height = LongestGene + 500, 500 # width depends on longest_gene. height depends on number of genes. Series of output images will be to scale.
for key, val in GeneFeatures.items():
	filename = './output/' + key + '.svg'
	surface = cairo.SVGSurface(filename, width, height)
	context = cairo.Context(surface)
	# Draw gene name 
	context.move_to(100, 20)
	context.select_font_face("Calibri", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL); context.set_font_size(20); context.set_source_rgb(0, 0, 0)
	context.show_text(key)
	# Draw introns
	context.set_line_width(2)
	context.move_to(100, (height / 2))
	context.line_to(100 + val[1], (height / 2))
	context.stroke()
	# Draw exons
	for i in val[3]: # iterate over exon start / end sites
		start = i[0]
		end = i[1]
		# print('start', start, 'end', end)
		context.set_source_rgb(0,0,0)
		context.stroke()
		context.rectangle(100 + start, height/2 - 50, end - start, height/2 - 150)
		context.fill()
	# Draw motifs - sorry about the extensive looping :< 
	for motif, start_sites in val[4].items(): # iterate over motif dictionary
		for site in start_sites: # iterate over sub-motifs list in motif dictionary
			for color in MotifColors: # color the motifs based on pre-determined colors 
				col1, col2, col3 = MotifColors[motif]; context.set_source_rgb(col1, col2, col3) # retrieve and set the colors associated with a motif
				context.rectangle(100 + site, height/2 - 50, len(motif), height/2 - 150)
				context.fill()
	# Draw legend
	context.set_source_rgb(0, 0, 0)
	context.set_line_width(5)
	### Box 
	context.rectangle(width - 300, 0 + 10, 290, 0 + 40 * NumberOfMotifs); context.stroke() # (canvas width - margin, 0 + height margin, ___, ___)
	sep = 1
	### Motif Colors + Text
	for motif in MotifColors:
		col1 = MotifColors[motif][0]; col2 = MotifColors[motif][1]; col3 = MotifColors[motif][2]
		context.set_source_rgb(col1, col2, col3)
		context.rectangle(width - 300 + 200, 0 + 10 + 15 * sep, 80, 20); context.fill() # (canvas width - canvas margin + legend width margin, 0 + legend box + legend height margin, height of image)
		context.select_font_face("Calibri", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL); context.set_font_size(15); context.set_source_rgb(0, 0, 0)
		context.set_source_rgb(0, 0, 0)
		context.move_to(width - 300 + 15, 0 + 10 + 15 + 15 * sep)
		context.show_text(motif)
		sep += 2
	surface.finish()

print("Done")