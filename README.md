# motif-mark
Draws genes and motifs given a FASTA and motif file

DESCRIPTION

Motif Mark (MM) draws genes and motifs (to scale) given a FASTA and motif file. 
Each gene in the FASTA file must contain a header, lower-case introns, and upper-case exons. The motif file must contain motifs (ambiguous nucleotides are allowed) separated by new lines. 
MM outputs a series of SVG images with introns, exons, and motifs clearly labeled. 
MM 1.0 takes in unlimited number of genes, but can only mark up to 16 unique motifs with the built-in color palette.

EXAMPLE

$ python motif_mark2.py -f MyGenes.fasta -m MyMotifs.txt

NOTE 

Have fun running MM 1.0! 
