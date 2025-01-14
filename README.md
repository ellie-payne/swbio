## README

This pipeline is written in Python and can be used to extract information from a DNA sequence, namely: open reading frames (ORFs), GC content (%), and restriction endonuclease recognition sites; including their count and position in the sequence.

## Packages required:
 - biopython
 - itertools

## How to use:
The program will ask how you would like to submit your query. In order to directly paste a sequence, you should select 'string' and paste your sequence to the command line. You will then be prompted to submit a sequence name. Once the analysis is complete, the program will generate a '.fasta' file of your sequence with your chosen name, and a '.txt' file named 'output_results.txt' containing the result of your analysis. If you would like to directly submit a '.fasta' file containing one or more sequences, you should select 'file' when asked how you would like to submit your query. Then you will be prompted to provide a path to your file. If you choose to submit a '.fasta' file, the only output from the program will be your analysis in an 'output_results.txt' file. 