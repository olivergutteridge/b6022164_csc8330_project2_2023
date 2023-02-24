from dna import DNA
import argparse

def cli():
    parser = argparse.ArgumentParser(prog = "dnaStat", description = "A simple program that produces statistical analysis of DNA sequences")
    parser.add_argument("file", help = "Input file for analysis", action = "store", type = str)
    parser.add_argument("odir", help = "Prefix for output directory", action = "store", type = str)
    parser.add_argument("-o", "--orf-cutoff", help = "Minimum length of ORF sequence (aa)", action = "store", type = int, required = False, default = 30)
    parser.add_argument("-b", "--basic-stats", help = "Run basic statistical analysis on input sequences", action = "store_true", required = False, default = False)
    parser.add_argument("-c", "--complex-stats", help = "Run complex statistical analysis on input sequences", action = "store_true", required = False, default = False)
    parser.add_argument("-t", "--translate", help = "Translate all sequences in input into all size reading frames", action = "store_true", required = False, default = False)
    args = parser.parse_args()

    dna = DNA(**vars(args))
    dna.core()