from dna import DNA
import argparse

def cli():
    parser = argparse.ArgumentParser(prog = "dnaStat", description = "A simple program that produces statistical analysis of DNA sequences")
    parser.add_argument("file", help = "Input file for analysis", action = "store", type = str)
    parser.add_argument("odir", help = "Prefix for output directory", action = "store", type = str)
    parser.add_argument("-o", "--orf-cutoff", help = "Minimum length of ORF sequence (aa)", action = "store", type = str)
    args = parser.parse_args()

    dna = DNA(**vars(args))
    dna.driver()