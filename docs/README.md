# :snake: dnaStat :snake:

## Overview

dnaStat is a python package for the analysis DNA sequence data. Tools can be selected to produce graphical analysis, locate ORFs and translate DNA sequences into all reading frames. 

Additional usage information can be found using the command 

```python3 dnaStat.py --help```

Dependencies can be install using the command

```pip3 install -r ./docs/requirements.txt```

The package was built as part of CSC8330 at Newcastle University.

## Arguments
| Argument | Description | Default | Required |
| -------- | ----------- | ------- | -------- |
| --help | dnaStat help information | n/a | n |
| file | input file for dnaStat analysis | n/a | y |
| out-dir | prefix for output directory | n/a | y |
| --min-orf | minimum length of ORF sequence | 90 | n |
| --basic-stat | basic statistal analysis of input sequences | False | n |
| --complex-stat | complex statistical analsysis of input sequences | False | n |
| --translate | translate each sequence from input into all reading frames and save to files | False | n |
| --save-orfs | locate all ORFs in each sequence from input and save to files | False | n |

## Examples

The file hoxC_sequences.fa is provided for example usage. The file contains 15 hoxC DNA sequences in .fa format from various organisms. At present, dnaStat is limited to .fa files. Please run all commands from the scripts directory or use relative paths.

Run basic statistical analysis 

```python3 dnaStat.py hoxC_sequences.fa test --basic-stat```
