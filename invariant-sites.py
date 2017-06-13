#!/usr/bin/env python3
# Script by JK

# Usage
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import re
from Bio import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import pandas as pd
from pandas import DataFrame
import numpy as np
from collections import Counter

# Functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Check file exists
def check_file(f):
	return os.path.isfile(f)

# Check if file is in FASTA format
def check_fasta(f):
	if not os.path.isfile(f) or os.path.getsize(f) < 1:
		return False
	with open(f, 'r') as fasta:
		if fasta.readline()[0] != '>':						# Check if header starts with ">"
			return False
		for line in fasta:
			line = line.strip()
			if not line or line[0] == '>':	
				continue
			if bool(re.search('[^ACTGactgNn\?\-]', line)):	# Check if there are non-nucleotide characters in sequence
				return False
	return True

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Counts and masks non-invariant sites from multi-FASTA alignment\n',
	usage='\n  %(prog)s [--out invariant.fa] FASTA')
parser.add_argument('fasta', metavar='FASTA', nargs=1, help='original multi-FASTA alignment file')
parser.add_argument('--out', metavar='FILE', nargs=1, help='specify output file with non-invariant sites masked')
parser.add_argument('--mask', metavar='X', default='X', nargs=1, help='symbol to use for masking (default=X)')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
args = parser.parse_args()

# Check input/output files
if not check_file(args.fasta[0]):
	err('ERROR: Cannot find "{}". Check file exists in the specified directory.'.format(args.fasta[0]))
if not check_fasta(args.fasta[0]):
	err('ERROR: Check "{}" is in FASTA format.'.format(args.fasta[0]))
if args.out:
	if check_file(args.out[0]):
		err('ERROR: "{}" already exists.'.format(args.out[0]))

# Read alignment into pandas dataframe
msg('Reading {} into dataframe ... please wait'.format(args.fasta[0]))
aln = AlignIO.read(args.fasta[0], 'fasta')
df = pd.DataFrame(np.array([list(record) for record in aln], np.character, order="F"), dtype=str)

# Iterate through each column to mask columns containing non-invariant sites including Ns and gaps
msg('Masking non-invariant sites ...')
for col in df:
	if bool(re.search('[^ACGT]', ''.join(set(df[col])))) or len(set(df[col])) > 1:
		df[col] = args.mask[0]

# Export 1st row as consensus (given they are all the same after masking variants)
newseq = df.iloc[0].tolist()

# Count nucleotides and variant sites in sequence
print('Column counts by nucleotide: (non-invariant = {})'.format(args.mask[0]))
bases = list(Counter(newseq))
bases.sort()
counts = {}
for nt in bases:
	counts[nt] = Counter(newseq)[nt]
	print(''.join([nt, ': ', str(counts[nt])]))
print('Total sites: {}'.format(len(newseq)))
print('Proportion invariant = {}'.format(round((counts['A']+counts['C']+counts['G']+counts['T'])/len(newseq),4)))

# Export masked sequences to FASTA output
newseqs = []
for record in SeqIO.parse(args.fasta[0], 'fasta'):
	newseqs.append(SeqRecord(Seq.Seq(''.join(newseq)), id=record.id, description=record.description))

if args.out:
	msg('Masked sequences saved to "{}" ... '.format(args.out[0]))
	SeqIO.write(newseqs, args.out[0], 'fasta')

sys.exit(0)
