###############################
# Name:   conservation.py     #
# Author: Santiago Gil Begue. #
###############################

from sys import argv
from os.path import exists, isdir
from math import log
from Bio import SeqIO
import operator

# Available formats.
formats = {'abi': 'Reads the ABI "Sanger" capillary sequence traces files, including the PHRED' +
                  'quality scores for the base calls. This allows ABI to FASTQ conversion. Note ' +
                  'each ABI file contains one and only one sequence (so there is no point in indexing the file).',
           'ace': 'Reads the contig sequences from an ACE assembly file. Uses Bio.Sequencing.Ace internally.',
           'clustal': 'The alignment format of Clustal X and Clustal W.',
           'embl': 'The EMBL flat file format. Uses Bio.GenBank internally.',
           'fasta': 'This refers to the input FASTA file format introduced for Bill Pearson\'s FASTA tool, where ' +
                    'each record starts with a ">" line. Resulting sequences have a generic alphabet by default.',
           'fastq': 'FASTQ files are a bit like FASTA files but also include sequencing qualities. In Biopython, ' +
                    '"fastq" (or the alias "fastq-sanger") refers to Sanger style FASTQ files which encode PHRED ' +
                    'qualities using an ASCII offset of 33. See also the incompatible "fastq-solexa" and ' +
                    '"fastq-illumina" variants used in early Solexa/Illumina pipelines, Illumina pipeline 1.8 ' +
                    'produces Sanger FASTQ.',
           'fastq-solexa': 'In Biopython, "fastq-solexa" refers to the original Solexa/Illumina style FASTQ files ' +
                           'which encode Solexa qualities using an ASCII offset of 64. See also what we call the ' +
                           '"fastq-illumina" format.',
           'fastq-illumina': 'In Biopython, "fastq-illumina" refers to early Solexa/Illumina style FASTQ files ' +
                             '(from pipeline version 1.3 to 1.7) which encode PHRED qualities using an ASCII offset ' +
                             'of 64. For good quality reads, PHRED and Solexa scores are approximately equal, so + '
                             'the "fastq-solexa" and "fastq-illumina" variants are almost equivalent.',
           'genbank': 'The GenBank or GenPept flat file format. Uses Bio.GenBank internally for parsing.', 
           'ig': 'This refers to the IntelliGenetics file format, apparently the same as the MASE alignment format.',
           'imgt': 'This refers to the IMGT variant of the EMBL plain text file format.',
           'nexus': 'The NEXUS multiple alignment format, also known as PAUP format. Uses Bio.Nexus internally.',
           'pdb-seqres': 'Reads a Protein Data Bank (PDB) file to determine the complete protein sequence as it ' +
                         'appears in the header (no dependency on Bio.PDB and NumPy).',
           'pdb-atom': 'Uses Bio.PDB to determine the (partial) protein sequence as it appears in the structure ' +
                       'based on the atom coordinate section of the file (requires NumPy).',
           'phd': 'PHD files are output from PHRED, used by PHRAP and CONSED for input. Uses Bio.Sequencing.Phd internally.',
           'phylip': 'PHYLIP files. Truncates names at 10 characters.',
           'pir': 'A "FASTA like" format introduced by the National Biomedical Research Foundation (NBRF) for the ' +
                  'Protein Information Resource (PIR) database, now part of UniProt.',
           'seqxml': 'Simple sequence XML file format.',
           'sff': 'Standard Flowgram Format (SFF) binary files produced by Roche 454 and IonTorrent/IonProton sequencing machines.',
           'stockholm': 'The Stockholm alignment format is also known as PFAM format.',
           'swiss': 'Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also the UniProt XML format.',
           'tab': 'Simple two column tab separated sequence files, where each line holds a record\'s identifier and ' +
                  'sequence. For example, this is used by Aligent\'s eArray software when saving microarray probes ' +
                  'in a minimal tab delimited text file.',
           'qual': 'Qual files are a bit like FASTA files but instead of the sequence, record space separated integer ' +
                   'sequencing values as PHRED quality scores. A matched pair of FASTA and QUAL files are often used ' +
                   'as an alternative to a single FASTQ file.',
           'uniprot-xml': 'UniProt XML format, successor to the plain text Swiss-Prot format.'}

# Check correct number of parameters.
error = False;
if len(argv) != 3 and len(argv) != 5:
	error = True
# Check if a valid format is given.
elif argv[2] not in formats:
	error = True
# Check if option -top is selected.
if len(argv) == 5:
	if argv[3] == '-top':
		isTop = True
		top = float(argv[4])
		if top < 0 or top > 1:
			print('Percentage must be in range [0-1]')
			exit(1)
	else:
		error = True
else: isTop = False

if error:
	print('Usage: python ' + argv[0] + ' <sequencies_file> <format> [ -top <percentage> ]')
	print('Tip: You can redirect it to a file by doing: ' + argv[0] + ' <sequencies_file> <format> [ -top <percentage> ] > output.txt\n')
	print('Available formats: ')
	print(formats.keys())
	print('\nUse: python ' + argv[0] + ' -h <format> for a format detail.')
	exit(1)

# Check if option is -h.
if argv[1] == '-h':
	print(formats[argv[2]])
	exit(1)

# Check that input file exists.
if not exists(argv[1]):
	print('> ERROR: File \'' + argv[1] + '\' doesn\'t exist.')
	print('  Exiting...')
	exit(1)
# Check that input file is not a directory.
if isdir(argv[1]):
	print('> ERROR: File \'' + argv[1] + '\' is a directory.')
	print('  Exiting...')
	exit(1)

# DNA molecule bases. IUPAC nucleotide code.
ADENINE 	= 'a'
CYTOSINE 	= 'c'
GUANINE 	= 'g'
THYMINE 	= 't'
URACIL		= 'u'
GAP1 		= '-'
GAP2		= '.'
DNA_BASES	= { ADENINE: [ADENINE],
                    CYTOSINE: [CYTOSINE],
                    GUANINE: [GUANINE],
                    THYMINE: [THYMINE],
                    URACIL: [THYMINE],
                    # Gaps.
                    GAP1 : [],
                    GAP2 : [],
                    # Combinations.
                    'r': [ADENINE, GUANINE],
                    'y': [CYTOSINE, THYMINE],
                    's': [GUANINE, CYTOSINE],
                    'w': [ADENINE, THYMINE],
                    'k': [GUANINE, THYMINE],
                    'm': [ADENINE, CYTOSINE],
                    'b': [CYTOSINE, GUANINE, THYMINE],
                    'd': [ADENINE, GUANINE, THYMINE],
                    'h': [ADENINE, CYTOSINE, THYMINE],
                    'v': [ADENINE, CYTOSINE, GUANINE],
                    'n': [ADENINE, CYTOSINE, GUANINE, THYMINE] }

# Amino acids. IUPAC amino acid code.
ALANINE           = 'a'
CYSTEINE          = 'c'
ASPARTIC_ACID     = 'd'
GLUTAMIC_ACID     = 'e'
PHENYLALANINE     = 'f'
GLYCINE           = 'g'
HISTIDINE         = 'h'
ISOLEUCINE        = 'i'
LYSINE            = 'k'
LEUCINE           = 'l'
METHIONINE        = 'm'
ASPARAGINE        = 'n'
PROLINE           = 'p'
GLUTAMINE         = 'q'
ARGININE          = 'r'
SERINE            = 's'
THREONINE         = 't'
VALINE            = 'v'
TRYPTOPHAN        = 'w'
TYROSINE          = 'y'
# Two additional amino acids: stop codons.
SELENOCYSTEINE    = 'u'
PYRROLYSINE       = 'o'
AMINO_ACIDS = { ALANINE: [ALANINE],
                CYSTEINE: [CYSTEINE],
                ASPARTIC_ACID: [ASPARTIC_ACID],
                GLUTAMIC_ACID: [GLUTAMIC_ACID],
                PHENYLALANINE: [PHENYLALANINE],
                GLYCINE: [GLYCINE],
                HISTIDINE: [HISTIDINE],
                ISOLEUCINE: [ISOLEUCINE],
                LYSINE: [LYSINE],
                LEUCINE: [LEUCINE],
                METHIONINE: [METHIONINE],
                ASPARAGINE: [ASPARAGINE],
                PROLINE: [PROLINE],
                GLUTAMINE: [GLUTAMINE],
                ARGININE: [ARGININE],
                SERINE: [SERINE],
                THREONINE: [THREONINE],
                VALINE: [VALINE],
                TRYPTOPHAN: [TRYPTOPHAN],
                TYROSINE: [TYROSINE],
                SELENOCYSTEINE: [SELENOCYSTEINE],
                PYRROLYSINE: [PYRROLYSINE],
                # Gaps.
                GAP1: [],
                GAP2: [],
                # Combinations
                'x': [ALANINE, CYSTEINE, ASPARTIC_ACID, GLUTAMIC_ACID,
                      PHENYLALANINE, GLYCINE, HISTIDINE, ISOLEUCINE,
                      LYSINE, LEUCINE, METHIONINE, ASPARAGINE, PROLINE,
                      GLUTAMINE, ARGININE, SERINE, THREONINE, VALINE,
                      TRYPTOPHAN, TYROSINE, SELENOCYSTEINE, PYRROLYSINE],
                'b': [ASPARTIC_ACID, ASPARAGINE],
                'z': [GLUTAMIC_ACID, GLUTAMINE],
                'j': [ISOLEUCINE, LEUCINE] }

# Array with the sequences to process.
try:
	sequences = list(SeqIO.parse(argv[1], argv[2]))
except Exception:
	print('> ERROR: Wrong file format.')
	print('  File \'' + argv[1] + '\' is not in \'' + argv[2] + '\' format.')
	print('  Exiting...')
	exit(1)

# Check that we have sequences to process. Format file may be incorrect.
if not sequences:
	print('> ERROR: No DNA sequences recognized.')
	print('  File \'' + argv[1] + '\' may not be in \'' + argv[2] + '\' format.')
	print('  Exiting...')
	exit(1)

# Check if the sequences are DNA or PROTEINS.
dna = 0; proteins = 0
for base in sequences[0].seq.lower():
	dna += base in DNA_BASES
	proteins += base in AMINO_ACIDS
isDNA = dna >= proteins

# Get maximum length of the sequences.
maxLength = max(len(sequence.seq) for sequence in sequences)

### DNA
if isDNA:
	# Dictionary with the occurrences of each base in each position.
	# bases[ADENINE][i] means how many sequences have an ADENINE base in the position i.
	bases = { ADENINE  : [0]*maxLength,
                  CYTOSINE : [0]*maxLength,
                  GUANINE  : [0]*maxLength,
                  THYMINE  : [0]*maxLength }
	target = DNA_BASES
### PROTEIN
else:
	# Dictionary with the occurrences of each amino acid in each position.
	# bases[ALANINE][i] means how many sequences have an ALANINE amino acid in the position i.
	bases = { ALANINE        : [0]*maxLength,
                  CYSTEINE       : [0]*maxLength,
                  ASPARTIC_ACID  : [0]*maxLength,
                  GLUTAMIC_ACID  : [0]*maxLength,
                  PHENYLALANINE  : [0]*maxLength,
                  GLYCINE        : [0]*maxLength,
                  HISTIDINE      : [0]*maxLength,
                  ISOLEUCINE     : [0]*maxLength,
                  LYSINE         : [0]*maxLength,
                  LEUCINE        : [0]*maxLength,
                  METHIONINE     : [0]*maxLength,
                  ASPARAGINE     : [0]*maxLength,
                  PROLINE        : [0]*maxLength,
                  GLUTAMINE      : [0]*maxLength,
                  ARGININE       : [0]*maxLength,
                  SERINE         : [0]*maxLength,
                  THREONINE      : [0]*maxLength,
                  VALINE         : [0]*maxLength,
                  TRYPTOPHAN     : [0]*maxLength,
                  TYROSINE       : [0]*maxLength,
                  SELENOCYSTEINE : [0]*maxLength,
                  PYRROLYSINE    : [0]*maxLength }
	target = AMINO_ACIDS

# Fill bases dictionary.
for sequence in sequences:
	for i in range(len(sequence.seq)):
		# No case sensitive.
		base = sequence.seq[i].lower()
		# Be robust to non DNA/PROTEIN format.
		if base in target:
			equivalentBases = target[base]
			for equivalentBase in equivalentBases:
				bases[equivalentBase][i] += 1.0 / len(equivalentBases)
		else:
			print('> ERROR: ' + sequence.id + ' may not be a DNA molecule or PROTEIN.')
			print('  Found \'' + base + '\' in position ' + str(i+1) + '.')
			print('  Exiting...')
			exit(1)

# Calculate conservation for each position.
conservations = {}
for i in range(maxLength):
	total = sum(bases[base][i] for base in bases)
	conservation = 0
	for base in bases:
		# When bases[base][i] == 0, we have 0*log(0), which is undefined.
		# So, avoid it, and add 0 (do nothing).
		if bases[base][i] != 0:
			freq = bases[base][i] / total
			conservation += freq * log(freq)
	# -top option selected.
	if isTop:
		if conservation < 0:
			conservations[i+1] = conservation
	# Normal output.
	else:
		print('Position ' + str(i+1) + ': ' + str(conservation))

# Return <top>% non-zero highest conservations.
if isTop:
	sorted_conservations = sorted(conservations.items(), key=operator.itemgetter(1), reverse=True)
	for i in range(int(len(sorted_conservations)*top)):
		print('Position ' + str(sorted_conservations[i][0]) + ': ' + str(sorted_conservations[i][1]))
        
