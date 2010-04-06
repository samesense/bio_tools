#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
import utils_graph, utils_motif, sys, random
"""
Functions for dealing with FASTA files.  The
format I'm using is:

>gene_name
seq
seq
seq
... .
"""

def getRandomDirtySeqs():
    """ Make random sequences.  Insert X's, *'s, and gaps. """

    residues = ['A', 'G', 'T', 'C', 'X', '*', '-']
    length = random.randint(10,100)
    num_seqs = random.randint(10,100)
    seqs = {}
    for seq_count in xrange(num_seqs):
        name = 'B.' + str(seq_count)
        seqs[name] = ''
        for pos in xrange(length):
            seqs[name] = seqs[name] + residues[random.randint(0,3)]
    return seqs

def getRandomSeqs():
    """ Make random sequences.  Insert X's, *'s, and gaps. """

    residues = ['A', 'G', 'T', 'C']
    length = random.randint(10,100)
    num_seqs = random.randint(10,100)
    seqs = {}
    for seq_count in xrange(num_seqs):
        name = 'B.' + str(seq_count)
        seqs[name] = ''
        for pos in xrange(length):
            seqs[name] = seqs[name] + residues[random.randint(0,3)]
    return seqs

def printFASTA_forGenes(gene_file, fasta_file, output_file):
    """ Given a file of genes and a FASTA file,
        print the fasta entries in the gene file.
    
    @param gene_file: one gene per line; you want seqs for these genes
    @param fasta_file: seqs
    @param output_file: where to stick the seqs
    """
    
    genes = {}
    gene_f = open(gene_file)
    for line in gene_f:
        genes[line.strip()] = True
    gene_f.close()

    fasta_f = open(fasta_file)
    fout = open(output_file, 'w')
    line = fasta_f.readline()
    while line != '':
        if line[0] == '>':
            if genes.has_key(line[1:].strip()):
                fout.write(line)
                line = fasta_f.readline()
                while line[0] != '>':
                    fout.write(line)
                    line = fasta_f.readline()
                    if line == '':
                        break
            else:
                line = fasta_f.readline()
                while line[0] != '>':
                    line = fasta_f.readline()
                    if line == '':
                        break
    fasta_f.close()
    fout.close()

def loadFASTA(fasta_file):
    """ Make a {} from a FASTA file.

    @param fasta_file: >name
                       seq
                       seq ...
    @return: fasta[gene] = seq (oneline)
    """

    fasta = {}
    fasta_f = open(fasta_file)
    name = ''
    seq = ''
    line = fasta_f.readline()
    while line != '':
        if line[0] == '>':
            if name != '':
                fasta[ name ] = seq
            name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = fasta_f.readline()
    fasta_f.close()

    if name != '':    
        fasta[name] = seq

    return fasta

def loadFASTA_forGenes(fasta_file, genes):
    """ Make a {} from a FASTA file.

    @param fasta_file: >name
                       seq
                       seq ...
    @param genes: {} of genes to get seqs for
    @return: fasta[gene] = seq (oneline)
    """

    fasta = {}
    fasta_f = open(fasta_file)
    name = ''
    seq = ''
    line = fasta_f.readline()
    while line != '':
        if line[0] == '>':
            if name != '' and genes.has_key(name):
                fasta[ name ] = seq
            name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = fasta_f.readline()
    fasta_f.close()

    if name != '' and genes.has_key(name):    
        fasta[name] = seq

    return fasta

def prettyPrint_runner(geneid, seq_string):
    """ Return FASTA formatted version of 
        this gene and its sequence.

    @param geneid: gene_name
    @param seq_string: FASTA seq on one line w/ no whitespace
    @return: formatted FASTA version
    """

    breakcount = 50
    line = '>' + geneid + '\n'
    while len(seq_string) > breakcount:
        line = line + seq_string[0:breakcount] + '\n'
        seq_string = seq_string[breakcount:]
    if seq_string != '':
        line = line + seq_string
    return line

def prettyPrint(geneid, seq_string):
    """ Print FASTA formatted version of this
        gene and its sequence.

    @param geneid: gene_name
    @param seq_string: FASTA seq on one line w/ no whitespace
    """

    print prettyPrint_runner(geneid, seq_string)

def prettyPrint_file(geneid, seq_string, afile):
    """ Write FASTA formatted gene/seq to this open file.

    @param geneid: gene_name
    @param seq_string: FASTA seq on one line w/ no whitespace
    @param afile: an open output file
    """
    
    afile.write( prettyPrint_runner(geneid, seq_string) + '\n')

def dumpFASTA(fasta, outfile):
    outf = open(outfile, 'w')
    for gene in fasta.keys():
        prettyPrint_file(gene, fasta[gene], outf)
    outf.close()

def getPos2Residue(seqs):
    """ Given a list of sequences, return a {} of positions, where each position has a list of the residues found at it. """

    pos2residue = {}
    for seq in seqs:
        for residue_index in xrange(len(seq)):
            key = residue_index
            if not pos2residue.has_key(key):
                pos2residue[key] = []
            pos2residue[key].append(seq[residue_index])
    return pos2residue

def countResidue(residueList):
    counts = {}
    for residue in residueList:
        if not counts.has_key(residue):
            counts[residue] = 0
        counts[residue] += 1
    return counts

def getAllGappedCols(keys, fasta):
    """ Residues are 0 indexed.  Find columns with all gaps. """

    seq_len = len(fasta[keys[0]])
    gapped = []
    for redisue in xrange(seq_len):
        gapped.append(True)
    for protein in keys:
        for residue in xrange(seq_len):
            call = fasta[protein][residue]
            if call == '-' or call == '*' or call == 'X':
                pass
            else:
                gapped[residue] = False
    gaps = []
    for residue in xrange(seq_len):
        if gapped[residue]:
            gaps.append(residue)
    return gaps

def getAllGappedColsCausedByX(keys, fasta):
    """ Residues are 0 indexed.  Find columns with all gaps. """

    seq_len = len(fasta[keys[0]])
    gapped = []
    for redisue in xrange(seq_len):
        gapped.append(True)
    residue2calls = {}
    for protein in keys:
        for residue in xrange(seq_len):
            call = fasta[protein][residue]
            if call == '-' or call == '*' or call == 'X':
                if not residue2calls.has_key(residue):
                    residue2calls[residue] = {}
                residue2calls[call] = True
            else:
                gapped[residue] = False
    gaps = []
    for residue in xrange(seq_len):
        if gapped[residue]:
            # only include the col if X has caused the gap
            if residue2calls.has_key('X') and residue2calls.has_key(' ') \
               or residue2calls.has_key('*') and residue2calls.has_key(' '):
                gaps.append(residue)
    return gaps

def getGapSplits(gap_ls):
    """ Given a list of columns that are all gaps. 
        Return a list of start/stops for new alignments split at these gaps. """

    merged_list = []
    for gap in gap_ls:
        utils_motif.mergeMotifs(gap, gap, merged_list)
    return merged_list

def assertSameLength(seq_names, fasta):
    """ Given {} of protein2seq, make sure all seqs are the same length. """

    length = len(fasta[seq_names[0]])
    for protein in seq_names:
        if len(fasta[protein]) != length:
            sys.stderr.write(protein + ' seqs are not same length')
            i=1/0

def mkZeroes(length):
    """ Return a str with all 0s """

    zeroes = ''
    for i in xrange(length):
        zeroes = zeroes + '0'
    return zeroes

def mkOnes(length):
    """ Return a str with all 1s """

    zeroes = ''
    for i in xrange(length):
        zeroes = zeroes + '1'
    return zeroes

def mkXline(length, X):
    """ Return a str with all 1s """

    zeroes = ''
    for i in xrange(length):
        zeroes = zeroes + str(X) + ' '
    return zeroes.strip()

def addSpacers(seq):
    """ Put a space btwn all residues. """

    new_seq = ''
    for char in seq:
        new_seq = new_seq + char + ' '
    return new_seq.strip()
