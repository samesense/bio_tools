#
# Copyright (C) 2004-2006 Victoria University of Wellington
#
# This file is part of the PFMFind module.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
# modified Feb'09 to run independent of GUI
# PE

from Bio import SubsMat
from Bio.SubsMat.FreqTable import FreqTable, FREQ
from Bio.Alphabet import IUPAC
from Bio.Align import Generic, AlignInfo
#from matrix import SCORE 

bg_dict = {'A': 0.075666295724018673, 'C': 0.016516288118360062,
           'E': 0.061870031303289649, 'D': 0.050642678152545494,
           'G': 0.069747892431969383, 'F': 0.040958473103989368,
           'I': 0.05779199724595098,  'H': 0.022834989231988598,
           'K': 0.056438990609013105, 'M': 0.023685622102844514,
           'L': 0.095945383365223735, 'N': 0.043386279391867545,
           'Q': 0.039429335425556614, 'P': 0.050869847166291719,
           'S': 0.073068362861969299, 'R': 0.053284643774968568,
           'T': 0.057493775368874442, 'W': 0.013726212301952652,
           'V': 0.065155475673901384, 'Y': 0.031487426645424192}

def henikoff_weights(seqs, alphabet, counts):
    order = {}
    for i, a in enumerate(alphabet):
        order[a] = i
    weights = [0.0]*len(seqs)
    w_sum = 0.0
    for i in range(len(counts)):
        r = 0
        for j in range(len(counts[i])):
            if counts[i][j] > 0.0: r += 1
        for k in range(len(seqs)):
            w = 1 / (r * counts[i][order[seqs[k][i]]])
            weights[k] += w
            w_sum += w
    c = len(seqs) / w_sum
    return [w*c for w in weights]

def freq_counts(seqs, alphabet, weights=None):
    if weights == None: weights = [1.0]*len(seqs)
    counts = [[0.0 for a in alphabet] for s in seqs[0]]
    order = {}
    for i, a in enumerate(alphabet):
        order[a] = i

    for i in range(len(seqs[0])):
        for j in range(len(seqs)):
            a = order[seqs[j][i]]
            counts[i][a] += weights[j]
    return counts

iteration = True
arg_list = [('Scale', 'real', 2.0),
            ('Weighting', ['None', 'Henikoff'], 'Henikoff'),
            ]

def get_matrix(seq_ls, scale, weight_type):

    seqs = seq_ls

    bcounts = freq_counts(seqs, "ACDEFGHIKLMNPQRSTVWY")
    if weight_type == 'None':
        weights = [1.0]*len(seqs)
    elif weight_type == 'Henikoff':
        weights = henikoff_weights(seqs, "ACDEFGHIKLMNPQRSTVWY",
                                   bcounts) 

    align = Generic.Alignment(IUPAC.protein)
    for i in range(len(seqs)):
        align.add_sequence("Seq #%d" % i, seqs[i], weight=weights[i])
    summary_align = AlignInfo.SummaryInfo(align)

    # Must get expected frequencies from our own (i.e. whole
    # database) background frequencies. Otherwise they would be
    # derived from the alignment which wouldn't be good if we
    # have a small sample
    ftab = FreqTable(bg_dict, FREQ)
    arm = SubsMat.SeqMat(summary_align.replacement_dictionary())
    lom = SubsMat.make_log_odds_matrix(arm, ftab, factor=scale,
                                       round_digit=0, keep_nd=0)

    PM = lom
    matrix_type = '?'#SCORE
    return PM, matrix_type, 0

def print_matrix(seq_ls_file):
    seq_ls = []
    f = open(seq_ls_file)
    for line in f:
        seq_ls.append(line.strip())
    f.close()
    [PM, matrix_type, i] = get_matrix(seq_ls, 10, 'None')
    print PM
