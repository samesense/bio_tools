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


from string import split, join, strip
import sys, os, copy
from transcendental import gamma, log, array, product, beta, lgamma, exp, zeros
from Bio.Align.AlignInfo import PSSM
from cStringIO import StringIO

# Background amino acid frequencies - as dictionary
## BKGRND_PROBS = {'A': 7.5610, 'C': 1.6504, 'D': 5.0605,
##                 'E': 6.1824, 'F': 4.0928, 'G': 6.9696,
##                 'H': 2.2818, 'I': 5.7749, 'K': 5.6397,
##                 'L': 9.5874, 'M': 2.3668, 'N': 4.3354,
##                 'P': 5.0832, 'Q': 3.9400, 'R': 5.3245,
##                 'S': 7.3014, 'T': 5.7451, 'V': 6.5107,
##                 'W': 1.3716, 'Y': 3.1464}

BKGRND_PROBS = {'A': 0.075666295724018673, 'C': 0.016516288118360062,
                'E': 0.061870031303289649, 'D': 0.050642678152545494,
                'G': 0.069747892431969383, 'F': 0.040958473103989368,
                'I': 0.05779199724595098,  'H': 0.022834989231988598,
                'K': 0.056438990609013105, 'M': 0.023685622102844514,
                'L': 0.095945383365223735, 'N': 0.043386279391867545,
                'Q': 0.039429335425556614, 'P': 0.050869847166291719,
                'S': 0.073068362861969299, 'R': 0.053284643774968568,
                'T': 0.057493775368874442, 'W': 0.013726212301952652,
                'V': 0.065155475673901384, 'Y': 0.031487426645424192}

def _print_vector(v, indent=0, floats_per_row=10, fp=sys.stdout,
                  indent_first_row=True):
    if indent_first_row: indent1 = indent
    else: indent1 = 0
    fp.write("%*s[" % (indent1, ""))
    k = range(len(v))
    for i in range(0, len(v), floats_per_row):
        for j in k[i: i+floats_per_row]:
            fp.write("%.7f, " % v[j])
        fp.write("\n%*s" % (indent+1, ""))
    fp.write("]")

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


def log_beta(x):
    return sum(lgamma(x)) - lgamma(sum(x))

class DirichletMix:
    def __init__(self, dic=None):
        if dic == None:
            self.mixture = []
            self.alpha = []
            self.order = {}
        else:
            self.__dict__.update(dic)
            self.order = {}
            for i, a in enumerate(self.alphabet):
                self.order[a] = i
            self.num_distr = len(self.mixture)
           
    def read(self, filename):
        fp = file(filename, 'r')
        actions = {'Name': "self.name = atoms[2]",
                   'Name=': "self.name = atoms[1]",
                   'Order': "self.alphabet = join(atoms[2:],'')",
                   'Order=': "self.alphabet = join(atoms[1:],'')",
                   'Mixture=': "self.mixture.append(float(atoms[1]))",
                   'Alpha=': "self.alpha.append(map(float, atoms[1:]))",
                   }

        for line in fp:
            atoms = split(strip(line))
            if not len(atoms): continue
            if atoms[0] == 'EndClassName': break
            if atoms[0] not in actions: continue
            else:
                exec(actions[atoms[0]])
        fp.close()
        for i, a in enumerate(self.alphabet):
            self.order[a] = i
        self.num_distr = len(self.mixture)

    def _print_as_PyDict(self, indent=0, floats_per_row=10, fp=sys.stdout,
                        indent_first_row=True):
        if indent_first_row: indent1 = indent
        else: indent1 = 0
        fp.write("%*s{'name': '%s',\n" % (indent1, "", self.name))
        fp.write("%*s 'alphabet': '%s',\n" % (indent, "", self.alphabet))
        fp.write("%*s 'mixture': " % (indent, ""))
        _print_vector(self.mixture, indent+12, floats_per_row, fp, False)
        fp.write(",\n")
        fp.write("%*s 'alpha': [" % (indent, ""))
        _print_vector(self.alpha[0], indent+11, floats_per_row, fp, False)
        fp.write(",\n")
        for i in range(1,len(self.alpha)):
            _print_vector(self.alpha[i], indent+11, floats_per_row, fp, True)
            fp.write(",\n")
        fp.write("%*s           ]}" % (indent, ""))

    # TO DO: FIX THE NOTATION SO THAT IT CORRESPONDS TO THE PAPER

    def _coeffs(self, k, n, sum_n):
        alpha = array(self.alpha[k][1:])
        sum_alpha = self.alpha[k][0]
        q = self.mixture[k]

        # Calculate the coefficient using logarithms
        log_coeff = log_beta(n+alpha) - log_beta(alpha)
        return q * exp(log_coeff)
        
##         P = 1.0
##         for i in range(len(alpha)):
##             if n[i]:
##                 P *= n[i] * beta(n[i], alpha[i])
##                 print n[i], alpha[i],  beta(n[i], alpha[i])
##         print P
##         B = sum_n * beta(sum_n, sum_alpha)
##         return q * B / P

    def _probs(self, i, coeff, n):
        a = array([self.alpha[k][i+1] for k in xrange(self.num_distr)])
        sum_a_j = array([sum(self.alpha[k][1:]) for k in xrange(self.num_distr)])
        return sum(coeff*(a+n[i])/(sum_a_j + sum(n)))

    def aa_vector(self, aa_dict):
        return [aa_dict[a] for a in self.alphabet]

    def block_counts(self, seqs, weights=None):
        return freq_counts(seqs, self.alphabet, weights)

        
    def _pos_probs(self, counts):
        n = array(counts)

        # *** NEW CODE FROM HERE ****
        sum_n = sum(n)
        p = array([0.0]*len(counts))
        log_coeffs = array([0.0]*self.num_distr)
        q = array(self.mixture)

        for j in xrange(self.num_distr):
            alpha_j = array(self.alpha[j][1:])
            log_coeffs[j] = log_beta(n+alpha_j) - log_beta(alpha_j)
            
        log_coeffs = log_coeffs - max(log_coeffs)
        coeffs = q * exp(log_coeffs)

        for i in xrange(len(counts)):
            a = array([self.alpha[j][i+1] for j in xrange(self.num_distr)])
            sum_a_j = array([sum(self.alpha[j][1:]) for j in xrange(self.num_distr)])
            N = (a+n[i])/(sum_a_j + sum(n))
            p[i] = sum(coeffs*N)

        return p / sum(p)    
        # **** END OF NEW CODE  **** 

##         coeff = array([self._coeffs(k, n, sum(n)) for k in range(self.num_distr)])
##         # coeff = coeff / sum(coeff)
##         X = array([self._probs(i, coeff, n) for i in xrange(len(counts))])
##         return X / sum(X)

    def _pos_log_odds(self, _pos_probs, bkgrnd, scale=1.0):
        n = len(_pos_probs)
        return [int(scale * log(_pos_probs[i]/bkgrnd[i]) / log(2.0)) for i in range(n)]

    def block_probs(self, block_counts):
        return [self._pos_probs(counts) for counts in block_counts]

    def block_log_odds(self, block_probs, bkgrnd, scale=1.0):
        return [self._pos_log_odds(probs, bkgrnd, scale) for probs in block_probs]
    
    def block2pssm(self, block_data, seq):
        pssm_info  = []
        for i in range(len(block_data)):
            score_dict = {}
            for a in self.alphabet:
                score_dict[a] = block_data[i][self.order[a]]
            pssm_info.append((seq[i], score_dict))
        return PSSM(pssm_info)

    def print_block_data(self, block_data, field_width=4, precision=4, dtype='int'):
        file_str = StringIO()
        if dtype == 'int':
            pdata = lambda x: "%*d " % (field_width, x)
        else:
            pdata = lambda x: "%*.*f " % (field_width, precision, x)
            
        file_str.write("   ")
        for i in range(len(block_data)):
            file_str.write("%*d " % (field_width, i))
        file_str.write("\n")
        for j, a in enumerate(self.alphabet):
            file_str.write(" %c " % a)
            for i in range(len(block_data)):
                file_str.write(pdata(block_data[i][j]))
            file_str.write("\n")
        return file_str.getvalue()
                      
def _Comp2PyScript(path, filename=None):
    """
    Converts all Dirichlet mixtures files (ending in comp)
    in the path to Python dictionaries and stores them
    in filename if provided (otherwise prints to stdout).
    """

    if filename == None:
        fp = sys.stdout
    else:
        fp = file(filename, 'w')

    names = os.listdir(path)
    ffunc = lambda s: s[-4:] == 'comp'
    names = filter(ffunc ,names)
    names.sort()

    fp.write("import DirichletMix\n\n")
    fp.write("get_mix = lambda name: DirichletMix.DirichletMix(DIR_MIX[name])\n") 
    fp.write("get_names = lambda : NAMES\n\n")

    k = range(len(names))
    fp.write("NAMES = [")
    for i in range(0, len(k), 3):
        for j in k[i: i+3]:
            fp.write("'%s', " % names[j])
        fp.write("\n%*s" % (9+1, ""))
    fp.write("]\n\n")
    
    fp.write("DIR_MIX = {\n")
    fp.write("    '%s': " % names[0])
    DM = DirichletMix()
    DM.read(names[0])
    DM._print_as_PyDict(4+len(names[0])+4, 4, fp, False)
    fp.write(",\n")
    
    for s in names[1:]:
        fp.write("%*s'%s': " % (4, "", s))
        DM = DirichletMix()
        DM.read(s)
        DM._print_as_PyDict(4+len(s)+4, 4, fp, False)                
        fp.write(",\n")
    fp.write("           }\n")
    
    if filename != None:
        fp.close()

