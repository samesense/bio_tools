#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions for p-values of enrichment of motifs
by permutation tests.
"""

import random
import utils_graph, utils_motif, utils_fasta

def scoreProteinLs_motifLs_variants(alist, annotations, motifs, motif_variants):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist.keys():
        if annotations.has_key(protein):
            for motif in motifs.keys():
                match = False
                if annotations[protein].has_key(motif):
                    for [st, stp, seq] in annotations[protein][motif]:
                        if motif_variants.has_key(seq):
                            match = True
                if match:
                    scores[motif] += 1
    return scores

def scoreProteinLs_motifLs_count_all(alist, annotations, motifs):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist.keys():
        if annotations.has_key(protein):
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    scores[motif] += len(annotations[protein][motif])
    return scores

def scoreProteinLsFreqs_motifLs(alist, annotations, motifs, fasta):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """
    total_len = 0
    for protein in alist.keys():
        total_len += len(fasta[protein])
    scores = scoreProteinLs_motifLs_count_all(alist, annotations, motifs)    
    freq_scores = {}
    for motif in scores.keys():
        freq_scores[motif] = float(scores[motif]) / float(total_len)
    return freq_scores

def scoreProteinLs_motifLs_binary(alist, annotations, motifs):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist.keys():
        if annotations.has_key(protein):
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    scores[motif] += 1#len(annotations[protein][motif])
    return scores

def scoreProteinLs_motifLs_notBinary(alist, annotations, motifs):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist.keys():
        if annotations.has_key(protein):
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    scores[motif] += len(annotations[protein][motif])
    return scores

def scoreProteinLs_motifLs_binary_proteinLs(alist, annotations, motifs):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist:
        if annotations.has_key(protein):
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    scores[motif] += 1#len(annotations[protein][motif])
    return scores

def scoreProteinLs_motifLs_notBinary_proteinLs(alist, annotations, motifs):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist:
        if annotations.has_key(protein):
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    scores[motif] += len(annotations[protein][motif])
    return scores

def scoreProteinLs_signatureLs_binary_proteinLs(alist, annotations, motifs, window):
    """ Given a list of genes, and annotations
        for these genes, count the # of times
        each given motif is seen.  Motifs can
        occur multiple times on a gene.

    @param alist: {} of genes
    @param annotations: {} annotations[protein][motif] = [motif_instances]
    @param motifs: {} of motifs
    @return: {} [motif] = # of times seen
    """

    scores = {}
    for motif in motifs.keys():
        scores[motif] = 0
    for protein in alist:
        if annotations.has_key(protein):
            for motif in motifs.keys():
                [motif1, motif2] = motif.split(':')
                found_match = False
                if annotations[protein].has_key(motif1) and annotations[protein].has_key(motif2):
                    for [st1, stp1, seq] in annotations[protein][motif1]:
                        for [st2, stp2, seq] in annotations[protein][motif2]:
                            if max([stp1, stp2]) - min([st1, st2]) + 1 <= window:
                                found_match = True
                if found_match:
                    scores[motif] += 1#len(annotations[protein][motif])
    return scores

def scoreProteinLs(alist, annotations, motifs):
    scores = 0
    for protein in alist.keys():
        if annotations.has_key(protein):
            motif_count = 0
            for motif in motifs.keys():
                if annotations[protein].has_key(motif):
                    motif_count += 1
            if motif_count == len(motifs.keys()):
                scores += 1
    return scores

def scoreBckgnd_randomDraw(size, alist, annotations, motifs):
    bck_ls = {}
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls[ keys[rand_index] ] = True
    return scoreProteinLs(bck_ls, annotations, motifs)

def scoreBckgnd_motifLs_binary(size, alist, annotations, motifs):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = {}
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls[ keys[rand_index] ] = True
    return scoreProteinLs_motifLs_binary(bck_ls, annotations, motifs) 

def scoreBckgnd_motifLs_binary_proteinLs(size, alist, annotations, motifs):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = []
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls.append(keys[rand_index])
    return scoreProteinLs_motifLs_binary_proteinLs(bck_ls, annotations, motifs)

def scoreBckgnd_motifLs_notBinary_proteinLs(size, alist, annotations, motifs):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = []
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls.append(keys[rand_index])
    return scoreProteinLs_motifLs_notBinary_proteinLs(bck_ls, annotations, motifs)

def scoreBckgnd_signatureLs_binary_proteinLs(size, alist, annotations, motifs, window):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = []
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls.append(keys[rand_index])
    return scoreProteinLs_signatureLs_binary_proteinLs(bck_ls, annotations, motifs, window)

def scoreBckgndFreqs_motifLs(size, alist, annotations, motifs, fasta):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = {}
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls[ keys[rand_index] ] = True
    return scoreProteinLsFreqs_motifLs(bck_ls, annotations, motifs, fasta) 

def scoreBckgnd_motifLs_variants(size, alist, annotations, motifs, motif_variants):
    """ Make a background permutation of this size by
        randomly drawing from the given list.  Return
        motif counts for this permutation.    
    
    @param size: # of background draws for this permutation
    @param alist: background pool
    @param annotations: {}[protein][motif] = [motif_instance]
    @param motifs: motifs to get background scores for
    @return: {} [motif] = # of times seen in this permutation
    """

    bck_ls = {}
    lenls = len(alist) - 1
    keys = alist.keys()
    while len(bck_ls) < size:
        rand_index = random.randint(0, lenls)
        bck_ls[ keys[rand_index] ] = True
    return scoreProteinLs_motifLs_variants(bck_ls, annotations, motifs, motif_variants)        

def findEnrichedAnnotations_old(fg_ls, bg_ls, annotations, motifs):
    score = scoreProteinLs(fg_ls, annotations, motifs)
    better_count = 0
    for i in xrange(1000):
        bck_score = scoreBckgnd_randomDraw( len(fg_ls),  bg_ls, annotations, motifs)
        if bck_score >= score:
            better_count += 1
    return [score, len(fg_ls), better_count]

def findEnrichedAnnotations_motifLs_binary(fg_ls, bg_ls, annotations, motifs, permutations):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do X background
         permutations.

    @param fg_ls: [] of foreground genes
    @param bg_ls: [] of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """

    scores = scoreProteinLs_motifLs_binary(fg_ls, annotations, motifs)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(permutations):
        bck_scores = scoreBckgnd_motifLs_binary( len(fg_ls),  bg_ls, annotations, motifs)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores

def findEnrichedAnnotations_motifLs_binary_proteinLs(fg_ls, bg_ls, annotations, motifs, permutations):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do X background
         permutations.

    @param fg_ls: [] of foreground genes
    @param bg_ls: {} of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """

    scores = scoreProteinLs_motifLs_binary_proteinLs(fg_ls, annotations, motifs)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(permutations):
        bck_scores = scoreBckgnd_motifLs_binary_proteinLs( len(fg_ls),  bg_ls, annotations, motifs)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores

def findEnrichedAnnotations_motifLs_notBinary_proteinLs(fg_ls, bg_ls, annotations, motifs, permutations):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do X background
         permutations.

    @param fg_ls: [] of foreground genes
    @param bg_ls: {} of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """
    scores = {}
    scores_fg = scoreProteinLs_motifLs_notBinary_proteinLs(fg_ls, annotations, motifs)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(permutations):
        bck_scores = scoreBckgnd_motifLs_notBinary_proteinLs( len(fg_ls),  bg_ls, annotations, motifs)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores_fg[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores_fg[motif], len(fg_ls), better_counts[motif] ]
    return scores

def findEnrichedSignatures_motifLs_binary_proteinLs(fg_ls, bg_ls, annotations, motifs, permutations, window):
    """  Check to see if the given signatures(motif1:motif2) are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do X background
         permutations.

    @param fg_ls: [] of foreground genes
    @param bg_ls: {} of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """

    scores = scoreProteinLs_signatureLs_binary_proteinLs(fg_ls, annotations, motifs, window)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(permutations):
        bck_scores = scoreBckgnd_signatureLs_binary_proteinLs(len(fg_ls),  bg_ls, 
                                                              annotations, motifs, window)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores

def findEnrichedAnnotationFreqs_motifLs(fg_ls, bg_ls, annotations, motifs, fasta_file, permutations):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do X background
         permutations.

    @param fg_ls: {} of foreground genes
    @param bg_ls: {} of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """
    fasta = utils_fasta.loadFASTA(fasta_file)
    scores = scoreProteinLsFreqs_motifLs(fg_ls, annotations, motifs, fasta)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(permutations):
        bck_scores = scoreBckgndFreqs_motifLs(len(fg_ls),  bg_ls, 
                                              annotations, motifs, fasta)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores

def findEnrichedAnnotations_motifLs_variants(fg_ls, bg_ls, annotations, motifs, permutations, motif_variants):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do X background
         permutations.

    @param fg_ls: [] of foreground genes
    @param bg_ls: [] of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """

    scores = scoreProteinLs_motifLs_variants(fg_ls, annotations, motifs, motif_variants)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(permutations):
        bck_scores = scoreBckgnd_motifLs_variants( len(fg_ls),  bg_ls, annotations, motifs, motif_variants)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    for motif in motifs.keys():
        scores[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores

def findEnrichedAnnotations(fg_ls, bg_ls, annotations):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do 1000 background
         permutations.

    @param fg_ls: [] of foreground genes
    @param bg_ls: [] of backgound genes
    @param annotations: protein 2 annotation; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """
    motifs = {}
    for protein in annotations.keys():
        for motif in annotations[protein].keys():
            motifs[motif] = True
    scores = scoreProteinLs_motifLs(fg_ls, annotations, motifs)
    better_counts = {}
    for motif in motifs.keys():
        better_counts[motif] = 0
    for i in xrange(1000):
        bck_scores = scoreBckgnd_motifLs( len(fg_ls),  bg_ls, annotations, motifs)
        for motif in motifs.keys():
            if bck_scores[motif] >= scores[motif]:
                better_counts[motif] += 1
    scores_ret = {}
    for motif in motifs.keys():
        scores_ret[motif] = [ scores[motif], len(fg_ls), better_counts[motif] ]
    return scores_ret

def findEnrichedAnnotations_file(fgnd_file, bgnd_file, annotation_file, motifs):
    """  Check to see if the given motifs are enriched in the fg,
         as compared to the bg.  Return motif {}, each entry
         is [#of proteins with motif, |fg|, # of times bg
         permutation is as good as fg].  Do 1000 background
         permutations.

    @param fgnd_file: file of foreground genes
    @param bgnd_file: file of backgound genes
    @param annotation_file: one annotation per line; 
    @param motifs: {} of motifs you are checking for enrichment
    @return: [#of proteins with motif, |fg|, # of times bg is as good as fg] for each motif
    """

    fg_ls = utils_graph.getNodes(fgnd_file)
    bg_ls = utils_graph.getNodes(bgnd_file)
    annotations = utils_motif.protein2annotation_forMotifs(annotation_file,
                                                           motifs)
    return findEnrichedAnnotations(fs_ls, bg_ls, annotations, motifs)   
    
