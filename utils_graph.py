#----------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#
#----------------------------------------
"""
Functions here help with network loading,
printing, and searching.
"""
import sys
from collections import defaultdict

def getNodes(afile):
    """ Make a {} from a file.  Each line
        of the file is one item in the {}.
    
    @param afile: file w/ one item per line
    @return: {} of nodes
    """

    nodes = dict()
    f_nodes = open(afile)
    for line in f_nodes:
        nodes[ line.strip() ] = True
    f_nodes.close()
    return nodes

def dumpNodes(afile, node_dict):
    """ Output {} to file.  One item
    per line.

    @param afile: write to this file
    @param node_dict: {} with items for writing
    """

    f_out = open(afile, 'w')
    for node in node_dict.keys():
        f_out.write(node + '\n')
    f_out.close()

def getEdges(afile):
    """ Make edge {} from file.
    {} format is d[n1][n2] = True

    @param afile: tab-delimited file; one edge per line; edges are bi-directional
    @return: {} of edges
    """

    edges = dict()
    f_edges = open(afile)
    for line in f_edges:
        sp = [x.strip() for x in line.split('\t')]
        [node1, node2] = sp[0:2]
        if node1 != node2:
            if not edges.has_key(node1): 
                edges[node1] = dict()
            if not edges.has_key(node2): 
                edges[node2] = dict()
            edges[node1][node2] = True
            edges[node2][node1] = True
    f_edges.close()
    return edges

def getEdges_pubmed(afile):
    """ Make edge {} from file.
    {} format is d[n1][n2] = True

    @param afile: tab-delimited file; one edge per line; edges are bi-directional
    @return: {} of edges
    """

    edges = dict()
    f_edges = open(afile)
    for line in f_edges:
        [node1, node2, pubmed] = [x.strip() for x in line.split('\t')]
        if node1 != node2:
            if not edges.has_key(node1): 
                edges[node1] = dict()
            if not edges.has_key(node2): 
                edges[node2] = dict()
            pubmed_ls = []
            for paper in pubmed.split(','):
                paper_id = paper.strip()
                if paper_id != '':
                    pubmed_ls.append(paper_id)
            edges[node1][node2] = pubmed_ls
            edges[node2][node1] = pubmed_ls
    f_edges.close()
    return edges

def getEdges_forProteins(afile, use_proteins):
    """ Make edge {} from file.
    {} format is d[n1][n2] = True

    @param afile: tab-delimited file; one edge per line; edges are bi-directional
    @return: {} of edges
    """

    edges = dict()
    f_edges = open(afile)
    for line in f_edges:
        [node1, node2] = [x.strip() for x in line.split('\t')][0:2]
        if use_proteins.has_key(node1) and use_proteins.has_key(node2):
            if node1 != node2:
                if not edges.has_key(node1): 
                    edges[node1] = dict()
                if not edges.has_key(node2): 
                    edges[node2] = dict()
                edges[node1][node2] = True
                edges[node2][node1] = True
    f_edges.close()
    return edges

def dumpEdges(afile, edge_dict):
    """ Write edges to file.
    
    @param afile: output is tab-delimited; bi-directional edges are not duplicated
    @param edge_dict: {} format is d[n1][n2] = True
    """
    
    f_out = open(afile, 'w')
    seen = {}
    for gene1 in edge_dict.keys():
        for gene2 in edge_dict[gene1].keys():
            key1 = gene1+':'+gene2
            key2 = gene2+':'+gene1
            if not seen.has_key(key1) and not seen.has_key(key2):
                seen[key1] = True
                seen[key2] = True
                f_out.write(gene1 + '\t' + gene2 + '\t' + edge_dict[gene1][gene2] + '\n')
    f_out.close()

def mkNodesFromEdges(afile):
    """ Given an edge file, make a {} of nodes.

    @param afile: tab-delimted file of edges; one per line; edges are bi-directional
    @return: {} of nodes
    """

    nodes = {}
    f_edges = open(afile)
    for line in f_edges:
        splitup = line.strip().split('\t')
        nodes[splitup[0]] = True
        nodes[splitup[1]] = True
    f_edges.close()
    return nodes

def getConnected(ls_file, net_file):
    """ Find all the genes in the network
    that are one away from these genes.

    @param ls_file: file of genes to find connections to; one gene per line
    @param net_file: tab-delimited file for the network; one edge per line; edges are bi-directional
    @return: {} of genes
    """

    seed_ls = getNodes(ls_file)
    net = getEdges(net_file)
    second_level_ls = {}
    for gene in seed_ls.keys():
        if net.has_key(gene):
            for gene2 in net[gene].keys():
                if gene != gene2:
                    second_level_ls[gene2] = True
    return second_level_ls

def intersectLists(list_of_lists):
    """ Given a list of {}, return the intersection
    as a {}.

    @param list_of_lists: [] of {}
    @return: {} of genes
    """

    if len(list_of_lists) == 0: return {}
    same = {}
    for item in list_of_lists[0].keys():
        count = 0
        for alist in list_of_lists[1:]:
            if alist.has_key(item): 
                count += 1
        if count == len(list_of_lists)-1:
            same[item] = True
    return same

def intersectFiles(list_of_files):
    """ Given a list of files, return the intersection
    as a {}.

    @param list_of_lists: [] of files
    @return: {} of genes
    """

    list_of_lists = []
    for afile in list_of_files:
        list_of_lists.append(getNodes(afile))
    return intersectLists(list_of_lists)
    
def unionLists(list_of_lists):
    """ Given a list of {}, return the union
    as a {}.

    @param list_of_lists: list of {}
    @return: {} of union
    """

    union = {}
    for alist in list_of_lists:
        for k in alist.keys():
            union[k] = True
    return union

def degree(gene, network):
    """ Find the degree of this node
    in this network.
    
    @param gene: node in the network
    @param network: edge {} format d[n1][n2] = True
    @return: deg of gene as int
    """

    return len(network[gene].keys())

def degree_ignoreSim(gene, network, sims, cutoff):
    """ Find the degree of this node
    in this network.
    
    @param gene: node in the network
    @param network: edge {} format d[n1][n2] = True
    @return: deg of gene as int
    """
    seen_genes = {}
    for neigh in network[gene]:
        found_match = False
        for seen in seen_genes:
            if neigh in sims[seen]:
                if sims[seen][neigh] < cutoff:
                    found_match = True
        if not found_match:
            seen_genes[neigh] = True
    return len(seen_genes.keys())

def avgConnectivity(gene_ls, network):
    """ Given a set of genes, return
    their average connectivity in the network.

    @param geneLs: {} of genes
    @param network: edge {} format d[n1][n2] = True
    @return: float for avg connectivity
    """

    [connectivity, total] = [0, 0]
    for gene in gene_ls.keys():
        if network.has_key(gene):
            connectivity += degree(gene, network)
            total += 1
    return float(connectivity) / float(total)

def getDegreeDistributionLs(gene_dict, network_edges):
    """ For these genes, return the degree distribution
        according to this network.

    @param gene_dict: {} w/ gene names
    @param network_edges {} of edges d[n1][n2] = True
    @return {} of degrees d[degree] = count
    """

    degreeDict = []
    for gene in gene_dict.keys():
        if network_edges.has_key(gene):
            degreeDict.append(degree(gene, network_edges))
            #if not degreeDict.has_key(deg):
            #    degreeDict[str(deg)] = 0
            #degreeDict[str(deg)] += 1
    return degreeDict

def getDegreeDistributionDict(gene_dict, network_edges):
    """ For these genes, return the degree distribution
        according to this network.

    @param gene_dict: {} w/ gene names
    @param network_edges {} of edges d[n1][n2] = True
    @return {} of degrees d[degree] = count
    """

    degreeDict = {}
    for gene in gene_dict.keys():
        if network_edges.has_key(gene):
            deg = degree(gene, network_edges)
            if not degreeDict.has_key(str(deg)):
                degreeDict[str(deg)] = 0
            degreeDict[str(deg)] += 1
    return degreeDict

def getProteinsWithXneighbors(network_edges, X):
    """ Return {} of proteins with X neighbors. """
    
    qualify = {}
    for protein in network_edges.keys():
        if degree(protein, network_edges) == X:
            qualify[protein] = True
    return qualify

def getProteinsWithAtLeastXneighbors(network_edges, X):
    """ Return {} of proteins with at least X neighbors. """
    
    qualify = {}
    for protein in network_edges.keys():
        if degree(protein, network_edges) >= X:
            qualify[protein] = True
    return qualify

def getTopXconnectedProteins(network_edges, X, bool_expand):
    """ Return {} of proteins in top X% connectivity.

        Rank by connectivity and take top X%. 

        If bool_expand == True, the hub list will be
        expanded to include all hubs with the degree
        cutoff.  If False, the degree cutoff will
        be the lowest cutoff that keeps the node
        count below X%."""
    
    degree_protein = []
    for protein in network_edges.keys():
        deg = degree(protein, network_edges)
        degree_protein.append([deg, protein])
    
    degree_protein.sort(reverse=True)
    total_proteins = len(network_edges.keys())
    top_x = int(float(X) * float(total_proteins) / float(100))
    deg = 0
    proteins = {}
    counter = 0
    while counter < top_x:
        [d, protein] = degree_protein[counter]
        proteins[protein] = True
        deg = d
        counter += 1
    if bool_expand:
        #extend list
        if counter < len(degree_protein):
            current_d = degree_protein[counter][0]
            while current_d == deg and counter < len(degree_protein):
                [d, protein] = degree_protein[counter]
                proteins[protein] = True
                current_d = d
                counter += 1
    else:
        # shrink list
        if counter < len(degree_protein):
            if degree_protein[counter][0] == deg:
                counter -= 1
            while degree_protein[counter][0] == deg:
                del proteins[degree_protein[counter][1]]
                counter -= 1

    return proteins

def getTopXconnectedProteins_deg(network_edges, X):
    """ Get the proteins corresponding to the top X degrees """
    
    degree2protein = defaultdict(dict)
    for protein in network_edges:
        deg = degree(protein, network_edges)
        degree2protein[deg][protein] = True
    degrees = degree2protein.keys()
    degrees.sort(reverse=True)
    total_degrees = len(degrees)
    top_x = int(float(X) * float(total_degrees) / float(100))

    proteins = {}
    counter = 0
    while counter < top_x:
        for protein in degree2protein[ degrees[counter] ]:
            proteins[protein] = True
        counter += 1

    return proteins

def split_into_topConnect_other(network_edges, X, expand):
    """ Return [top connected, all others] """

    hubs = getTopXconnectedProteins(network_edges, 
                                    X, expand)
    nonHubs = {}
    for gene in network_edges:
        if not gene in hubs:
            nonHubs[gene] = True
    return [hubs, nonHubs]    

def getTopXconnectedProteins_ignoreSim(network_edges, sims, X, simCut):
    """ Return {} of proteins in top X% connectivity.

        Rank by connectivity and take top X%. """
    
    degree_protein = []
    for protein in network_edges.keys():
        deg = degree_ignoreSim(protein, network_edges, sims, simCut)
        degree_protein.append([deg, protein])
    
    degree_protein.sort(reverse=True)
    total_proteins = len(network_edges.keys())
    top_x = int(float(X) * float(total_proteins) / float(100))
    deg = 0
    proteins = {}
    counter = 0
    while counter < top_x:
        [d, protein] = degree_protein[counter]
        proteins[protein] = True
        deg = d
        counter += 1
    #for [d, protein] in degree_protein:
    #    print d, protein
    current_d = degree_protein[counter][0]
#    print '//', current_d
    while current_d == deg:
        [d, protein] = degree_protein[counter]
        proteins[protein] = True
        current_d = d
        counter += 1
    #print '//', current_d
    return proteins

def getNeighbors(genes, edges):
   neighbors = {}
   for gene in genes.keys():
       if edges.has_key(gene):
           for neighbor in edges[gene].keys():
               neighbors[neighbor] = True
   return neighbors

def countEdges(network):
    edges = {}
    for p1 in network.keys():
        for p2 in network[p1].keys():
            k1 = p1+':'+p2
            k2 = p2+':'+p1
            if not edges.has_key(k1) and not edges.has_key(k2):
                edges[k1] = True
    return len(edges.keys())
