#--------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
import MySQLdb, sys, string
import utils_graph

gocats = ['Process', 'Component', 'Function']

gocatAlias = dict()
gocatAlias['cellular_component'] = 'Component'
gocatAlias['molecular_function'] = 'Function'
gocatAlias['biological_process'] = 'Process'
gocatAlias['universal'] = 'Universal'

def findWorkingTerm(termLs, d):
    for t in termLs:
        if d.has_key(t):
            return t
    print 'your fucked: utils_go.py'
    print 'cannot find your GO term'
    sys.exit(0)

# connect to the GO database
def getGOConnection_Cursor():
    connection = MySQLdb.connect(host="mysql.ebi.ac.uk", \
                                 user='go_select', \
                                 db='go_latest', \
                                 port=4085, \
                                 passwd='amigo')
    cursor = connection.cursor()
    return [connection, cursor]

# for a given species,
#  for each GO term for that species
#   return a dict of the term's proteins and
#          an empty dict for the term's children
#          an empty dict for the term's parents
def getTermProperties(tax_id):
    ls = dict()
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT term_id, gene_product_id from association a inner join gene_product b on a.gene_product_id=b.id where b.species_id=(select id from species where ncbi_taxa_id=" + tax_id + ");")
        data = cursor.fetchall()
        for item in data:
            if not ls.has_key(item[0]):
                ls[ item[0] ] = [dict(), dict(), dict()] #proteins, children
            ls[ item[0] ][0][ item[1] ] = True
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()

    cursor.close()
    connection.close()
    return ls

# get the parent/child relations
# for all terms
# return 2 lists
#  child2Parent and parent2Child
def getTermRelations():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        # in term table
        #term2 (pair[1]) is child
        #term1 (pair[0]) is parent
        cursor.execute("select term1_id, term2_id, relationship_type_id from term2term;")
        data = cursor.fetchall()
        child2Parent = dict()
        parent2Child = dict()
        for pair in data:
            if not child2Parent.has_key(pair[1]):
                child2Parent[ pair[1] ] = dict()
            child2Parent[pair[1]][ pair[0] ] = True

            if not parent2Child.has_key(pair[0]):
                parent2Child[ pair[0] ] = dict()
            parent2Child[ pair[0] ][ pair[1] ] = True
            
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()

    cursor.close()
    connection.close()
    return [child2Parent, parent2Child]

# return a dictionary and a maximum distance
#  from the root (dis from root to root is 0)
# indexed by distances from the root terms
# values are the GO terms found at a
#  specific distance
# a term can be at multiptle distances
#  from the root b/c there can be multiple paths
#  to a term
def getDistanceFromRoot():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT term.id, p.distance FROM term INNER JOIN graph_path AS p ON (p.term2_id=term.id) INNER JOIN term AS root ON (p.term1_id=root.id) WHERE root.is_root=1;")
        data = cursor.fetchall()
        disDict = dict()
        maxDist = 0
        for item in data:
            if not disDict.has_key(str(item[1])):
                if int(item[1]) > maxDist:
                    maxDist = int(item[1])
                disDict[str(item[1])] = dict()
            disDict[str(item[1])][item[0]] = True
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
        
    cursor.close()
    connection.close()
    return [disDict, maxDist]

def getMaxDistanceFromRoot():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT term.id, p.distance FROM term INNER JOIN graph_path AS p ON (p.term2_id=term.id) INNER JOIN term AS root ON (p.term1_id=root.id) WHERE root.is_root=1;")
        data = cursor.fetchall()
        distDict = dict()
        for item in data:
            distDict[ item[0] ] = item[1]
            
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
        
    cursor.close()
    connection.close()
    return distDict

# assigns children to all terms in termProperties
#  and assigns proteins and children to terms
#  that can be reached using child2Parent and a term
#  in termProperties
# GO doesn't give proteins for all its terms
#  so I must find this information using the
#  hierarchy
def fillInTermProperties(termProperties, child2Parent):
    distInfo = getDistanceFromRoot()
    distDict = distInfo[0]
    maxDist = distInfo[1]
    for dis in xrange(maxDist+1):
        for term in distDict[str(maxDist-dis)].keys():
            if termProperties.has_key(term):
                proteins = termProperties[term][0].keys()
                if child2Parent.has_key(term):
                    for p in child2Parent[term].keys():
                        if not termProperties.has_key(p):
                            termProperties[p] = [dict(), dict(), dict()]
                        for k in proteins:
                            termProperties[p][0][k] = True
                        termProperties[p][1][term] = True
                        for c in termProperties[term][1].keys():
                            termProperties[p][1][c] = True
            else:
                termProperties[term] = [dict(),dict(), dict()]
                if child2Parent.has_key(term):
                    for p in child2Parent[term].keys():
                        if not termProperties.has_key(p):
                            termProperties[p] = [dict(), dict(), dict()]
                        termProperties[p][1][term] = True    

# assigns children and parents to all terms in termProperties
#  and assigns proteins and children to terms
#  that can be reached using child2Parent or parent2Child
#  and a term in termProperties
# GO doesn't give proteins for all its terms
#  so I must find this information using the
#  hierarchy
def fillInTermPropertiesWParents(termProperties, child2Parent, parent2Child):
    fillInTermProperties(termProperties, child2Parent)
    for term in termProperties.keys():
        children = termProperties[term][1]
        for c in children.keys():
            termProperties[c][2][term] = True
                
# for each GO term id
# get its GOterm and cocat
def getTermIDtoGOterm():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT id, acc, term_type, name FROM term WHERE term_type='biological_process' or term_type='cellular_component' or term_type='molecular_function' or term_type='universal';")
        data = cursor.fetchall()        
        termInfo = dict()
        for item in data:
            termInfo[ item[0] ] = [item[1], gocatAlias[item[2]], item[3]]
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
    cursor.close()
    connection.close()
    return termInfo

def term2term_id():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT id, acc, term_type, name FROM term WHERE term_type='biological_process' or term_type='cellular_component' or term_type='molecular_function' or term_type='universal';")
        data = cursor.fetchall()        
        termInfo = dict()
        for item in data:
            if not termInfo.has_key(item[1]): termInfo[item[1]]={}
            termInfo[ item[1] ][ item[0] ] = True
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
    cursor.close()
    connection.close()
    return termInfo

def term_id2term():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT id, acc, term_type, name FROM term WHERE term_type='biological_process' or term_type='cellular_component' or term_type='molecular_function' or term_type='universal';")
        data = cursor.fetchall()        
        termInfo = dict()
        for item in data:
            if not termInfo.has_key(item[0]): termInfo[item[0]]={}
            termInfo[ item[0] ][item[1]] = True
    except MySQLdb.OperationalError, message:
        print message[0]
        print message[1]
        sys.exit()
    cursor.close()
    connection.close()
    return termInfo

# given 2 dicts of terms,
#  find the intersection,
#  split the intersecion into 3 gocats and print
def mkGOStandardFile(terms1, terms2, afile):
    ls = dict()
    for y in terms1.keys():
        if terms2.has_key(y):
            ls[y]=True
    termInfo = getTermIDtoGOterm()
    files = dict()
    for gocat in gocats:
        files[gocat] = open(afile+'.'+gocat, 'w')
    for item in ls.keys():
        itemName = termInfo[item][0]
        itemCat = termInfo[item][1]
        itemDesc = termInfo[item][2]
        files[itemCat].write(itemName + '\t' + itemDesc + '\n')
    for gocat in gocats:
        files[gocat].close()

def getTerm_proteins_children_parents(tax_id):
    [child2Parent, parent2Child] = getTermRelations()
    termProperties = getTermProperties(tax_id)
    fillInTermPropertiesWParents(termProperties, child2Parent, parent2Child)
    return termProperties

def getTermDesc(term):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT name from term WHERE id='" + \
                       term + "';")
        desc = cursor.fetchall()[0][0]
        return desc
    except:
        print 'error in getTermDesc: term =', term
        sys.exit(0)

def getTermSynonym(term):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute("SELECT acc FROM term INNER JOIN term_synonym ON (term.id=term_synonym.term_id) WHERE acc_synonym='" + \
                       term + "';")
        data = cursor.fetchall()
        ls = []
        for item in data:
            ls.append(item[0])
        return ls
    except:
        print 'error in getTermSynonym: term =', term
        sys.exit(0)
        

def getGeneSymbolsPerSpecies(species_commonName):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute('SELECT symbol, id from gene_product where species_id=(select id from species where common_name=\"' + species_commonName + '\");')
        data = cursor.fetchall()
        d = dict()
        for list in data:
            d[ list[0] ] = list[1]
        return d
    except:
        print 'error in getAnnotatedGeneProductsBySpecies'
        sys.exit(0)

def getGeneSynonyms():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute('SELECT gene_product_id, product_synonym from gene_product_synonym;')
        data = cursor.fetchall()
        d = dict()
        for list in data:
            if not d.has_key( list[0] ): d[ list[0] ] = dict()
            d[ list[0] ][ list[1] ] = True
        return d
    except:
        print 'error in getGeneSynonyms'
        sys.exit(0)

def getTermsPerGeneProduct():
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute('SELECT term_id, gene_product_id from association;')
        data = cursor.fetchall()
        d = dict()
        for list in data:
            if not d.has_key( list[1] ): d[ list[1] ] = dict()
            d[ list[1] ][ list[0] ] = True
        return d
    except:
        print 'error in getTermsPerGeneProduct'
        sys.exit(0)

def getTermsPerGeneProductPerSpecies(species_commonName):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute('SELECT term_id, gene_product_id from association a inner join gene_product b on a.gene_product_id=b.id where b.species_id=(select id from species where common_name=\"' + species_commonName + '\");')
        data = cursor.fetchall()
        d = dict()
        for list in data:
            if not d.has_key( list[1] ): d[ list[1] ] = dict()
            d[ list[1] ][ list[0] ] = True
        return d
    except:
        print 'error in getTermsPerGeneProductPerProtein'
        sys.exit(0)

def getAncestors(goTerm):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute('SELECT ancestor.*, graph_path.distance, graph_path.term1_id AS ancestor_id FROM term child, graph_path, term ancestor WHERE child.id=graph_path.term2_id AND ancestor.id=graph_path.term1_id AND child.id=' + goTerm + ';')
        data = cursor.fetchall()
        d = dict()
        for list in data:
            d[ list[0] ] = True
        return d
    except:
        print 'getAncestors'
        sys.exit(0)

def getDescendants(goTerm):
    try:
        [connection, cursor] = getGOConnection_Cursor()
        cursor.execute('SELECT DISTINCT descendant.id, descendant.acc, descendant.name, descendant.term_type FROM term INNER JOIN graph_path ON (term.id=graph_path.term1_id) INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id) WHERE term.name='+ goTerm + ' AND distance <> 0 ;')
        data = cursor.fetchall()
        d = dict()
        for list in data:
            d[ list[0] ] = True
        return d
    except:
        print 'getDescendants'
        sys.exit(0)

def expandTerm(goTerm):
    ancestors = getAncestors(goTerm)
    descendants = getDescendants(goTerm)
    d = dict()
    for a in ancestors.keys():
        d[a] = True
    for c in descendants.keys():
        d[c] = True
    return d

# GO_terms will be converted to GO_ids
def parseGene2GO(gene2go_file, genes):
    GO_ID2database_id = term2term_id()
    gene2go = {}
    f = open(gene2go_file)
    f.readline()
    for line in f.xreadlines():
        [tax_id, GeneID, GO_ID,
         Evidence, Qualifier,
         GO_term, PubMed, Category] = map(string.strip,
                                          line.split('\t'))
        if genes.has_key(GeneID):
            if not gene2go.has_key(GeneID): gene2go[GeneID] = {}
            if not gene2go[GeneID].has_key(Category):
                gene2go[GeneID][Category] = {}
            if GO_ID2database_id.has_key(GO_ID):
                gene2go[GeneID][Category][ GO_ID2database_id[GO_ID].keys()[0] ] = True
    f.close()
    return gene2go

def getTaxID(geneLs, gene2go_file):
    taxa = 'missing'
    f = open(gene2go_file)
    f.readline()
    for line in f.xreadlines():
        [tax_id, GeneID, GO_ID,
         Evidence, Qualifier, GO_term,
         PubMed, Category] = map(string.strip,
                                 line.split('\t'))
        if geneLs.has_key(GeneID):
            taxa = tax_id
            break
    f.close()
    return taxa

def threshold(termProperties, thresh):
    for k in termProperties.keys():
        if len(termProperties[k][0].keys()) < thresh:
            del termProperties[k]
    meetsThreshold = dict()
    for k in termProperties.keys():
        addIt = True
        children = termProperties[k][1]
        for c in children.keys():
            if termProperties.has_key(c):
                addIt = False
        if addIt: meetsThreshold[k] = True
    return meetsThreshold

def getInformativeTerms(tax_id, thresh):
    [child2Parent, parent2Child] = getTermRelations()
    termProperties = getTermProperties(tax_id)
    fillInTermProperties(termProperties, child2Parent)
    return threshold(termProperties, thresh)

def getGOAnnotations_old(taxID, genes2go, child2Parents, ignoreISS):
    workingTermDict = dict()
    for gene in genes2go.keys():
        for category in genes2go[gene].keys():
            toAdd = dict()
            for term in genes2go[gene][category].keys():
                if not child2Parents.has_key(term):
                    if not workingTermDict.has_key(term):
                        termLs = getTermSynonym(term)
                        workingTermDict[term] = termLs
                    else:
                        termLs = workingTermDict[term]
                    term = findWorkingTerm(termLs, child2Parents)
                parents = child2Parents[term]
                for p in parents:
                    if not unknownGOterms.has_key(p):
                        toAdd[p] = True
            for a in toAdd.keys():
                genes2go[gene][category][a] = True

def getGOAnnotations(taxID, genes2go, child2Parents):
    missing = 0
    workingTermDict = dict()
    for gene in genes2go.keys():
        for category in genes2go[gene].keys():
            for term in genes2go[gene][category].keys():
                for a in child2Parents[term].keys():
                    genes2go[gene][category][a] = True

def getChild2Parents(tax_id, termID2termName):
    term2Parents = dict()
    termProperties = getTerm_proteins_children_parents(tax_id)
    for term in termProperties.keys():
        #termName = termID2termName[term][0]
        term2Parents[term] = dict()
        for pID in termProperties[term][2]:
            term2Parents[term][ pID ] = True
    return term2Parents

# call me
# just use a file w/ a list of genes
# and the gene2go file from ncbi
#
# you get each gene and all it's GO terms, their distances from the root, and their cats
def annotate(entrez_gene_ls_file, gene2go_file):
    genes = utils_graph.getNodes(entrez_gene_ls_file)
    tax_id = getTaxID(genes, gene2go_file)
    informative_terms = getInformativeTerms(tax_id, 0)
    gene2go = parseGene2GO(gene2go_file, genes)   
    termId2termName = getTermIDtoGOterm()
    child2Parents = getChild2Parents(tax_id, termId2termName)
    getGOAnnotations(tax_id, gene2go, child2Parents)
    termID2term = term_id2term()
    distances = getMaxDistanceFromRoot()
    for protein in gene2go.keys():
        for category in gene2go[protein].keys():
            for goTerm in gene2go[protein][category].keys():
                print protein + '\t' + termID2term[goTerm].keys()[0] + '\t' + str(distances[goTerm]) + '\t' + category

def getEntrezGenesForTerm(gene2go_file, term, tax_id):
    #immune system process GO:0002376
    #9606
    proteins_with_term = {}
    genes = {}
    gene2go_open = open(gene2go_file)
    for line in gene2go_open:
        genes[ line.split()[1].strip() ] = True
    gene2go_open.close()
    tax_id = getTaxID(genes, gene2go_file)
    informative_terms = getInformativeTerms(tax_id, 0)
    gene2go = parseGene2GO(gene2go_file, genes)   
    termId2termName = getTermIDtoGOterm()
    child2Parents = getChild2Parents(tax_id, termId2termName)
    getGOAnnotations(tax_id, gene2go, child2Parents)
    termID2term = term_id2term()
    distances = getMaxDistanceFromRoot()
    for protein in gene2go.keys():
        for category in gene2go[protein].keys():
            for goTerm in gene2go[protein][category].keys():
                term_ids = termID2term[goTerm]
                if term_ids.has_key(term):
                    proteins_with_term[protein] = True
    return proteins_with_term

def loadAnnotations(afile, levels, go_cats):
    gene2go = {}
    f = open(afile)
    for line in f.xreadlines():
        [gene, go, level, cat] = map(string.strip, line.split('\t'))
        if levels.has_key(level) and go_cats.has_key(cat):
            if not gene2go.has_key(gene):
                gene2go[gene] = {}
            gene2go[gene][go] = True
    f.close()
    return gene2go
    
#def getEntrezGenesForTerm(term_name, tax_id):
    #immune system process GO:0002376
    #9606
#    d = {}
#    try:
#        [connection, cursor] = getGOConnection_Cursor()
#        cursor.execute("SELECT term.name AS superterm_name, gene_product.symbol AS gp_symbol, gene_product.id AS gp_id, species.* FROM term INNER JOIN graph_path ON (term.id=graph_path.term1_id) INNER JOIN association ON (graph_path.term2_id=association.term_id)INNER JOIN gene_product ON (association.gene_product_id=gene_product.id) INNER JOIN species ON (gene_product.species_id=species.id) INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id) WHERE term.name='" + term_name + "' AND species.id=(select id from species where ncbi_taxa_id=" + tax_id + ");")
#        data = cursor.fetchall()
#        d = dict()
#        for list in data:
#            print list
            #d[ list[0] ] = True
#        return d
#    except:
#        print 'getEntrezGenesForTerm'
#        sys.exit(0)
#    return d
    
