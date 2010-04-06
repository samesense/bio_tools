#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
import sys, os, utils_graph
from SOAPpy import WSDL

def parseKO_name(description_st):
    sp = description_st.split('\n')
    ko_index = 1
    name_index = 2
    if len(sp[0]) > 0:
        ko_index = 0
        name_index = 1
    ko  = 'ko:' + sp[ko_index].split()[1].lstrip().strip()
    name= sp[name_index].split()[1].lstrip().strip()
    return [ko, name]

def parseKO_EC(description_st):    
    sp = description_st.split('\n')
    ko_index = 1
    if len(sp[0]) > 0:
        ko_index = 0
    ko  = 'ko:' + sp[ko_index].split()[1].lstrip().strip()
    ec = ''
    for line in sp[ko_index+1:]:
        if line != '':
            sp_break = line.split()
            if sp_break[0].strip() == 'DBLINKS':
                ec = sp_break[2].strip().lstrip()
                break
    return [ko, ec] 

def parseGene_name(description_st):
    sp = description_st.split('\n')
    gene_index = 1
    name_index = 2
    if len(sp[0]) > 0:
        gene_index = 0
        name_index = 1
    gene  = sp[gene_index].split()[1].lstrip().strip()
    name= sp[name_index].split()[1].lstrip().strip()
    return [gene, name]

def dict2str(d):
    s = ''
    for item in d.keys():
        s = s + item + ' '
    return s.strip()        

def getKOnames(ko_dict):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    ls = serv.bget( dict2str(ko_dict) ).split('///')
    if ls[0]=='':
        ls = ls[1:]
    if ls[-1] == '\n':
        ls = ls[:-1]   
    map_ls = map(parseKO_name, ls)
    ret_d = dict()
    for item in map_ls:
        ret_d[ item[0] ] = item[1]
    return ret_d

def getKO2EC_EC2KO(ko_dict):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    ls = serv.bget( dict2str(ko_dict) ).split('///')
    if ls[0]=='':
        ls = ls[1:]
    if ls[-1] == '\n':
        ls = ls[:-1]   
    map_ls = map(parseKO_EC, ls)
    ret_d = dict()
    ec2ko = dict()
    for item in map_ls:
        if item[1] != '':
            if not ec2ko.has_key(item[1]): ec2ko[ item[1] ] = dict()
            ec2ko[ item[1] ][ item[0] ] = True
            ret_d[ item[0] ] = item[1]
    return [ret_d, ec2ko]

def getGeneNames(gene_dict):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    #print serv.bget( dict2str(gene_dict) )
    ls = serv.bget( dict2str(gene_dict) ).split('///')
    if ls[0]=='':
        ls = ls[1:]
    if ls[-1] == '\n':
        ls = ls[:-1]   
    map_ls = map(parseGene_name, ls)
    ret_d = dict()
    prefix = gene_dict.keys()[0].split(':')[0]
    for item in map_ls:
        ret_d[ prefix + ':' + item[0] ] = item[1]
    return ret_d

def getGene2KO(path):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    kos = serv.get_kos_by_pathway(path)
    gene2ko = dict()
    for ko in kos:
        if ko=='ko:K07363': print 'have ko in gene2KO'
        genes = serv.get_genes_by_ko(ko)
        #ko = serv.bget(ko).split('NAME')[1].split('DEFINITION')[0].strip()
        for gene in genes:
            gene = gene['entry_id']
            if not gene2ko.has_key(gene):
                gene2ko[gene] = dict()
            gene2ko[gene][ko] = True
            if ko=='ko:K07363':
                print 'genes', gene
    return gene2ko

def junk(path):
   wsdl = 'http://soap.genome.jp/KEGG.wsdl'
   serv = WSDL.Proxy(wsdl)
   elements = serv.get_elements_by_pathway('path:hsa04660')
   print elements
   
def printEdges(species_code):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    for path in pathways:
        path = path['entry_id']
        elements = serv.get_elements_by_pathway(path)
        relations = serv.get_element_relations_by_pathway(path)
        elementid2gene = dict()
        for element in elements:
            if element['type'] == 'gene':
                elementid = element['element_id']
                names = element['names']
                elementid2gene[elementid] = dict()
                for n in names:
                    elementid2gene[elementid][n] = True

        for relation in relations:
            k1 = relation['element_id1']
            k2 = relation['element_id2']
            if elementid2gene.has_key(k1) and elementid2gene.has_key(k2): 
                elements1 = elementid2gene[ k1 ]
                elements2 = elementid2gene[ k2 ]
                for e1 in elements1.keys():
                    for e2 in elements2.keys():
                        for s in relation['subtypes']:
                            print e1.split(':')[1] + '\t' + \
                                  e2.split(':')[1] + '\t' + \
                                  relation['type'] + '\t' + \
                                  s['relation'] + '\t' + path

def printPathwayDesc(species_code):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    for path in pathways:
        print path['entry_id'] + '\t' + path['definition']

def printGenes2Path(species_code):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    for path in pathways:
        genes = serv.get_genes_by_pathway(path['entry_id'])
        for g in genes:
            print g + '\t' + path['entry_id']

def getPathwayBoxes(path):
    edges = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    elements = serv.get_elements_by_pathway(path)
    relations = serv.get_element_relations_by_pathway(path)
    elementid2gene = dict()
    for element in elements:
        if element['type'] == 'gene':            
            elementid = element['element_id']
            names = element['names']
            elementid2gene[elementid] = dict()
            for n in names:
                elementid2gene[elementid][n.split(':')[-1]] = True
    return elementid2gene

def getHitBoxes(path, geneLs):
    element2box = getPathwayBoxes(path)
    hit_boxes = {}
    for e in element2box.keys():
        for gene in element2box[e].keys():
            if geneLs.has_key(gene):
                hit_boxes[e] = True
                break
    return hit_boxes

def getPathwayEdges(path):
    edges = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    elements = serv.get_elements_by_pathway(path)
    relations = serv.get_element_relations_by_pathway(path)
    elementid2gene = dict()
    for element in elements:
        if element['type'] == 'gene':
            elementid = element['element_id']
            names = element['names']
            elementid2gene[elementid] = dict()
            for n in names:
                elementid2gene[elementid][n] = True

    for relation in relations:
        k1 = relation['element_id1']
        k2 = relation['element_id2']
        if elementid2gene.has_key(k1) and elementid2gene.has_key(k2): 
            elements1 = elementid2gene[ k1 ]
            elements2 = elementid2gene[ k2 ]
            for e1 in elements1.keys():
                for e2 in elements2.keys():
                    for s in relation['subtypes']:
                        gene1 = e1.split(':')[1]
                        gene2 = e2.split(':')[1]
                        if not edges.has_key(gene1):
                            edges[gene1] = dict()
                        if not edges.has_key(gene2):
                            edges[gene2] = dict()
                        if not edges[gene1].has_key(gene2):
                            edges[gene1][gene2] = dict()
                        if not edges[gene2].has_key(gene1):
                            edges[gene2][gene1] = dict()
                        edges[gene1][gene2][ relation['type'] + '/' + s['relation'] ] = True
                        edges[gene2][gene1][ relation['type'] + '/' + s['relation'] ] = True
    return edges

def getPathways(species_code):
    """ Ex path:hsa04010 """

    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    pathways = serv.list_pathways(species_code)
    paths = {}
    for path in pathways:
        paths[ path['entry_id'] ] = True
    return paths

# these are stripped of the prefix
def getPathwayGenes(path):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    ls = {}
    for gene in serv.get_genes_by_pathway(path):
        ls[gene.split(':')[1]] = True
    return ls

def getKODictForPathway(path):
    complex2gene = dict()
    gene2complex = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    elements = serv.get_elements_by_pathway(path)
    for element in elements:
        if element['type'] == 'gene':            
            names = element['names']
            elementID = str(element['element_id'])
            complex2gene[elementID] = dict()
            for n in names:
                n=n.split(':')[1]
                complex2gene[elementID][n] = True
                if not gene2complex.has_key(n): gene2complex[n] = dict()
                gene2complex[n][elementID] = True
    return [complex2gene, gene2complex]

def getKOPathDict(path):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    use_KOs_ls = serv.get_kos_by_pathway(path)
    use_KOs = dict()
    for item in use_KOs_ls:
        use_KOs[item] = True
    return use_KOs

def getGraphicLinks(path):
    use_KOs = dict()
    gene2ko = getGene2KO(path)
    edges = dict()
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    use_KOs_ls = serv.get_kos_by_pathway(path)
    for item in use_KOs_ls:
        use_KOs[item] = True
    #print use_KOs.keys()
    elements = serv.get_elements_by_pathway(path)
    relations = serv.get_element_relations_by_pathway(path)
    elementid2gene = dict()
    for element in elements:
        if element['type'] == 'gene':
            elementid = element['element_id']
            names = element['names']
            elementid2gene[elementid] = dict()
            for n in names:
                elementid2gene[elementid][n] = True

    for relation in relations:
        k1 = relation['element_id1']
        k2 = relation['element_id2']
        if elementid2gene.has_key(k1) and elementid2gene.has_key(k2): 
            elements1 = elementid2gene[ k1 ]
            elements2 = elementid2gene[ k2 ]
            for gene1 in elements1.keys():
                ko1 = ''
                if gene2ko.has_key(gene1):
                    ko1s = gene2ko[gene1]
                    for k in ko1s.keys():
                        if use_KOs.has_key(k):
                            if ko1 != '':
                                print 'WARNING, multiple assignment', ko1s
                            ko1 = k
                if ko1 == 'ko:K07363': print 'found'
                for gene2 in elements2.keys():
                    for s in relation['subtypes']:
                        ko2 = ''                        
                        if gene2ko.has_key(gene2):
                            ko2s = gene2ko[gene2]
                            for k in ko2s.keys():
                                if use_KOs.has_key(k):
                                    if ko2 != '':
                                        print 'WARNING, multiple assignment', ko2s
                                    ko2 = k
                        if ko2 == 'ko:K07363': print 'found'
                        if ko1 != '' and ko2 != '':
                            #print ko1, ko2
                            if not edges.has_key(ko1):
                                edges[ko1] = dict()
                            if not edges.has_key(ko2):
                                edges[ko2] = dict()
                            if not edges[ko1].has_key(ko2): edges[ko1][ko2] = dict()
                            if not edges[ko2].has_key(ko1): edges[ko2][ko1] = dict()
                            edges[ko1][ko2][ relation['type'] + '/' + s['relation']  ] = True
                            edges[ko2][ko1][ relation['type'] + '/' + s['relation']  ] = True
                        elif ko2 != '':
                            #print gene1, ko2
                            if not edges.has_key(gene1):
                                edges[gene1] = dict()
                            if not edges.has_key(ko2):
                                edges[ko2] = dict()
                            if not edges[gene1].has_key(ko2): edges[gene1][ko2] = dict()
                            if not edges[ko2].has_key(gene1): edges[ko2][gene1] = dict()
                            edges[gene1][ko2][ relation['type'] + '/' + s['relation']  ] = True
                            edges[ko2][gene1][ relation['type'] + '/' + s['relation']  ] = True
                        elif ko1 != '':
                            #print ko1, gene2
                            if not edges.has_key(ko1):
                                edges[ko1] = dict()
                            if not edges.has_key(gene2):
                                edges[gene2] = dict()
                            if not edges[ko1].has_key(gene2): edges[ko1][gene2] = dict()
                            if not edges[gene2].has_key(ko1): edges[gene2][ko1] = dict()
                            edges[ko1][gene2][ relation['type'] + '/' + s['relation']  ] = True
                            edges[gene2][ko1][ relation['type'] + '/' + s['relation']  ] = True
                        else:
                            #print gene1, gene2
                            if not edges.has_key(gene1):
                                edges[gene1] = dict()
                            if not edges.has_key(gene2):
                                edges[gene2] = dict()
                            if not edges[gene1].has_key(gene2): edges[gene1][gene2] = dict()
                            if not edges[gene2].has_key(gene1): edges[gene2][gene1] = dict()
                            edges[gene1][gene2][ relation['type'] + '/' + s['relation']  ] = True
                            edges[gene2][gene1][ relation['type'] + '/' + s['relation']  ] = True
    return edges

def mkGene2GI(gene_dict):
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    #serv.get_genes_by_organism('sce')
    # this must happen in groups of 100
    # b/c of the limit allowed by the query
    count = 1
    new_d = {}
    for gene in gene_dict.keys():
        new_d[gene]=True
        if count % 100 == 0:            
            print serv.bget( dict2str(new_d) ).strip()
            new_d = {}
        count += 1
    print serv.bget( dict2str(new_d) ).strip()#.split('///')
    #for item in ls:
    #    print item

def prepAndColor(pathway, gene_file, html_color):
    """ Color these genes for this pathway.
        Does not color Cell Communication (path:hsa01430).

    @param pathway: KEGG pathway, format path:hsa04010
    @param gene_file: one gene per line
    @param html_color_1: color, html format
    """
    gene_dict = utils_graph.getNodes(gene_file)
    obj_ls = []
    fore_gnd = []
    back_gnd = []
    for gene in gene_dict.keys():
        obj_ls.append(gene)
        fore_gnd.append('black')
        back_gnd.append(html_color)

    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    url = serv.color_pathway_by_objects(pathway, obj_ls, fore_gnd, back_gnd)
    name = pathway.split(':')[-1] + '_' + html_color
    os.system('wget --output-document=' + name + ' ' + url)
    return name

def prepAndColor_gene_dict(pathway, gene_dict, html_color):
    """ Color these genes for this pathway.
        Does not color Cell Communication (path:hsa01430).

    @param pathway: KEGG pathway, format path:hsa04010
    @param gene_dict: {} of genes
    @param html_color_1: color, html format
    """
    obj_ls = []
    fore_gnd = []
    back_gnd = []
    for gene in gene_dict.keys():
        obj_ls.append(gene)
        fore_gnd.append('black')
        back_gnd.append(html_color)

    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    url = serv.color_pathway_by_objects(pathway, obj_ls, fore_gnd, back_gnd)
    name = pathway.split(':')[-1] + '_' + html_color + '.gif'
    os.system('wget --output-document=' + name + ' ' + url)
    return name

def colorPathwayWith2GeneSets(pathway, 
                              gene_file_1, gene_file_2, 
                              html_color_1, html_color_2,
                              blending_percent,
                              output_dir):
    """ Color these genes, and make the overlap a mix of the colors.
        Does not color Cell Communication (path:hsa01430).

    @param pathway: KEGG pathway, format path:hsa04010
    @param gene_file_1: one gene per line
    @param gene_file_2: one gene per line
    @param html_color_1: color for first set
    @parma html_color_2: color for second
    @param blending_percent: int 0..100 for composite --dissolve input
    @param output_dir: directory to put figure
    """

    name1 = prepAndColor(pathway, gene_file_1, html_color_1)
    name2 = prepAndColor(pathway, gene_file_2, html_color_2)
    name1_alpha = name1 + '_alpha'
    name2_alpha = name2 + '_alpha'
    gimp_prefix = "gimp --no-interface --batch '(python-fu-color-kegg "
    os.system(gimp_prefix
              + "RUN-NONINTERACTIVE \""
              + name2 + "\" \"" + name1 + "\" \""
              + name1_alpha + "\" (list 191 255 191))' --batch '(gimp-quit 1)'")
    os.system(gimp_prefix
              + "RUN-NONINTERACTIVE \""
              + name1 + "\" \"" + name2 + "\" \""
              + name2_alpha + "\" (list 191 255 191))' --batch '(gimp-quit 1)'")
    os.system("composite "
              + name1_alpha + ' '
              + name2_alpha + " -dissolve \""
              + str(blending_percent) + "%\" "
              + output_dir + pathway.split(':')[1] + '.png')
    os.system('rm ' + name1)
    os.system('rm ' + name2)
    os.system('rm ' + name1_alpha)
    os.system('rm ' + name2_alpha)
    
def getGenesBySpecies(species):
    """ Return {} of all KEGG genes for this speices. 
        Genes do not have the species prefix.  """
    
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    gene_ls = serv.get_genes_by_organism(species, 1, 50000)
    gene_dict = {}
    for gene in gene_ls:
        gene_dict[gene.split(':')[1]] = True
    return gene_dict

def getPathwayGenesBySpecies(species):
    """ Return {} of all KEGG pathway genes for this species. 
        Genes do not have the species prefix.  """
    kegg_genes = {}
    pathways = getPathways('hsa')
    for p in pathways.keys():
        pathway_genes = getPathwayGenes(p)
        for g in pathway_genes.keys():
            kegg_genes[g] = True
    return kegg_genes    

