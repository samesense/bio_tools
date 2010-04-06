#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions for automating calls to DAVID.
"""
import PyMozilla, utils_kegg

def fullMonty(astring):
    """ Remove hsa, commas, and whitespace from
        this term.  This only works for
        human proteins.

    @param astring: string to be formatted
    @return: formatted string
    """

    if astring == '':
        return astring
    if astring[-1] == 'O':
        astring = astring[0:-4]
    elif astring[-8:-5] == 'hsa':
        astring = astring[0:-10]
    return astring.strip().lstrip().strip(',')

def cat2term2protein_sig(iterable):
    """ Take the output from DAVID chart
        and organize by category, term, protein.
    
    @param iterable: open file or list that has DAVID output
    @return: [Category][Term]['count' | 'pvalue' | 'genes'={} ]
    """

    cat2term2protein = {}
    for line in iterable:
        splitup = [x.strip() for x in line.split('\t')]
        category = splitup[0]
        term = splitup[1]
        count = splitup[2]
        percent = splitup[3].strip('%')
        pvalue = splitup[-3]
        genes = splitup[5]         
        if not cat2term2protein.has_key(category):
            cat2term2protein[category] = {}
        genes = genes.strip(',').split(',')
        gene_d = {}
        for gene in genes: 
            gene_d[gene.strip().lstrip()] = True        
        cat2term2protein[category][term] = {}
        cat2term2protein[category][term]['count'] = int(count)
        cat2term2protein[category][term]['pathway_total'] = str(percent) + '/' + str(int( float(count) * float(100) / float(percent) ))
        cat2term2protein[category][term]['pvalue'] = float(pvalue)
        cat2term2protein[category][term]['genes'] = gene_d
    return cat2term2protein

def cat2term2protein_sigFile(afile):
    """ Parse a DAVID chart file into
        category, term, proteins.
    @param afile: DAVID annotation output
    @return: [Category][Term]['count' | 'pvalue' | 'genes'={} ]
    """

    f_david = open(afile)
    f_david.readline()
    cat2term2protein = cat2term2protein_sig(f_david)
    f_david.close()
    return cat2term2protein

def cat2term2protein_all(iterable):
    """ Take the output from DAVID annotation
        and organize by category, term, protein.
    
    @param iterable: open file or list that has DAVID output
    @return: [Category][Term]['count' | 'genes'={} ]
    """

    cat2term2protein = {}
    categories = [x.strip() for x in iterable[0].split('\t')[3:7]]

    for cat in categories:
        cat2term2protein[cat] = {}
    
    for line in iterable[1:]:
        splitup = [x.strip() for x in line.split('\t')]
        for i in xrange(len(categories)):
            cat = categories[i]
            terms = [x.strip() for x in splitup[i+3].strip(',').split(':')]
            for index in xrange(len(terms)-1):
                if terms[index] != '':
                    term = terms[index].split()[-1] + ':' \
                           + fullMonty(terms[index+1])
                    if not cat2term2protein[cat].has_key(term):
                        cat2term2protein[cat][term] = {}
                        #cat2term2protein[cat][term]['count'] = 0
                        cat2term2protein[cat][term]['genes'] = {}                  
                    cat2term2protein[cat][term]['genes'][splitup[0]] = True
    for cat in cat2term2protein.keys():
        for term in cat2term2protein[cat].keys():
            cat2term2protein[cat][term]['count'] = len(cat2term2protein[cat][term]['genes'].keys())      
    return cat2term2protein

def cat2term2protein_allFile(afile):
    """ Parse a DAVID annotation file into
        category, term, proteins.
    @param afile: DAVID annotation output
    @return: [Category][Term]['count' | 'genes'={} ]
    """

    f_david = open(afile)
    cat2term2protein = cat2term2protein_all(f_david.readlines())
    f_david.close()
    return cat2term2protein

def getSigKEGGgenes_sigFile(afile, cutoff):
    """ Return {} of genes in significant KEGG
        pathways with p-values less than the
        cutoff, as described in the DAVID file.

    @param afile: DAVID chart file
    @param cutoff: p-value cut; less than #
    @return: {} of genes in sig KEGG pathways
    """

    genes = {}
    cat2term2protein = cat2term2protein_sigFile(afile)
    for path in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][path]['pvalue'] < cutoff:
            for gene in cat2term2protein['KEGG_PATHWAY'][path]['genes'].keys():
                genes[gene] = True
    return genes

def getAllTerms_web(gene_file):
    """ Given a file with one gene per line, return
        the parsed DAVID annotation for the genes in
        the form [Category][Term]['count' | 'genes'={} ].
    
    @param gene_file: one gene per line in this file
    @return string representation of downloaded webpage
    """
    david_prefix = 'http://david.abcc.ncifcrf.gov/'
    moz_emu = PyMozilla.MozillaEmulator(cacher=None, trycount=0)
    files = [['fileBrowser', 'testData', '']]
    fields = []
    moz_emu.download(david_prefix + 'tools.jsp')
    fdata = file(gene_file).read()
    fields = []
    fields.append(['idType', 'ENTREZ_GENE_ID'])
    fields.append(['Mode', 'file'])
    fields.append(['uploadType', 'list'])

    files = []
    files.append(['fileBrowser', 'testData', fdata])
    
    moz_emu.post_multipart(david_prefix + 'tools.jsp',
                           fields, files)
    david_cats = 'GOTERM_MF_5,GOTERM_MF_2,KEGG_PATHWAY'
    moz_emu.download(david_prefix + 'template.css')
    moz_emu.download(david_prefix + 'scripts/sidebar.js')
    moz_emu.download(david_prefix + 'scripts/summary.js')
    user_download_page = moz_emu.download(david_prefix
                                          + 'annotationReport.jsp?annot=,'
                                          + david_cats + '&currentList=0')
    index = user_download_page.find('UserDownload')
    addr = user_download_page[index+13:].split('.')[0]
    return moz_emu.download(david_prefix + 'UserDownload/' + addr + '.txt')

def getAllTerms(gene_file):
    """ Given a file with one gene per line, return
        the parsed DAVID annotation for the genes in
        the form [Category][Term]['count' | 'genes'={} ].
    
    @param gene_file: one gene per line in this file
    @return [Category][Term]['count' | 'genes'={} ]
    """

    web_page = getAllTerms_web(gene_file)
    # page ends in a newline
    return cat2term2protein_all(web_page.split('\n')[:-1])

def getSigTerms_web(fg_file, bg_file):
    """ Given a file with one gene per line, return
        the parsed DAVID significance chart for the genes in
        the form [Category][Term]['count' | 'genes'={} ].
    
    @param fg_file: one gene per line in this file
    @param bg_file: one gene per line; backround gene set
    @return string representing DAVID webpage
    """

    david_prefix = 'http://david.abcc.ncifcrf.gov/'
    #david_cats = 'GOTERM_MF_5,GOTERM_MF_2,KEGG_PATHWAY,GOTERM_BP_5,GOTERM_BP_2'
    david_cats = '22,25,36,39,48'
    moz_emu = PyMozilla.MozillaEmulator(cacher=None, trycount=0)
    files = [['fileBrowser', 'fgData', '']]
    fields = []
    moz_emu.download(david_prefix + 'tools.jsp')
    fdata = file(fg_file).read()
    fields = []
    fields.append(['idType', 'ENTREZ_GENE_ID'])
    fields.append(['Mode', 'file'])
    fields.append(['uploadType', 'list'])

    files = []
    files.append(['fileBrowser', 'testData', fdata])
    
    moz_emu.post_multipart(david_prefix + 'tools.jsp',
                           fields, files)
    
    moz_emu.download(david_prefix +  'template.css')
    moz_emu.download(david_prefix + 'scripts/sidebar.js')
    moz_emu.download(david_prefix + 'scripts/summary.js')

    # upload the background
    fdata = file(bg_file).read()
    fields = []
    fields.append(['idType', 'ENTREZ_GENE_ID'])
    fields.append(['Mode', 'file'])
    fields.append(['uploadType', 'population'])

    files = []
    files.append(['fileBrowser', 'bgData', fdata])
    
    moz_emu.post_multipart('http://david.abcc.ncifcrf.gov/tools.jsp',
                           fields, files)
    
    moz_emu.download(david_prefix + 'template.css')
    moz_emu.download(david_prefix + 'scripts/sidebar.js')
    moz_emu.download(david_prefix + 'summary.jsp')

    user_download_page = moz_emu.download(david_prefix
                                          + 'chartReport.jsp?annot='
                                          + david_cats + '&currentList=0')
    index = user_download_page.find('download/')
    #print user_download_page
    #print user_download_page[index+13:][0:50]
    addr = user_download_page[index+9:].split('.txt')[0]

    #print 'http://david.abcc.ncifcrf.gov/UserDownload/' + addr + '.txt'
    return moz_emu.download('http://david.abcc.ncifcrf.gov/data/download/'
                            + addr  + '.txt')

def getSigTerms(fg_file, bg_file):
    """ Given a file with one gene per line, return
        the parsed DAVID significance chart for the genes in
        the form [Category][Term]['count' | 'genes'={} ].
    
    @param fg_file: one gene per line in this file
    @param bg_file: one gene per line; backround gene set
    @return [Category][Term]['count' | 'p-value' | 'pathway_total' | 'genes'={} ]
    """

    web_page = getSigTerms_web(fg_file, bg_file)
    # page ends in a newline
    return cat2term2protein_sig(web_page.split('\n')[1:-1])

def getSigKEGGgenes(motif_search_results, background, cutoff):
    """ Given a foreground and background {}, find
        significant KEGG pathways with DAVID and return
        {} of foreground genes that are in these pathways
        for the given p-value cutoff.
        
    @param motif_search_results: foreground {}
    @param backgground: background gene {}
    @param cuttoff: float p-value cutoff
    @return: {} of foreground genes in sig KEGG pathways
    """

    genes = {}
    cat2term2protein = getSigTerms(motif_search_results, background)
    for path in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][path]['pvalue'] < cutoff:
            for gene in cat2term2protein['KEGG_PATHWAY'][path]['genes'].keys():
                genes[gene] = True
    return genes

def getSigKEGGpathways(motif_search_results, background, cutOff):
    """ Given a foreground and background {}, find
        significant KEGG pathways with DAVID 
        for the given p-value cutoff.
        
    @param motif_search_results: foreground {}
    @param backgground: background gene {}
    @param cuttoff: float p-value cutoff
    @return: {} of sig KEGG pathways
    """

    pathways = {}
    cat2term2protein = cat2term2protein_sigFile(afile)
    paths = {}
    for kegg in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][kegg]['pvalue'] < cutOff:
            paths['path:' + kegg.split(':')[0]] = True
    for kegg in paths.keys():
        pathways[kegg] = True
    return pathways

def getSigKEGGpathwayGenes(motif_search_results, background, cutOff):
    """ Given a foreground and background {}, find genes in
        significant KEGG pathways with DAVID 
        for the given p-value cutoff.
        
    @param motif_search_results: foreground {}
    @param backgground: background gene {}
    @param cuttoff: float p-value cutoff
    @return: {} genes in sig KEGG pathways
    """

    sig_genes = {}
    cat2term2protein = getSigTerms(motif_search_results, background)
    paths = {}
    for kegg in cat2term2protein['KEGG_PATHWAY'].keys():
        if cat2term2protein['KEGG_PATHWAY'][kegg]['pvalue'] < cutOff:
            paths['path:' + kegg.split(':')[0]] = True
    for kegg in paths.keys():
        path_genes = utils_kegg.getPathwayGenes(kegg)
        for gene in path_genes.keys():
            sig_genes[gene] = True
    return sig_genes
