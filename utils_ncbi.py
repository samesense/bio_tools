""" Run retrieval scripts on weekends or between 9 pm and 5 am Eastern Time weekdays for any series of more than 100 requests. """
import utils_graph, Bio.Entrez
import os


# page grabber for mkFASTA
def wget(query, outputName):
    os.system("wget -tries=10 -O " + outputName
              + " 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
              + query.strip(',') + "&rettype=gp'")

# page grabber for mkFASTA
def wget_fasta(query, outputName):
    os.system("wget -O " + outputName
              + " 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
              + query.strip(',') + "&rettype=fasta'")

# parser for mkFASTA
# strip away the sequence
# from the wget file
def parseWget(wget_file, fout):
    f = open(wget_file)
    for line in f:
        print line.strip()
    f.close()
    os.system('rm ' + wget_file)

# input
#       file w/ a list of GIs
#
# output
#       fasta file @ the specified location
def mkFASTAfromFile(geneLsFile):
    genes = utils_graph.getNodes(geneLsFile)
    query = ''
    count = 0
    for gene in genes:
        query = query + gene + ','
        count += 1
        if count % 500 == 0:
            wget_name = 'ncbi.query_' + str(count/500)
            wget_files.append(wget_name)
            wget_fasta(query, wget_name)
        query = ''
        wget_fasta(query, wget_name)
    wget_name = 'ncbi.query_' + str(count/500 + 1)
    wget_fasta(query, wget_name)
    wget_files.append(wget_name)
    for f in wget_files:
        parseWget(f, fout)

def mkFASTAfromDict(genes):
    query = ''
    count = 0
    for gene in genes:
        query = query + gene + ','        
    wget_name = 'ncbi.query_temp'
    wget_fasta(query, wget_name)
    f = open(wget_name)
    for line in f:
        print line.strip()
    f.close()
    os.system('rm ' + wget_name)
        
def getChromosomeLocation(entrez_gene_ls):
    handle = Bio.Entrez.efetch(db="unigene", id="10298", rettype="genome")
    print handle.read()

