import GenomeDiagram
from Bio import SeqFeature as SeqFeatureClass
# hack to remove aliasing problem in BioPython
from Bio.SeqFeature import SeqFeature

def renderDiagram(protein_dict, output_file):
    """ Given [] of proteins, each w/ 
        [] of [start, stop name], make
        annotation diagram.
    """

    genome = GenomeDiagram.GDDiagram('test') 
    for protein in protein_dict.keys():
        a_new = genome.new_track(1, scale=0, name=protein).new_set('feature')
        for [start, stop, name] in protein_dict[protein]:
            if name=='LIG_PDZ_3':
                color = 'grey'
            else:
                color = 'black'
            seq_loc = SeqFeatureClass.FeatureLocation(start, stop)
            a_new.add_feature(SeqFeature(location=seq_loc, strand=1, id=name),
                              colour=color)    
    
    genome.draw(format='linear', fragments=1)
    genome.write(output_file, 'PDF')
