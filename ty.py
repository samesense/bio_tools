import utils_graph
edges = utils_graph.getEdges('/home/perry/bioperry/Projects/Thesis/Data/Network/Human/HPRD/hprd.intr')
p=utils_graph.getTopXconnectedProteins(edges, 10)
