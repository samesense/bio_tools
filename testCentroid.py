import utils_cluster

genes = {'1':True, '2':True}
gene2feature = {'1':{'f1':1}, '2':{'f2':1}}
print utils_cluster.getCentroid(genes, gene2feature)
