import sys, os, numpy, math, utils_graph, utils_stats

def loadClusters(afile):
    clusters = []
    f_open = open(afile)
    for line in f_open:
        [cluster, gene] = line.strip().split('\t')
        cluster_index = int(cluster)
        if len(clusters) < cluster_index + 1:
            clusters.append({})
        clusters[cluster_index][gene] = True
    f_open.close()

    return clusters

def getClusterGenes(clusters):
    genes = {}
    for cluster in clusters:
        for gene in cluster.keys():
            genes[gene] = True
    return genes

def getBestSim(neighbor, neighbors_choose, blast_sims):
    if neighbor in neighbors_choose:
        return float(1)
    best_sim = float(0)
    for choice in neighbors_choose:
        if blast_sims[neighbor].has_key(choice):
            new_score = blast_sims[neighbor][choice]
            if new_score > best_sim:
                best_sim = new_score
    return best_sim

def computeSim(neighbors1, neighbors2, blast_sims):
    blast_sum1 = float(0)
    blast_sum2 = float(0)
    for neighbor in neighbors1:
        blast_sum1 += getBestSim(neighbor, neighbors2, blast_sims)
    for neighbor in neighbors2:
        blast_sum2 += getBestSim(neighbor, neighbors1, blast_sims)
    return max([blast_sum1 / float(len(neighbors1)),
               blast_sum2 / float(len(neighbors2))])

def getBlast(blast_file):
    blast_sims = {}
    blast_f = open(blast_file)
    for line in blast_f:
        [ver1, ver2, val] = line.strip().split('\t')
        if not blast_sims.has_key(ver1):
            blast_sims[ver1] = {}
        if not blast_sims.has_key(ver2):
            blast_sims[ver2] = {}
        blast_sims[ver1][ver2] = float(val)
    blast_f.close()
    return blast_sims

def getEntrez2Version(translation_file):
    entrez2version = {}
    f_translation = open(translation_file)
    for line in f_translation:
        [version, entrez] = line.strip().split('\t')
        if not entrez2version.has_key(entrez):
            entrez2version[entrez]= {}
        entrez2version[entrez][version] = True
    f_translation.close()
    return entrez2version

def gmeans(feature_file, cluster_count, algorithm):
    """ algorithm is 
           s: spherical k-means algorithm (default)
           e: euclidean k-means algorithm
           b: information bottleneck algorithm
	   k: kullback_leibler k-means algorithm
	   d: diametric k-means algorithm

           return {} of clusterid 2 gene
    """

    os.system('cut -f 1 ' + feature_file
              + ' > xs')
    os.system('cut -f 2 ' + feature_file
              + ' > ys')
    os.system('words2ids.pl xs > xs.dict')
    os.system('words2ids.pl ys > ys.dict')
    os.system('rm xs ys')
    
    os.system('encodepairs.pl ' + feature_file
              + ' xs.dict ys.dict > temp.coo')
    os.system('sim2ccs.pl -t temp.coo')
    os.system('/opt/gmeans/gmeans- temp.coo -c '
              + str(cluster_count) + ' -a ' 
              + algorithm)

    #print 'GENE TO CLUSTER'
    cluster2gene= {}
    x_dict_file = open('xs.dict')
    results_file = open('temp.coo_tfn_doctoclus.'
                        + str(cluster_count))
    results_file.readline()
    line = x_dict_file.readline()
    while line != '':
        gene = line.split()[1].strip()
        cluster_assignment = results_file.readline().strip()
        if not cluster2gene.has_key(cluster_assignment):
            cluster2gene[cluster_assignment] = {}
        cluster2gene[cluster_assignment][gene] = True
        line = x_dict_file.readline()
    results_file.close()
    x_dict_file.close()
    os.system('rm xs.dict ys.dict temp.coo_dim temp.coo_row_ccs temp.coo_col_ccs temp.coo temp.coo_tfn_nz temp.coo_tfn_doctoclus.' + str(cluster_count))
    return cluster2gene

def getGmeansClusters(cluster_file):
    clusters = {}
    fopen = open(cluster_file)
    line = fopen.readline()    
    while line != '':
        [gene, cluster] = line.strip().split('\t')
        if not clusters.has_key(cluster):
            clusters[cluster] = {}
        clusters[cluster][gene] = True
        line = fopen.readline()
    fopen.close()
    return clusters

def getGmeansClusters_entrez(cluster_file):
    clusters = {}
    fopen = open(cluster_file)
    line = fopen.readline()   
    while line != '':
        [gene, cluster] = line.strip().split('\t')
        if not clusters.has_key(cluster):
            clusters[cluster] = {}
        clusters[cluster][gene] = True
        line = fopen.readline()
    fopen.close()
    return clusters

def loadFeatures(feature_file):
    pass

def dumpFeatures(feature_mat, afile):
    pass

def getFeatureFreq(feature_mat):
    pass

def within_sum_sqs(cluster_file, feature_file):
    pass

def mkNullFeatures(freq):
    pass

def gmeans_gap_stat(feature_file, sample_count, regex_for_clusters):
    pass

def mkNullData(feature_file, sample_count, dump_dir):
    feature_mat = loadFeatures(feature_file)
    freq = getFeatureFreq(feature_mat)
    for random in xrange(sample_count):
        rand_data = mkNullFeatures(freq)
        dumpFeatures(rand_data,
                     dump_dir + 'sample_' + str(random) + '.features')

def getDistance(gene, centroid_features, gene2feature):
    #dotprod distance abs value
    features = {}
    for feature in gene2feature[gene].keys():
        features[feature] = True
    for feature in centroid_features.keys():
        features[feature] = True
    feature_keys = features.keys()
    gene_vec = []
    centroid_vec = []
    normalized_gene_vec_feature = float(1) / \
                                  math.sqrt(len(gene2feature[gene].keys()))
    for feature in feature_keys:
        if gene2feature[gene].has_key(feature):
            gene_vec.append(normalized_gene_vec_feature)
        else:
            gene_vec.append(float(0))
        if centroid_features.has_key(feature):
            centroid_vec.append(centroid_features[feature])
        else:
            centroid_vec.append(float(0))
    cosine = numpy.dot(numpy.array(gene_vec), numpy.array(centroid_vec))
    if float(cosine) > float(1):
        print 'cosine is', cosine
        sys.exit(0)
    return numpy.dot(numpy.array(gene_vec), numpy.array(centroid_vec))

def getDotProd(gene1, gene2, gene2feature):
    features = {}
    for feature in gene2feature[gene1].keys():
        features[feature] = True
    for feature in gene2feature[gene2].keys():
        features[feature] = True
    feature_keys = features.keys()
    gene1_vec = []
    gene2_vec = []
    normalized_gene1_vec_feature = float(1) / \
                                   math.sqrt(len(gene2feature[gene1].keys()))
    normalized_gene2_vec_feature = float(1) / \
                                   math.sqrt(len(gene2feature[gene2].keys()))
    for feature in feature_keys:
        if gene2feature[gene1].has_key(feature):
            gene1_vec.append(normalized_gene1_vec_feature)
        else:
            gene1_vec.append(float(0))
        if gene2feature[gene2].has_key(feature):
            gene2_vec.append(normalized_gene2_vec_feature)
        else:
            gene2_vec.append(float(0))
    cosine = numpy.dot(numpy.array(gene1_vec), numpy.array(gene2_vec))
    #if float(cosine) > float(1):
    #    print 'cosine is', cosine
    #    sys.exit(0)
    return cosine

def l2_normalize(features):
    sq_sum = float(0)
    for feature in features.keys():
        sq_sum += features[feature]**2
    norm = math.sqrt(sq_sum)
    norm_features = {}
    for feature in features.keys():
        norm_features[feature] = features[feature] / norm
    return norm_features

def getCentroid(genes, gene2feature):
    #l2 norm of sum of l2 norms of genes
    features = {}
    norms = {}
    for gene in genes:
        norms[gene] = float(1) / math.sqrt(len(gene2feature[gene].keys()))
        for feature in gene2feature[gene].keys():
            features[feature] = float(0)
    for feature in features.keys():
        for gene in genes:
            if gene2feature[gene].has_key(feature):
                features[feature] += norms[gene]
    return l2_normalize(features)  

def getDotProdSum(dotprods, gene, cluster_genes):
    sum = float(0)
    for cluster_gene in cluster_genes.keys():
        k1 = cluster_gene + ':' + gene
        k2 = gene + cluster_gene
        if dotprods.has_key(k1):
            sum += dotprods[k1]
        else:
            sum += dotprods[k2]
    return sum

def getDotProds(test_genes, train_genes, gene2feature):
    dotprods = {}
    for gene1 in test_genes.keys():
        for gene2 in train_genes.keys():
            dotprods[gene1 + ':' + gene2] = getDotProd(gene1, gene2, gene2feature)
    return dotprods            

#def getCentroid(genes, gene2feature):
#    #l2 norm of sum of l2 norms of gene vecs
#    feature_totals = {}
#    for gene in genes:
#        for feature in gene2feature[gene].keys():
#            if not feature_totals.has_key(feature):
#                feature_totals[feature] = 0
#            feature_totals[feature] += 1
#    sum_of_sqrs = 0
#    for feature in feature_totals.keys():
#        sum_of_sqrs += feature_totals[feature]**2
#    divisor = math.sqrt(sum_of_sqrs)
#    for feature in feature_totals.keys():
#        feature_totals[feature] = float(feature_totals[feature]) / divisor
#    return feature_totals

def findBestCluster(gene, sig_clusters, clusters, gene2feature):
    best_clusters = []
    best_distance = float(0)
    #best_precision = float(-1)
    #best_virus_protein = 'na'
    for cluster in clusters.keys():
        centroid_features = getCentroid(clusters[cluster].keys(),
                                        gene2feature)
        distance = getDistance(gene, centroid_features, gene2feature)
        if distance < float(0):
            print 'DISTANCE NEGATIVE'
        if distance > float(0):
            if distance > best_distance:
                best_distance = distance
                best_clusters = [cluster]
                #best_virus_protein = virus_protein
                #best_precision = precision
            elif distance == best_distance:
                best_clusters.append(cluster)                
    if len(best_clusters) > 1:
        print 'CLUSTER TIE'
        sys.stderr.write('CLUSTER TIE\n')
        sys.exit(0)
    best_cluster = best_clusters[0]
    virus_precision = []
    if sig_clusters.has_key(best_cluster):
        for virus in sig_clusters[best_cluster].keys():
            virus_precision.append([virus, 
                                     sig_clusters[best_cluster][virus]])
    else:
        virus_precision = []
    return virus_precision
        
    #return ['NEF', '90']
    return [best_vir, best_precision] 

def findClosestCluster(gene, clusters, gene2feature):
    closest_cluster = []
    best_sim = float(1000)
    for cluster in clusters.keys():
        centroid_features = getCentroid(clusters[cluster].keys(),
                                        gene2feature)
        sim = getDistance(gene, centroid_features, gene2feature)
        
        if sim < float(0):
            print 'DISTANCE NEGATIVE'
        if sim > best_sim:
            best_sim = sim
            closest_cluster = [cluster]
        else:
            closest_cluster.append(cluster)
    if len(closest_cluster) > 1:
        return 'tie'
    return closest_cluster[0] 

def findClosestClusterGivenCentroids(gene, clusters, centroids, gene2feature):
    closest_cluster = []
    best_sim = float(-2)
    for cluster in clusters.keys():
        sim = getDistance(gene, centroids[cluster], gene2feature)
        #if sim != float(0):
        #    print 'non zero', sim
        if sim < float(0):
            print 'DISTANCE NEGATIVE'
        if sim > best_sim:
            best_sim = sim
            closest_cluster = [cluster]
        elif sim == best_sim:
            closest_cluster.append(cluster)

    if len(closest_cluster) > 1:
        #if len(closest_cluster) != 680:
        #    print 'len', len(closest_cluster), best_sim
        #    print 'gene', gene2feature[gene]
        #    for cluster in closest_cluster:
        #        feature_ls = []
        #        for feature in centroids[cluster].keys():
        #            if gene2feature[gene].has_key(feature):
        #                feature_ls.append([feature, centroids[cluster][feature]])
        #        print 'cluster', feature_ls, getDistance(gene, 
        #                                                 centroids[cluster], 
        #                                                 gene2feature)
        #    sys.exit(0)
        #print 'tie', len(closest_cluster)
        return 'tie'
    return closest_cluster[0]

def findClosestClusterGivenDotProds(gene, clusters, dotProds):
    closest_cluster = []
    best_sim = float(-2)
    for cluster in clusters.keys():
        sim = getDotProdSum(dotProds, gene, clusters[cluster])        
        if sim < float(0):
            print 'DISTANCE NEGATIVE'
        if sim > best_sim:
            best_sim = sim
            closest_cluster = [cluster]
        else:
            closest_cluster.append(cluster)
    if len(closest_cluster) > 1:
        print 'len', len(closest_cluster), best_sim
        return 'tie'
    return closest_cluster[0]  

def loadFeatureVecs(features_file):
    gene2feature = {}
    afile = open(features_file)
    for line in afile:
        [gene, feature] = line.strip().split('\t')
        if not gene2feature.has_key(gene):
            gene2feature[gene] = {}
        gene2feature[gene][feature] = True
    afile.close()
    return gene2feature

def getEnrichedClusters(gmeans_clusters, target_sets, background_genes, pval_cut):
    """ Return {} of cluster_id to [enriched_cat, enriched_genes {}, pval]
        target_sets is {} of cat to genes.
    """

    # check to make sure everything is in background
    for cluster in gmeans_clusters.keys():
        for gene in gmeans_clusters[cluster].keys():
            if not background_genes.has_key(gene):
                print 'cluster gene not in background genes'
                sys.exit(0)
    for target in target_sets.keys():
        for gene in target_sets[target].keys():
             if not background_genes.has_key(gene):
                print 'target set gene not in background genes'
                sys.exit(0)

    ret_clusters = {}
    background_len = len(background_genes.keys())
    for a_set in target_sets.keys():
        target_genes = target_sets[a_set]
        target_len = len(target_genes.keys())
        for cluster in gmeans_clusters.keys():
            cluster_genes = gmeans_clusters[cluster]
            cluster_len = len(cluster_genes.keys())
            match_genes = utils_graph.intersectLists([cluster_genes,
                                                       target_genes])
            match_len = len(match_genes.keys())
            pval = utils_stats.prob3(background_len,
                                     cluster_len,
                                     target_len,
                                     match_len)
            if pval < pval_cut:
                if not ret_clusters.has_key(cluster):
                    ret_clusters[cluster] = {}
                cluster_per = int(float(100)*float(match_len)/float(cluster_len))
                target_per = int(float(100)*float(match_len)/float(target_len))
                ret_clusters[cluster][a_set] = [match_genes, cluster_per, target_per, pval]
    return ret_clusters
            
