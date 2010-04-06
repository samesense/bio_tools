import scipy.stats, sys, decimal, inspect
import statlib.stats
import rpy2.robjects as robjects
r = robjects.r
r('library(tseries)')
r('library(diptest)')
r('library(nortest)')

wil_func = r('function(x, y) wilcox.test(x,y)')
wil_func_gtr = r("function(x, y) wilcox.test(x,y,alternative='greater')")
wil_func_less = r("function(x, y) wilcox.test(x,y,alternative='less')")
fisher = r("function(o1, o2, e1, e2) fisher.test(matrix(c(o1,o2,e1,e2),nrow=2))")
fisher_positive = r("function(o1, o2, e1, e2) fisher.test(matrix(c(o1,o2,e1,e2),nrow=2),alternative='greater')")
jarque = r('function(x) jarque.bera.test(x)$p.value')
unimodal_permute = r('function(x) dip(runif(x[1]))')
unimodal = r('function(x) dip(x)')
pcc_func = r('function(x, y) cor.test(x,y)')

anderson = r('function(x) ad.test(x)')

def normality_test_anderson(ls):
    """ p-values greater than .05 means normal """

    result = anderson(robjects.FloatVector(ls))
    #print 'result', result
    return result

def unimodal_test(ls):
    """ The hartigan dip test stat measures the depearture from unimodal,
        so smaller is better (more unimodal)
    """

    sample_size = 100
    ls.sort()
    n = robjects.FloatVector([len(ls)])
    test_stat = unimodal(robjects.FloatVector(ls))[0]
    btr_eql = 0
    for i in xrange(sample_size):
        test_stat_permute = unimodal_permute(n)[0]
        if test_stat_permute <= test_stat:
            btr_eql += 1
    return [btr_eql, sample_size]

def normality_test_jarque(ls):
    rFeature1 = robjects.FloatVector(ls)
    result = jarque(rFeature1)
    
    return result[0]

def pcc(features1, features2):
    """ pearson's correlation coefficient from R
        features are lists of floats """

    rFeature1 = robjects.FloatVector(features1)
    rFeature2 = robjects.FloatVector(features2)
    return pcc_func(rFeature1, rFeature2)[3][0]

def wilcox(ls1, ls2):
    rLs1 = robjects.FloatVector(ls1)
    rLs2 = robjects.FloatVector(ls2)
    result = wil_func(rLs1, rLs2)
#    print result, result[1], result[2], result[3]
    return result[2][0]

def wilcox_gtr(ls1, ls2):
    rLs1 = robjects.FloatVector(ls1)
    rLs2 = robjects.FloatVector(ls2)
    result = wil_func_gtr(rLs1, rLs2)
#    print result, result[1], result[2], result[3]
    return result[2][0]

def wilcox_less(ls1, ls2):
    rLs1 = robjects.FloatVector(ls1)
    rLs2 = robjects.FloatVector(ls2)
    result = wil_func_less(rLs1, rLs2)
#    print result, result[1], result[2], result[3]
    return result[2][0]

def mannwhitney(ls1, ls2):
    return scipy.stats.mannwhitneyu(ls1, ls2)[1]

def wilcoxon_pval(ls1, ls2):
    return scipy.stats.wilcoxon(ls1, ls2)[1]

def fisher_pval(row1, row2):
    """ alternative is two.sided, greater, less
        alternative is not implemented yet,
        all tests are 2 sided. 
    """

    result = fisher(row1[0], row1[1],
                    row2[0], row2[1])
    #print result
    return result[0][0]

def fisher_positive_pval(row1, row2):
    """ alternative is two.sided, greater, less
        alternative is not implemented yet,
        all tests are 2 sided. 

        Are most values in the upper left and lower right?
        http://www.langsrud.com/fisher.htm
    """

    result = fisher_positive(row1[0], row1[1],
                             row2[0], row2[1])
    #print result
    return result[0][0]

#def chi2_pval(observed_ls, expected_ls):
#    dof = len(observed_ls) - 1
#    sum = float(0)
#    for index in xrange(len(observed_ls)):
#        sum += float((observed_ls[index] - expected_ls[index])**2) / \
#               float(expected_ls[index])
#    print 'stat', sum
#    return float(1) - float(scipy.stats.chi2.cdf(sum, dof))

def combo3_decimal(m, n, mult, num_start):
    t = max([n,m-n])
    t2 = min([n,m-n])
    #num = 1
    i_n = num_start
    i_d = t2
    while i_n > t or i_d > 1:
        if i_n > t:
            mult *= float(i_n)
            i_n -= 1
        if i_d > 1:
            mult /= float(i_d)
            i_d -= 1
    return mult

def combo3(m, n, mult, num_start):
    t = max([n,m-n])
    t2 = min([n,m-n])
    #num = 1
    i_n = num_start
    i_d = t2
    while i_n > t or i_d > 1:
        if i_n > t:
            mult *= float(i_n)
            i_n -= 1
        if i_d > 1:
            mult /= float(i_d)
            i_d -= 1
    return mult

def combo5(m, n, p1, num_start):
    t = max([n,m-n])
    t2 = min([n,m-n])
    i_n = num_start
    i_d = t2
    mult = float(1)
    used = False
    while i_n > t or i_d > 1:
        #print i_n, i_d, mult
        if i_n > t:
            temp = mult * float(i_n)
            while temp == float('inf') and (i_d > 1  or not used):
                if not used:
                    temp = mult * p1
                    used = True
                elif id > 1:
                    mult /= float(i_d)
                    temp = mult * float(i_n)
                    i_d -= 1
            mult = temp
            i_n -= 1
        if i_d > 1:
            temp = mult / float(i_d)
            while temp == float(0) and i_n > t:
                mult *= float(i_n)
                temp = mult / float(i_d)
                i_n -= 1
            mult = temp
            i_d -= 1
    if mult == float('inf'): sys.exit(0)
    return mult

def prob3(g, d, o, m):
    """ g is the background set you are drawing from
        d is the total foreground (all predictions)
        o is the size of the target set (HHE) in the background set
        m is the intersection of o and d

        prob of matching m ore more in o when choosing d from g
        
        implemented with R's hypergeometric test
        r.phyper(HHP_HHE_match, HHE, total_genes - HHE, HHP)
    """

    if m == 0:
        return float(1)
    kwargs = {'lower.tail':False}
    # m is observation, match, i
    # o is number of black balls, proteins with domain
    # n is number of red balls, proteins w/o domain (total{g}-o)
    # d is HHP, draws, N
    return r.phyper(m, o, g-o, d, **kwargs)[0]

def prob3_old(g, d, o, m):
    """ g is the background set you are drawing from
        d is the total foreground (all predictions)
        o is the size of the target set (HHE) in the background set
        m is the intersection of o and d

        prob of matching m ore more in o when choosing d from g
    """
    
    s = float(0)
    for x in xrange(m,min([d+1,o+1])):
        num_start = d
        #c = binomial(d, x)
        #c = combo(d, x)
        #c1 = combo(o, x)
        #c2 = combo(g-o, d-x)
        p1 = float(1)
        for i in xrange(0,x):
            #print 'i', o-i, g-i
            temp = p1 * float(o-i)/float(g-i)
            prod = float(num_start)
            while temp == float(0) and num_start > x:
                temp = prod * p1 * float(o-i)/float(g-i)
                if temp == float('inf'): print 'infinite float'
                num_start -= 1
                prod *= float(num_start)
            p1 = temp
        p2 = float(1)
        #print 'p1', p1
        for j in xrange(0,d-x):
            #print 'x', x, 'j', j, ':::::', g-o-j, g-j-x
            p2 *= float(g-o-j)/float(g-j-x)
        #if p2 == float('inf'): print 'infinite float'
        #print c, p1, p2
        c = combo3(d, x, p1, num_start)
        if c == float('inf'):
            #print p1, p2, d, x, num_start
            #print d, x, p2, num_start
            s += p2 * combo5(d, x, p1, num_start)
            #if s == float('inf'): print d, x, p2, p1, num_start
        else:
            s += p2 * combo3(d, x, p1, num_start)
        
    return s

def fac(n):
    value = reduce(lambda i, j : i * j, range(1, n + 1))
    return value

def com(x, y):
    if y == 0:
        return 1
    return fac(x) / (fac(y) * fac(x-y))

def probOfMatchingXinTwoSets(set1_size, set2_size, match_size, total):
    #x = decimal.Decimal(set1_size)
    #y = decimal.Decimal(set2_size)
    #m = decimal.Decimal(match_size)
    #t = decimal.Decimal(total)
    prob = decimal.Decimal(0)
    denom = decimal.Decimal(com(total, set1_size) * com(total, set2_size))
    for m in xrange(match_size, min([set1_size, set2_size])+1):
        num = decimal.Decimal(com(total, m) * com(total-m, set1_size-m) \
              * com(total-set1_size, set2_size-m))        
        if denom - num < 0:
            sys.stderr.write('ERROR: probOfMatchXinTwoSets has denom less than num\n')
            sys.exit(0)
        prob += num/denom
    return prob
#fac(total)**2/(fac(total-set1_size)*fac(set1_size)*fac(total-set2_size)*fac(set2_size))

def normalize(ls):
    return statlib.stats.zs(ls)
