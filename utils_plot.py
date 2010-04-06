import pylab, time, os
import numarray, Gnuplot, random

def plotPtLs(pt_ls, titles, g):
    r = random.randint(1, 100)
    plot_command = "plot 'tmp"
    rm_command = 'rm tmp'
    for index in xrange(len(pt_ls)):
        tmp_suffix = str(r + index)
        f = open('tmp' + tmp_suffix, 'w')
        #print pt_ls[index][0]
        for [x,y] in pt_ls[index]:
            f.write(str(x) + '\t' + str(y) + '\n')
        f.close()    
        plot_command = plot_command + tmp_suffix + "' title '" \
                       + titles[index] + "', 'tmp"
        rm_command = rm_command + tmp_suffix + ' tmp'
    time.sleep(5)     
    g(plot_command.rstrip(", 'tmp"))
    time.sleep(5)
    os.system(rm_command.rstrip(' tmp'))

def plotPtLs_lines(pt_ls, titles, g):
    r = random.randint(1, 100)
    plot_command = "plot 'tmp"
    rm_command = 'rm tmp'
    for index in xrange(len(pt_ls)):
        tmp_suffix = str(r + index)
        f = open('tmp' + tmp_suffix, 'w')
        #print pt_ls[index][0]
        for [x,y] in pt_ls[index]:
            f.write(str(x) + '\t' + str(y) + '\n')
        f.close()    
        plot_command = plot_command + tmp_suffix + "' title '" \
                       + titles[index] + "' with linespoints, 'tmp"
        rm_command = rm_command + tmp_suffix + ' tmp'
    time.sleep(5)     
    g(plot_command.rstrip(", 'tmp"))
    time.sleep(5)
    os.system(rm_command.rstrip(' tmp'))

def heatPlot(matrix_file, out_file, title):
    g = Gnuplot.Gnuplot(debug=0)
    g("set terminal pdf")
    g("set output '" + out_file + "'")
    g('set tic scale 0')
    g("set clabel 'P-value'")
    g("set cbtics ('1' 0, '.8' 20, '.6' 40, '.4' 60, '.2' 80, '0' 100)")
    g("set title '" + title + "'")
    g.xlabel('UVXY')
    g.ylabel('Scope')
    g("set ytics ('1' 0, '2' 1, '3' 2, '4' 3, '5' 4, '6' 5, '7' 6, '8' 7)")
    g("set xtics ('AAAA' 0, 'ATAG' 50, 'CGCA' 100, 'GCCG' 150, 'TAGA' 200, 'TTGG' 250)") 
    #g('set palette rgbformula -7,2,-7')
    g('set xrange [-0.5:255.5]')
    g('set yrange [-0.5:7.5]')
    g('unset key')
    #g('unset cbtics')
    g('set view map')
    g("plot '" + matrix_file + "' matrix with image")
    g('set terminal x11')
    g('set output')    

def boxPlot(dataset_ls, title, filename):
    r = pylab.boxplot(dataset_ls)
    pylab.setp(r['medians'], color='black')
    F = pylab.gcf()
    DefaultSize = F.get_size_inches()
    F.set_figsize_inches( (DefaultSize[0]*10, DefaultSize[1]*2) )
    pylab.title(title)
    pylab.savefig(filename)

# A helper function to make histograms
def histOutline(dataIn, binsIn=None):
    """
    Make a histogram that can be plotted with plot() so that
    the histogram just has the outline rather than bars as it
    usually does.
    """
    if (binsIn == None):
        (en, eb) = pylab.matplotlib.mlab.hist(dataIn, bins=50, normed=True)
        binsIn = eb
    else:
        (en, eb) = pylab.matplotlib.mlab.hist(dataIn, bins=binsIn)

    stepSize = binsIn[1] - binsIn[0]

    bins = numarray.zeros(len(eb)*2 + 2, type=numarray.Float)
    data = numarray.zeros(len(eb)*2 + 2, type=numarray.Float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        data[2*bb + 1] = en[bb]
        data[2*bb + 2] = en[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)

def pdfHist(data, b, filename):
    (n, bins_res, patches) = pylab.hist(data, bins=b)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', .75)
    pylab.plot(bins_res)
    pylab.savefig(filename, dpi=72)

def pdfHistFromFile(input_file, b, filename):
    data = []
    f = open(input_file)
    for line in f:
        data.append(float(line.strip()))
    f.close()
    (n, bins_res, patches) = pylab.hist(data, bins=b)
    pylab.setp(patches, 'facecolor', 'g', 'alpha', .75)
    pylab.plot(bins_res)
    pylab.savefig(filename, dpi=72)


def pdfHist4(data1, data2, data3, data4, b, filename):
    (n, bins_res1, patches1) = pylab.hist(data1, bins=b)
    (n, bins_res2, patches2) = pylab.hist(data2, bins=b)
    (n, bins_res3, patches3) = pylab.hist(data3, bins=b)
    (n, bins_res4, patches4) = pylab.hist(data4, bins=b)
    pylab.setp(patches1, 'facecolor', 'g', 'alpha', .75)
    pylab.setp(patches2, 'facecolor', 'r', 'alpha', .75)
    pylab.setp(patches3, 'facecolor', 'b', 'alpha', .75)
    pylab.setp(patches4, 'facecolor', 'y', 'alpha', .75)
    pylab.plot(bins_res1, bins_res2, bins_res3, bins_res4)
    pylab.savefig(filename, dpi=72)

def pdfHist4_outline(data1, data2, data3, data4, labels, title, filename):
    #(n, bins_res1, patches1) = pylab.hist(data1, bins=b)
    (bins1, n1) = histOutline(data1)
    #(n, bins_res2, patches2) = pylab.hist(data2, bins=b)
    (bins2, n2) = histOutline(data2)
    (bins3, n3) = histOutline(data3)
    (bins4, n4) = histOutline(data4)
    #(n, bins_res3, patches3) = pylab.hist(data3, bins=b)
    #(n, bins_res4, patches4) = pylab.hist(data4, bins=b)
    #pylab.setp(patches1, 'facecolor', 'g', 'alpha', .75)
    #pylab.setp(patches2, 'facecolor', 'r', 'alpha', .75)
    #pylab.setp(patches3, 'facecolor', 'b', 'alpha', .75)
    #pylab.setp(patches4, 'facecolor', 'y', 'alpha', .75)
    #pylab.plot(bins_res1, bins_res2, bins_res3, bins_res4)
    pylab.plot(bins1, n1, 'g-', label=labels[0])
    pylab.plot(bins2, n2, 'r--', label=labels[1])
    pylab.plot(bins3, n3, 'b-', label=labels[2])
    pylab.plot(bins4, n4, 'k--', label=labels[3])
    pylab.legend()
    pylab.title(title)
    pylab.savefig(filename, dpi=72)
    pylab.clf()

def pdfHist2_outline(data1, data4, labels, title, filename):
    (bins1, n1) = histOutline(data1)    
    (bins4, n4) = histOutline(data4)
    pylab.plot(bins1, n1, 'g-', label=labels[0])   
    pylab.plot(bins4, n4, 'k--', label=labels[1])
    pylab.legend()
    pylab.title(title)
    pylab.savefig(filename, dpi=72)
    pylab.clf()

def pdfHist1_outline(data1, title, filename):
    (bins1, n1) = histOutline(data1)
    pylab.plot(bins1, n1, 'g-')   
    pylab.title(title)
    pylab.savefig(filename, dpi=72)
    pylab.clf()
