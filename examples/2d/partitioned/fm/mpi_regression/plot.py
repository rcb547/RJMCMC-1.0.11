
import numpy

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

def plot_surface(filename, xmin, xmax, ymin, ymax, xsamples, ysamples):

    xaxis = numpy.linspace(xmin, xmax, xsamples)
    yaxis = numpy.linspace(ymin, ymax, ysamples)
    X, Y = numpy.meshgrid(xaxis, yaxis)

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    if len(lines) != ysamples:
        print 'error: expected 100 lines'
        return 
        
    Z = numpy.zeros((xsamples, ysamples))
    for i, line in enumerate(lines):
        for j, z in enumerate(map(float, line.split())):
            Z[i, j] = z

    #
    # Original data
    #
    f = open("data.txt", "r")
    lines = f.readlines()
    f.close()
    datax = []
    datay = []
    dataz = []
    for line in lines:
        x, y, z = map(float, line.split())
        datax.append(x)
        datay.append(y)
        dataz.append(z)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_surface(X, Y, Z, rstride=5, cstride=5, alpha=0.3)
    ax.scatter(datax, datay, dataz)

    ax.set_xlabel('X')
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel('Y')
    ax.set_ylim(ymin, ymax)
    ax.set_zlabel('Z')

    plt.show()

def plot_comparison_surface(filename1, filename2, xmin, xmax, ymin, ymax, xsamples, ysamples):

    xaxis = numpy.linspace(xmin, xmax, xsamples)
    yaxis = numpy.linspace(ymin, ymax, ysamples)
    X, Y = numpy.meshgrid(xaxis, yaxis)

    f = open(filename1, 'r')
    lines = f.readlines()
    f.close()

    if len(lines) != ysamples:
        print 'error: expected %d lines' % ysamples
        return 
        
    Z1 = numpy.zeros((xsamples, ysamples))
    for i, line in enumerate(lines):
        for j, z in enumerate(map(float, line.split())):
            Z1[i, j] = z

    f = open(filename2, 'r')
    lines = f.readlines()
    f.close()

    if len(lines) != ysamples:
        print 'error: expected %d lines' % ysamples
        return 
        
    Z2 = numpy.zeros((xsamples, ysamples))
    for i, line in enumerate(lines):
        for j, z in enumerate(map(float, line.split())):
            Z2[i, j] = z

    #
    # Original data
    #
    f = open("data.txt", "r")
    lines = f.readlines()
    f.close()
    datax = []
    datay = []
    dataz = []
    for line in lines:
        x, y, z = map(float, line.split())
        datax.append(x)
        datay.append(y)
        dataz.append(z)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_surface(X, Y, Z1, rstride=5, cstride=5, alpha=0.3, color='r')
    ax.plot_surface(X, Y, Z2, rstride=6, cstride=6, alpha=0.3, color='b')
    ax.scatter(datax, datay, dataz)

    ax.set_xlabel('X')
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel('Y')
    ax.set_ylim(ymin, ymax)
    ax.set_zlabel('Z')

    plt.show()

def plot_partition_count_histogram():

    f = open('regression.partitions', 'r')
    history = map(int, f.readlines())
    f.close()

    f = open('regression.partition_hist', 'r')
    counts = map(int, f.readlines())
    f.close()

    fig = plt.figure(1)

    a = fig.add_subplot(2,1,1)
    
    a.plot(history)
    a.set_title('History of No. Partitions')

    a = fig.add_subplot(2,1,2)

    npart = map(lambda x: x + 2, range(len(counts)))
    a.bar(npart, counts, align='center')
    a.set_xticks(npart)

    a.set_title('No. Partitions Histogram')
    a.set_xlabel('No. Partitions')
    a.set_ylabel('Counts')

    plt.show()

