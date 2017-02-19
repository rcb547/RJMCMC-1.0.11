
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
    ax.set_xlim(-40, 40)
    ax.set_ylabel('Y')
    ax.set_ylim(-40, 40)
    ax.set_zlabel('Z')
    ax.set_zlim(-20, 50)

    plt.show()

def plot_image(filename, xmin, xmax, ymin, ymax, xsamples, ysamples):

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

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.imshow(Z, extent=[xmin, xmax, ymin, ymax])

    plt.show()

def plot_partition_count_histogram():

    f = open('gaussian.partitions', 'r')
    history = map(int, f.readlines())
    f.close()

    f = open('gaussian.partition_hist', 'r')
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
