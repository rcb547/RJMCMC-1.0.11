
import matplotlib
import matplotlib.pyplot

def plot_history(variable):
    
    f = open('synthdata.%s' % variable, 'r')
    history = map(float, f.readlines())
    f.close()
    
    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(2,1,1)

    cvar = variable.capitalize()

    a.plot(history)
    a.set_title('%s History' % cvar)
    a.set_xlabel('Iteration')
    a.set_ylabel(cvar)

    a = fig.add_subplot(2,1,2)

    a.hist(history, bins = 50)
    a.set_title('%s Histogram' % cvar)
    a.set_xlabel(cvar)
    a.set_ylabel('Count')
                
    matplotlib.pyplot.show()

def load_xy(filename):
    f = open(filename, 'r')
    data = map(lambda x: map(float, x.split()), f.readlines())
    f.close()
    return data

def plot_curve(variable):

    data = load_xy('data.txt')
    var = load_xy('synthdata.%s' % variable)

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(1,1,1)
    
    dx, dy, dn = zip(*data)
    vx, vy = zip(*var)
    cvar = variable.capitalize()

    a.plot(dx, dy, 'k.', label = 'Data')
    a.plot(vx, vy, 'r-', label = cvar)
    
    a.set_title('%s vs Data' % cvar)
    a.legend(loc = 'lower right')

    matplotlib.pyplot.show()

def plot_partition_x_histogram():

    data = load_xy('data.txt')
    hist = load_xy('synthdata.partition_x_hist')
    
    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(2,1,1)
    
    dx, dy, dn = zip(*data)
    vx, vy = zip(*hist)

    a.plot(dx, dy, 'k.')

    a = fig.add_subplot(2,1,2)
    a.bar(vx, vy, align='center')
    
    matplotlib.pyplot.show()
    
def plot_partition_count_histogram():

    f = open('synthdata.partitions', 'r')
    history = map(int, f.readlines())
    f.close()

    f = open('synthdata.partition_hist', 'r')
    counts = map(int, f.readlines())
    f.close()

    fig = matplotlib.pyplot.figure(1)

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

    matplotlib.pyplot.show()



