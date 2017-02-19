
import matplotlib
import matplotlib.pyplot

def plot_history(variable):
    
    f = open('%s.txt' % variable, 'r')
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

    data = load_xy('actual.txt')
    var = load_xy('%s.txt' % variable)

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(1,1,1)
    
    dx, cdy, dy = zip(*data)
    vx, vy = zip(*var)
    cvar = variable.capitalize()

    a.plot(dx, dy, 'k.', label = 'Data')
    a.plot(vx, vy, 'r-', label = cvar)
    
    a.set_title('%s vs Data' % cvar)
    a.legend(loc = 'lower right')

    matplotlib.pyplot.show()

def plot_partition_x_histogram():

    data = load_xy('actual.txt')
    hist = load_xy('partition_x_hist.txt')
    
    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(2,1,1)
    
    dx, cdy, dy = zip(*data)
    vx, vy = zip(*hist)

    xmin = min([dx[0], vx[0]])
    xmax = max([dx[len(dx) - 1], vx[len(vx) - 1]])

    a.plot(dx, dy, 'k.')
    a.set_title('Data')
    a.set_xlim(xmin, xmax)
    
    a = fig.add_subplot(2,1,2)
    a.bar(vx, vy, align='center', width=vx[1] - vx[0])
    a.set_title('Partition Location Histogram')
    a.set_ylabel('Count')
    a.set_xlabel('x')
    a.set_xlim(xmin, xmax)
    
    matplotlib.pyplot.show()
    
def plot_partition_count_histogram():

    f = open('partitions.txt', 'r')
    history = map(int, f.readlines())
    f.close()

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(2,1,1)
    
    a.plot(history)
    a.set_title('History of No. Partitions')

    a = fig.add_subplot(2,1,2)

    minp = min(history)
    maxp = max(history)

    print minp, maxp, range(minp, maxp)

    bins = map(lambda x: float(x) - 0.5, range(minp, maxp + 2))
    print bins
    a.hist(history, bins = bins)
    a.set_xticks(range(minp, maxp + 1))

    a.set_title('No. Partitions Histogram')
    a.set_xlabel('No. Partitions')
    a.set_ylabel('Counts')

    matplotlib.pyplot.show()



