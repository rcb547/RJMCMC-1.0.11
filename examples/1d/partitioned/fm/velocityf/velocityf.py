
import matplotlib
import matplotlib.pyplot

def plot_history(variable):
    
    f = open('velocityf.%s' % variable, 'r')
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
    
    dx, dy, dn = zip(*data)
    vx, vy = zip(*var)
    cvar = variable.capitalize()

    a.plot(dx, dy, 'k.', label = 'Data')
    a.plot(vx, vy, 'r-', label = cvar)
    
    a.set_title('%s vs Data' % cvar)
    a.legend(loc = 'lower right')

    matplotlib.pyplot.show()

def plot_credible():

    data = load_xy('actual.txt')
    cred = load_xy('credible.txt')

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(1,1,1)
    
    dx, dy, dn = zip(*data)
    vx, vmin, vmax = zip(*cred)

    a.plot(dx, dy, 'k.', label = 'Data')
    a.plot(vx, vmin, 'r-', label = 'Credible Min')
    a.plot(vx, vmax, 'b-', label = 'Credible Max')
    
    a.set_title('Credible Range vs Data')
    a.legend(loc = 'lower right')

    matplotlib.pyplot.show()

def plot_partition_x_histogram():

    data = load_xy('actual.txt')
    hist = load_xy('partition_x_hist.txt')
    
    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(2,1,1)
    
    dx, dy, dn = zip(*data)
    vx, vy = zip(*hist)

    a.plot(dx, dy, 'k.')

    a = fig.add_subplot(2,1,2)
    a.bar(vx, vy, align='center', width = vx[1] - vx[0])
    
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



