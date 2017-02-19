
import matplotlib
import matplotlib.pyplot

def plot_history(variable):
    
    f = open('cubic.%s' % variable, 'r')
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
    var = load_xy('cubic.%s' % variable)

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(1,1,1)
    
    dx, dy, n = zip(*data)
    vx, vy = zip(*var)
    cvar = variable.capitalize()

    a.plot(dx, dy, 'k.', label = 'Data')
    a.plot(vx, vy, 'r-', label = cvar)
    
    a.set_title('%s vs Data' % cvar)
    a.legend(loc = 'upper left')

    matplotlib.pyplot.show()
    
def plot_credible():

    data = load_xy('data.txt')
    lvar = load_xy('cubic.credible_min')
    hvar = load_xy('cubic.credible_max')

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(1,1,1)
    
    dx, dy, n = zip(*data)
    lx, ly = zip(*lvar)
    hx, hy = zip(*hvar)

    a.plot(dx, dy, 'k.', label = 'Data')
    a.plot(lx, ly, 'g-', label = 'Credible Min.')
    a.plot(hx, hy, 'b-', label = 'Credible Max.')
    
    a.set_title('Credible Range vs Data')
    a.legend(loc = 'upper left')

    matplotlib.pyplot.show()

def plot_order():

    f = open('cubic.order', 'r')
    counts = map(int, f.readlines())
    f.close()

    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(1,1,1)
    
    order = range(len(counts))
    a.bar(order, counts, align='center')
    a.set_xticks(order)

    a.set_title('Order Histogram')
    a.set_xlabel('Order')
    a.set_ylabel('Counts')

    matplotlib.pyplot.show()



