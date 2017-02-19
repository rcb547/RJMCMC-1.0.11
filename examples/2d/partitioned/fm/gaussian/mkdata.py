
import math
import random

xmin = -50.0
xmax = 50.0

ymin = -50.0
ymax = 50.0

sigmax = 25.0
sigmay = 15.0

mag = 25.0

def g(x, y):
    #
    # A Gaussian
    #
    global xmin, xmax, ymin, ymax, mag

    R = min([(xmax - xmin), (ymax - ymin)])/4.0
    x0 = (xmax + xmin)/2.0
    y0 = (ymax + ymin)/2.0

    dx = x - x0
    dy = y - y0

    return mag * math.exp(-((dx*dx)/(sigmax*sigmax) + (dy*dy)/(sigmay*sigmay))/2.0)

if __name__ == '__main__':

    noise = 2.0
    nsamplepoints = 500

    xt = [random.uniform(xmin, xmax) for i in range(nsamplepoints)]
    yt = [random.uniform(ymin, ymax) for i in range(nsamplepoints)]

    data = [(x, y, g(x, y) + random.normalvariate(0.0, noise)) for (x, y) in zip(xt, yt)]
    
    f = open('data.txt', 'w')
    f.write('\n'.join(map(lambda x: ' '.join(map(str, x)), data)))
    f.write('\n')
    f.close()

