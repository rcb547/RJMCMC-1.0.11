
import math
import random

xmin = -50.0
xmax = 50.0

ymin = -50.0
ymax = 50.0

disc_value = 25.0

def g(x, y):
    #
    # A Square
    #
    global xmin, xmax, ymin, ymax, disc_value

    R = min([(xmax - xmin), (ymax - ymin)])/4.0
    x0 = (xmax + xmin)/2.0
    y0 = (ymax + ymin)/2.0

    dx = math.fabs(x - x0)
    dy = math.fabs(y - y0)

    if (dx < R) and (dy < R):
        return disc_value
    else:
        return 0.0

if __name__ == '__main__':

    noise = 2.0
    nsamplepoints = 500

    xt = [random.uniform(xmin, xmax) for i in range(nsamplepoints)]
    yt = [random.uniform(ymin, ymax) for i in range(nsamplepoints)]

    data = [(x, y, g(x, y) + random.normalvariate(0.0, noise)) for (x, y) in zip(xt, yt)]
    
    f = open('data.txt', 'w')
    f.write('\n'.join(map(lambda x: ' '.join(map(str, x)), data)))
    f.close()

