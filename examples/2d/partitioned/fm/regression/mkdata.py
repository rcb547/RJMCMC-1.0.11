
import math
import random

def linspace(minv, maxv, n):
    return map(lambda x: float(x)/float(n - 1)*(maxv - minv) + minv, range(n))

def g(x, y):
    r = math.sqrt(x*x + y*y)
    if (r < 25.0):
        return 1.0
    else:
        return 0.0

if __name__ == '__main__':

    xmin = -50.0
    xmax = 50.0
    ymin = -50.0
    ymax = 50.0

    npoints = 51
    
    disc_value = 25.0
    sigma = 5.0

    x = linspace(xmin, xmax, npoints)
    y = linspace(ymin, ymax, npoints)

    f = open('data.txt', 'w')
    for yp in y:
        for xp in x:
            f.write('%f %f %f\n' % (xp, yp, disc_value * g(xp, yp) + random.normalvariate(0.0, sigma)))
    f.close()
    
