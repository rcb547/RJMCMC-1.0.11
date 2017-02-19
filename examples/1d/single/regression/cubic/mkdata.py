import random

def f(x):

    a = 0.00015
    b = -0.005
    c = -0.550
    d = 3.50

    return a*x*x*x + b*x*x + c*x + d

def linspace(minv, maxv, n):
    return map(lambda x: float(x)/float(n - 1) * (maxv - minv) + minv, range(n))

if __name__ == '__main__':

    xmin = 0.0
    xmax = 100.0
    noise = 10.0
    npoints = 100

    x = linspace(xmin, xmax, npoints);
    y = map(f, x)
    yn = map(lambda x: random.normalvariate(x, noise), y)
    n = [noise] * len(yn)

    f = open('data.txt', 'w')
    f.write('\n'.join(map(lambda x: ' '.join(map(str, x)), zip(x, yn, n))))
    f.close()

