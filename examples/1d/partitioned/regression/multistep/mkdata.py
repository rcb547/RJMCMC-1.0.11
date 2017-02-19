import random

xmin = 0.0
xmax = 100.0

def f(x):
    global xmin, xmax

    cx1 = xmin + (xmax - xmin)/3.0
    cx2 = xmin + 2.0*(xmax - xmin)/3.0
    if (x < cx1):
        return 0.0
    elif (x < cx2):
        b = 10.0
        return 10.0
    else:
        return -15.0

def linspace(minv, maxv, n):
    return map(lambda x: float(x)/float(n-1) * (maxv - minv) + minv, range(n))

if __name__ == '__main__':
    
    npoints = 100
    noise = 5.0

    x = linspace(xmin, xmax, npoints);
    y = map(f, x)
    yn = map(lambda x: random.normalvariate(x, noise), y)
    n = [noise] * len(yn)

    f = open('data.txt', 'w')
    f.write('\n'.join(map(lambda x: ' '.join(map(str, x)), zip(x, yn, n))))
    f.close()

