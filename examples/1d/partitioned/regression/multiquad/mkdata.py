import random

xmin = 0.0
xmax = 100.0

def f(x):
    global xmin, xmax

    cx = (xmax + xmin)/2.0
    if (x < cx):
        aa = 0.05
        ab = -0.5
        ac = 35.0

        return aa*x*x + ab*x + ac
    else:
        
        ba = -0.05
        bb = 8.5
        bc = -300.0

        return ba*x*x + bb*x + bc

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

