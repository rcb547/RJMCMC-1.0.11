import random

def real_function(x):

    if x < 5.0:
        return x
    else:
        return 10.0 - x

xmin = 0.0
xmax = 10.0
npoints = 100
x = map(lambda x: float(x)/float(npoints - 1) * (xmax - xmin), range(npoints))
y = map(real_function, x)

f = open('realdata.txt', 'w')
for (x, y) in zip(x, y):
    f.write('%f %f\n' % (x, y))
f.close()

npoints = 20
xr = map(lambda x: float(x)/float(npoints - 1) * (xmax - xmin), range(npoints))
yr = map(real_function, xr)

sigma = 1.0

y = map(lambda x: x + random.normalvariate(0.0, sigma), yr)

f = open('data.txt', 'w')
for (x, y) in zip(xr, y):
    f.write('%f %f\n' % (x, y))
f.close()



