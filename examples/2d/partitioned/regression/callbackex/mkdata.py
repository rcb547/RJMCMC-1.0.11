
import math
import random

xmin = -50.0
xmax = 50.0

ymin = -50.0
ymax = 50.0

def g(x, y):
    #
    # A 2D Gaussian function
    #
    global xmin, xmax, ymin, ymax

    x0 = (xmin + xmax)/2.0;
    y0 = (ymin + ymax)/2.0;
    theta = math.pi/6;

    sigma_x = (xmax - xmin)/8
    sigma_x2 = sigma_x * sigma_x
    sigma_y = (ymax - ymin)/4
    sigma_y2 = sigma_y*sigma_y

    ctheta = math.cos(theta)
    stheta = math.sin(theta)

    a = (ctheta*ctheta)/(2.0*sigma_x2) + (stheta*stheta)/(2.0*sigma_y2)
    b = -math.sin(2.0*theta)/(4.0*sigma_x2) + math.sin(2.0*theta)/(4.0*sigma_y2)
    c = (stheta*stheta)/(2.0*sigma_x2) + (ctheta*ctheta)/(2.0*sigma_y2)
    A = 50.0;
    
    dx = x - x0
    dy = y - y0
    return A * math.exp(-(a*dx*dx + 2*b*dx*dy + c*dy*dy));

if __name__ == '__main__':

    noise = 2.0
    nsamplepoints = 500

    xt = [random.uniform(xmin, xmax) for i in range(nsamplepoints)]
    yt = [random.uniform(ymin, ymax) for i in range(nsamplepoints)]

    data = [(x, y, g(x, y) + random.normalvariate(0.0, noise)) for (x, y) in zip(xt, yt)]
    
    f = open('data.txt', 'w')
    f.write('\n'.join(map(lambda x: ' '.join(map(str, x)), data)))
    f.close()
