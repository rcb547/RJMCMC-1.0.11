
import os
import sys

import matplotlib
import matplotlib.pyplot

if __name__ == '__main__':

    f = open('natural.misfit', 'r')
    misfit = map(float, f.readlines())
    f.close()
    
    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(111)

    a.plot(misfit)
    a.set_title('Misfit History')
    a.set_xlabel('Iteration')
    a.set_ylabel('Misfit')
                
    matplotlib.pyplot.show()
        
