
import os
import sys

import matplotlib
import matplotlib.pyplot

if __name__ == '__main__':

    f = open('history.txt', 'r')
    history = map(float, f.readlines())
    f.close()
    
    fig = matplotlib.pyplot.figure(1)

    a = fig.add_subplot(2,1,1)

    a.plot(history)
    a.set_title('Value History')
    a.set_xlabel('Iteration')
    a.set_ylabel('Value')

    a = fig.add_subplot(2,1,2)

    a.hist(history)
    a.set_title('Value Histogram')
    a.set_xlabel('Value')
    a.set_ylabel('Count')
                
    matplotlib.pyplot.show()
        
