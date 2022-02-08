import numpy as np
import matplotlib.pyplot as plt

#5#
#a)Random numbers

def uniform_deviate(number, lower, upper):
    n = np.random.uniform(lower, upper, number)
    plt.hist(n, bins='auto')
    
    return n

uniform_deviate(6,1, 46)

#b#
def pdf_deviate(number, lower, upper):
    x = uniform_deviate(number, lower, upper)
    y = np.arccos(1-2*x)
    #for p(y)=1/2 sin(y) => F^-1=cos^-1(1-2x)
    
    
    hist = plt.hist(y, bins='auto')
    plt.show(hist)
    return y
    
#pdf_deviate(100000, 0, 1)

#c#
#pick sin(x) as comparison function

samples = list()
def rejection_method(number, lower, upper):
    #x = uniform_deviate(number, lower, upper)
    x = np.linspace(lower, upper, number)
    p = np.random.uniform(0, np.sin(x),number)
    y = (2/np.pi)* np.sin(x)**2
    
    if len(samples) == number:
        print('samples',samples)
        return samples
    
    elif y < p:
        pass
    
    elif y > p:
        samples.append(x)
        return samples
    
    
    
#samples=rejection_method(10000, 0, np.pi)
#print(samples)