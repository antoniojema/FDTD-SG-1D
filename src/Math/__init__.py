import numpy as np

def gaussian(x, mean, sigma, amplitude=1):
    return amplitude * np.exp(-(x-mean)*(x-mean)/(2*sigma*sigma))

def gaussian_derivative(x, mean, sigma, amplitude=1):
    return (mean-x)*gaussian(x, mean, sigma, amplitude)
