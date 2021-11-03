import numpy as np

def FT_one(t0, delta_t, values, w):
    values_calc = 0.5*(values[1:]+values[:-1])
    ti_calc = np.arange(t0, t0+delta_t*len(values_calc)+0.5*delta_t, delta_t)
    
    return delta_t * np.sum(values * np.exp(-2*np.pi*1j*w*(ti_calc+0.5*delta_t)))


def FT(t0, delta_t, values, w):
    def f(freq):
        return FT_one(t0, delta_t, values, freq)
    
    return list(map(f, w))


def FT_one_inhomogeneous(t, values, w):
    delta_t = t[1:]-t[:-1]
    t_calc = 0.5*(t[1:]+t[:-1])
    values_calc = 0.5*(values[1:]+values[:-1])
    
    return np.sum(delta_t * values_calc * np.exp(-2*np.pi*1j*w*t_calc))


def FT_inhomogeneous(t, values, w):
    def f(freq):
        return FT_one_inhomogeneous(t, values, freq)
    
    return list(map(f, w))


def gaussian(t, amplitude, mean, sigma):
    return amplitude * np.exp(-(t-mean)**2 / (2*sigma**2))


def ppw_to_sigma(ppw, delta, c0=299792458, dbDecay=3):
    return np.sqrt(dbDecay/40*np.log(10)) * delta*ppw/(np.pi*c0)


def freq_to_sigma(freq, dbDecay=3):
    return np.sqrt(dbDecay/40*np.log(10)) * 1/(np.pi*freq)


def freq_to_ppw(freq, delta=0.1, c=299792458):
    return c/(delta*freq)


def ppw_to_freq(ppw, delta=0.1, c=299792458):
    return c/(delta*ppw)
