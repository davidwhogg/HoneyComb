import numpy as np
import matplotlib.pyplot as plt

def displacement(t, x0, dx0dt, n=1/(3.7*3600*24), w0=3e-3, A=.1):
    '''
    From Chaplin et al. (1997).
    x0, dx0/dt = initial displacement and velocity (m, m/s).
    w0 = natural angular frequency (Hz).
    n = damping constant (s-1).
    A = amplitude of the (white) forcing function, d(t-t0).
    Default values copied from the paper.
    '''
    wn2 = w0**2 - n**2
    C = (A + dx0dt + (n * x0)) / wn2**.5
    ent = np.exp(-n * t)
    return C * ent * np.sin(wn2)**.5 * t + x0 * ent * np.cos(wn2)**.5 * t

def velocity(t, x0, dx0dt, n=1/(3.7*3600*24), w0=3e-3, A=.1):
    wn2 = w0**2 - n**2
    C = (A + dx0dt + (n * x0)) / wn2**.5
    ent = np.exp(-n * t)
    return (C * ent * wn2**.5 - (x0 * n * ent)) * np.cos(wn2)**.5 * t \
            - (x0  * ent * wn2**.5 + (C * n * ent)) * np.sin(wn2)**.5 * t

def velocity2Lum(v_osc, Teff):
    '''
    Equation 5 of Kjeldsen & Bedding (1995).
    v_osc in m/s, Teff in K. Returns dL/L (ppm).
    '''
    return v_osc / Teff / 5777 * 17.7

if __name__ == "__main__":
    t = np.linspace(0, 10*24*3600, 1000)
    d = displacement(t)
    v1 = velocity(t)
    v2 = velocity(t, .1, .1)
#     dLL = velocity2Lum(v, 5000)

    plt.clf()
    plt.plot(t/24/3600, v1)
    plt.plot(t/24/3600, v2)
    plt.show()
#     plt.savefig("test")
