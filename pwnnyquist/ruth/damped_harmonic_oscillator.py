import numpy as np
import matplotlib.pyplot as plt

def velocity(t, n, w0, A, x0, dx0dt):
    wn2 = w0**2 - n**2
    C = (A + dx0dt + n*x0) / wn2**.5
    ent = np.exp(-n*t)
    return ( C*ent*wn2**.5 - x0*n*ent ) * np.cos(wn2**.5*t) \
            - ( x0*ent*wn2**.5 + C*n*ent ) * np.sin(wn2**.5*t)

def position(t, n, w0, A, x0, dx0dt):
    wn2 = w0**2 - n**2
    C = (A + dx0dt + n*x0) / wn2**.5
    ent = np.exp(-n*t)
    return C*ent*np.sin(wn2**.5*t) + x0*ent*np.cos(wn2**.5*t)

def general_solution(t, w0, n, A):
    return A * np.exp(-n*t/2) * np.cos(w0*t)

if __name__ == "__main__":

    t = np.arange(0, 30, .1)
    dx0dt = 0.
    w0 = 3e-3 * 24 * 3600
    print w0
    tau = 3.7
    n = 1./tau
    n = .02
    A = 1.
    x0 = 0.

    x = position(t, n, w0, A, x0, dx0dt)
    dxdt = velocity(t, n, w0, A, x0, dx0dt)

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(t, x)
    plt.subplot(2, 1, 2)
    plt.plot(t, dxdt)
    plt.show()
