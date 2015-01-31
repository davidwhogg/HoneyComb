import numpy as np
import superpgram

starts = np.arange(0., 10.1, 1.)
stops = np.arange(0.99, 11.2, 1.)
print starts.shape, stops.shape
data = np.random.normal(size=starts.shape)
ivars = np.ones_like(starts)
wavenumbers = np.arange(10., 15.1, 0.5)

print superpgram.superpgram(starts, stops, data, ivars, wavenumbers)
