import kplr
client = kplr.API()
import numpy as np
import matplotlib.pyplot as plt
import astero
data = astero.astero()

kids = data.iKID
for kid in kids:
    print str(int(kid))
    star = client.star(str(int(kid)))
    star.get_light_curves(fetch=True)
