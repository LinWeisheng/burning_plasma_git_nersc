import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
import numpy as np
from Settings import *

# collect data
n_average = collect('n_average',path = path_a)
Ne = collect('Ne',path = path_a)

V=0

for px in range(0,nx):
    for py in range(0,ny):
        V = V+dx[px,py]*dy[px,py]*jacobi[px,py]
print(V*2*3.14)
print('total particle is'+str(V*2.314*n_average[0]*1e19))





