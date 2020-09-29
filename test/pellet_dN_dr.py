import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from matplotlib import animation
from Settings import *

pellet_dN_dr= collect('pellet_dN_dr',path=path_a)
plt.plot(x,pellet_dN_dr[:,mid,0,0]*1e19*200000)

plt.show()