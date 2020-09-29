from boutpy.boututils.datafile import DataFile
import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
import numpy as np
from matplotlib import animation
from Settings import *

# collect data
P_fusion_total = collect('P_fusion_total',path = path_a)/3.5*17.6


# create figure
plt.figure()
plt.title(r'$P_{fusion}$')
plt.plot(t,P_fusion_total[:])
# plt.plot(t,pellet_timestep[3000:])

# plt.show()
# print(len(n_average))
plt.savefig('fusion_power.png',dpi=300)
with open('data.txt','a') as file_object:
    file_object.write("fusion power is " + str(np.mean(P_fusion_total[-100:])) + "\n")