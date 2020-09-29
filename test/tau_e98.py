import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from Settings import *


# collect data
tau_e98 = collect('tau_e98',path = path_a)


# create figure
plt.figure()
plt.title(r'$\tau_{e98}$')
plt.plot(t,tau_e98[:])
# plt.plot(t,pellet_timestep[3000:])

#plt.show()
# print(len(n_average))
plt.savefig('tau_e98.png',dpi=300)
