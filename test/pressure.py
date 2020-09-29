import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from Settings import *

# collect data
pressure = collect('pressure',path = path_a)
Ni = collect('Ni',path = path_a)

# create figure
plt.figure()
plt.xticks(fontsize=10)
plt.yticks([65,70,75,80,80.99],[r'$65$',r'$70$',r'$75$',r'$80$',r'$80.99$'],fontsize=10)

plt.title(r'time evolution of pressure at $\psi$=0.94',fontsize=10)
plt.plot(t,pressure[63,46,0,:])
plt.axhline(y=80.99,ls=":",c="b") #添加水平直线
plt.xlabel('$time /s$',fontsize=10)
plt.ylabel('$pressure /kPa$',fontsize=10)
# plt.plot(t,pellet_timestep[3000:])

# print(x[63])
# print(x[64])
# plt.show()
# print(len(n_average))
plt.savefig('pressure.png',dpi=300)
