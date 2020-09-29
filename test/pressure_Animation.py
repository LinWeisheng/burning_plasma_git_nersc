import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from matplotlib import animation
from Settings import *

# collect data
Ni = collect('Ni', path=path_a)
Ni_t = collect('Ni_t', path=path_a)
Ni_he = collect('Ni_he', path=path_a)
Ne = collect('Ne', path=path_a)
Ti = collect('Ti', path=path_a)
Ti_t = collect('Ti_t', path=path_a)
Ti_he = collect('Ti_he', path=path_a)
Te = collect('Te', path=path_a)
D_pressure = Ni*Ti*T0*1.602/1000
T_pressure = Ni_t*Ti_t*T0*1.602/1000
He_pressure = Ni_he*Ti_he*T0*1.602/1000
e_pressure = Ne*Te*T0*1.602/1000
alpha_pressure = collect('alpha_pressure', path=path_a)
pressure = collect('pressure', path=path_a)

# create figure
fig, ax = plt.subplots()
ax.set_xlim(0, 1)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

line1, = ax.plot(x, D_pressure[:, mid, 0, 0])
line2, = ax.plot(x, T_pressure[:, mid, 0, 0])
line3, = ax.plot(x, He_pressure[:, mid, 0, 0])
line4, = ax.plot(x, e_pressure[:, mid, 0, 0])
line5, = ax.plot(x, alpha_pressure[:, mid, 0, 0])
line6, = ax.plot(x, pressure[:, mid, 0, 0])

plt.xlabel('$\psi$',fontsize=10)
plt.ylabel('$pressure /kPa$',fontsize=10)
plt.title('pressure of different species',fontsize=10)
lengend_line = plt.legend([line1, line2, line3, line4, line5, line6], ['$Deuterium pressure$', '$Tritium pressure$', 
                                                                       '$Helium pressure$', '$electron pressure$', 
                                                                       r'$\alpha$ pressure', 'total pressure'],
                          loc='upper right', fontsize = 10)
text_pt = plt.text(0.6, 500, '', fontsize=16)


# evolve function
def animate(i):
    line1.set_ydata(D_pressure[:, mid, 0, 10*i])
    line2.set_ydata(T_pressure[:, mid, 0, 10*i])
    line3.set_ydata(He_pressure[:, mid, 0, 10*i])
    line4.set_ydata(e_pressure[:, mid, 0, 10*i])
    line5.set_ydata(alpha_pressure[:, mid, 0, 10*i])
    line6.set_ydata(pressure[:, mid, 0, 10*i])

    f = i * tbar * 20 * 10
    text_pt.set_text("t=%f s" % f)
    return line1, line2, line3, line4, line5, line6, text_pt,


# init
def init():
    line1.set_ydata(D_pressure[:, mid, 0, 0])
    line2.set_ydata(T_pressure[:, mid, 0, 0])
    line3.set_ydata(He_pressure[:, mid, 0, 0])
    line4.set_ydata(e_pressure[:, mid, 0, 0])
    line5.set_ydata(alpha_pressure[:, mid, 0, 0])
    line6.set_ydata(pressure[:, mid, 0, 0])

    return line1, line2, line3, line4, line5, line6


ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=300,
                              init_func=init,
                              interval=1,
                              blit=False)
# plt.show()
ani.save('pressure.gif',writer='imagemagick',dpi=300)
