import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from matplotlib import animation
from Settings import *


# collect data
Ni = collect('Ni', path=path_a)
Ni_t = collect('Ni_t', path=path_a)
Ni_he = collect('Ni_he', path=path_a)
N_alpha = collect('N_alpha', path=path_a)
Ne = collect('Ne', path=path_a)

# create figure
fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 12)

line1, = ax.plot(x, Ni[:, mid, 0, 0])
line2, = ax.plot(x, Ni_t[:, mid, 0, 0])
line3, = ax.plot(x, Ni_he[:, mid, 0, 0])
line4, = ax.plot(x, N_alpha[:, mid, 0, 0])
line5, = ax.plot(x, Ne[:, mid, 0, 0])

plt.xlabel('$\psi$')
plt.ylabel('density/$10^{19}m^{-3}$')
lengend_line = plt.legend([line1, line2, line3, line4, line5], ['$n_D$', '$n_T$', '$n_{He}$', r'$n_{\alpha}$', '$n_e$'],
                          loc='upper right')
text_pt = plt.text(0.1, 0.8, '', fontsize=16)


# evolve function
def animate(i):
    line1.set_ydata(Ni[:, mid, 0, 10 * i])
    line2.set_ydata(Ni_t[ :, mid, 0,10 * i])
    line3.set_ydata(Ni_he[ :, mid, 0,10 * i])
    line4.set_ydata(N_alpha[ :, mid, 0,10 * i])
    line5.set_ydata(Ne[:, mid, 0,10*i])
    f = i * tbar * 20 * 10
    text_pt.set_text("t=%f s" % f)
    return line1, line2, line3, line4, line5, text_pt,


# init
def init():
    line1.set_ydata(Ni[ :, mid, 0,0])
    line2.set_ydata(Ni_t[ :, mid,0, 0])
    line3.set_ydata(Ni_he[ :, mid, 0,0])
    line4.set_ydata(N_alpha[:, mid, 0,0])
    line5.set_ydata(Ne[:, mid, 0,0])

    return line1, line2, line3, line4, line5,


ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=300,
                              init_func=init,
                              interval=1,
                              blit=False)
# plt.show()
ani.save('Ni_Animation.gif',writer='imagemagick',dpi=300)
