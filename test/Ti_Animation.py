import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from matplotlib import animation
from Settings import *

# collect data
Ti = collect('Ti', path=path_a)*T0/1000
Ti_t = collect('Ti_t', path=path_a)*T0/1000
Ti_he = collect('Ti_he', path=path_a)*T0/1000
Te = collect('Te', path=path_a)*T0/1000

# create figure
fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 35)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

line1, = ax.plot(x, Ti[:, mid, 0, 0])
line2, = ax.plot(x, Ti_t[:, mid, 0, 0])
line3, = ax.plot(x, Ti_he[:, mid, 0, 0])
line4, = ax.plot(x, Te[:, mid, 0, 0])

plt.xlabel('$\psi$',fontsize=10)
plt.ylabel('temperature /keV',fontsize=10)
lengend_line = plt.legend([line1, line2, line3 , line4], ['$n_D$', '$n_T$', '$n_{He}$', r'$n_{\alpha}$', '$n_e$'],
                          loc='upper right',fontsize=10)
text_pt = plt.text(0.1, 0.8, '', fontsize=16)


# evolve function
def animate(i):
    line1.set_ydata(Ti[:, mid, 0, 10 * i])
    line2.set_ydata(Ti_t[ :, mid, 0,10 * i])
    line3.set_ydata(Ti_he[ :, mid, 0,10 * i])
    line4.set_ydata(Te[:, mid, 0,10*i])
    f = i * tbar * 20 * 10
    text_pt.set_text("t=%f s" % f)
    return line1, line2, line3, line4, text_pt,


# init
def init():
    line1.set_ydata(Ti[ :, mid, 0,0])
    line2.set_ydata(Ti_t[ :, mid,0, 0])
    line3.set_ydata(Ti_he[ :, mid, 0,0])
    line4.set_ydata(Te[:, mid, 0,0])

    return line1, line2, line3, line4,


ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=300,
                              init_func=init,
                              interval=1,
                              blit=False)
# plt.show()
ani.save('Ti_Animation.gif',writer='imagemagick',dpi =300)
