import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from matplotlib import animation
from Settings import *


# collect data
S_pellet = collect('S_pellet', path=path_a)

# create figure
fig, ax = plt.subplots()
ax.set_xlim(0, 1)


line1, = ax.plot(x, S_pellet[:, mid, 0, 0],linestyle='--')


plt.xlabel('$\psi$')
plt.ylabel('quantity/$10^{19}$')
text_pt = plt.text(0.1, 0.8, '', fontsize=16)


# evolve function
def animate(i):
    line1.set_ydata(S_pellet[:, mid, 0, 10 * i])
    f = i * tbar * 20 * 10
    text_pt.set_text("t=%f s" % f)
    return line1,  text_pt,


# init
def init():
    line1.set_ydata(S_pellet[ :, mid, 0,0])


    return line1,


ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=300,
                              init_func=init,
                              interval=1,
                              blit=False)
# plt.show()
ani.save('S_pellet_Animation.gif',writer='imagemagick',dpi=300)
