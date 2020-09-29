import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
from Settings import *
from matplotlib import animation

# collect data
Diff_psi = collect('Diff_psi', path = path_a)*Lbar*Lbar/tbar

# create figure
fig, ax = plt.subplots()
ax.set_xlim(0, 1)

line1, = ax.plot(x, Diff_psi[:,mid,0,0])

plt.xlabel('$\psi$')
plt.ylabel('Diffusion coef $m^{2}/s$')
# lengend_line = plt.legend(line1,'Diff',loc='upper right')
text_pt = plt.text(0.1, 0.2, '', fontsize=16)

# evolve function
def animate(i):
    line1.set_ydata(Diff_psi[:,mid,0,10*i])

    f=i*tbar*20*10
    text_pt.set_text("t=%f s" % f)
    return line1,text_pt,

#init
def init():
    line1.set_ydata(Diff_psi[ :, mid, 0,0])


    return line1,

ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=300,
                              init_func=init,
                              interval=1,
                              blit=False)

# plt.show()
ani.save('Diff_Animation.gif',writer='imagemagick',dpi=300)