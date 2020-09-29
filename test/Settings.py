from boutpy.boututils.datafile import DataFile
from boutpy.boutdata.collect import collect
import matplotlib.pylab as plt
import numpy as np

# grid file
g = DataFile('../CFETRHybridit3530_rechange.nc')
psi = (g.read('psixy') - g.read('psi_axis'))/(g.read('psi_bndry')-g.read('psi_axis'))

# Normailization for CFETRHybridit3530
tbar = 4.317739e-04
Lbar = 9.434243
Ni0 = 1e19
T0 = 10
x = psi[:,46]
PI=3.14159
rp = 0.2
np_pellet = 5.96e22
dx = g.read('dx')
dy = g.read('dy')
Bpxy = g.read('Bpxy')
hthe = g.read('hthe')
jacobi = hthe/Bpxy
nx = g.read('nx')
ny = g.read('ny')

t = [value*tbar*20 for value in range(0,3001)]

# specify the path
path_a =\
    '/global/homes/p/pkulinws/bout/examples/burning_plasma_pellet_ELM_for_paper_add_displacement'\
    '/data'


# middle plane position y direaction
mid = 46

r_cloud = collect('r_cloud', path=path_a)
T_pellet=collect('T_pellet', path=path_a)
pellet_rp_profile=collect('pellet_rp_profile',path=path_a)
L_c=collect('L_c',path=path_a)
pellet_displacement=collect('pellet_displacement',path=path_a)
add_pellet_displacement=collect('add_pellet_displacement',path=path_a)
pellet_dN_dr=collect('pellet_dN_dr',path=path_a)
R_new=collect('R_new',path=path_a)
print(r_cloud[:,mid,0,0])
print(L_c[:,mid,0,0])
print(pellet_displacement[:,mid,0,0])
print(x-pellet_displacement[:,mid,0,0])
print(R_new)
plt.figure()
# plt.plot(x,g.read('Rxy')[:,mid])
plt.plot(x-add_pellet_displacement[:,mid,0,0],R_new[:,mid,0,0])
# plt.plot(x,pellet_displacement[:,mid,0,0])
# plt.plot(x,add_pellet_displacement[:,mid,0,0])

plt.show()
# print(g.read('psi_bndry'))
# print(g.read('psi_axis'))
# print(g.list())
# print(T_pellet[:,mid,0,0])
# print(pellet_rp_profile[:,mid,0,0])


