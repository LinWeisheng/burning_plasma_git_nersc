from boutpy.boutdata.collect import collect
import numpy as np
from Settings import *

# collect data
pellet_timestep = collect('pellet_timestep',path = path_a)
P_fu_total = collect('P_fusion_total',path = path_a)
S_fu_total = P_fu_total/3.5/1.602

# caculate burn fraction
Burn_particle = 0
time = pellet_timestep[-1]*10/20
for i in range(0, int(time)):
    Burn_particle += S_fu_total[-i]*tbar*20*1e19
Total_particle = 4/3*np_pellet*PI*(rp**3)/2*10
burn_fraction = Burn_particle/Total_particle


print('Burn fraction is', burn_fraction)
with open('data.txt','a') as file_object:
    file_object.write("burn fraction is " + str(burn_fraction) + "\n")
