import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
import numpy as np
from Settings import *

# collect data
Ti = collect('Ti',path = path_a)*T0/1000

T_average = np.zeros(3001)
for t in range (0,3001):
    a=0
    b=0
    for px in range(0,nx):
        for py in range(0,ny):
            a = a+Ti[px,py,0,t]*dx[px,py]*dy[px,py]*jacobi[px,py]
            b = b+dx[px,py]*dy[px,py]*jacobi[px,py]
    T_average[t] = a/b

t = [value*tbar*20 for value in range(0,3001)]




# middle plane position y direaction
mid = 46


# pellet_timestep = collect('pellet_timestep',path = path_a)

# create figure
plt.figure(0)
plt.title(r'$Ti_{avearge}$',fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel(r'time/s',fontsize=10)
plt.ylabel(r'temperature /keV',fontsize=10)
plt.plot(t,T_average[:])
# plt.plot(t,pellet_timestep[3000:])
plt.savefig('T_average_t.png',dpi=300)

# plt.show()
T_average_time=np.mean(T_average[-100:])
diff_T_average = np.zeros(len(T_average[-100:]))

for i in range(0,len(T_average[-100:])):
    diff_T_average[i] = T_average[-100+i]-T_average_time
print(diff_T_average)
print(np.max(diff_T_average))
print(np.min(diff_T_average))
print(np.max(diff_T_average)/T_average_time)
print(np.min(diff_T_average)/T_average_time)

T_max = np.zeros(nx)
T_mean = np.zeros(nx)
T_min = np.zeros(nx)
for i in range(0,nx):
    T_max[i] = np.max(Ti[i,mid,0,-100:])
    T_mean[i] = np.mean(Ti[i,mid,0,-100:])
    T_min[i] = np.min(Ti[i,mid,0,-100:])
print(T_max)
print(T_mean)
print(T_min)

T_max_diff = (T_max-T_mean)/T_mean
T_min_diff = (T_min-T_mean)/T_mean

T_max_diff_abs = (T_max-T_mean)
T_min_diff_abs = (T_min-T_mean)

print(x)
plt.figure(1)
line1,=plt.plot(x,T_max)
line2,=plt.plot(x,T_mean)
line3,=plt.plot(x,T_min)
plt.xlabel('$\psi$')
plt.ylabel('temperature/keV$')
lengend_line = plt.legend([line1, line2, line3], ['$T_max$', '$T_mean$', '$T_min$'],
                          loc='upper right')
plt.savefig('T_average.png',dpi=300)

plt.figure(2)
plt.plot(x,T_max_diff)
plt.plot(x,T_min_diff)

plt.figure(3)
plt.plot(x,T_max_diff_abs)
plt.plot(x,T_min_diff_abs)
# print(len(n_average))
# plt.savefig('n_average.png')
# plt.show()
a=np.max(T_max_diff_abs)
b=np.min(T_min_diff_abs)
c=np.max(T_max_diff)
d=np.min(T_min_diff)
e=np.mean(T_average[-100:])

with open('data.txt','a') as file_object:
    file_object.write("positive temperature fluctuation at fuel position is " + str(a) + "\n")
    file_object.write("negative temperature fluctuation at fuel position is " + str(b) + "\n")
    file_object.write("relative temperature fluctuation at fuel position is " + str(c) + "\n")
    file_object.write("relative temperature fluctuation at fuel position is " + str(d) + "\n")
    file_object.write("average volume temperature is" + str(e) + "\n")