import matplotlib.pylab as plt
from boutpy.boutdata.collect import collect
import numpy as np
from Settings import *

# collect data
n_average = collect('n_average',path = path_a)
Ne = collect('Ne',path = path_a)

# create figure
plt.figure(0)
plt.title(r'$Ne_{avearge}$',fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel(r'time/s',fontsize=10)
plt.ylabel(r'electron average density $10^{19} m^{3}$',fontsize=10)
plt.plot(t,n_average[:])
plt.savefig('n_average_t.png',dpi=300)

n_average_time=np.mean(n_average[-100:])
diff_n_average = np.zeros(len(n_average[-100:]))

for i in range(0,len(n_average[-100:])):
    diff_n_average[i] = n_average[-100+i]-n_average_time
print(diff_n_average)
print(np.max(diff_n_average))
print(np.min(diff_n_average))
print(np.max(diff_n_average)/n_average_time)
print(np.min(diff_n_average)/n_average_time)

n_max = np.zeros(nx)
n_mean = np.zeros(nx)
n_min = np.zeros(nx)
for i in range(0,nx):
    n_max[i] = np.max(Ne[i,mid,0,-100:])
    n_mean[i] = np.mean(Ne[i,mid,0,-100:])
    n_min[i] = np.min(Ne[i,mid,0,-100:])
print(n_max)
print(n_mean)
print(n_min)

n_max_diff = (n_max-n_mean)/n_mean
n_min_diff = (n_min-n_mean)/n_mean

n_max_diff_abs = (n_max-n_mean)
n_min_diff_abs = (n_min-n_mean)

print(x)
plt.figure(1)
line1,=plt.plot(x,n_max)
line2,=plt.plot(x,n_mean)
line3,=plt.plot(x,n_min)
plt.xlabel('$\psi$')
plt.ylabel('density/$10^{19}m^{-3}$')
lengend_line = plt.legend([line1, line2, line3], ['$n_max$', '$n_mean$', '$n_min$'],
                          loc='upper right')
plt.savefig('n_average.png',dpi=300)

plt.figure(2)
plt.plot(x,n_max_diff)
plt.plot(x,n_min_diff)

plt.figure(3)
plt.plot(x,n_max_diff_abs)
plt.plot(x,n_min_diff_abs)
# print(len(n_average))
# plt.savefig('n_average.png',dpi=300)
# plt.show()
a=np.max(n_max_diff_abs)
b=np.min(n_min_diff_abs)
c=np.max(n_max_diff)
d=np.min(n_min_diff)
e=np.mean(np.mean(n_average[-100:]))

with open('data.txt','a') as file_object:
    file_object.write("positive density fluctuation at fuel position is " + str(a) + "\n")
    file_object.write("negative density fluctuation at fuel position is " + str(b) + "\n")
    file_object.write("relative positive density fluctuation at fuel position is " + str(c) + "\n")
    file_object.write("relative negative density fluctuation at fuel position is " + str(d) + "\n")
    file_object.write("volume average density is  " + str(e) + "\n")