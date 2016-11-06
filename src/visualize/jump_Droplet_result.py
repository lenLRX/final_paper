import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math
import struct
#Droplet_Dynamics_on_Surface
data = open('../JumpDroplet_example/JumpDroplet_example/result_5600.dat','r')
dummy = data.readline()
dummy = data.readline()
dummy = data.readline()
size = (101,101,101)

fieldU = np.zeros(size,dtype=float)
fieldV = np.zeros(size,dtype=float)

speedfield = np.zeros(size,dtype=float)

def speed(x,y):
    return speedfield[x][y]


for i in xrange(size[0]):
    for j in xrange(size[1]):
        for k in xrange(size[2]):
            l = data.readline().strip()
            split_l = l.split("   ")
            rho = split_l[-2]
            speedfield[i][j][k] = rho

fig = plt.figure()

'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''
print np.max(speedfield)
print np.min(speedfield)
print speedfield
#plt.quiver(fieldU,fieldV,scale = 3)
C = plt.contour(speedfield[50],10)
plt.clabel(C)

#Y, X = np.mgrid[0:257, 0:257]
#plt.streamplot(X,Y,fieldU,fieldV,density=3)

plt.show()