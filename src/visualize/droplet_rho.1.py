import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math
import struct
#Droplet_Dynamics_on_Surface
data = open('../droplet_example/droplet_example/cavity_100.data','rb')

size = (101,101)

fieldU = np.zeros(size,dtype=float)
fieldV = np.zeros(size,dtype=float)

speedfield = np.zeros(size,dtype=float)

def speed(x,y):
    return speedfield[x][y]


for i in xrange(size[0]):
    for j in xrange(size[1]):
        p = data.read(8)
        if(p == ''):
            print p
            print (i,j)
        pair = struct.unpack('d',p)
        speedfield[i][j] = pair[0]

fig = plt.figure()
plt.xlim = (0,80)
plt.ylim = (0,80)

print np.max(speedfield)
'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''

#plt.quiver(fieldU,fieldV,scale = 3)
#C = plt.contour(speedfield,[0.4,0.8,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.5])
C = plt.contour(speedfield)
plt.clabel(C)

#Y, X = np.mgrid[0:257, 0:257]
#plt.streamplot(X,Y,fieldU,fieldV,density=3)

plt.show()