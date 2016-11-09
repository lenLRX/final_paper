import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math
import struct
#Droplet_Dynamics_on_Surface
data = open('../SRCone/SRCone/gs-0.2_200.data','rb')

size = (40,50)

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

'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''
print np.max(speedfield)
#plt.quiver(fieldU,fieldV,scale = 3)
C = plt.contour(speedfield,levels = [3.4])
plt.clabel(C)

#Y, X = np.mgrid[0:257, 0:257]
#plt.streamplot(X,Y,fieldU,fieldV,density=3)

plt.show()