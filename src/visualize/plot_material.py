import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math
import struct

data = open('../SRCone/SRCone/lattice_attr','rb')

size = (80,160)

fieldU = np.zeros(size,dtype=float)
fieldV = np.zeros(size,dtype=float)

rho = np.zeros(size,dtype=float)

solid = np.ones(size,dtype=float)


for i in xrange(size[0]):
    for j in xrange(size[1]):
        if i == size[0] - 1 or i == 0 or abs(j - size[1]/2)**2 <= i*i /64.0:
            solid[i][j] = 0.0
        p = data.read(4)
        if(p == ''):
            print p
            print (i,j)
        pair = struct.unpack('i',p)
        rho[i][j] = pair[0]


axes = plt.gca()
axes.set_xlim(0,size[1] - 1)
axes.set_ylim(0,size[0] - 1)
#axes.set_xticklabels([])
#axes.set_yticklabels([])


'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''
print np.max(rho)

plt.imshow(rho,cmap=plt.get_cmap("gray"),interpolation="none")

#Y, X = np.mgrid[0:257, 0:257]
#plt.streamplot(X,Y,fieldU,fieldV,density=3)

plt.show()