import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math
import struct
import sys

#data = open('../SRCone/SRCone/speed_field-1_10.data','rb')
data = open(sys.argv[1],'rb')

size = (80,160)

fieldU = np.zeros(size,dtype=float)
fieldV = np.zeros(size,dtype=float)

solid = np.ones(size,dtype=float)

speedfield = np.zeros(size,dtype=float)

def speed(x,y):
    return speedfield[x][y]

p_max = 0

for i in xrange(size[0]):
    for j in xrange(size[1]):
        if i == size[0] - 1 or i == 0 or abs(j - size[1]/2)**2 <= i*i /64.0:
            solid[i][j] = 0.0
        p = data.read(16)
        if(p == ''):
            print p
            print (i,j)
        pair = struct.unpack('dd',p)
        fieldU[i][j] = pair[0]
        fieldV[i][j] = pair[1]
        if abs(pair[1] > p_max):
            p_max = pair[1]
        speedfield[i][j] = math.sqrt(pair[0] * pair[0] + pair[1] * pair[1])

print p_max

plt.imshow(solid,cmap=plt.get_cmap("gray"),interpolation="none")

'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''

#plt.quiver(fieldU,fieldV,scale = 3)
#plt.contour(speedfield)

Y, X = np.mgrid[0:size[0], 0:size[1]]
plt.streamplot(X,Y,fieldU,fieldV)

plt.show()