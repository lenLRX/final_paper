import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math
import struct

data = open('../lid_driven_flow/lid_driven_flow/cavity_80000.data','rb')

size = (257,257)

fieldU = np.zeros(size,dtype=float)
fieldV = np.zeros(size,dtype=float)

speedfield = np.zeros(size,dtype=float)

def speed(x,y):
    return speedfield[x][y]


for i in xrange(size[0]):
    for j in xrange(size[1]):
        p = data.read(16)
        if(p == ''):
            print p
            print (i,j)
        pair = struct.unpack('dd',p)
        fieldU[i][j] = pair[0]
        fieldV[i][j] = pair[1]
        speedfield[i][j] = math.sqrt(pair[0] * pair[0] + pair[1] * pair[1])
print data.read(16)
fig = plt.figure()
plt.xlim = (0,256)
plt.ylim = (0,256)

'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''

#plt.quiver(fieldU,fieldV,scale = 3)
plt.contour(speedfield)

plt.show()