import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import math
import struct
#Droplet_Dynamics_on_Surface
data = open('../JumpDroplet_example/JumpDroplet_example/cavity_5500.data','rb')

size = (101,101)

x,y=np.mgrid[-50:51,-50:51]

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
        speedfield[i][j] = math.log(pair[0])

fig = plt.figure()

ax = p3.Axes3D(fig)

ax.plot_surface(x,y,speedfield,rstride=3,cstride=3)

'''
fig_quiver = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
fig_contour = fig.add_subplot(2,1,1,xlim = (0,256),ylim = (0,256))
'''
print np.max(speedfield)



#Y, X = np.mgrid[0:257, 0:257]
#plt.streamplot(X,Y,fieldU,fieldV,density=3)

plt.show()