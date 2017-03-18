import math
import numpy as np

e = [[] for i in xrange(19)]
print e

xdim = 160
ydim = 160
zdim = 80
Zco = 64.0

e[0] = [0.0, 0.0, 0.0]
e[1] = [1.0, 0.0, 0.0]
e[2] = [-1.0, 0.0, 0.0]
e[3] = [0.0, 1.0, 0.0]
e[4] = [0.0, -1.0, 0.0]
e[5] = [0.0, 0.0, 1.0]
e[6] = [0.0, 0.0, -1.0]
e[7] = [1.0, 1.0, 0.0]
e[8] = [-1.0, -1.0, 0.0]
e[9] = [1.0, -1.0, 0.0]
e[10] = [-1.0, 1.0, 0.0]
e[11] = [1.0, 0.0, 1.0]
e[12] = [-1.0, 0.0, -1.0]
e[13] = [1.0, 0.0, -1.0]
e[14] = [-1.0, 0.0, 1.0]
e[15] = [0.0, 1.0, 1.0]
e[16] = [0.0, -1.0, -1.0]
e[17] = [0.0, 1.0, -1.0]
e[18] = [0.0, -1.0, 1.0]


def calc_dist(x, y, z, k):
    local_x = (x - xdim / 2)
    local_y = (y - ydim / 2)
    print local_x, local_y, z
    q = 0
    a = abs(e[k][0]) + abs(e[k][1]) - abs(e[k][2]) / Zco
    b = 2 * e[k][0] * local_x + 2 * e[k][1] * local_y - 2.0 / Zco * e[k][2] * z
    print
    c = local_x * local_x + local_y * local_y - z * z / Zco
    print a,b,c
    delta = (b * b - 4 * a * c)
    sqrt_delta = math.sqrt(delta)
    t1 = (-b + sqrt_delta) / (2 * a)
    t2 = (-b - sqrt_delta) / (2 * a)
    print t1, t2

if __name__ == "__main__":
    calc_dist(80, 70, 71, 16)
