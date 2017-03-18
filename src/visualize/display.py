import os
import matplotlib.pyplot as plt

import math

from energy import *

g_x = []

g_x_velo = []

g_pair = []

g_z_pair = []

img_dir = "img"

def visualize_file(fname,radius):

    first = True

    rho_prev = -1

    z_critical = -1

    with open(fname) as fp:
        x_max = -1.0
        z_max = -1.0
        x_arr = []
        x_velo = []
        z_velo = []
        btm_density = []
        for line_no,line in enumerate(fp):
            x_arr.append(line_no)
            ll = line.strip().split("   ")
            z_ = float(ll[1])
            z_velo.append(z_)
            x_velo.append(float(ll[5]))
            if float(ll[5]) > x_max:
                x_max = float(ll[5])
            
            
            btm_density.append(float(ll[4]))

            if float(ll[4])  < 3.5 and rho_prev > 3.5:
                z_critical = z_
            
            rho_prev = float(ll[4])

            '''
            if float(ll[1]) > z_max:
                first = False
                z_max = float(ll[1])
            '''
                
    
    g_x_velo.append(x_max)

    g_pair.append((radius,x_max))

    if(z_critical > 0):
        g_z_pair.append((radius,z_critical))
    
    fig = plt.figure()
    line_x_velo, = plt.plot(x_arr,x_velo,label = "x_velo")
    line_z_velo, = plt.plot(x_arr,z_velo,label = "z_velo")
    plt.legend([line_x_velo,line_z_velo],["x_velocity","z_velocity"])
    plt.title("R = %d velocity"%radius)

    plt.savefig(os.path.join(img_dir,"R%d_velocity.jpg"%radius))

    fig = plt.figure()
    bottom_density, = plt.plot(x_arr,btm_density)
    plt.title("R = %d density"%radius)
    plt.savefig(os.path.join(img_dir,"R%d_density.jpg"%radius))

def visualize_dir(dir_name,radius):
    files = os.listdir(dir_name)
    max_num = 0
    max_file = ""
    
    for fname in files:
        flist = fname.split("_")
        if len(flist) > 2 and int(flist[2]) > max_num:
            max_num = int(flist[2])
            max_file = fname
    
    visualize_file(os.path.join(dir_name,max_file),radius)

def visualize(base_dir):
    base_dir_files = os.listdir(base_dir)
    for file_name in base_dir_files:
        if file_name.startswith("R1"):
            visualize_dir(os.path.join(base_dir,file_name),int(file_name.split("_")[1]))
    
    fig = plt.figure()

    g_pair.sort(key=lambda x:x[0])

    plt.title("x_max")
    plt.plot([x[0] for x in g_pair],[x[1] for x in g_pair])
    plt.savefig("x_max.jpg")

    fig = plt.figure()

    g_z_pair.sort(key=lambda x:x[0])

    plt.title("z_velo m/s")
    plt.plot([x[0] for x in g_z_pair],[get_real_velocity(x[1],0.53,x[0]) for x in g_z_pair])
    plt.plot([x[0] for x in g_z_pair],[calc_velocity(calc_Ek(get_real_radius(30,0.53),get_real_radius(x[0],0.53)),(get_real_radius(30,0.53)**3 + get_real_radius(x[0],0.53)**3) ** (1.0/3.0)) for x in g_z_pair])
    plt.savefig("z_velo.jpg")

if __name__ == "__main__":
    visualize(".\\")
    plt.show()