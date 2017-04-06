import os
import matplotlib.pyplot as plt

import math

from energy import *

#tau_list = [0.53,0.57,0.61,0.65,0.69,0.75,0.78,0.8,0.83,0.85,1,1.2,1.35,1.45]#
tau_list = [0.85]

tau_velo = {}

#g_pair = []

#g_z_pair = []

img_dir = "img"

def visualize_file(base_dir,fname,radius,tau):

    first = True

    rho_prev = -1

    z_critical = -1
    x_critical = -1

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
                x_critical = float(ll[5])
            
            rho_prev = float(ll[4])

            '''
            if float(ll[1]) > z_max:
                first = False
                z_max = float(ll[1])
            '''


    #g_pair.append((radius,x_max))

    if(z_critical > 0):
        tau_velo[str(tau)].append((radius,z_critical,x_critical))
    
    
    fig = plt.figure()
    line_x_velo, = plt.plot(x_arr,x_velo,label = "x_velo")
    line_z_velo, = plt.plot(x_arr,z_velo,label = "z_velo")
    plt.legend([line_x_velo,line_z_velo],["x_velocity","z_velocity"])
    plt.title("R = %d velocity"%radius)

    plt.savefig(os.path.join(os.path.join(base_dir,img_dir),"R%d_velocity.jpg"%radius))

    fig = plt.figure()
    bottom_density, = plt.plot(x_arr,btm_density)
    plt.title("R = %d density"%radius)
    plt.savefig(os.path.join(os.path.join(base_dir,img_dir),"R%d_density.jpg"%radius))
    

def visualize_dir(base_dir,dir_name,radius,tau):
    files = os.listdir(dir_name)
    max_num = 0
    max_file = ""
    
    for fname in files:
        flist = fname.split("_")
        #print flist
        if len(flist) > 2 and int(flist[2]) > max_num:
            max_num = int(flist[2])
            max_file = fname
    
    visualize_file(base_dir,os.path.join(dir_name,max_file),radius,tau)

def visualize(base_dir,tau):
    try:
        os.mkdir(os.path.join(base_dir,img_dir))
    except:
        pass
    base_dir_files = os.listdir(base_dir)
    for file_name in base_dir_files:
        if file_name.startswith("R1"):
            visualize_dir(base_dir,os.path.join(base_dir,file_name),int(file_name.split("_")[1]),tau)
    
    '''
    fig = plt.figure()

    g_pair.sort(key=lambda x:x[0])

    plt.title("x_max")
    plt.plot([x[0] for x in g_pair],[x[1] for x in g_pair])
    plt.savefig(os.path.join(os.path.join(base_dir,img_dir),"x_max.jpg"))
    '''

    fig = plt.figure()

    tau_velo[str(tau)].sort(key=lambda x:x[0])

    #print tau, tau_velo

    plt.title("velocity m/s")
    z_velo, = plt.plot([get_real_radius(x[0],tau)*10**6 for x in tau_velo[str(tau)]],[get_real_velocity(x[1],tau,x[0]) for x in tau_velo[str(tau)]])
    x_velo, = plt.plot([get_real_radius(x[0],tau)*10**6 for x in tau_velo[str(tau)]],[get_real_velocity(x[2],tau,x[0]) for x in tau_velo[str(tau)]])
    theoretical_velo, = plt.plot([get_real_radius(x[0],tau)*10**6 for x in tau_velo[str(tau)]],[calc_velocity(calc_Ek(get_real_radius(30,tau),get_real_radius(x[0],tau)),(get_real_radius(30,tau)**3 + get_real_radius(x[0],tau)**3) ** (1.0/3.0)) for x in tau_velo[str(tau)]])
    plt.legend([z_velo,x_velo,theoretical_velo],["z_velo","x_velo","theoretical_velo"],loc= 0)
    plt.savefig(os.path.join(os.path.join(base_dir,img_dir),"velocity.jpg"))

def summary():
    fig = plt.figure()
    plt.title("dimensionless speed vs ratio")

    legend_handles = []
    legend_names = []

    for tau in tau_list:
        z_velo, = plt.plot([x[0] / 30.0 for x in tau_velo[str(tau)]],[x[1] for x in tau_velo[str(tau)]])
        x_velo, = plt.plot([x[0] / 30.0 for x in tau_velo[str(tau)]],[x[2] for x in tau_velo[str(tau)]])
        Ek_R = [(calc_Ek(get_real_radius(30,tau),get_real_radius(x[0],tau)),(get_real_radius(30,tau)**3 + get_real_radius(x[0],tau)**3) ** (1.0/3.0)) for x in tau_velo[str(tau)]]
        #theoretical_velo, = plt.plot([x[0] / 30.0 for x in tau_velo[str(tau)]],[get_dimensionless_velocity(calc_velocity(*x),tau,x[0]) for x in Ek_R])
        legend_handles.append(z_velo)
        legend_handles.append(x_velo)
        #legend_handles.append(theoretical_velo)

        legend_names.append("z_velo tau = %s"%str(tau))
        legend_names.append("x_velo tau = %s"%str(tau))
        #legend_names.append("theoretical_velo tau = %s"%str(tau))
    
    #plt.legend(legend_handles,legend_names,loc = 0)
    plt.savefig("dimesionless.jpg")


def visualize_all(tau_list):
    for tau in tau_list:
        tau_velo[str(tau)] = []
        #print tau_velo
        visualize(os.path.join(".","varR.tau" + str(tau)),tau)
    summary()

if __name__ == "__main__":
    visualize_all(tau_list)
    #visualize(".\\")
    #plt.show()