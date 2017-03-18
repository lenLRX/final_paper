import math

def get_real_radius(r,tau):
    coef2 = ((tau-0.5)/3*6.698)/(math.sqrt(6.698*0.18299))
    coef1 = 1.0087 * 10 ** -3 / math.sqrt(998.23 * 72.75 * 10 ** -3)
    return r * (coef1 / coef2) ** 2

sigma_real = 72.75 * 10 ** -3
miu_real = 1.0087 * 10 ** -3
rho_real = 998.23

def get_dimensionless_velocity(speed,tau,R2):
    R_new = (30**3 + R2**3) ** (1.0/3.0)
    U = math.sqrt(sigma_real/(rho_real * get_real_radius(R_new,tau)))
    return speed /U

def get_real_velocity(speed,tau,R2):
    Ur = 332.532 * math.sqrt(3)

    R_new = (30**3 + R2**3) ** (1.0/3.0)

    U = math.sqrt(sigma_real/(rho_real * get_real_radius(R_new,tau)))

    real_speed = U * speed
    return real_speed





def calc_Ek(R1,R2):
    Vsum_old = 4.0/3.0 * math.pi * (R1**3 + R2**3)
    Asum_old = 4 * math.pi * (R1**2 + R2**2)
    R_new = (R1**3 + R2**3) ** (1.0/3.0)
    A_new =  4 * math.pi * R_new ** 2
    dAlv = Asum_old - A_new 
    dEs = sigma_real * dAlv

    Evis = 36 * math.pi * miu_real * (math.sqrt(sigma_real * R1 ** 3/rho_real) + math.sqrt(sigma_real * R2 ** 3/rho_real))

    Evis = 0

    Eh = R_new * Vsum_old * rho_real - rho_real * R1 * 4.0/3.0 * math.pi * R1**3\
    - rho_real * R2 * 4.0/3.0 * math.pi * R2**3

    R_btm = math.sin(20.0/180.0 * math.pi) * R_new

    Asl = math.pi * R_btm ** 2

    Ew = sigma_real * (1 + math.cos(160.0/180.0 * math.pi)) * Asl

    

    Ek = dEs - Evis - Eh - Ew
    #print Ew,Ek,Eh,Evis,dEs
    return Ek

def calc_velocity(Ek,R):
    V = 4.0/3.0 * math.pi * R ** 3
    M = V * rho_real
    if(Ek < 0):
        return 0
    velocity = math.sqrt(2 * Ek / M)
    print velocity
    return velocity

if __name__ == "__main__":
    tau = 0.57
    R1_real = get_real_radius(30,tau)
    print R1_real
    for i in range(1,31):
        R2_real = get_real_radius(i,tau)

        R_new_real = (R1_real**3 + R2_real**3) ** (1.0/3.0)

        Ek = calc_Ek(R1_real,R2_real)
        print i,Ek

        #calc_velocity(Ek,R_new_real)