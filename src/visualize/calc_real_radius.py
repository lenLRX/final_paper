import math

def get_real_radius(r,tau):
    coef2 = ((tau-0.5)/3*6.698)/(math.sqrt(6.698*0.18299))
    coef1 = 1.0087 * 10 ** -3 / math.sqrt(998.23 * 72.75 * 10 ** -3)
    return r * (coef1 / coef2) ** 2