import os
import sys
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd)
from godunov import mega_foo as mega_foo_cython
from godunov import get_ray_URP as get_ray_URP_cython
from godunov import border_wall_URPDD as border_wall_URP_cython
import numpy as np
from riman_python import get_F_13_15, get_ray_URP, mega_foo, border_wall_URP



def Godunov_fluxes_python(p1, ro1, u1, c1, p2, ro2, u2, c2, vbi, gamma, e_foo, **kwargs):
    suc, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P = \
        mega_foo(p_1=p1, ro_1=ro1, u_1=u1, c_1=c1, 
                 p_2=p2, ro_2=ro2, u_2=u2, c_2=c2, p_0=0, gamma=gamma, eps_F=1e-5, n_iter_max=100)
    
    if not suc:
        return 0,0,0
    rayU,rayR,rayP = get_ray_URP(vbi, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P,
                         p_1=p1, ro_1=ro1, u_1=u1, c_1=c1, p_2=p2, ro_2=ro2, u_2=u2, c_2=c2, gamma=gamma)
    
    rayE = e_foo(p=rayP, ro=rayR)
    # формула 14.4
    M_j = rayR * (rayU - vbi)
    J_j = rayP + M_j * rayU
    E_j = (rayE + 0.5*rayU**2)*M_j + rayP*rayU

    return M_j, J_j, E_j

def flux_wall_border(left_border, p, ro, u, c, vbi, gamma, e_foo, **kwargs):
    rayU,rayR,rayP = \
        border_wall_URP(
            p=p, ro=ro, u=u, c=c, 
            vbi=vbi, left_border=left_border, 
            p_0=0, gamma=gamma, eps_F=1e-5, n_iter_max=100)
    
    rayE = e_foo(p=rayP, ro=rayR)
    # формула 14.4
    M_j = rayR * (rayU - vbi)
    J_j = rayP + M_j * rayU
    E_j = (rayE + 0.5*rayU**2)*M_j + rayP*rayU

    return M_j, J_j, E_j    

def flux_wall_border_cython(left_border, p, ro, u, c, vbi, gamma, e_foo, **kwargs):
    rayU,rayR,rayP,*rest = \
        border_wall_URP_cython(
            p=p, ro=ro, u=u, c=c, 
            vbi=vbi, left_border=left_border, 
            p_0=0, gamma=gamma, eps_F=1e-5, n_iter_max=100)
    
    rayE = e_foo(p=rayP, ro=rayR)
    # формула 14.4
    M_j = rayR * (rayU - vbi)
    J_j = rayP + M_j * rayU
    E_j = (rayE + 0.5*rayU**2)*M_j + rayP*rayU

    return M_j, J_j, E_j 
    


def Godunov_fluxes_cython(p1, ro1, u1, c1, p2, ro2, u2, c2, vbi, gamma, e_foo, **kwargs):
    suc, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P = \
        mega_foo_cython(p_1=p1, ro_1=ro1, u_1=u1, c_1=c1, 
                 p_2=p2, ro_2=ro2, u_2=u2, c_2=c2, p_0=0, gamma=gamma, eps_F=1e-5, n_iter_max=100)
    
    if not suc:
        return 0,0,0
    rayU,rayR,rayP = get_ray_URP_cython(vbi, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P,
                         p_1=p1, ro_1=ro1, u_1=u1, c_1=c1, p_2=p2, ro_2=ro2, u_2=u2, c_2=c2, gamma=gamma)
    
    rayE = e_foo(p=rayP, ro=rayR)
    # формула 14.4
    M_j = rayR * (rayU - vbi)
    J_j = rayP + M_j * rayU
    E_j = (rayE + 0.5*rayU**2)*M_j + rayP*rayU

    return M_j, J_j, E_j