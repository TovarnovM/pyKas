# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False
from tube cimport InterpXY

cdef class DirectBall(object):
    cdef public dict opts  
    cdef double I_k,alpha_k,ro,omega,S,d,fi,q,F_0,delta,f,k,R_g,v_T,sigma_T,p_f,p0,l_max,t_max,dt
    cdef InterpXY dpsi_dz

    cpdef double[:] euler_step(self, double dt, double[:] y, double[:] dy) 
    cpdef void init_consts(self)
    cdef void get_dydt(self, double t, double[:] y, double[:] res)
    cdef int p_f_foo(self, double p, double V)
    cdef double get_dpsi(self, double z, double dzdt)


cdef class DirectBallMany(object):
    cdef public dict opts
    cdef public list dpsi_dz_s
    cdef double[:] I_ks,alpha_ks,ros,omegas, fs, T_1s, ks
    cdef public object stop_foo
    
    cdef double S,d,fi,q,F_0,k,R_g,v_T,sigma_T,p_f,p0,l_max,t_max,dt, omega, W_kam
    cdef int n_powders

    cpdef double[:] euler_step(self, double dt, double[:] y, double[:] dy)
    cpdef void init_consts(self)
    cdef void get_dydt(self, double t, double[:] y, double[:] res)
    cdef int p_f_foo(self, double p, double V)
    cdef double get_dpsi(self, double z, double dzdt)