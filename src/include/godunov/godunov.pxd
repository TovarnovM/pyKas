# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False, wraparound=False

cpdef inline double abs(double x) nogil

cpdef  double get_e_13_1(double p, double ro, double p_0, double c_0, double gamma) nogil

cpdef  double get_p_13_1(double e, double ro, double p_0, double c_0, double gamma) nogil

cpdef  double get_p_0(double ro_0, double c_0, double gamma) nogil

cpdef  double get_R_13_2(double ro, double p, double P, double p_0, double gamma) nogil

cpdef  double get_a_13_4_shock(double ro, double p, double P, double p_0, double gamma) nogil

cpdef  double get_c_13_8(double p, double ro, double p_0, double gamma) nogil

cpdef double get_a_13_11_discharge(double ro, double p, double P, double c, double p_0, double gamma) nogil

cpdef double get_f_13_16(double P, double p_k, double ro_k, double c_k, double p_0, double gamma) nogil

cpdef double get_df_13_17(double P, double p_k, double ro_k, double c_k, double p_0, double gamma) nogil

cpdef double get_ddf_13_18(double P, double p_k, double ro_k, double c_k, double p_0, double gamma) nogil

cpdef  double get_F_13_15(double P, double p_1, double ro_1, double c_1, double p_2, double ro_2, double c_2, double p_0, double gamma) nogil

cpdef  double get_dF_13_15(double P, double p_1, double ro_1, double c_1, double p_2, double ro_2, double c_2, double p_0, double gamma) nogil

cpdef  double get_ddF_13_15(double P, double p_1, double ro_1, double c_1, double p_2, double ro_2, double c_2, double p_0, double gamma) nogil

cpdef (double, double, double) get_Us_13_22(double p_1, double p_2, double ro_1, double c_1, double c_2, double p_0, double gamma) nogil

cdef class MegaFooResult:
    cdef public bint suc, UD_left, UD_right
    cdef public double D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P

cdef class Border_URPDD_Result:
    cdef public double rU,rR,rP,D_1,D_2,U_kr

cpdef MegaFooResult mega_foo_cython(double p_1, double ro_1, double u_1, double c_1, \
             double p_2, double ro_2, double u_2, double c_2, \
             double p_0, double gamma, double eps_F=*, int n_iter_max=*)

cpdef Border_URPDD_Result border_wall_URPDD_result(bint left_border, double vbi, double p, double ro, double u, double c, \
             double p_0, double gamma, double eps_F=*, int n_iter_max=*)

cpdef (double, double, double, double, double,double) border_wall_URPDD(bint left_border, double vbi, double p, double ro, double u, double c, \
             double p_0, double gamma, double eps_F=*, int n_iter_max=*) nogil

cpdef (bint, bint, bint, double, double, double, double, double, double, double, double) \
    mega_foo(double p_1, double ro_1, double u_1, double c_1, \
             double p_2, double ro_2, double u_2, double c_2, \
             double p_0, double gamma, double eps_F=*, int n_iter_max=*) nogil

cpdef (double, double, double) get_ray_URP( \
    double ray_W, bint UD_left, bint UD_right, double D_1, double D_star_1, double U, \
    double D_star_2, double D_2, double R_1, double R_2, double P, \
    double p_1, double ro_1, double u_1, double c_1, \
    double p_2, double ro_2, double u_2, double c_2, double gamma) nogil