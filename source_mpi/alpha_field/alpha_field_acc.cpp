
double alpha_field_tt()
{

double f_tt  = 

pow[a, -2.0]*


(3.0*(1.0 - 2.0*alpha)*(a_t)*pow[a, 5.0]*pow[f_t, 5.0] + 2.0*(phi_t)*pow[a, 6.0]*(-1.0 + alpha + 2.0*pow[alpha, 2.0])*
pow[f_t, 5.0] - 2.0*(1.0 + 8.0*phi*(-2.0 + alpha)*(-1.0 + alpha) + 2.0*(-3.0 + alpha)*alpha)*
(a_t)*pow[a, 3.0]*pow[f_t, 3.0]*(pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0]) + pow[a, 4.0]*pow[f_t, 3.0]*
(-((-1.0 + 2.0*alpha)*(-4.0*(1.0 + 4.0*phi)*(-1.0 + alpha)*((f_tx)*(f_x) + (f_ty)*(f_y) + (f_tz)*(f_z)) + (f_t)*(-(f_xx) - (f_yy) - 4.0*phi*((f_xx) + (f_yy)) - (1.0 + 4.0*phi)*(f_zz) + 2.0*(-1.0 + alpha)*(f_x)*(phi_x) + 2.0*(-1.0 + alpha)*(f_y)*(phi_y) + 2.0*(-1.0 + alpha)*(f_z)*(phi_z)))) + 4.0*(1.0 + (-4.0 + alpha)*alpha)*(phi_t)*(pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0]))+
(pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0])*
(-4.0*(f_x)*(f_xy)*(f_y) - 16*phi*(f_x)*(f_xy)*(f_y) + 4.0*alpha*(f_x)*(f_xy)*(f_y) + 16*phi*alpha*(f_x)*(f_xy)*(f_y) + (f_yy)*pow[f_x, 2.0] + 4.0*phi*(f_yy)*pow[f_x, 2.0] - 2.0*(f_y)*(phi_y)*pow[f_x, 2.0] + 2.0*alpha*(f_y)*(phi_y)*pow[f_x, 2.0] - 2.0*(phi_x)*pow[f_x, 3.0] + 2.0*alpha*(phi_x)*pow[f_x, 3.0] - (f_yy)*pow[f_y, 2.0] - 4.0*phi*(f_yy)*pow[f_y, 2.0] + 2.0*alpha*(f_yy)*pow[f_y, 2.0] + 8.0*phi*alpha*(f_yy)*pow[f_y, 2.0] - 2.0*(f_x)*(phi_x)*pow[f_y, 2.0] + 2.0*alpha*(f_x)*(phi_x)*pow[f_y, 2.0] + (1.0 + 4.0*phi)*(f_zz)*(pow[f_x, 2.0] + pow[f_y, 2.0]) + (1.0 + 4.0*phi)*(f_xx)*((-1.0 + 2.0*alpha)*pow[f_x, 2.0] + pow[f_y, 2.0]) + 2.0*(-1.0 + alpha)*(f_z)*(2.0*(1.0 + 4.0*phi)*((f_x)*(f_xz) + (f_y)*(f_yz)) + (phi_z)*(pow[f_x, 2.0] + pow[f_y, 2.0])) - 2.0*(phi_y)*pow[f_y, 3.0] + 2.0*alpha*(phi_y)*pow[f_y, 3.0] + ((f_xx) + 4.0*phi*(f_xx) + (1.0 + 4.0*phi)*(f_yy) + (1.0 + 4.0*phi)*(-1.0 + 2.0*alpha)*(f_zz) + 2.0*(-1.0 + alpha)*(f_x)*(phi_x) + 2.0*(-1.0 + alpha)*(f_y)*(phi_y))*pow[f_z, 2.0] + 2.0*(-1.0 + alpha)*(phi_z)*pow[f_z, 3.0]) +
 a*(-5.0 + 2.0*alpha)*(a_t)*(f_t)*pow[pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0], 2.0] - 2.0*(f_t)*pow[a, 2.0]*(2.0*(-1.0 + alpha)*((f_tx)*(f_x) + (f_ty)*(f_y) + (f_tz)*(f_z))*(pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0]) + (f_t)*(2.0*(f_x)*(f_xy)*(f_y) + 16*phi*(f_x)*(f_xy)*(f_y) - 6.0*alpha*(f_x)*(f_xy)*(f_y) - 48*phi*alpha*(f_x)*(f_xy)*(f_y) + 4.0*(f_x)*(f_xy)*(f_y)*pow[alpha, 2.0] + 32*phi*(f_x)*(f_xy)*(f_y)*pow[alpha, 2.0] - 4.0*phi*(f_yy)*pow[f_x, 2.0] + alpha*(f_yy)*pow[f_x, 2.0] + 8.0*phi*alpha*(f_yy)*pow[f_x, 2.0] + 2.0*(f_y)*(phi_y)*pow[f_x, 2.0] - 4.0*alpha*(f_y)*(phi_y)*pow[f_x, 2.0] + 2.0*(f_y)*(phi_y)*pow[alpha, 2.0]*pow[f_x, 2.0] + 2.0*(phi_x)*pow[f_x, 3.0] - 4.0*alpha*(phi_x)*pow[f_x, 3.0] + 2.0*(phi_x)*pow[alpha, 2.0]*pow[f_x, 3.0] + (f_yy)*pow[f_y, 2.0] + 4.0*phi*(f_yy)*pow[f_y, 2.0] - 2.0*alpha*(f_yy)*pow[f_y, 2.0] - 16*phi*alpha*(f_yy)*pow[f_y, 2.0] + 2.0*(f_x)*(phi_x)*pow[f_y, 2.0] - 4.0*alpha*(f_x)*(phi_x)*pow[f_y, 2.0] + 2.0*(f_yy)*pow[alpha, 2.0]*pow[f_y, 2.0] + 16*phi*(f_yy)*pow[alpha, 2.0]*pow[f_y, 2.0] + 2.0*(f_x)*(phi_x)*pow[alpha, 2.0]*pow[f_y, 2.0] + (alpha + phi*(-4.0 + 8.0*alpha))*(f_zz)*(pow[f_x, 2.0] + pow[f_y, 2.0]) + (f_xx)*((1.0 + 2.0*(-1.0 + alpha)*alpha + 4.0*phi*pow[1.0 - 2.0*alpha, 2.0])*pow[f_x, 2.0] + (alpha + phi*(-4.0 + 8.0*alpha))*pow[f_y, 2.0]) + 2.0*(-1.0 + alpha)*(f_z)*((1.0 + 8.0*phi)*(-1.0 + 2.0*alpha)*((f_x)*(f_xz) + (f_y)*(f_yz)) + (-1.0 + alpha)*(phi_z)*(pow[f_x, 2.0] + pow[f_y, 2.0])) + 2.0*(phi_y)*pow[f_y, 3.0] - 4.0*alpha*(phi_y)*pow[f_y, 3.0] + 2.0*(phi_y)*pow[alpha, 2.0]*pow[f_y, 3.0] + ((alpha + phi*(-4.0 + 8.0*alpha))*(f_xx) - 4.0*phi*(f_yy) + alpha*(f_yy) + 8.0*phi*alpha*(f_yy) + 2.0*(f_x)*(phi_x) - 4.0*alpha*(f_x)*(phi_x) + (f_zz)*(1.0 + 2.0*(-1.0 + alpha)*alpha + 4.0*phi*pow[1.0 - 2.0*alpha, 2.0]) + 2.0*(f_y)*(phi_y)*pow[-1.0 + alpha, 2.0] + 2.0*(f_x)*(phi_x)*pow[alpha, 2.0])*pow[f_z, 2.0] + 2.0*(phi_z)*pow[-1.0 + alpha, 2.0]*pow[f_z, 3.0]) + (-3.0 + alpha)*(phi_t)*pow[pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0], 2.0]))

*pow[(1.0 - 2.0*alpha)*pow[a, 2.0]*pow[f_t, 2.0] + pow[f_x, 2.0] + pow[f_y, 2.0] + pow[f_z, 2.0], -2.0];


}

