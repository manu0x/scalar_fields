
double field_acc_approx(double f_t,double f_sx[3],double f_t_x[3],double phi,double phi_t,double phi_sx[3])
{
	double acc,phi_x,phi_x,phi_y,phi_z,phi_x,f_x,f_y,f_z,f_tx,f_ty,f_tz;    
	
	phi_x = phi_sx[0];	phi_y = phi_sx[1];	phi_z = phi_sx[2];
	f_x = f_sx[0];	f_y = f_sx[1];	f_z = f_sx[2];
	f_tx = f_t_x[0];	f_ty = f_t_x[1];	f_tz = f_t_x[2];

	acc = -(pow(a, -3.0)*pow(-1.0 + 2.0*alpha, -1.0)*pow(-1.0 + 2.0*alpha*phi, -1.0)*pow(f_t, -1.0)*
	  	(3.0*(-1.0 + 2.0*alpha*phi)*(a_t)*pow(a, 2.0)*pow(f_t, 2.0) + 2.0*(1.0 + alpha)*(phi_t)*pow(a, 3.0)*pow(f_t, 2.0)
			 + (-1.0 + 2.0*(-2.0 + alpha)*phi)*(a_t)*(pow(f_x, 2.0) + pow(f_y, 2.0) + pow(f_z, 2.0)) +
			 a*(-((-1.0 + 2.0*alpha)*(-1.0 + 2.0*(-2.0 + alpha)*phi)*((f_tx)*(f_x) + (f_ty)*(f_y) + (f_tz)*(f_z))) 
			+ (f_t)*(  ((4.0-2.0*alpha)*phi+1.0)*lap_f  
			+ 2.0*(f_x)*(phi_x) - 2.0*(alpha*((f_x)*(phi_x)) + (-1.0 + alpha)*(f_y)*(phi_y)) 
			- 2.0*(-1.0 + alpha)*(f_z)*(phi_z)) + (phi_t)*(pow(f_x, 2.0) + pow(f_y, 2.0) + pow(f_z, 2.0)))));

	return acc;

}

