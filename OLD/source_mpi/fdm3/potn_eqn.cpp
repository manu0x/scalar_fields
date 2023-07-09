

double potn_vel_eqn(double a,double a_t,double phi,double psi_amp2,double omega_dm_0,double psi_amp2_b=0.0,double a0=1.0)
{     
     double potn_a,lambda=1.0;

	
	potn_a = -phi/(a) - (c_box*100000.0/a_t)*(c_box*100000.0/a_t)*(1.0/(6.0*a*a))*(psi_amp2 - 3.0*H0*H0*omega_dm_0*a0*a0*a0) ;

	return(potn_a);
}

