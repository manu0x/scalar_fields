

double potn_approx_vel(double a,double a_t,double phi,double lap_phi,double f_t,double omega_de_0)
{     
     double potn_t,lambda=1.0,x_power_approx;//  x_power_approx (approx.)= pow(0.5*f_t*H0*f_t*H0, alpha)
	x_power_approx = 3.0*H0*(1.0-omega_de_0)*pow(a,-3.0*(1.0+w))/(4.0*twopie*(2.0*alpha-1.0));
     potn_t  =   (pow(a, -1.0)*pow(Mfield, -1.0)*pow(a_t, -1.0)
			*(-(Mfield*(-2.0*(lap_phi) + 3.0*omega_de_0*pow(a, 2.0) + (-3.0 + 6.0*phi)*pow(a_t, 2.0))) 
				+ 2.0*(-1.0 + 2.0*alpha)*G*(-1.0 + 2.0*alpha*phi)*twopie*pow(a, 2.0)*x_power_approx))/3.0;

	//printf("lap_potn %lf %lf\n",potn_t, lap_phi);

	return(potn_t);
}

