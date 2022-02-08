

double potn_vel_eqn(double a,double a_t,double phi,double lap_phi,double f_t,double omega_dm_0,double Xb_0,double a0=1.0)
{     
     double potn_t,lambda=1.0,x_power_approx;//  x_power_approx (approx.)= pow(0.5*f_t*f_t/Mfield, alpha-1)
	x_power_approx = 3.0*H0*H0*(omega_dm_0)*pow(a,-3.0*(1.0+w))/(4.0*twopie*Xb_0*(2.0*alpha-1.0));


     potn_t  =   omega_dm_0*H0*H0*a0*a0*a0/(2.0*a_t*a*a) + 2.0*twopie*G*x_power_approx*0.5*f_t*f_t*(2.0*alpha*phi-1.0)/(3.0*a_t) + lap_phi/(3.0*a*a_t);

	//printf("lap_potn %lf %lf\n",potn_t, lap_phi);

	return(potn_t);
}

