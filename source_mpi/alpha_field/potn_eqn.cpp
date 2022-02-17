

double potn_vel_eqn(double a,double a_t,double phi,double lap_phi,double f_t,double omega_dm_0,double Xb_0,double a0=1.0)
{     
     double potn_t,lambda=1.0,x_power_approx;//  x_power_approx (approx.)= pow(0.5*f_t*f_t/Mfield, alpha-1)
	double Xb = Xb_0*pow(a0/a,6.0/(2.0*alpha-1.0));
	x_power_approx = 3.0*H0*H0*(omega_dm_0)*pow(a,-3.0*(1.0+w))/(4.0*twopie*Xb*(2.0*alpha-1.0));
	
	potn_t = -phi*(a_t/a) + lap_phi/(3.0*a_t*a) + 2.0*twopie*G*a*(2.0*alpha-1.0)*x_power_approx*((Xb/(H0*H0))-0.5*f_t*f_t/(1.0+2.0*phi))/(3.0*a_t);


// potn_t = -phi*(a_t/a) + lap_phi/(3.0*a_t*a) + 2.0*twopie*G*a*(2.0*alpha-1.0)*x_power_approx*((Xb/(H0*H0))-0.5*f_t*f_t/(1.0+2.0*phi))/(3.0*a_t); (x->dimless(x))

    
//	printf("lap_potn %.10lf %.10lf %.10lf\n",potn_t, (2.0*alpha-1.0)*2.0*twopie*G*a*x_power_approx*((Xb/(H0*H0))-0.5*f_t*f_t)/(3.0*a_t) ,
//	2.0*(2.0*alpha-1.0)*twopie*G*a*alpha*x_power_approx*f_t*f_t*phi/(3.0*a_t)  );

	/*
		 potn_t  = // (2.0*alpha-1.0)*2.0*twopie*G*a*x_power_approx*((Xb/(H0*H0))-0.5*f_t*f_t)/(3.0*a_t) 
		//+ 2.0*(2.0*alpha-1.0)*twopie*G*a*alpha*x_power_approx*f_t*f_t*phi/(3.0*a_t) 
		+ lap_phi/(3.0*a*a_t)- phi*(a_t/a);


	*/


	return(potn_t);
}

