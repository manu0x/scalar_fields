

double potn_vel_eqn(double a,double a_t,double phi,double f_a,double omega_dm_0,double Xb,double *bin_coef,double a0=1.0)
{     
     double potn_t,lambda=1.0,epsilon,x_power_approx;//  x_power_approx (approx.)= pow(0.5*f_t*f_t/Mfield, alpha-1)
	
	int i;

	if(X_POWER)
	{
	  
	  epsilon = ((0.5*f_a*f_a*a_t*a_t)/Xb - 1.0);
	  x_power_approx = 1.0;

	   for(i=1;i<binomial_n;++i)
	   {
		x_power_approx+= bin_coef[i]*pow(epsilon,(double)(i));
		

	   }
	

	}

	else
	x_power_approx = 0.5*f_a*f_a*a_t*a_t/Xb;

	
	potn_t = -phi/(a) - 0.5*(H0*H0)*(omega_dm_0)*pow(a/a0,-3.0*(1.0+w))*(a/a_t)*(x_power_approx - 1.0)/a_t ;


	return(potn_t);
}

