
double field_acc_eqn(double f_a,double phi,double phi_a,double a,double a_t,double a_tt)
{
	double numer,denom;	
	double acc,phi_x,phi_y,phi_z,f_x,f_y,f_z,f_tx,f_ty,f_tz,f_ssqr;    
	
/*	phi_x = phi_sx[0];	phi_y = phi_sx[1];	phi_z = phi_sx[2];
	f_x = f_sx[0];	f_y = f_sx[1];	f_z = f_sx[2];
	f_tx = f_t_x[0];	f_ty = f_t_x[1];	f_tz = f_t_x[2];
	f_ssqr = f_x*f_x + f_y*f_y + f_z*f_z;

*/	
	numer = 3.0*(-1.0+2.0*alpha*phi)*f_a/a + 2.0*(1.0+alpha)*f_a*phi_a  ;
			// + (-1.0+2.0*(-2.0+alpha)*phi)*(a_t/a)*f_ssqr/(f_t*a*a);
	//numer+= ( phi_t*f_ssqr/(f_t*a*a) - (-1.0+2.0*alpha)*(-1.0+2.0*(-2.0+alpha)*phi)*( f_x*f_tx + f_y*f_ty + f_z*f_tz )/(f_t*a*a) 
	//	   +( -2.0*(-1.0+alpha)*( f_x*phi_x + f_y*phi_y + f_z*phi_z )/(a*a) + (1.0-2.0*(-2.0+alpha)*phi)*lap_f/(a*a)  )     )/(H0*H0);

	denom = (1.0-2.0*alpha)*(-1.0+2.0*alpha*phi);	



	acc = (numer/denom) - a_tt*f_a/(a_t*a_t) + 4.0*phi*lap_f/((-1.0+2.0*alpha)*a*a*a_t*a_t);

	//printf("acc %lf %.15lf\n",  pow(f_t, -1.0),f_t);

	return acc;


	/*
	
		
	numer = 3.0*(a_t/a)*f_t ;
	

	denom = (-1.0+2.0*alpha);
	

	acc = numer/denom;*/

}

