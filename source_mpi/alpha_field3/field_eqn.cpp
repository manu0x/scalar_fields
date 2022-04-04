
double field_eqn(double f_a,double phi,double phi_a,double a,double da,double a_t,double a_tt)
{
	double numer,denom;	
	double acc,phi_x,phi_y,phi_z,f_x,f_y,f_z,f_tx,f_ty,f_tz,f_ssqr;    
	
/*	phi_x = phi_sx[0];	phi_y = phi_sx[1];	phi_z = phi_sx[2];
	f_x = f_sx[0];	f_y = f_sx[1];	f_z = f_sx[2];
	f_tx = f_t_x[0];	f_ty = f_t_x[1];	f_tz = f_t_x[2];
	f_ssqr = f_x*f_x + f_y*f_y + f_z*f_z;

*/	

	
	denom = (2.0*alpha-1.0)*(2.0*alpha*phi-1.0);	

	numer = f_a*(1.0 +  0.5*da*(3.0/a - 6.0*alpha*phi - 2.0*phi_a*(1.0+alpha))/denom - 0.5*da*a_tt/(a_t*a_t)   );
			

	



	acc = numer;

	return acc;




}





double field_acc_eqn(double f_a,double phi,double phi_a,double lap_f,double a,double a_t,double a_tt)
{
	double numer,denom;	
	double acc;    
	

	denom = (2.0*alpha-1.0)*(2.0*alpha*phi-1.0);	

	numer =  f_a*((3.0/a - 6.0*alpha*phi - 2.0*phi_a*(1.0+alpha))/denom - a_tt/(a_t*a_t)   ) - lap_f/(a_t*a*a_t*a*(2.0*alpha-1.0));
			

	



	acc = numer;
	return acc;



}

