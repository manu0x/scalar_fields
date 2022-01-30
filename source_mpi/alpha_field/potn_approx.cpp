

double potn_approx_vel(double a,double a_t,double phi,double lap_phi,double f_t,double omega_de_0)
{     
     double potn_t,lambda=1.0;
     potn_t  =   (pow(2.0, -1.0 - alpha)*pow(a, -1.0)*pow(Mfield, -1.0)*pow(a_t, -1.0)
			*(-(Mfield*pow(2.0, alpha)*(-2.0*(lap_phi) + 3.0*omega_de_0*pow(a, 2.0) + (-3.0 + 6.0*phi)*pow(a_t, 2.0))) 
				+ 4.0*(-1.0 + 2.0*alpha)*G*(-1.0 + 2.0*alpha*phi)*twopie*pow(a, 2.0)*pow(pow(f_t, 2.0), alpha)))/3.0;

	return(potn_t);
}

