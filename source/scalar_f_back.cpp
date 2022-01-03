void set_back_cosmo(double &a0,double &ai,double &Hi,double &omega_dm_ini)
{
	double z = 99.0;
	double alpha = 1.0;// Mass in 10^(-22) eV;	

	c_box = 2.99;
	pc_box = 3.0857;
	hbar_box = 6.582119569;//eVs 
	lenfac = 1.0;
	omega_dm_ini = 0.29;
	h = 0.7;
	hbar_by_m = hbar_box*c_box*(1e-8)/(alpha*pc_box);	
		
	
	
	H0 = lenfac*(h/c_box)*0.001;
	printf("H0 %lf hbar_by_m  %.15lf\n",H0,hbar_by_m);
	
	
	a0 = 1.0;
	ai = a0/(1.0+z);
	Hi =   H0*sqrt(omega_dm_ini*pow(a0/ai,3.0)+ (1.0-omega_dm_ini));


}


double dlogD_dloga(double a)
{
	return (1.0);

}
