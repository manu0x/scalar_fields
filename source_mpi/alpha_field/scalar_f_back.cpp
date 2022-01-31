void set_back_cosmo(double &a0,double &ai,double &Hi,double &omega_dm_0,param_alpha p)
{
	double z = p.z_ini;
	alpha = p.loc_alpha;
	w = 1.0/(2.0*alpha-1.0);

	c_box = p.loc_c_box;
	pc_box = p.loc_pc_box;
	hbar_box = p.loc_hbar_box;//eVs 
	lenfac = p.loc_lenfac;
	omega_dm_0 = p.omega_dm_0;
	h = p.loc_h;
	//hbar_by_m = hbar_box*c_box*(1e-8)/(alpha*pc_box);	
	space_mpc_to_dimless = p.loc_space_mpc_to_dimless; ////	\tilde{x} (dimensionless) = physical{x (In Mpc)}*space_mpc_to_dimless  
		
	
	
	H0 = lenfac*(h/c_box)*0.001;
	//printf("H0 %lf h  %e\n",H0,h);
	
	
	a0 = p.a0;
	ai = a0/(1.0+z);
	Hi =   H0*sqrt(omega_dm_0*pow(a0/ai,3.0*(1.0+w))+ (1.0-omega_dm_0));

	
	p.print_to_file(fp_sim_info);
	

}


double dlogD_dloga(double a)
{
	return (1.0);

}
