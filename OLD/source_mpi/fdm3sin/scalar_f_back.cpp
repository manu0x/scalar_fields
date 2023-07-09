void set_back_cosmo(double &a0,double &ai,double &Hi,double &omega_dm_0,param_fdm p)
{
	double z = p.z_ini;
	double kji,kj0;
	
	w = 0.0;//1.0/(2.0*alpha-1.0);
	cs2 = 1.0;///(2.0*alpha-1.0);

	c_box = p.loc_c_box;
	pc_box = p.loc_pc_box;
	hbar_box = p.loc_hbar_box;//eVs 
	lenfac = p.loc_lenfac;
	omega_dm_0 = p.omega_dm_0;
	h = p.loc_h;
	alpha = p.loc_alpha_m22;
	hbar_by_m = 0.001*hbar_box*h*c_box*c_box/(alpha*pc_box);	
	space_mpc_to_dimless = p.loc_space_mpc_to_dimless; ////	\tilde{x} (dimensionless) = physical{x (In Mpc)}*space_mpc_to_dimless  

	method = p.loc_method;
	tot_steps = p.tot_steps;

	nwave[0] = p.nx; nwave[1] = p.ny; nwave[2] = p.nz;
	nwave_amp = sqrt(nwave[0]*nwave[0] + nwave[1]*nwave[1] + nwave[2]*nwave[2]);

	T = twopie/alpha; 	
	
	
	H0 = 100.0;
	
	
	
	a0 = 2.0*T;
	ai = 0.0;
	Hi =   H0*sqrt(omega_dm_0*pow(a0/ai,3.0*(1.0+w))+ (1.0-omega_dm_0));
	kji = pow(6.0*omega_dm_0/(1.0+p.z_ini),0.25)*sqrt(H0/hbar_by_m);
	kj0 = pow(6.0*omega_dm_0,0.25)*sqrt(H0/hbar_by_m);
	
	printf("H0 %lf Hi  %e\n",H0,Hi);
	printf("Jeans length at starting z is %lf kj is %lf\n",1.0/kji,kji);
	printf("Jeans length at z = 0 is %lf kj is %lf\n",1.0/kj0,kj0);

	fprintf(fp_sim_info,"Jeans length at starting z is %lf kj is %lf\n",1.0/kji,kji);
	fprintf(fp_sim_info,"Jeans length at z = 0 is %lf kj is %lf\n",1.0/kj0,kj0);
	fprintf(fp_sim_info,"\nnwave is %lf  %lf  %lf\n\n",nwave[0],nwave[1],nwave[2]);
	fprintf(fp_sim_info,"\nalpha is %lf\n\n",alpha);
	fprintf(fp_sim_info,"\nT is %lf  a0 is %lf  ai is %lf\n\n",T,a0,ai);
	fprintf(fp_sim_info,"\ntot_steps are %lf\n\n",tot_steps);
	fprintf(fp_sim_info,"\nhbar_by_m is %lf\n\n",hbar_by_m);

	

	
	p.print_to_file(fp_sim_info);
	

}


double dlogD_dloga(double a)
{
	return (1.0);

}
