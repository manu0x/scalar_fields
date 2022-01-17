double ini_power_spec(double ksqr)
{
	return (1e-5);
}

double res_limits(double max_potn,double vmax,double dx,double a,double &dt_limit,double &len_res )
{
	//############ Uses constraints from May & Springel  2101.01828
	double t_constrain_1,t_constrain_2,t_smaller_constrain,vmax_cap;
	len_res  = M_PI*hbar_by_m/vmax;
	vmax_cap  = M_PI*hbar_by_m/dx;

	t_constrain_1 = (4.0/(3.0*M_PI))*H0*a*a*len_res*len_res/hbar_by_m;
	t_constrain_2 = 2.0*M_PI*hbar_by_m*H0/fabs(max_potn);

	if (t_constrain_1<t_constrain_2)
		t_smaller_constrain = t_constrain_1;
	else
		t_smaller_constrain = t_constrain_2;

	
	
	dt_limit = t_smaller_constrain;

	printf("\nResolution Constraints @ z = %lf\n",(1.0/a)-1.0);
	printf("\tt_constrain_1 %e\n\tt_constrain_2 %e\n\tdt_limit %e\n\tlen_res req. %e\n\tdx is %e\n\tvmax_cap is %e\n",
							t_constrain_1,t_constrain_2,dt_limit,len_res,dx,vmax_cap);
	
	return vmax_cap;

	


}


void initialise_mpi(int * ind,int *ind_loc,fdm_psi_mpi &psi,metric_potential_mpi &phi,double k_grid[][3],int kbin_grid[],double a0,double ai,double Hi,double omega_dm_ini,double *dx,double &dk,int & kbins,ini_power_generator pk,gauss_rand_field_gen_mpi grf,bool use_hdf5_format,int cum_lin_ind)
{
      

    
      int ci,i,j,k,ret;
      int xcntr[3]={-1,-1,-1};
      
      double ktmp,maxkmagsqr = 0.0,minkmagsqr = 1e10;
     
      double a = ai;
      double a_t = a*Hi;
      double a_ti = ai*Hi;

      FILE *fpstoreini = fopen("initial.txt","w");	
      hid_t file;
     

      double L[3],boxlength;
      
           
      int tN = ind_loc[0]*ind_loc[1]*ind_loc[2]; 
      int kbin_count[tN];
	double x_grid[tN][3];
   //   double kmag_grid[tN];
      int n[3]{ind_loc[0],ind_loc[1],ind_loc[2]};
      int loc_ind[3],err_hold;

      double ini_dc[tN],ini_theta[tN],pwr_spec[tN],vel[tN][3],dc[tN];
      double psi_amp,psi_r_val,psi_i_val,poisson_rhs;
      double f_ini;
      double potn;
      double k_nyq;

      double max_potn=0.0,vmax=0.0,vmax_cap;
      double dt_limit,len_res;

      double a3a03omega = pow(a/a0,3.0)/omega_dm_ini;
	
	int my_corank;
      
	MPI_Comm_rank(cart_comm,&my_corank);
	     
	 
	f_ini=dlogD_dloga(a);

	//double kf = twopie*lenfac/(64.0);
	boxlength = 1.0/h;
	
        dx[0] = boxlength/((double)(ind[0]-1));	dx[1] = boxlength/((double)(ind[1]-1));	dx[2] = boxlength/((double)(ind[2]-1));
	L[0] = boxlength;	L[1] = boxlength;	L[2] = boxlength;
	dk = 1.0/boxlength;
	kbins = 0;
	printf("dk %lf\n",dk);
	k_nyq = 0.5/(dx[0]);
	
////////////////////// For MPI ////////////////////////////////
	xcntr[0] = cum_lin_ind -1;
////////////////////////////////////////////////////////////////
	
	//ini_rand_field();
	//  read_ini_rand_field();
        
	for(ci = 0;ci <tN; ++ci)
	{
		kbin_count[ci]=0;
	}
	


	for(ci = 0;ci <tN; ++ci)
	{
		
		
		if((ci%(n[2]*n[1]))==0)
		 ++xcntr[0];
		if((ci%(n[2]))==0)
		 ++xcntr[1];
		 ++xcntr[2];
		ktmp=0.0;
		
		for(j=0;j<3;++j)
		{	
			

			if((xcntr[j]%n[j])<=(n[j]/2))
				{
					
				 k_grid[ci][j] = ((double)(xcntr[j]%n[j]))/L[j];

				  ktmp+= k_grid[ci][j]*k_grid[ci][j];

					
				}
			else
				{ 
				 k_grid[ci][j] = ((double)((xcntr[j]%n[j])-n[j]))/L[j];

				 ktmp+= k_grid[ci][j]*k_grid[ci][j];

				
				  
				}
		
			 
			x_grid[ci][j] = ((double)xcntr[j])*dx[j];
			
		
		}
		
		
		


		
	
			
		if(ktmp>maxkmagsqr)
		maxkmagsqr = (ktmp);
		if((ktmp>=0.0)&&(minkmagsqr>ktmp))
		minkmagsqr = ktmp;
		

		kbin_grid[ci] = (int)(sqrt(ktmp)/(dk));
		//kmag_grid[ci] = sqrt(ktmp);
		 //printf("yo  %d  %lf\n",kmag_grid[ci],sqrt(ktmp));
		kbin_count[kbin_grid[ci]] = kbin_count[kbin_grid[ci]]+1;

		if(kbin_grid[ci]>kbins)
		kbins=kbin_grid[ci];
		
			

		
		
		
      	}	


	//ini_rand_field_2(ind,k_grid,kbin_grid,kbins, dk,
	//					ini_dc,ini_theta,a,a0,a_t,pk);
	//ini_rand_field(ind,kmag_grid,ini_dc,ini_theta,a,a0,a_t,pk);

	grf.gen(k_grid,ini_dc,ini_theta, pk,a_t,a,a0,f_ini);
	
	double ggg;
	for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			loc_ind[0] = i;  loc_ind[1] = j;  loc_ind[2] = k;
			psi_amp = sqrt(omega_dm_ini*pow(a0/ai,3.0)*(1.0+ini_dc[ci]));
			psi_r_val = psi_amp*cos(ini_theta[ci]);
			psi_i_val = psi_amp*sin(ini_theta[ci]);

			err_hold =  psi.update(loc_ind, psi_r_val, psi_i_val);

			//ggg = 1.5*H0*H0*a*a*omega_dm_ini*pow(a0/ai,3.0)*(ini_dc[ci]);

			

			poisson_rhs = 1.5*H0*H0*a*a*(psi_amp*psi_amp - omega_dm_ini*pow(a0/ai,3.0));

			//printf("P rhs %.15lf %.15lf %.10lf %.10lf %.10lf\n",poisson_rhs,ggg,1.5*H0*H0,psi_i_val,ini_theta[ci]);
			phi.update_4pieGpsi(ci,poisson_rhs);
			


		    }
		  }

		}

	
	ret = calculate_vel_from_psi(ind_loc,dx,psi, vel,vmax,ai);
	
	for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			loc_ind[0] = i;  loc_ind[1] = j;  loc_ind[2] = k;
			potn = phi.get_potential(ci);
			psi_amp = sqrt(omega_dm_ini*pow(a0/ai,3.0)*(1.0+ini_dc[ci]));
			psi_r_val = psi_amp*cos(ini_theta[ci]);
			psi_i_val = psi_amp*sin(ini_theta[ci]);
			
			if(fabs(potn)>=max_potn)
				max_potn = fabs(potn);

			

			//fprintf(fpstoreini,"%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n",
			//				a,dx[0]*i,dx[1]*j,dx[2]*k,psi_r_val,psi_i_val,potn,ini_dc[ci],ini_theta[ci],vel[ci][0],vel[ci][1],vel[ci][2]);
			//printf("ci %d   %.15lf\n",ci,potn);

	            }
		  }

		}
	//fclose(fpstoreini);
/*
	FILE *fpwr_spec = fopen("spec_test.txt","w");
	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	cal_spectrum(ini_dc,kbin_grid, kbins,n,pwr_spec, dk,a/ai,fpwr_spec);
	fclose(fpwr_spec);

	file = H5Fcreate("test_initial.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	initial_hdf5_write(n, psi, phi,file,dc,k_grid,x_grid,a3a03omega, ai,true);
	H5Fclose(file);
	
	vmax_cap=res_limits(max_potn, vmax, dx[0],ai,dt_limit, len_res);
	
*/
	printf("\nInitialization Complete from rank %d.\n",my_corank);
	/*printf("\n Length details\n");
	printf("	L is %lf\n",L[0]);
	printf("	dx is %lf req len res %lf\n",dx[0],len_res);
	printf("\nK details:\n	dk is %lf  per MPc\n",dk/lenfac);
	printf("	kmin %lf kmax %lf\n",sqrt(minkmagsqr),sqrt(maxkmagsqr));
	printf("	k_nyquist is %.5lf\n",k_nyq);

	printf("	kbins is %d   %d\n",kbins,(int)((sqrt(maxkmagsqr)-sqrt(minkmagsqr))/dk));
	printf("	Max vel is %e  /c\n",vmax);
	printf("	Max vel cap is %e  /c\n",vmax_cap);
	printf("	de Broglie wave is %e\n",twopie*hbar_by_m/vmax);
	printf("	Max potential is %.10lf  /c^2\n",max_potn);


	/*printf("\n Nyquist Wavenumber is %lf",M_PI/dx[0]);
	printf("\n	Min k_mag is %lf per MPc:: corr lmbda is %.16lf MPc",1.0/(dx[0]*lenfac*((double) n)),dx[0]*lenfac*((double) n));
	printf("\n	Max k_mag is %lf per MPc:: corr lmbda is %.16lf Mpc",sqrt(maxkmagsqr)/lenfac,lenfac/sqrt(maxkmagsqr));
	printf("\n	kbins is %d\n",kbins);

	printf("\nLengthscales:");
	printf("\n	Grid Length is %.6lf MPc",dx[0]*lenfac*((double) n));
	printf("\n	dx is %.16lf MPc\n",dx[0]*lenfac);

	*/
	
	
	  

}


void initial_hdf5_write(int *ind,fdm_psi_mpi psi,metric_potential_mpi phi,hid_t filename,double *dc,double k_grid[][3],double x_grid[][3],
													double a3a03omega,double a,bool get_dc=true)
{	
	herr_t status_psi,status_phi,status;	
	hid_t file,dtype,dspace,dataset;
	hsize_t dim[3],odim[2];
	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];
	odim[0] = ind[0]*ind[1]*ind[2];
	odim[1] = 3;

	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);
	dspace = H5Screate_simple(3, dim, NULL);

	status_psi=psi.write_hdf5_psi(filename, dtype, dspace,dc,a3a03omega,a,get_dc);
	
	status_phi=phi.write_hdf5_potn(filename, dtype);

	H5Sclose(dspace);
	H5Tclose(dtype);

	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);
	dspace = H5Screate_simple(2, odim, NULL);
	dataset = H5Dcreate(filename, "k_grid", dtype, dspace,
			H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, dtype, H5S_ALL, H5S_ALL,
		      				H5P_DEFAULT,k_grid);
	H5Dclose(dataset);

	

	H5Sclose(dspace);
	H5Tclose(dtype);

	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);
	dspace = H5Screate_simple(2, odim, NULL);
	dataset = H5Dcreate(filename, "x_grid", dtype, dspace,
			H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, dtype, H5S_ALL, H5S_ALL,
		      				H5P_DEFAULT,x_grid);
	H5Dclose(dataset);

	

	H5Sclose(dspace);
	H5Tclose(dtype);




}



////////////// Alternative codes for ini field generation...Currently (Jan 3rd, 2022) not being used....######################################################

void ini_rand_field(int * ind,double *kmag_grid,double * ini_dc,double * ini_theta_return,double a,double a0,double a_t,ini_power_generator gen)
{	init_genrand(time(0));
	int i,cnt,tN; 
	double ksqr,muk,sigk;
	double a1,a2,b1,b2,a_rand,b_rand;
	double f_ini;
	


	FILE *fpinirand = fopen("initial_rand_field.txt","w");

	
	
	tN = ind[0]*ind[1]*ind[2];

	fftw_plan ini_del_plan;
	fftw_plan ini_theta_plan;

	fftw_complex *F_ini_del;
	fftw_complex *ini_del;
	fftw_complex *F_ini_theta;
	fftw_complex *ini_theta;

	F_ini_del = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	ini_del = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	F_ini_theta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	ini_theta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);

	init_genrand(time(0));

	f_ini  = dlogD_dloga(a);
	
	for(cnt=0;cnt<tN;++cnt)
	{	
		 	    ksqr = kmag_grid[cnt];
			   // sigk  = sqrt(ini_power_spec(sqrt(ksqr)));
			if(ksqr>0.0)	
			    {sigk  = gen.test_spec(sqrt(ksqr));}// printf("jhfhfbsfbn %lf %lf\n",h*sqrt(ksqr),h);}
			     muk = sqrt(sigk)/pow(twopie,1.5);
		 	     a1 = genrand_res53();
 			     a2 = genrand_res53(); 
			   // b1 = genrand_res53();
 			  //  b2 = genrand_res53();
			     a_rand = (muk*(sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2)));
			     b_rand = (muk*(sqrt(-2.0*log(a1))*sin(2.0*M_PI*a2)));
				
			    
			    if(ksqr>0.0)
			    { 
			      F_ini_del[cnt][0] = a_rand;	F_ini_del[cnt][1] = b_rand;
			      F_ini_theta[cnt][0] =  (a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][0]/(ksqr*hbar_by_m);
			      F_ini_theta[cnt][1] =  (a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][1]/(ksqr*hbar_by_m);
			     
			    }	
			    else
			    { F_ini_theta[cnt][0] =  0.0;
			      F_ini_theta[cnt][1]  = 0.0;
			      F_ini_del[cnt][0] = 0.0;
			      F_ini_del[cnt][1] = 0.0;
			    }	



	}




	ini_del_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2], F_ini_del, ini_del, FFTW_BACKWARD, FFTW_ESTIMATE);
	ini_theta_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2], F_ini_theta, ini_theta, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	

	fftw_execute(ini_del_plan);
	fftw_execute(ini_theta_plan);
	

	
	for(cnt=0;cnt<tN;++cnt)
	{
		
		ini_del[cnt][0] = ini_del[cnt][0]/(tN); ini_del[cnt][1] = ini_del[cnt][1]/(tN); 
		ini_theta[cnt][0] = ini_theta[cnt][0]/(tN); ini_theta[cnt][1] = ini_theta[cnt][1]/(tN); 
		ini_dc[cnt] = ini_del[cnt][0];
		ini_theta_return[cnt] = ini_theta[cnt][0];
		

		fprintf(fpinirand,"%d\t%.16lf\t%.16lf\t%.16lf\n",
					cnt,kmag_grid[cnt]*sqrt(kmag_grid[cnt]),ini_del[cnt][0],ini_theta[cnt][0]);


		

	}
    

	 fftw_free(F_ini_del);
	 fftw_free(ini_del);
	 fftw_destroy_plan(ini_del_plan);
	
	 fftw_free(F_ini_theta);
	 fftw_free(ini_theta);
	 fftw_destroy_plan(ini_theta_plan);
	

	
}




double unit_normal()
{

		double a1,a2,r;
		a1 = genrand_res53();
 		a2 = genrand_res53(); 
		r = (sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2));

		return(r);


}



void ini_rand_field_2(int * ind,double k_grid[][3],int *kbin_grid,int kbins,double dk,
						double * ini_dc,double * ini_theta_return,double a,double a0,double a_t,ini_power_generator pk)
{	init_genrand(time(0));
	int i,cnt,tN; 
	double ksqr,pk_val;
	double a1,a2,b1,b2,a_rand,b_rand;
	double f_ini;
	double pwr_spec[kbins+1];	
	double dtN;
	






	FILE *fpinirand = fopen("initial_rand_field.txt","w");
	FILE *fppwr = fopen("ini_power.txt","w");
	FILE *fppwr_check = fopen("ini_power_check.txt","w");

	
	
	tN = ind[0]*ind[1]*ind[2];
	dtN = (double)tN;

	fftw_plan ini_del_f_plan;
	fftw_plan ini_del_b_plan;
	fftw_plan ini_theta_b_plan;

	fftw_complex *F_ini_del;
	fftw_complex *ini_del;
	fftw_complex *F_ini_theta;
	fftw_complex *ini_theta;

	F_ini_del = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	ini_del = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	F_ini_theta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	ini_theta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);

	ini_del_f_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2], ini_del, F_ini_del, FFTW_FORWARD, FFTW_ESTIMATE);
	ini_del_b_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2],  F_ini_del,ini_del, FFTW_BACKWARD, FFTW_ESTIMATE);
	ini_theta_b_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2],  F_ini_theta,ini_theta, FFTW_BACKWARD, FFTW_ESTIMATE);

	init_genrand(time(0));

	for(cnt=0;cnt<tN;++cnt)
	{
		ini_del[cnt][0] = unit_normal();
		ini_del[cnt][1] = 0.0;
	}
	
	fftw_execute(ini_del_f_plan);

	f_ini  = dlogD_dloga(a);
	
	for(cnt=0;cnt<tN;++cnt)
	{	
		 	    

			ksqr = (k_grid[cnt][0]*k_grid[cnt][0]+k_grid[cnt][1]*k_grid[cnt][1]+k_grid[cnt][2]*k_grid[cnt][2]);

			pk_val = pk.test_spec(sqrt(ksqr));	
		
			if(ksqr<=0.0)
			{
			  F_ini_del[cnt][0] = 0.0;
			  F_ini_del[cnt][1] = 0.0;

			  F_ini_theta[cnt][0] =  0.0;
			  F_ini_theta[cnt][1]  = 0.0;
			}
			else
			{

				F_ini_del[cnt][0] = F_ini_del[cnt][0]/sqrt(pow(twopie,3.0));
				F_ini_del[cnt][1] = F_ini_del[cnt][1]/sqrt(pow(twopie,3.0));

				F_ini_theta[cnt][0] =  (a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][0]/(ksqr*hbar_by_m);
			        F_ini_theta[cnt][1] =  (a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][1]/(ksqr*hbar_by_m);

			}



	}



	fftw_execute(ini_del_b_plan);
	fftw_execute(ini_theta_b_plan);
	

	
	for(cnt=0;cnt<tN;++cnt)
	{
		
		ini_del[cnt][0] = ini_del[cnt][0]/(dtN); ini_del[cnt][1] = ini_del[cnt][1]/(dtN); 
		ini_theta[cnt][0] = ini_theta[cnt][0]/(dtN); ini_theta[cnt][1] = ini_theta[cnt][1]/(dtN); 
		ini_dc[cnt] = ini_del[cnt][0];
		ini_theta_return[cnt] = ini_theta[cnt][0];

		ksqr = (k_grid[cnt][0]*k_grid[cnt][0]+k_grid[cnt][1]*k_grid[cnt][1]+k_grid[cnt][2]*k_grid[cnt][2]);
		pk_val = pk.test_spec(sqrt(ksqr));

		

		fprintf(fpinirand,"%d\t%.16lf\t%.16lf\t%.16lf\n",
					cnt,ksqr,ini_del[cnt][0],ini_theta[cnt][0]);


		

	}


	cal_spectrum(ini_dc,kbin_grid, kbins,ind,pwr_spec, dk,1.0,fppwr);

	for(cnt=0;cnt<=kbins;++cnt)
	{
		
		

		ksqr = ((double)(cnt+1))*dk*(0.5);
		pk_val = pk.test_spec(ksqr);

		

		fprintf(fppwr_check,"%d\t%.16lf\t%.16lf\t%.16lf\n",
					cnt,ksqr,pwr_spec[cnt],pk_val);


		

	}


	fclose(fpinirand);
	fclose(fppwr);
	fclose(fppwr_check);
    

	 fftw_free(F_ini_del);
	 fftw_free(ini_del);
	 fftw_free(F_ini_theta);
	 fftw_free(ini_theta);

	 fftw_destroy_plan(ini_del_f_plan);
	 fftw_destroy_plan(ini_del_b_plan);
	 fftw_destroy_plan(ini_theta_b_plan);

	

	
}



