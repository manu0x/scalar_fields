double ini_power_spec(double ksqr)
{
	return (1e-5);
}
/*
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

*/

void read_dc_from_hdf5(string fname,double *dc,double *theta,int *ind, int cum_lin_ind)
{

	herr_t status;	
	hid_t file,dtype,dspace_glbl,dspace,dset_glbl,plist_id,g_id;
	MPI_Info info  = MPI_INFO_NULL;

	hsize_t dim[3],odim[2];

	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];
	odim[0] = ind[0]*ind[1]*ind[2];
	odim[1] = 1;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, cart_comm, info);

	
	char fchar_name[fname.length()];
	fname.copy(fchar_name,fname.length()-1,0);
	printf("FINI_dc %s %d\n",fchar_name,fname.length());
	
	file = H5Fopen(fchar_name, H5F_ACC_RDONLY, plist_id);
	H5Pclose(plist_id);

	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);


	//dspace_glbl = H5Screate_simple(2, odim, NULL);

	//g_id = H5G.open(file, "/");

	dset_glbl = H5Dopen(file, "/dc",H5P_DEFAULT);

	dspace_glbl = H5Dget_space(dset_glbl);

	hsize_t count[2],offset[2];
	count[0] = ind[0]*ind[1]*ind[2]; count[1] = odim[1]; offset[0] = cum_lin_ind*ind[1]*ind[2]; offset[1] = 0;

	

	dspace = H5Screate_simple(2, count, NULL);

	

	H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	

	status = H5Dread(dset_glbl, dtype, dspace, dspace_glbl, plist_id, dc);

	H5Dclose(dset_glbl);
	H5Pclose(plist_id);
	H5Sclose(dspace);


	dset_glbl = H5Dopen(file, "/theta",H5P_DEFAULT);

	dspace_glbl = H5Dget_space(dset_glbl);


	count[0] = ind[0]*ind[1]*ind[2]; count[1] = odim[1]; offset[0] = cum_lin_ind*ind[1]*ind[2]; offset[1] = 0;

	

	dspace = H5Screate_simple(2, count, NULL);

	

	H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	

	status = H5Dread(dset_glbl, dtype, dspace, dspace_glbl, plist_id, theta);

	H5Dclose(dset_glbl);
	H5Pclose(plist_id);
	H5Sclose(dspace);




	int my_corank,mpi_check;
/*      
	MPI_Comm_rank(cart_comm,&my_corank);

	FILE *fpstoreini;

	if(my_corank==0)
	fpstoreini = fopen("readtest1.txt","w");
	else
	if(my_corank==1)
	fpstoreini = fopen("readtest2.txt","w");
	else
	if(my_corank==2)
	fpstoreini = fopen("readtest3.txt","w");
	else
	if(my_corank==3)
	fpstoreini = fopen("readtest4.txt","w");

	int i,j,k,ci;

	for(i=0;i<ind[0];++i)
	{
		for(j=0;j<ind[1];++j)
		{
		   for(k=0;k<ind[2];++k)
		   {
			ci = (ind[2]*ind[1])*i + ind[2]*j + k;
			

			

			fprintf(fpstoreini,"%d\t%d\t%d\t%.10lf\n",
							cum_lin_ind+i,j,k,dc[ci]);

	           }
		}

	}

	fclose(fpstoreini);
*/

}

void initialise_mpi(int * ind,int *ind_loc,fdm_poisson_mpi &psi,metric_potential_poisson_mpi &phi,
					metric_potential_poisson_mpi_ini &poisson_phi,
				double k_grid[][3],int kbin_grid[],double a0,double ai,double Hi,double omega_dm_0,double & Xb_0,double *dx,double &dk,int & kbins,
				ini_power_generator pk,gauss_rand_field_gen_mpi grf,bool use_hdf5_format,double boxlength,
											double da,int cum_lin_ind,string fini_name="None")
{
      

    
      int ci,i,j,k,ret;
      int xcntr[3]={cum_lin_ind-1,-1,-1};
      
      double ktmp,maxkmagsqr = 0.0,minkmagsqr = 1e10;
     
      double a = ai;
      double a_t = a*Hi;
      double a_ti = ai*Hi;

      FILE *fpstoreini;
     	
      
    

      double L[3];
      
           
      int tN = ind_loc[0]*ind_loc[1]*ind_loc[2]; 
      int * kbin_count = new int[tN];
      double * ini_dc = new double[tN];
      double *ini_theta= new double[tN];
      double *pwr_spec=new double[tN];
      double *dc = new double[tN];
	
	double x_grid[tN][3];
   //   double kmag_grid[tN];
      int n[3]{ind_loc[0],ind_loc[1],ind_loc[2]};
      int loc_ind[3],err_hold;

     

      double poisson_rhs;
      double psi_amp, psi_val[2];
      double f_ini;
      double potn;
      double k_nyq;

      double max_potn=0.0,vmax=0.0,vmax_cap;
      double dt_limit,len_res;

      double a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;

/////////////////////////	MPI and hdf5 variables	///////////////////////////////
	
	int my_corank,mpi_check;
      
	MPI_Comm_rank(cart_comm,&my_corank);

	if(my_corank==0)
	fpstoreini = fopen("initial1.txt","w");
	else
	if(my_corank==1)
	fpstoreini = fopen("initial2.txt","w");
	else
	if(my_corank==2)
	fpstoreini = fopen("initial3.txt","w");
	else
	if(my_corank==3)
	fpstoreini = fopen("initial4.txt","w");

	MPI_Info info  = MPI_INFO_NULL;

	hid_t file;
	hid_t plist_id;
	     
	 
	f_ini=dlogD_dloga(a);

	//double kf = twopie*lenfac/(64.0);
	
	space_mpc_to_dimless=1.0;
	if(my_corank==0)
	 printf("\nbox_length %lf\n",boxlength);
        dx[0] = boxlength/((double)(ind[0]-1));	dx[1] = boxlength/((double)(ind[1]-1));	
	dx[2] = boxlength/((double)(ind[2]-1));
	L[0] = boxlength*space_mpc_to_dimless;	L[1] = boxlength*space_mpc_to_dimless;	L[2] = boxlength*space_mpc_to_dimless;
	dk = 1.0/(boxlength*space_mpc_to_dimless);
	kbins = 0;
	
	k_nyq = 0.5/(dx[0]);
	
////////////////////// For MPI ////////////////////////////////
	xcntr[0] = cum_lin_ind -1;
/////////////////////////////////////////////////////////////////////////////////////////////////
	
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
			

			if((xcntr[j]%ind[j])<=(ind[j]/2))
				{
					
				 k_grid[ci][j] = ((double)(xcntr[j]%ind[j]))/L[j];

				  ktmp+= k_grid[ci][j]*k_grid[ci][j];

					
				}
			else
				{ 
				 k_grid[ci][j] = ((double)((xcntr[j]%ind[j])-ind[j]))/L[j];

				 ktmp+= k_grid[ci][j]*k_grid[ci][j];

				
				  
				}
		
			 
			x_grid[ci][j] = ((double)(xcntr[j]%ind[j]))*dx[j];
			
		
		}
		
		
		

		
		
	
			
		if(ktmp>maxkmagsqr)
		maxkmagsqr = (ktmp);
		if((ktmp>=0.0)&&(minkmagsqr>ktmp))
		minkmagsqr = ktmp;
		

		kbin_grid[ci] = (int)(sqrt(ktmp)/(dk));
		//kmag_grid[ci] = sqrt(ktmp);
		// printf("yo  %d  %lf\n",kbin_grid[ci],sqrt(ktmp));
		kbin_count[kbin_grid[ci]] = kbin_count[kbin_grid[ci]]+1;

		if(kbin_grid[ci]>kbins)
		kbins=kbin_grid[ci];
		
			

		
		
		
      	}	

	
	//ini_rand_field_2(ind,k_grid,kbin_grid,kbins, dk,
	//					ini_dc,ini_theta,a,a0,a_t,pk);
	//ini_rand_field(ind,kmag_grid,ini_dc,ini_theta,a,a0,a_t,pk);
	
	if(fini_name=="None")
	grf.gen(k_grid,ini_dc, pk,a_t,a,a0,f_ini);
	else
	read_dc_from_hdf5(fini_name,ini_dc,ini_theta,ind_loc, cum_lin_ind);
	
	
	for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			

		
			ci = (n[2]*n[1])*i + n[2]*j + k;
			loc_ind[0] = i;  loc_ind[1] = j;  loc_ind[2] = k;

			 ini_dc[ci] = 0.0;

			 poisson_rhs = 1.5*omega_dm_0*H0*H0*pow(a0/ai,3.0*(1.0+w))*ini_dc[ci];

	
			poisson_phi.update_4pieGpsi(ci,poisson_rhs);
			


		    }
		  }

		}

	
	
	poisson_phi.solve_poisson(k_grid, ai, Hi);
	int chk; double fa[2];
	int cnt2;
	for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			cnt2 = (n[2]*n[1])*(i+cum_lin_ind) + n[2]*j + k;
			loc_ind[0] = i;  loc_ind[1] = j;  loc_ind[2] = k;
			potn = poisson_phi.get_potential(ci);
			

			
			
			psi_amp = sqrt(3.0*H0*H0*omega_dm_0*a0*a0*a0*(1.0+ini_dc[ci]));
			psi_val[0] = psi_amp*cos(ini_theta[ci]);
			psi_val[1] = psi_amp*sin(ini_theta[ci]);
			
			
			psi.update_fdm(ci,psi_val);
			psi.update_amp2_value(loc_ind);

			phi.update_value(loc_ind, potn);
			
			

			
			//potn = phi.get_potential(loc_ind);
			
			if(fabs(potn)>=max_potn)
				max_potn = fabs(potn);

			

			fprintf(fpstoreini,"%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.22lf\t%.15lf\n",
							a,dx[0]*i,dx[1]*j,dx[2]*k,fa[0],fa[1],potn,ini_dc[ci]);

	            }
		  }

		}

	

	
	
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, cart_comm, info);

	file = H5Fcreate("test_initial.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);
	//initial_hdf5_write_mpi(ind,ind_loc, psi, phi,file,dc,k_grid,x_grid,a3a03omega, ai,cum_lin_ind,true);

	H5Fclose(file);

	
	
	

	
	printf("df dk %lf\n",dk);
/*
	FILE *fpwr_spec = fopen("spec_test.txt","w");
	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	cal_spectrum(ini_dc,kbin_grid, kbins,n,pwr_spec, dk,a/ai,fpwr_spec);
	fclose(fpwr_spec);

	
	
	
	vmax_cap=res_limits(max_potn, vmax, dx[0],ai,dt_limit, len_res);
	
*/
	


	fprintf(fp_sim_info,"\n######## Ini Info   ##########\n");
	fprintf(fp_sim_info,"Info from cart rank %d\n",my_corank);
	fprintf(fp_sim_info,"	L is %lf\n",L[0]);
	fprintf(fp_sim_info,"	dx is %lf req len res %lf\n",dx[0],len_res);
	fprintf(fp_sim_info,"\nK details:\n	dk is %lf  per MPc\n",dk/lenfac);
	fprintf(fp_sim_info,"	kmin %lf kmax %lf\n",dk,sqrt(maxkmagsqr));
	fprintf(fp_sim_info,"	k_nyquist is %.5lf\n",k_nyq);

	fprintf(fp_sim_info,"	kbins is %d   %d\n",kbins,(int)((sqrt(maxkmagsqr)-sqrt(minkmagsqr))/dk));

	
	
	/*printf("\n Length details\n");
	
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
	printf("\nInitialization Complete from rank %d.\n",my_corank);
	MPI_Barrier(cart_comm);

      delete[] kbin_count ;
      delete[] ini_dc;
      delete[] ini_theta;
      delete[] pwr_spec;
      delete[] dc ;

	fclose(fpstoreini);
	
	  

}


void initial_hdf5_write_mpi(int *ind,int *ind_loc,fdm_poisson_mpi psi,metric_potential_poisson_mpi phi,hid_t filename,
								double *dc,double k_grid[][3],double x_grid[][3],
													double a3a03omega,double a,int cum_lin_ind,bool get_dc=true)
{	
	herr_t status_psi,status_phi,status;	
	hid_t file,dtype,dspace_glbl_psi,dspace_glbl_potn,dspace_glbl_grid,dspace,dset_glbl_k,dset_glbl_x,plist_id;
	hsize_t dim[3],odim[2],gdim[2];
	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];
	odim[0] = ind[0]*ind[1]*ind[2];
	odim[1] = 2;

	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);
/*	dspace_glbl_psi = H5Screate_simple(3, dim, NULL);
	dspace_glbl_potn = H5Screate_simple(2, odim, NULL);

	

	status_psi = psi.write_hdf5_psi_mpi(filename, dtype, dspace_glbl_psi,dc,a3a03omega,a,cum_lin_ind,get_dc);

	H5Sclose(dspace_glbl_psi);

	status_phi = phi.write_hdf5_potn_mpi(filename, dtype, dspace_glbl_potn);

	H5Sclose(dspace_glbl_potn);
*/

	gdim[0] = ind[0]*ind[1]*ind[2];
	gdim[1] = 3;

	dspace_glbl_grid = H5Screate_simple(2, gdim, NULL);

	

	dset_glbl_k = H5Dcreate(filename, "k_grid", dtype, dspace_glbl_grid,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t count[2],offset[2];
	count[0] = ind_loc[0]*ind_loc[1]*ind_loc[2]; count[1] = gdim[1]; offset[0] = cum_lin_ind*ind_loc[1]*ind_loc[2]; offset[1] = 0;

	dspace = H5Screate_simple(2, count, NULL);

	

	H5Sselect_hyperslab(dspace_glbl_grid, H5S_SELECT_SET, offset, NULL, count, NULL);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	

	status = H5Dwrite(dset_glbl_k, dtype, dspace, dspace_glbl_grid,
		      				plist_id, k_grid);

	H5Dclose(dset_glbl_k);
	H5Pclose(plist_id);
	H5Sclose(dspace);
	

	dset_glbl_x = H5Dcreate(filename, "x_grid", dtype, dspace_glbl_grid,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dspace = H5Screate_simple(2, count, NULL);

	H5Sselect_hyperslab(dspace_glbl_grid, H5S_SELECT_SET, offset, NULL, count, NULL);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dwrite(dset_glbl_x, dtype, dspace, dspace_glbl_grid,
		      				plist_id,x_grid);

	H5Dclose(dset_glbl_x);
	H5Pclose(plist_id);
	H5Sclose(dspace);
	H5Sclose(dspace_glbl_grid);



	H5Tclose(dtype);




}





