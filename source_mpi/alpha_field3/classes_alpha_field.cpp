


class scalar_field_3d_mpi
{

	protected:
	 double ***f,***fup,***f_t,****f_x,***f_lap;
         int n[3],cum_lin_ind;
         char ini_status{'N'};


	public:

	scalar_field_3d_mpi(int *n_arr,int cum_lin_ind_arr,bool need_lap=false,bool need_space_grads=false)
	{
	   int i,j,nx;
	   cum_lin_ind = cum_lin_ind_arr;
	   n[0] = n_arr[0];	n[1] = n_arr[1];	n[2] = n_arr[2];
	   ini_status = 'I';
	   nx = n[0]+4;
	   f = new double** [nx] ;
	   fup = new double** [nx] ;
	   //f_t = new double** [nx] ;
	   if(need_lap)
	   f_lap = new double** [n[0]] ;
	   if(need_space_grads)
	   { f_x = new double *** [3];
	     f_x[0] = new double** [n[0]];
	     f_x[1] = new double** [n[0]];
	     f_x[2] = new double** [n[0]];  
	   }




	double *f_pool;
	double *fup_pool;
	//double *f_t_pool;
	double *f_lap_pool;
	double *f_x_pool, *f_y_pool, *f_z_pool;



	f_pool = new double [nx*n[1]*n[2]];
	fup_pool = new double [nx*n[1]*n[2]];
	//f_t_pool = new double [n[1]*n[2]];

			
	if(need_space_grads)
	 { f_x_pool = new double [n[0]*n[1]*n[2]];
	   f_y_pool = new double [n[0]*n[1]*n[2]];
	   f_z_pool = new double [n[0]*n[1]*n[2]];
	}

	if(need_lap)
	 { f_lap_pool = new double [n[0]*n[1]*n[2]];
		    
	 }



	for(i=0;i<(nx);++i)
	   {
		f[i] = new double* [n[1]];
		fup[i] = new double* [n[1]];
		//f_t[i] = new double* [n[1]];

		

		
		 

		if((need_lap)&&(i<n[0]))
		{ f_lap[i] = new double* [n[1]];
		   
		}

		if((need_space_grads)&&(i<n[0]))
		{ f_x[0][i] = new double* [n[1]];
		  f_x[1][i] = new double* [n[1]];
		  f_x[2][i] = new double* [n[1]];

		  
		}


		
		



		for(int j=0;j<n[1];++j)
	     	{
		  
		 
		  f[i][j] = f_pool;
		  fup[i][j] = fup_pool;
		  //f_t[i][j] = f_t_pool;

		  f_pool+=n[2];
		  fup_pool+=n[2];
		  //f_t_pool+=n[2];
		  

		   if((need_lap)&&(i<n[0]))
		   { f_lap[i][j] = f_lap_pool;
	   	     f_lap_pool+=n[2];
		    }
		   if((need_space_grads)&&(i<n[0]))
		   { f_x[0][i][j] = f_x_pool ;
		     f_x[1][i][j] = f_y_pool ;
		     f_x[2][i][j] = f_z_pool ;

		     f_x_pool+=n[2];
		     f_y_pool+=n[2];
		     f_z_pool+=n[2];
		   }
	

		 }

	    }



	if((need_lap)||(need_space_grads))
	  {	
		if((need_lap)&&(need_space_grads))
 		 fprintf(fp_sim_info,"Field allocated with arrays for space der and laplacian %d %d %d\n",n[0],n[1],n[2]);
		else
		  if(need_lap)
			fprintf(fp_sim_info,"Field allocated with array for laplacian %d %d %d %d\n",n[0],n[1],n[2],nx);	
		  else
			fprintf(fp_sim_info,"Field allocated with space ders\n");

	 }
	else
	   fprintf(fp_sim_info,"Field allocated withOUT arrays for space der and laplacian\n");
		

	}


	double get_field(int *ind,code1 c)
	{
		double h;
		switch(c){
				case give_f:
					h = f[ind[0]+2][ind[1]][ind[2]];
					return(h);
				case give_f_t:
					cout<<"Reddd\n";
					return(f_t[ind[0]+2][ind[1]][ind[2]]);
				case give_f_x:
					cout<<"Green\n";
					return(f_x[0][ind[0]][ind[1]][ind[2]]);
				case give_f_y:
					return(f_x[1][ind[0]][ind[1]][ind[2]]);
				case give_f_z:
					return(f_x[2][ind[0]][ind[1]][ind[2]]);
				case give_f_lap:
					return(f_lap[ind[0]][ind[1]][ind[2]]);
				default:
					return(f[ind[0]+2][ind[1]][ind[2]]);

			}
		
	}

	int get_field_spt_der(int *ind,double der[3])
	{
		
		der[0] = f_x[0][ind[0]][ind[1]][ind[2]];
		der[1] = f_x[1][ind[0]][ind[1]][ind[2]];
		der[2] = f_x[2][ind[0]][ind[1]][ind[2]];

		if(isnan(der[0]+der[2]+der[2]))
			return(-1);
		else
			return(1);
		
		
	}
		

	void test_ind()
	{
		int i,j,k,l[3];
		double h;
		printf("ff %d %d %d\n",n[0],n[1],n[2]);
		FILE *fp = fopen("ind_check.txt","w");
		
		for(i=-2;i<(n[0]+2);++i)
		{
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<(n[2]);++k)
				{
					l[0] = i; l[1] = j;l[2] = k;
					h = get_field(l,give_f);//f[l[0]+2][l[1]][l[2]];
					fprintf(fp,"h ijk %.10lf %d %d %d\n",h,i,j,k);

				}

			}



		}
		printf("\ntest_ind SUCCESS\n");
		fclose(fp);

	}

	int cal_spt_grads(int *ind,double *dx,bool laplacian=false,bool spt_grad=false)
	{
	   int i,j;
	   double m[3],lapsum=0.0;
	   int ind_lw[3]{ind[0]+2,ind[1],ind[2]};
	   int ind_l1[3]{ind_lw[0],ind_lw[1],ind_lw[2]},ind_l2[3]{ind_lw[0],ind_lw[1],ind_lw[2]},ind_r1[3]{ind_lw[0],ind_lw[1],ind_lw[2]},
												ind_r2[3]{ind_lw[0],ind_lw[1],ind_lw[2]};
	   for(i=0;i<3;++i)
	   {
	     
	  	if(i!=0)
	     	{ind_l1[i] = (n[i]+ind_lw[i]-1)%n[i];
	     	 ind_l2[i] = (n[i]+ind_lw[i]-2)%n[i];

	     	 ind_r1[i] = (ind_lw[i]+1)%n[i];
	     	 ind_r2[i] = (ind_lw[i]+2)%n[i];
	      	}
	        else
		{ind_l1[i] = (ind_lw[i]-1);
	     	 ind_l2[i] = (ind_lw[i]-2);

	     	 ind_r1[i] = (ind_lw[i]+1);
	     	 ind_r2[i] = (ind_lw[i]+2);
	      	}
	    
	     if(spt_grad==true)
	     {
	       m[0] = (-8.0*f[ind_l1[0]][ind_l1[1]][ind_l1[2]]+8.0*f[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	       m[1] = (f[ind_l2[0]][ind_l2[1]][ind_l2[2]]-f[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	       m[2] = m[0] + m[1] ;
		
		
	       f_x[i][ind[0]][ind[1]][ind[2]] = m[2]/(12.0*dx[i]);
	      }

	     if(laplacian==true)
	      {	m[0] = (16.0*f[ind_l1[0]][ind_l1[1]][ind_l1[2]]+16.0*f[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	      	m[1] = (-f[ind_l2[0]][ind_l2[1]][ind_l2[2]]-f[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	      	m[2] = m[0] + m[1] -30.0*f[ind_lw[0]][ind_lw[1]][ind_lw[2]];
	     	lapsum+= (m[2]/(12.0*dx[i]*dx[i]));
		//printf("\ni %d j %d k %d potn_val %.10lf potn_lap  %.15lf\n",ind_lw[0],ind_lw[1],ind_lw[2],m[1],f[ind_r1[0]][ind_r1[1]][ind_r1[2]]);
		//printf("left1 i %d j %d k %d  %.10lf\n",ind_l1[0],ind_l1[1],ind_l1[2],f[65][62][0]);
	     	
	      }
		ind_l1[i] = ind_lw[i];
	    	ind_l2[i] = ind_lw[i];

	     	ind_r1[i] = ind_lw[i];
	     	ind_r2[i] = ind_lw[i];

	   } 

	if(laplacian==true)
	  f_lap[ind[0]][ind[1]][ind[2]]=lapsum;

	  if(isnan(lapsum))
		return(0);
	  else
		return (1);
	
	}


	void switch_fs()
	{
		int i,j,k;
		for(i=0;i<(n[0]+4);++i)
		{
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<(n[2]);++k)
				{
					f[i][j][k] =  fup[i][j][k] ;

				}

			}



		}



	}

	int update_field(int * ind,double fu,int upd=1,double f_tu=0.0)
	{
		
		if(upd)
		fup[ind[0]+2][ind[1]][ind[2]] = fu;
		else
		f[ind[0]+2][ind[1]][ind[2]] = fu;
		//f_t[ind[0]][ind[1]][ind[2]] = f_tu;

		if(isnan(fu+f_tu))
			return(0);
		else
			return(1);

	}

	

	herr_t write_hdf5_mpi(hid_t filename,hid_t dtype,hid_t dset_glbl)
	{
		hid_t plist_id;
		hid_t dspace,dspace_glbl;
		herr_t status ;
		hsize_t count[3],offset[3];
		count[0] = n[0]; count[1] = n[1]; count[2] = n[2];
		offset[0] = cum_lin_ind; offset[1] = 0; offset[2] = 0;

		dspace = H5Screate_simple(3, count, NULL);
		

		dspace_glbl = H5Dget_space(dset_glbl);


		H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		plist_id = H5Pcreate(H5P_DATASET_XFER);
    		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		status = H5Dwrite(dset_glbl, dtype, dspace, dspace_glbl,
		      				plist_id, &f[2][0][0]);

	
    		H5Sclose(dspace);
    		H5Sclose(dspace_glbl);
    		H5Pclose(plist_id);
  
 

		
		return status;

	}

	int mpi_send_recv()
	{
		int mpi_check = 0;	
		MPI_Request send_req,recv_req;
	
	
		int left_rank,right_rank,my_cart_rank;
		MPI_Status status;
		mpi_check=MPI_Cart_shift(cart_comm,0,+1,&my_cart_rank,&right_rank);
		mpi_check=MPI_Cart_shift(cart_comm,0,-1,&my_cart_rank,&left_rank);

		int sendtag = 1,recvtag = 1;


		mpi_check = MPI_Isend(&f[2][0][0],2, c_x_plain, left_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&f[n[0]+2][0][0],2, c_x_plain,right_rank , recvtag,cart_comm, &recv_req);

		//mpi_check = MPI_Isend(&f[3][0][0],1, c_x_plain, left_rank,sendtag, cart_comm,&send_req);
		//mpi_check = MPI_Irecv(&f[n[0]+3][0][0],1, c_x_plain,right_rank , recvtag,cart_comm, &recv_req);

		

		


		mpi_check = MPI_Isend(&f[n[0]][0][0],2, c_x_plain, right_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&f[0][0][0], 2, c_x_plain,left_rank , recvtag,cart_comm, &recv_req);

		//mpi_check = MPI_Isend(&f[n[0]+1][0][0],1, c_x_plain, right_rank,sendtag, cart_comm,&send_req);
		//mpi_check = MPI_Irecv(&f[1][0][0], 1, c_x_plain,left_rank , recvtag,cart_comm, &recv_req);
	




	/*	mpi_check = MPI_Sendrecv_replace(&f[0][0][0], 2, c_x_plain,
                         left_rank, sendtag, left_rank, recvtag,
                         cart_comm, &status);

		

		mpi_check = MPI_Sendrecv_replace(&f[n[0]+2][0][0], 2, c_x_plain,
                         right_rank, sendtag, right_rank, recvtag,
                        cart_comm, &status);
*/
		//printf("sending done for cart rank %d %d\n",my_cart_rank,n[0]);
		
		return(mpi_check);
	}


};


class param_cosmo_sim
{
	
	public:
	double omega_dm_0,loc_h,z_ini,a0;
	double box_length;
	int box_n;

	double loc_c_box,
		loc_pc_box,
		loc_hbar_box , 
		loc_lenfac,loc_space_mpc_to_dimless;


	
	void load_default_unit_const()	
	{
		loc_c_box = 2.99;
		loc_pc_box = 3.0857;
		loc_hbar_box = 6.582119569;//eVs 
		loc_lenfac = 1.0;




	}



	

	void load_default_cosmo_sim()
	{
		load_default_unit_const();
		a0 = 1.0;
		omega_dm_0 = 0.29;
		loc_h = 0.7;
		//hbar_by_m = hbar_box*c_box*(1e-8)/(alpha*pc_box);	
		loc_space_mpc_to_dimless = 0.001/loc_c_box; ////	\tilde{x} (dimensionless) = physical{x (In Mpc)}*space_mpc_to_dimless  

		box_length = 1.0;
		box_n = 128;
		z_ini = 99.0;



	}


	void print_to_file_cosmo_sim(FILE *fp)
	{

		fprintf(fp,"\n######## Cosmo info: ##########\n");
		fprintf(fp,"h = %lf\n",loc_h);
		fprintf(fp,"omega_dm_0 = %lf\n",omega_dm_0);
		
		fprintf(fp,"\n######## Unit boxes: ##########\n");
		fprintf(fp,"c_box = %lf\n pc_box = %lf\nhbar_box = %lf\nlenfac = %lf\n",loc_c_box,loc_pc_box,loc_hbar_box,loc_lenfac);
		fprintf(fp,"\n######## Box details: #########\n");
		fprintf(fp,"box_length = %lf\n",box_length);
		fprintf(fp,"box_N = %d\n",box_n);



	}
	



};


class param_alpha: public param_cosmo_sim
{
	public:
	double loc_hbar_by_m22, loc_alpha;
	string fini_dc;
	void load_defaults()
	{
		load_default_cosmo_sim();

		loc_alpha = 1e8;
		fini_dc = "None";	
	

	}

	void print_to_file(FILE *fp)
	{


		print_to_file_cosmo_sim(fp);
		
		fprintf(fp,"\n######## alpha info: ##########\n");
		fprintf(fp,"alpha = %lf\n",loc_alpha);
		fprintf(fp,"Mfield = %lf\n",Mfield);
		


	}

	
};





class linear_poisson_field_mpi
{

	private:
	double k[3],delta_k[3],delta_a_k[3];
	public:
	linear_poisson_field_mpi(double ki[3],double delta_ki,double delta_a_ki)
	{

		k[0] = ki[0]; k[1] = ki[1]; k[2] = ki[2];
		//k[0] = 0.5*h;
		delta_k[0] = delta_ki; delta_k[1] = delta_ki; delta_k[2] = delta_ki;
		delta_a_k[0] = delta_a_ki; delta_a_k[1] = delta_a_ki; delta_a_k[2] = delta_a_ki;
	}
	
	void evolve(double ai,double da,double omega_dm_0,double H0,double a0=1.0)
	{
		int i;
		double a,z,a_t,ak,kdelta_k[3],kdelta_a_k[3], acc1[3], acc2[3],alpha_lin,beta_lin,omega,A;
		double ini_delta_k[3],ini_delta_a_k[3];
		FILE *fp_lin[3];
		//da=0.01*da;
		for(i=0;i<3;++i)
		{
			ini_delta_k[i] = delta_k[i];
			ini_delta_a_k[i] = delta_a_k[i];
		}

		fp_lin[0] = fopen("linear_min.txt","w");
		fp_lin[1] = fopen("linear_mid.txt","w");
		fp_lin[2] = fopen("linear_max.txt","w");

		fprintf(fp_sim_info,"\n//#### LINEAR Calc...####\n");
		fprintf(fp_sim_info,"k_min %lf (h/Mpc)\n",k[0]*h);
		fprintf(fp_sim_info,"k_mid %lf (h/Mpc)\n",k[1]*h);
		fprintf(fp_sim_info,"k_max %lf (h/Mpc)\n",k[2]*h);
		

		for(a=ai;a<=a0;a+=da)
		{
			a_t = a*H0*sqrt(omega_dm_0*pow(a0/a,3.0*(1.0+w))+ (1.0-omega_dm_0));
			omega = omega_dm_0*pow(a/a0,-3.0*(1.0+w))*H0*H0*a*a/(a_t*a_t);
			A = -1.5*(a/a0)*(1.0+w)*omega_dm_0*pow(a/a0,-3.0*(1.0+w)-1.0)/(omega_dm_0*pow(a/a0,-3.0*(1.0+w))+(1.0-omega_dm_0));
			z = a0/a -1.0;
			
			beta_lin = -3.0/a - (A - 3.0*(2.0*w-cs2))/a;
			ak = a + da;
	
			for(i=0;i<3;++i)
			{
			  
			  fprintf(fp_lin[i],"%lf\t%lf\t%lf\t%lf\t%lf\n",a,a/ai,z,delta_k[i],delta_k[i]/ini_delta_k[i]);
				
			  alpha_lin = 1.5*omega*(1.0-6.0*cs2+8.0*w-3.0*w*w)/(a*a) - k[i]*k[i]*cs2/(a_t*a_t*a*a) ;
			
			  acc1[i] = alpha_lin*delta_k[i] + beta_lin*delta_a_k[i];
			 if(a==ai)
			   printf("for i %d k is %lf  (h/Mpc)\n",i,k[i]/h);
			// printf("i %d acc1 %.10lf  %.10lf  %.10lf\n",i,k[i],cs2,beta_lin);

			  kdelta_a_k[i] = delta_a_k[i] + da*acc1[i];

			  delta_k[i] = (delta_k[i]*(1.0 + 0.25*alpha_lin*da*da ) + (da + 0.5*beta_lin*da*da )*delta_a_k[i])/(1.0 - 0.25*alpha_lin*da*da);

			}

			a_t = ak*H0*sqrt(omega_dm_0*pow(a0/ak,3.0*(1.0+w))+ (1.0-omega_dm_0));
			omega = omega_dm_0*pow(a/a0,-3.0*(1.0+w))*H0*H0*ak*ak/(a_t*a_t);
			A = -1.5*(a/a0)*(1.0+w)*omega_dm_0*pow(ak/a0,-3.0*(1.0+w)-1.0)/(omega_dm_0*pow(ak/a0,-3.0*(1.0+w))+(1.0-omega_dm_0));
			
			beta_lin = -3.0/ak - (A - 3.0*(2.0*w-cs2))/ak;
	
			for(i=0;i<3;++i)
			{
				
			  alpha_lin = 1.5*omega*(1.0-6.0*cs2+8.0*w-3.0*w*w)/(ak*ak) - k[i]*k[i]*cs2/(a_t*a_t*ak*ak) ;
			
			  acc2[i] = alpha_lin*delta_k[i] + beta_lin*kdelta_a_k[i];

			  delta_a_k[i] = delta_a_k[i] + 0.5*da*(acc1[i]+acc1[i]);

			  

			}
			

		}


	}


};




class metric_potential_poisson_mpi
{

	private:
	int n[3];
	int n_loc[3];
	int cum_lin_ind;
	int potential;
	scalar_field_3d_mpi phi_val;

	double *pot_p;

	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	
	ptrdiff_t alloc_local, local_n0, local_0_start;
	
	public:
	
	metric_potential_poisson_mpi(int *ind,int *ind_loc,int cum_lin_ind_ar,int ptn,bool lb=false,bool sgb=false):phi_val(ind_loc,cum_lin_ind_ar,lb,sgb)
	{
		int l = ind[0]*ind[1]*ind[2];
		n[0]=ind[0];n[1]=ind[1];n[2]=ind[2];
		n_loc[0]=ind_loc[0];n_loc[1]=ind_loc[1];n_loc[2]=ind_loc[2];

		cum_lin_ind = cum_lin_ind_ar;


		pot_p =  new double [ind_loc[0]*ind_loc[1]*ind_loc[2]];
		potential = ptn;

		const ptrdiff_t n0 = n[0];
		const ptrdiff_t n1 = n[1];
		const ptrdiff_t n2 = n[2];
	
		alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
                                 cart_comm,
                                 &local_n0, &local_0_start);
		
		fpGpsi = fftw_alloc_complex(alloc_local);
		fpGpsi_ft = fftw_alloc_complex(alloc_local);


		plan_pois_f = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               fpGpsi, fpGpsi_ft,
                              cart_comm, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_pois_b = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               fpGpsi_ft, fpGpsi,
                              cart_comm, FFTW_BACKWARD, FFTW_ESTIMATE);

	}


	void solve_poisson(double k_grid[][3],double Xb,double a,double a_t,double da)
	{
		int i,j,k,ci,ind[3]{0,0,0},r;
		double k2fac;
		fftw_execute(plan_pois_f);
		double sqrt_tN = sqrt((double)(n[0]*n[1]*n[2])); 
		double dtN = (double)(n[0]*n[1]*n[2]);

		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
			k2fac = twopie*twopie*(k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);
			
			//if(k2fac>0.0)
			{
			  if(potential)	
			  {	fpGpsi_ft[ci][0] = fpGpsi_ft[ci][0]/(1.0+da*k2fac/(3.0*a_t*a_t*a));
			 	fpGpsi_ft[ci][1] = fpGpsi_ft[ci][1]/(1.0+da*k2fac/(3.0*a_t*a_t*a));
			  }

			  else
			  {
				fpGpsi_ft[ci][0] = fpGpsi_ft[ci][0]/(1.0+da*k2fac/((2.0*alpha-1.0)*a_t*a*a));
			 	fpGpsi_ft[ci][1] = fpGpsi_ft[ci][1]/(1.0+da*k2fac/((2.0*alpha-1.0)*a_t*a*a));

			  }
		
			 fpGpsi_ft[ci][0] = fpGpsi_ft[ci][0]/(dtN);
			 fpGpsi_ft[ci][1] = fpGpsi_ft[ci][1]/(dtN);

			 
				
			}

			if(k2fac<=0.0)
			{

				 fpGpsi_ft[ci][0] = Xb;///(dtN);
			 	 fpGpsi_ft[ci][1] = 0.0;///(dtN);


				if(potential)
			      {
				fpGpsi_ft[ci][0] = 0.0;
		        	fpGpsi_ft[ci][1] = 0.0;


			      }

			}

	
		/*	else
			{
                   	    if(potential)
			    {
				fpGpsi_ft[ci][0] = 0.0;
		        	fpGpsi_ft[ci][1] = 0.0;


			    }
			    else
			    {   fpGpsi_ft[ci][0] = Xb;
		        	fpGpsi_ft[ci][1] = 0.0;

			    }
			 	
			}
		*/
			

		    }

		  }

		}
		
		fftw_execute(plan_pois_b);
	

	}


	
	int calc_vel(int * ind,double &potn_vel,double f_t,double potn,double potn_a,double a,double a_t,double a_tt,double *dx,double omega_dm_0,double Xb)
	{
		int ci;		
		
		
		
		ci =  (n_loc[2]*n_loc[1])*ind[0] + n_loc[2]*ind[1] + ind[2];
		
		//potn_val = fpGpsi[ci][0] ;
		
		
			
		//potn_vel = potn_approx_vel(a,a_t,potn_val,lap_potn,f_t,omega_de_0);
		//printf("lap_potn %lf %lf %lf\n",lap_potn,potn_val,f_t);
	
	
		if(potential)
		potn_vel = potn_vel_eqn(a,a_t,potn,f_t,omega_dm_0,Xb);
		else
		{ potn_vel = field_acc_eqn(f_t,potn,potn_a,a,a_t,a_tt);
		  //if(ci==134)
			//printf("potn_vel %.10lf\n",potn_vel);	
		}
		
		

		
		
		if(isnan(potn_vel))
			return (-1);
		else
			return(1);
	}

	double get_potential(int ci,int get_imag=0)
	{
		if(!get_imag)
		return (fpGpsi[ci][0]);	
		else
		return(fpGpsi_ft[ci][0]);

	}

	
	double cal_X_4vel(int * ind,double a,double a_t,double phi)
	{
		double fa_val,fa_t_val,s_der[3],X;	
		int c1;	
		
		c1 =  (n_loc[2]*n_loc[1])*ind[0] + n_loc[2]*ind[1] + ind[2];
		fa_t_val = get_potential(c1);
		fa_t_val = fa_t_val*a_t;

		//c1 = f_alpha.get_field_spt_der(ind,s_der);
		
		
		X = fa_t_val*fa_t_val;///(1.0+2.0*phi);//  - (s_der[0]*s_der[0]+s_der[1]*s_der[1]+s_der[2]*s_der[2])/(a*a*(1.0-2.0*phi));
		X = 0.5*X;

		return(X);

	}
	
	void update_4pieGpsi(int ci,double val)
	{
		
		fpGpsi[ci][0] = val;
		fpGpsi[ci][1] = 0.0;
	
	}


	


	void update_value(int *indi,double phiv)
	{

		phi_val.update_field(indi,phiv,0);


	}

	double get_value(int *indi)
	{
		double val;
		val = phi_val.get_field(indi,give_f);

		return val;


	}


	
	



	herr_t write_hdf5_values_mpi(hid_t filename,hid_t dtype,hid_t dspace_glbl,hid_t dspace_dc_glbl,double *dc,double a0,double a,double a_t,double Xb,
							metric_potential_poisson_mpi phi,int cum_lin_ind,bool get_dc=false)
	{

		hid_t dataset,dset_glbl,dspace,plist_id;
		herr_t status ;
		int tN  = n[0]*n[1]*n[3];
		int i,j,k,locind[3],ci;
		double fa_val,fa_t_val,x4val,rho_fa,phival;

		

		if(potential)
		{ 
		  
		  dset_glbl = H5Dcreate(filename, "potential", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		  status = phi_val.write_hdf5_mpi( filename,dtype,dset_glbl);

		  H5Dclose(dset_glbl);


		}

	
		else
		{ dset_glbl = H5Dcreate(filename, "field_a", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		  status = phi_val.write_hdf5_mpi( filename,dtype,dset_glbl);

		  H5Dclose(dset_glbl);

		 if(get_dc)
		 {	

		 
			
		  for(i=0;i<n_loc[0];++i)
		   {
			for(j=0;j<n_loc[1];++j)
			{
				for(k=0;k<n_loc[2];++k)
				{	ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
					
					phival = phi.get_value(locind);	
					

					x4val = cal_X_4vel(locind,a,a_t,phival);	
							

					//rho_fa = x4val*(3.0*H0*H0/(4.0*a3a03omega*twopie*Xb_0));
					//Xb = Xb_0*pow(a0/a,6.0/(2.0*alpha-1.0));
					
					
					dc[ci] = (x4val/Xb)-1.0;
					//if(ci<10)
					//printf("ci %d dc %.10lf  %.10lf %.10lf  %.10lf\n",ci,dc[ci],x4val,Xb,phival);
					 

				}
	
			}

		  }
		
		

		  dataset = H5Dcreate(filename, "dc_alpha_field", dtype, dspace_dc_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		

		  hsize_t count[1],offset[1];
	  	  count[0] = n_loc[0]*n_loc[1]*n_loc[2];
		  offset[0] = cum_lin_ind*n[1]*n[2]; 

		  dspace = H5Screate_simple(1, count, NULL);
		 // printf("POTENTIAL %d %d %d\n",potential,count[0],n[0]);
		  
		  H5Sselect_hyperslab(dspace_dc_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		  
		  plist_id = H5Pcreate(H5P_DATASET_XFER);
    		  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		  status = H5Dwrite(dataset, dtype, dspace, dspace_dc_glbl,
		      				plist_id,dc);

		  H5Sclose(dspace);
		  H5Dclose(dataset);

	         }




	
		 
	    }
		  
		
		return(status);

	}



	





};




class metric_potential_poisson_mpi_ini
{

	private:
	int n[3];
	int n_loc[3];
	int cum_lin_ind;
	//scalar_field_3d phi;

	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	
	ptrdiff_t alloc_local, local_n0, local_0_start;
	
	public:
	
	metric_potential_poisson_mpi_ini(int *ind,int *ind_loc,int cum_lin_ind_ar,bool lb=false,bool sgb=false)//:phi(ind,lb,sgb)
	{
		int l = ind[0]*ind[1]*ind[2];
		n[0]=ind[0];n[1]=ind[1];n[2]=ind[2];
		n_loc[0]=ind_loc[0];n_loc[1]=ind_loc[1];n_loc[2]=ind_loc[2];

		cum_lin_ind = cum_lin_ind_ar;


		const ptrdiff_t n0 = n[0];
		const ptrdiff_t n1 = n[1];
		const ptrdiff_t n2 = n[2];
	
		alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
                                 cart_comm,
                                 &local_n0, &local_0_start);
		
		fpGpsi = fftw_alloc_complex(alloc_local);
		fpGpsi_ft = fftw_alloc_complex(alloc_local);


		plan_pois_f = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               fpGpsi, fpGpsi_ft,
                              cart_comm, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_pois_b = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               fpGpsi_ft, fpGpsi,
                              cart_comm, FFTW_BACKWARD, FFTW_ESTIMATE);

	}


	void solve_poisson(double k_grid[][3],double a,double Hc)
	{
		int i,j,k,ci,ind[3]{0,0,0},r;
		double k2fac;
		fftw_execute(plan_pois_f);
		double sqrt_tN = sqrt((double)(n[0]*n[1]*n[2])); 
		double dtN = (double)(n[0]*n[1]*n[2]);

		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
			k2fac = twopie*twopie*(k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);
			
			if(k2fac>0.0)
			{fpGpsi_ft[ci][0] = -fpGpsi_ft[ci][0]/((k2fac/(a*a)) +3.0*Hc*Hc);
			 fpGpsi_ft[ci][1] = -fpGpsi_ft[ci][1]/((k2fac/(a*a)) +3.0*Hc*Hc);
		
			 fpGpsi_ft[ci][0] = -fpGpsi_ft[ci][0]/(dtN);
			 fpGpsi_ft[ci][1] = -fpGpsi_ft[ci][1]/(dtN);
				
			}	
			else
			{fpGpsi_ft[ci][0] = 0.0;
			 fpGpsi_ft[ci][1] = 0.0;
			}

			

		    }

		  }

		}
		
		fftw_execute(plan_pois_b);
	

	}


	void update_4pieGpsi(int ci,double val)
	{
		
		fpGpsi[ci][0] = val;
		fpGpsi[ci][1] = 0.0;
	
	}


	double get_potential(int ci,int get_imag=0)
	{
		if(!get_imag)
		return (fpGpsi[ci][0]);	
		else
		return(fpGpsi[ci][1]/sqrt(fpGpsi[ci][0]*fpGpsi[ci][0]+fpGpsi[ci][1]*fpGpsi[ci][1]));

	}

	void write_potential(FILE *fp_ptn,double *dx,double a3a03omega,double a)
	{	
		int i,j,k,ci;
		
		for(i=0;i<n_loc[0];++i)
		{
			for(j=0;j<n_loc[1];++j)
			{
				for(k=0;k<n_loc[2];++k)
				{
					 ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;
				
					  fprintf(fp_ptn,"%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n",a,dx[0]*i,dx[1]*j,dx[2]*k,fpGpsi[ci][0]);
				



				}
	
			}

		}
		
		fprintf(fp_ptn,"\n\n\n\n");

	}

/*
	herr_t write_hdf5_potn_mpi(hid_t filename,hid_t dtype,hid_t dspace_glbl)
	{

		hid_t dset_glbl;
		herr_t status;
		hid_t  dspace,plist_id;
		hsize_t     dim[2];
		dim[0] = n_loc[0]*n_loc[1]*n_loc[2]; dim[1]=2;

		
    		status = H5Tset_order(dtype, H5T_ORDER_LE);	

		dset_glbl = H5Dcreate(filename, "potential_poisson", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


		hsize_t count[2],offset[2];if(get_dc)
		{	

		 
			
		 for(i=0;i<n[0];++i)
		   {
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<n[2];++k)
				{	ci = (n[2]*n[1])*i + n[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
					
					phival = phi.get_potential(ci);	
					

					x4val = cal_X_4vel(locind,a,a_t,phival);	
							

					//rho_fa = x4val*(3.0*H0*H0/(4.0*a3a03omega*twopie*Xb_0));
					//Xb = Xb_0*pow(a0/a,6.0/(2.0*alpha-1.0));
					
					//if(ci==10)
					//dc[ci] = (x4val/Xb)-1.0;
					// printf("ci %d dc %.10lf  %.10lf %.10lf\n",ci,dc[ci],x4val,Xb);

				}
	
			}

		  }
		


		dataset = H5Dcreate(filename, "dc_alpha_field", dtype, dspace_dc_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		

		hsize_t count[1],offset[1];
		count[0] = n[0]*n[1]*n[2];
		offset[0] = cum_lin_ind*n[1]*n[2]; 

		dspace = H5Screate_simple(1, count, NULL);


		H5Sselect_hyperslab(dspace_dc_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		plist_id = H5Pcreate(H5P_DATASET_XFER);
    		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		status = H5Dwrite(dataset, dtype, dspace, dspace_dc_glbl,
		      				plist_id,dc);

		H5Sclose(dspace);
		H5Dclose(dataset);

	      }

		count[0] = dim[0]; count[1] = 2; offset[0] = cum_lin_ind; offset[1] = 0;

		dspace = H5Screate_simple(2, count, NULL);



		H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		plist_id = H5Pcreate(H5P_DATASET_XFER);
    		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		status = H5Dwrite(dset_glbl, dtype, dspace, dspace_glbl,
		      				plist_id, fpGpsi);

	
    		H5Sclose(dspace);
    		H5Pclose(plist_id);



		return status;

	}

*/

/*	int update(int * ind,double phi_val)
	{
		int c1;
		c1 = phi.update_field(ind,phi_val);
		
		return (c1);
	}

*/


};


///##################################################################################################################################################################






class ini_power_generator
{

	private:
	double *p,*k;
	const char *fname;
	double **c,max_dx=0.0,x_min;
	int point_cnt,checkspl;
	
	
	public:

	ini_power_generator(const char *name)
	{
		double a,b,pre_x=0.0;
		int i = 0,fre=1;
		fname =  name;		

		FILE *fp = fopen(fname,"r");
		while(fre!=EOF)
		{
			fre=fscanf(fp,"       %lf       %lf\n",&a,&b);

			if(i==0)
			{
			
				max_dx = a;	
				x_min = a;
				pre_x = a;	


			}

			else
			if(fre!=EOF)
			{
			   if(fabs(a-pre_x)>max_dx)
				max_dx = fabs(a-pre_x);
			   if(a<x_min)
				x_min = a;
			   pre_x = a;

			}
			
			if(fre!=EOF)
			{
			//printf("Reading %d %d %lf %.8lf\n",i,fre,a,b);
			 ++i;
			}

			
		
			
		}

		fclose(fp);
		point_cnt = i;
		p = new double[point_cnt];	
		k = new double[point_cnt];
		c = new double* [point_cnt];	
		for(i=0;i<point_cnt;++i)
			c[i] =  new double[3];	
		fp = fopen(fname,"r");
		i=0;
		fre=1;
		double cload[point_cnt][3];
		while(fre!=EOF)
		{
			fre=fscanf(fp,"%lf\t%lf\n",&a,&b);
			if(fre!=EOF)
			{*(k+i)=a;
			 *(p+i)=b;
			 //printf("ff %d %lf %.10lf\n",i,*(k+i),*(p+i));
			}
			++i;
		}
		
		fclose(fp);
		fprintf(fp_sim_info,"	doing spline..\n");
		
		checkspl = spline(k,p,point_cnt,cload); 
		if(checkspl)
		 fprintf(fp_sim_info,"\nALERT !!!! Spline fit failed....checkspl %d  %.10lf\n",checkspl,max_dx);
		else
		 fprintf(fp_sim_info,"	spline done successfully..\n");
	
		for(i=0;i<point_cnt;++i)
		{
			c[i][0] = cload[i][0];
			c[i][1] = cload[i][1];
			c[i][2] = cload[i][2];
			 

		}

		
		
	
		

	}

	double get_ini_spectrum(double x)
	{

	  double Val;
	  int cur_ind = (int)((x-x_min)/max_dx);
	  int ind_match=0;//printf("Yo %d %lf %lf\n",cur_ind,x,max_dx);
	while(cur_ind<point_cnt)
	 { if((x>=k[cur_ind])&&(x<=k[cur_ind+1]))
    	  {//printf("cur_ind is %d\n",cur_ind);
            ind_match = 1;
            Val = (p[cur_ind]+c[cur_ind][0]*(x-k[cur_ind])+c[cur_ind][1]*(x-k[cur_ind])*(x-k[cur_ind])
                                                  +c[cur_ind][2]*(x-k[cur_ind])*(x-k[cur_ind])*(x-k[cur_ind]));
          	break;
          }
	  ++cur_ind;
	  
	 }

	if(!ind_match)
	fprintf(fp_sim_info,"\nALERT index overrun for spline...for x %lf\n",x);

	return Val;
		

	}

	double test_spec(double x)
	{
		return(1.0/pow(x,2.0));
	}

	void check(double x)
	{
		int i;
		double v;
		FILE *fpcheck = fopen("spline_check.txt","w");
		for(i=0;i<point_cnt;++i)
		{
			v = get_ini_spectrum(*(k+i));
			fprintf(fpcheck,"%.10lf\t%.10lf\t%.10lf\n",*(k+i),v,*(p+i));	
			
			//printf("Check result is v %lf @ x %lf\n",v,x);
		}
		
		
		

	}



};



void cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
{	int i,j;
	int tN = s[0]*s[1]*s[2];
	double dtN = (double)(tN);
	double delta_pw;
	double kbincnt[kbins+1];

	fftw_complex *dens_cntrst; fftw_complex *Fdens_cntrst;
	fftw_plan spec_plan;

	dens_cntrst = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tN );
	Fdens_cntrst = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tN);
	
	for(i=0;i<tN;++i)
	{
		dens_cntrst[i][0] = f[i]; 
		dens_cntrst[i][1] = 0.0;
		if(i<=kbins)
		pwspctrm[i] = 0.0;
		kbincnt[kbingrid[i]] = 0.0;
	}

	
	spec_plan = fftw_plan_dft_3d(s[0],s[1],s[2], dens_cntrst, Fdens_cntrst, FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(spec_plan);
	fftw_free(dens_cntrst);

	

	for(i=0;i<tN;++i)
	{
		
		pwspctrm[kbingrid[i]]+=  (Fdens_cntrst[i][1]*Fdens_cntrst[i][1] + Fdens_cntrst[i][0]*Fdens_cntrst[i][0])/(dtN);
		++kbincnt[kbingrid[i]];
		

	}
	
	for(i=0;i<=kbins;++i)
	{

		if(kbincnt[i]!=0)
	        {  //delta_pw = sqrt(pwspctrm[i]*i*i*i*dk*dk*dk/(2.0*M_PI*M_PI*kbincnt[i]));  
		    delta_pw = 1.0/(i*dk);

		   fprintf(fspwrite,"%lf\t%lf\t%.13lf\t%.13lf\t%.13lf\n",
							abyai,i*dk,pwspctrm[i]/(kbincnt[i]),delta_pw,pwspctrm[i]/(kbincnt[i]*abyai*abyai));
		   pwspctrm[i]=pwspctrm[i]/(kbincnt[i]);

		}

	

	}

	

	
	
	fftw_free(Fdens_cntrst);
	fftw_destroy_plan(spec_plan);
	fprintf(fspwrite,"\n\n\n\n");
}


class gauss_rand_field_gen_mpi
{
	private:
	int n[3],n_loc[3];
	//scalar_field_3d phi;

	fftw_complex *field;
	fftw_complex *field_ft;



	fftw_plan plan_grf_f;
	fftw_plan plan_grf_b;
	

	ptrdiff_t alloc_local, local_n0, local_0_start;

	

	public:
	int ntN; 
	gauss_rand_field_gen_mpi(int *ind,int *ind_loc)
	{
		int l = ind[0]*ind[1]*ind[2];
		ntN = l;
		n[0]=ind[0];n[1]=ind[1];n[2]=ind[2];
		n_loc[0]=ind_loc[0];n_loc[1]=ind_loc[1];n_loc[2]=ind_loc[2];


		const ptrdiff_t n0 = n[0];
		const ptrdiff_t n1 = n[1];
		const ptrdiff_t n2 = n[2];

		alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
                                 cart_comm,
                                 &local_n0, &local_0_start);
		
		field = fftw_alloc_complex(alloc_local);
		field_ft = fftw_alloc_complex(alloc_local);

		

		plan_grf_f = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               field, field_ft, 
                              cart_comm, FFTW_FORWARD, FFTW_ESTIMATE);

		plan_grf_b = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               field_ft, field,  cart_comm, FFTW_BACKWARD, FFTW_ESTIMATE); 

		



	}
	

	double unit_normal()
	{

		double a1,a2,r;
		a1 = genrand_res53();
 		a2 = genrand_res53(); 
		r = (sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2));

		return(r);


	}

	void test_normal()
	{
		int i; double r;
		FILE *fpnorm = fopen("norm.txt","w");
		for(i=0;i<20000;++i)
		{	
			r = unit_normal();			
			fprintf(fpnorm,"%lf\n",r);
			

		}
		fclose(fpnorm);

	}


	void gen(double k_grid[][3],double *ini_dc, ini_power_generator p_k,double a_t,double a,double a0,double f_ini)
	{	int i,j,k,ci;
		double ksqr,pk_val,dtN;

		
		dtN = (double)(n[0]*n[1]*n[2]);
		

		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
						
			
			field[ci][0] = sqrt(twopie)*unit_normal();
			field[ci][1] = 0.0;
			

		    }

		  }

		}
		
		fftw_execute(plan_grf_f);

		
		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
			ksqr = (k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);

						
			if(ksqr==0.0)
			{field_ft[ci][0] = 0.0;
			 field_ft[ci][1] = 0.0;

			
			}
			else
			{
				pk_val = 0.001;//p_k.test_spec(sqrt(ksqr));				
					
				field_ft[ci][0] = sqrt(pk_val)*field_ft[ci][0]/dtN;
				field_ft[ci][1] = sqrt(pk_val)*field_ft[ci][1]/dtN;

				


			}

			 

			
			

		    }

		  }

		}
		
		fftw_execute(plan_grf_b);
		


		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
			ini_dc[ci] = field[ci][0];
			
			

		    }

		  }

		}
		
	

	fftw_free(field);
	fftw_free(field_ft);
	
	
	fftw_destroy_plan(plan_grf_f);
	fftw_destroy_plan(plan_grf_b);
	
	}







};

////##################################################

/*   if(0)
		   {	

		 
			
		   for(i=0;i<n[0];++i)
		    {
			 for(j=0;j<n[1];++j)
			 {
				for(k=0;k<n[2];++k)
				{	ci = (n[2]*n[1])*i + n[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
					
					phival = phi.get_value(locind);	
					

					x4val = cal_X_4vel(locind,a,a_t,phival);	
							

					//rho_fa = x4val*(3.0*H0*H0/(4.0*a3a03omega*twopie*Xb_0));
					//Xb = Xb_0*pow(a0/a,6.0/(2.0*alpha-1.0));
					
					//if(ci==10)
					dc[ci] = (x4val/Xb)-1.0;
					// printf("ci %d dc %.10lf  %.10lf %.10lf\n",ci,dc[ci],x4val,Xb);

				}
	
			}

		    }
		


		  dataset = H5Dcreate(filename, "dc_alpha_field", dtype, dspace_dc_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		

		hsize_t count[1],offset[1];
		count[0] = n[0]*n[1]*n[2];
		offset[0] = cum_lin_ind*n[1]*n[2]; 

		dspace = H5Screate_simple(1, count, NULL);


		H5Sselect_hyperslab(dspace_dc_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		plist_id = H5Pcreate(H5P_DATASET_XFER);
    		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		//status = H5Dwrite(dataset, dtype, dspace, dspace_dc_glbl,
		  //    				plist_id,dc);

		H5Sclose(dspace);
		H5Dclose(dataset);

	      }*/
