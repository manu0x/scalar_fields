


class scalar_field_3d_mpi
{

	protected:
	 double ***f,***f_t,****f_x,***f_lap;
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
	   //f_t = new double** [nx] ;
	   if(need_lap)
	   f_lap = new double** [n[0]] ;
	   if(need_space_grads)
	   { f_x = new double *** [3];
	     f_x[0] = new double** [n[0]];
	     f_x[1] = new double** [n[0]];
	     f_x[2] = new double** [n[0]];  
	   }

/*	   for(i=0;i<(nx);++i)
	   {
		f[i] = new double* [n[1]];
		f_t[i] = new double* [n[1]];
		if((need_lap)&&(i<n[0]))
		f_lap[i] = new double* [n[1]];
		if((need_space_grads)&&(i<n[0]))
		{ f_x[0][i] = new double* [n[1]];
		  f_x[1][i] = new double* [n[1]];
		 f_x[2][i] = new double* [n[1]];
		}
		for(int j=0;j<n[1];++j)
	     	{
		  f[i][j] = new  double[n[2]] ;
		  f_t[i][j] = new  double[n[2]] ;
		  if((need_lap)&&(i<n[0]))
		   f_lap[i][j] = new  double[n[2]] ;
		  if((need_space_grads)&&(i<n[0]))
		  { f_x[0][i][j] = new  double[n[2]] ;
		     f_x[1][i][j] = new  double[n[2]] ;
		     f_x[2][i][j] = new  double[n[2]] ;
		  }
	

		 }

	    }
*/

/*	for(i=0;i<(nx);++i)
	   {
		f[i] = new double* [n[1]];
		//f_t[i] = new double* [n[1]];

		double *f_pool;
		//double *f_t_pool;
		double *f_lap_pool;
		double *f_x_pool, *f_y_pool, *f_z_pool;

		
		 

		if((need_lap)&&(i<n[0]))
		{ f_lap[i] = new double* [n[1]];
		   
		}

		if((need_space_grads)&&(i<n[0]))
		{ f_x[0][i] = new double* [n[1]];
		  f_x[1][i] = new double* [n[1]];
		  f_x[2][i] = new double* [n[1]];

		  
		}


		f_pool = new double [n[1]*n[2]];
		  //f_t_pool = new double [n[1]*n[2]];

			
		if((need_space_grads)&&(i<n[0]))
		  { f_x_pool = new double [n[1]*n[2]];
		    f_y_pool = new double [n[1]*n[2]];
		    f_z_pool = new double [n[1]*n[2]];
		  }

		if((need_lap)&&(i<n[0]))
		  { f_lap_pool = new double [n[1]*n[2]];
		    
		  }
		



		for(int j=0;j<n[1];++j)
	     	{
		  
		 
		  f[i][j] = f_pool;
		  //f_t[i][j] = f_t_pool;

		  f_pool+=n[2];
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
*/


	double *f_pool;
	//double *f_t_pool;
	double *f_lap_pool;
	double *f_x_pool, *f_y_pool, *f_z_pool;



	f_pool = new double [nx*n[1]*n[2]];
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
		  //f_t[i][j] = f_t_pool;

		  f_pool+=n[2];
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
		for(i=0;i<(n[0]);++i)
		{
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<(n[2]);++k)
				{
					l[0] = i; l[1] = j;l[2] = k;
					h = get_field(l,give_f);//f[l[0]+2][l[1]][l[2]];
					//printf("ijk %lf %d %d %d\n",h,i,j,k);

				}

			}



		}
		printf("\ntest_ind SUCCESS\n");

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
	     
	  
	     	ind_l1[i] = (n[i]+ind_lw[i]-1)%n[i];
	     	ind_l2[i] = (n[i]+ind_lw[i]-2)%n[i];

	     	ind_r1[i] = (ind_lw[i]+1)%n[i];
	     	ind_r2[i] = (ind_lw[i]+2)%n[i];
	      
	    
	    
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

	     	
	      }
		ind_l1[i] = ind_lw[i];
	    	ind_l2[i] = ind_lw[i];

	     	ind_r1[i] = ind_lw[i];
	     	ind_r2[i] = ind_lw[i];

	   } 

	if(laplacian==true)
	  f_lap[ind[0]][ind_lw[1]][ind_lw[2]]=lapsum;

	  if(isnan(lapsum))
		return(0);
	  else
		return (1);
	
	}

	int update_field(int * ind,double fu,double f_tu=0.0)
	{
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


class param_fdm: public param_cosmo_sim
{
	public:
	double loc_hbar_by_m22, loc_alpha;
	void load_defaults()
	{
		load_default_cosmo_sim();

		loc_hbar_by_m22 = 1.0;	
	

	}

	void print_to_file(FILE *fp)
	{


		print_to_file_cosmo_sim(fp);
		
		fprintf(fp,"\n######## fdm info: ##########\n");
		fprintf(fp,"(in units of 10^-22 eV ) h_bar/m = %lf\n",loc_hbar_by_m22);

		


	}

	
};


class fdm_psi_mpi
{
	private:
	scalar_field_3d_mpi psi_r;
	scalar_field_3d_mpi psi_i;
	int n[3];
	
	public:  
	//int n[3];
	fdm_psi_mpi(int *ind,int cum_lin_ind,bool lb=false,bool sgb=false):psi_r(ind,cum_lin_ind,lb,sgb),psi_i(ind,cum_lin_ind,lb,sgb)
	{
		n[0] = ind[0];  n[1] = ind[1];  n[2] = ind[2];
	}
	
	void test_ind()
	{
		psi_r.test_ind();

	}


	void test_ind2()
	{   int locind[3],i,j,k,ci;
		double psi_r_val;
		for(i=0;i<n[0];++i)
		{
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<n[2];++k)
				{	ci = (n[2]*n[1])*i + n[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
						psi_r_val=psi_r.get_field(locind,give_f);

				}
			}
		}

		printf("\nTTTTTTTTTESSST 222 Done\n");
	

	}

	int calc_vel(int * ind,double *v,double potn,double a,double a_t,double *dx)
	{
		int c1;		
		double psi_r_lap,psi_i_lap,psi_r_val,psi_i_val;
		psi_r_val = psi_r.get_field(ind,give_f);
		psi_i_val = psi_i.get_field(ind,give_f);
		c1 = psi_r.cal_spt_grads(ind,dx,true);
		c1 = psi_i.cal_spt_grads(ind,dx,true);
		psi_r_lap = psi_r.get_field(ind,give_f_lap);
		psi_i_lap = psi_i.get_field(ind,give_f_lap);
		v[0] = -1.5*(a_t/a)*psi_r_val;
		v[0]+= ((-0.5*hbar_by_m*psi_i_lap/(a*a) + potn*psi_i_val/hbar_by_m))/H0;
		//printf("back  %lf lap %.10lf  potn  %.10lf\n", -1.5*(a_t/a)*psi_r_val,-0.5*hbar_by_m*psi_i_lap/(a*a) ,potn*psi_i_val/hbar_by_m);
		v[1] = -1.5*(a_t/a)*psi_i_val;
		v[1]+= ((0.5*hbar_by_m*psi_r_lap/(a*a) - potn*psi_r_val/hbar_by_m))/H0;
		//printf("back  %lf lap %.10lf  potn  %.10lf\n\n", -1.5*(a_t/a)*psi_i_val,0.5*hbar_by_m*psi_r_lap/(a*a),potn*psi_r_val/hbar_by_m);

		if(isnan(v[0]+v[1]))
			return (-1);
		else
			return(1);
	}

	double cal_mod(int * ind)
	{
		double psi_r_val,psi_i_val,mag;		
		
		psi_r_val = psi_r.get_field(ind,give_f);
		psi_i_val = psi_i.get_field(ind,give_f);

		
		mag = sqrt(psi_r_val*psi_r_val+psi_r_val*psi_r_val);

		return(mag);

	}

	int update(int * ind,double psir,double psii)
	{
		int c1,c2;
		c1 = psi_r.update_field(ind,psir);
		c2 = psi_i.update_field(ind,psii);
		
		return (c1*c2);
	}

	int get_psi(int *ind,double *psi_ret,code1 c = give_f)
	{
		int c1=1;		
		
		psi_ret[0]= psi_r.get_field(ind,c);
		psi_ret[1] = psi_i.get_field(ind,c);

		return c1;




	}

	int mpi_send_recv()
	{
		int mpi_check;
		mpi_check=psi_r.mpi_send_recv();
		mpi_check=psi_i.mpi_send_recv();
		return(mpi_check);


	}


	void write_psi(FILE *fp_psi,double *dc,double *dx,double a3a03omega,double a,bool get_dc=false, bool get_psi=true)
	{	//printf("a3  %.10lf\n",a3a03omega);
		int i,j,k,locind[3],ci;
		double psi_r_val,psi_i_val,psi_amp2;
		for(i=0;i<n[0];++i)
		{
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<n[2];++k)
				{	ci = (n[2]*n[1])*i + n[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
					psi_r_val=psi_r.get_field(locind,give_f);
					psi_i_val=psi_i.get_field(locind,give_f);	
							

					if(get_dc)
					{
					  psi_amp2 = psi_r_val*psi_r_val + psi_i_val*psi_i_val;		
					  dc[ci]= a3a03omega*psi_amp2 - 1.0;
					  fprintf(fp_psi,
						"%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n",
							a,dx[0]*i,dx[1]*j,dx[2]*k,psi_r_val,psi_i_val,psi_amp2,dc[ci],a3a03omega);
					}

					else
					{
					  fprintf(fp_psi,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,dx[0]*i,dx[1]*j,dx[2]*k,psi_r_val,psi_i_val);


					}


				}
	
			}

		}
		
		fprintf(fp_psi,"\n\n\n\n");

	}

	herr_t write_hdf5_psi_mpi(hid_t filename,hid_t dtype,hid_t dspace_glbl,double *dc,double a3a03omega,double a,int cum_lin_ind,bool get_dc=false)
	{

		hid_t dataset,dset_glbl_r,dset_glbl_i,dspace,plist_id;
		herr_t status ;
		int tN  = n[0]*n[1]*n[3];
		int i,j,k,locind[3],ci;
		double psi_r_val,psi_i_val,psi_amp2,h;

		

		

    /*
     * Create the dataset with default properties and close filespace.
     */
   		dset_glbl_r = H5Dcreate(filename, "psi_r", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dset_glbl_i = H5Dcreate(filename, "psi_i", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    		
		

		status = psi_r.write_hdf5_mpi( filename,dtype,dset_glbl_r);
		
		
		status = psi_i.write_hdf5_mpi( filename,dtype,dset_glbl_i);

		 FILE *test0; FILE *test1;
		if(get_dc)
		{	

		 
			
		 for(i=0;i<n[0];++i)
		   {
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<n[2];++k)
				{	ci = (n[2]*n[1])*i + n[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
					
						
					psi_r_val=psi_r.get_field(locind,give_f);
					psi_i_val=psi_i.get_field(locind,give_f);	
							

					
					 psi_amp2 = psi_r_val*psi_r_val + psi_i_val*psi_i_val;		
					  dc[ci]= a3a03omega*psi_amp2 - 1.0;
					 

				}
	
			}

		  }
		


		dataset = H5Dcreate(filename, "dc", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		

		hsize_t count[1],offset[1];
		count[0] = n[0]*n[1]*n[2];
		offset[0] = cum_lin_ind*n[1]*n[2]; 

		dspace = H5Screate_simple(1, count, NULL);


		H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		plist_id = H5Pcreate(H5P_DATASET_XFER);
    		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		//status = H5Dwrite(dataset, dtype, dspace, dspace_glbl,
		      		//		plist_id,dc);

		H5Sclose(dspace);
		H5Dclose(dataset);

	      }

		H5Dclose(dset_glbl_r);
		H5Dclose(dset_glbl_i);
		return(status);

		

	}



	void write_hdf5_psi_test()
	{

		
		int i,j,k,locind[3],ci,my_corank;
		double psi_r_val,psi_i_val,psi_amp2,h;


		//printf("gg %d %d %d\n",n[0],n[1],n[2]);	

		  MPI_Comm_rank(cart_comm,&my_corank);
			
		 for(i=0;i<n[0];++i)
		   {
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<n[2];++k)
				{	ci = (n[2]*n[1])*i + n[2]*j + k;
					locind[0]=i; locind[1]=j; locind[2]=k;
					//if(cum_lin_ind==0)
					//fprintf(test0,"%d %d %d %d\n",i,locind[0],j,k);
					//else
					//fprintf(test1,"%d %d %d %d\n",i,locind[0],j,k);
						
					psi_r_val=psi_r.get_field(locind,give_f);
					psi_i_val=psi_i.get_field(locind,give_f);	
							

					
					 psi_amp2 = psi_r_val*psi_r_val + psi_i_val*psi_i_val;		
					 h= psi_amp2 - 1.0;
					 

				}
	
			}

		  }

		//printf("bfjkbkbw TESTING DONE\n");

	}


	



};






class metric_potential_mpi
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
	
	metric_potential_mpi(int *ind,int *ind_loc,int cum_lin_ind_ar,bool lb=false,bool sgb=false)//:phi(ind,lb,sgb)
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


	void solve_poisson(fdm_psi_mpi psi,double k_grid[][3])
	{
		int i,j,k,ci,ind[3]{0,0,0},r;
		double k2fac;
		fftw_execute(plan_pois_f);
		double sqrt_tN = sqrt((double)(n[0]*n[1]*n[2])); 

		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
			k2fac = twopie*twopie*(k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);
			
			if(k2fac>0.0)
			{fpGpsi_ft[ci][0] = -fpGpsi_ft[ci][0]/(k2fac*sqrt_tN*sqrt_tN);
			 fpGpsi_ft[ci][1] = -fpGpsi_ft[ci][1]/(k2fac*sqrt_tN*sqrt_tN);
				
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


	double get_potential(int ci)
	{

		return (fpGpsi[ci][0]);	

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


	herr_t write_hdf5_potn_mpi(hid_t filename,hid_t dtype,hid_t dspace_glbl)
	{

		hid_t dset_glbl;
		herr_t status;
		hid_t  dspace,plist_id;
		hsize_t     dim[2];
		dim[0] = n_loc[0]*n_loc[1]*n_loc[2]; dim[1]=2;

		
    		status = H5Tset_order(dtype, H5T_ORDER_LE);	

		dset_glbl = H5Dcreate(filename, "potential", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


		hsize_t count[2],offset[2];
		count[0] = dim[0]; count[1] = 2; offset[0] = cum_lin_ind; offset[1] = 0;

		dspace = H5Screate_simple(2, count, NULL);



		H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
		plist_id = H5Pcreate(H5P_DATASET_XFER);
    		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		
		
		status = H5Dwrite(dset_glbl, dtype, dspace, dspace_glbl,
		      				plist_id, fpGpsi);

	
    		H5Sclose(dspace);
    		H5Pclose(plist_id);
		H5Dclose(dset_glbl);



		return status;

	}



/*	int update(int * ind,double phi_val)
	{
		int c1;
		c1 = phi.update_field(ind,phi_val);
		
		return (c1);
	}

*/


};


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
		return(1.0/pow(x,3.0));
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

	fftw_complex *theta;
	fftw_complex *theta_ft;

	fftw_plan plan_grf_f;
	fftw_plan plan_grf_b;
	fftw_plan plan_grf_theta_b;

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

		theta = fftw_alloc_complex(alloc_local);
		theta_ft = fftw_alloc_complex(alloc_local);

		plan_grf_f = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               field, field_ft, 
                              cart_comm, FFTW_FORWARD, FFTW_ESTIMATE);

		plan_grf_b = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               field_ft, field,  cart_comm, FFTW_BACKWARD, FFTW_ESTIMATE); 

		plan_grf_theta_b = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               	   theta_ft, theta, cart_comm, FFTW_BACKWARD, FFTW_ESTIMATE);



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


	void gen(double k_grid[][3],double *ini_dc,double *ini_theta, ini_power_generator p_k,double a_t,double a,double a0,double f_ini)
	{	int i,j,k,ci;
		double ksqr,pk_val,dtN;

		//printf("ini_dc %lf %d %d %d\n",ini_dc[ci],n_loc[0],n_loc[1],n_loc[2]);
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

			 theta_ft[ci][0] = 0.0;
			 theta_ft[ci][1] = 0.0;
			}
			else
			{
				pk_val = p_k.test_spec(sqrt(ksqr));				
					
				field_ft[ci][0] = sqrt(pk_val)*field_ft[ci][0]/dtN;
				field_ft[ci][1] = sqrt(pk_val)*field_ft[ci][1]/dtN;

				theta[ci][0] =  (a_t*Hi/a)*f_ini*(a/a0)*(a/a0)*field_ft[ci][0]/(ksqr*hbar_by_m);
				theta_ft[ci][1] =  (a_t*Hi/a)*f_ini*(a/a0)*(a/a0)*field_ft[ci][1]/(ksqr*hbar_by_m);
				//printf("Hi %lf\n",Hi);


			}

			 

			
			

		    }

		  }

		}
		
		fftw_execute(plan_grf_b);
		fftw_execute(plan_grf_theta_b);
		

		for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;			
			ini_dc[ci] = field[ci][0];
			ini_theta[ci] = theta[ci][0];
			
			

		    }

		  }

		}
		
	

	fftw_free(field);
	fftw_free(field_ft);
	fftw_free(theta);
	fftw_free(theta_ft);
	
	fftw_destroy_plan(plan_grf_f);
	fftw_destroy_plan(plan_grf_b);
	fftw_destroy_plan(plan_grf_theta_b);	
	}


///////////////////##################### Codes just for testing...############################//////////////////////////////////////////////////////////////////////////////
	void gen_check(double k_grid[][3],double f[][2],double ft[][2], ini_power_generator p_k)
	{	int i,j,k,ci;
		double ksqr,pk_val,dtN;

		dtN = (double)(n[0]*n[1]*n[2]);

		for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;			
						
			
			field[ci][0] = sqrt(twopie)*unit_normal();
			field[ci][1] = 0.0;
			

		    }

		  }

		}

		fftw_execute(plan_grf_f);


		for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;			
			ksqr = (k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);

			pk_val = p_k.get_ini_spectrum(sqrt(ksqr));			
			if(ksqr==0.0)
			{field_ft[ci][0] = 0.0;
			 field_ft[ci][1] = 0.0;
			}
			else
			{

				field_ft[ci][0] = sqrt(pk_val)*field_ft[ci][0]/dtN;
				field_ft[ci][1] = sqrt(pk_val)*field_ft[ci][1]/dtN;


			}

			ft[ci][0] = field_ft[ci][0];
			ft[ci][1] = field_ft[ci][1];
			

		    }

		  }

		}
		
		fftw_execute(plan_grf_b);


		for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;			
			f[ci][0] = field[ci][0];
			f[ci][1] = field[ci][1];
			

		    }

		  }

		}
		
	

	//fftw_free(field);
	//fftw_free(field_ft);
	
	//fftw_destroy_plan(plan_grf_f);
	//fftw_destroy_plan(plan_grf_b);	
	}
		
	void check(double k_grid[][3],ini_power_generator p_k,int *kbin_grid,int kbins,double dk,double abyai)
	{
		
		FILE *fp = fopen("grf2.txt","w");
		FILE *fp3 = fopen("3.txt","w");
		FILE *fppwr = fopen("pwr_grf.txt","w");

		double f[ntN][2],ft[ntN][2],pwr_spec[ntN],fld[ntN],pk_val,ksqr,kv,delta_pw;
		double kbincnt[kbins+1];
		double pwspctrm[kbins+1];
		int i,j,k,ci,maxi=0;

		gen_check(k_grid,f,ft,p_k);

		for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;

			fld[ci] = f[ci][0];
			fprintf(fp3,"%lf\t%lf\t%lf\n",fld[ci],f[ci][0],f[ci][1]);
				
			pwr_spec[ci] = 0.0;
			kbincnt[kbin_grid[ci]] = 0.0;
			if(kbin_grid[ci]>maxi)
				maxi = kbin_grid[ci];

		    }
		  }
		}


		for(i=0;i<ntN;++i)
		{ ++kbincnt[kbin_grid[i]];
			//fprintf(fp,"%d\t%lf\n",i,fld[i]);
		
		}

		cal_spectrum(fld,kbin_grid, kbins,n,pwr_spec, dk,abyai,fppwr);
		printf("hsjk %d\t%lf\t%d\n",kbins,kbincnt[0],maxi);
		for(i=0;i<=kbins;++i)
		{

		   //if(kbincnt[i]!=0)
	           {  delta_pw = sqrt(pwspctrm[i]*i*i*i*dk*dk*dk/(2.0*M_PI*M_PI*kbincnt[i]));  
			kv = (((double)(i+1))*dk*0.5);
			
			pk_val = p_k.get_ini_spectrum(kv);
			
		      fprintf(fp,"%.16lf\t%.16lf\t%.13lf\t%.13lf\t%.13lf\n",
							kv,pwr_spec[i],pk_val,f[i][0],f[i][1]);

		   }

	

		}


	}


	void stats_check(double k_grid[][3],ini_power_generator p_k,int *kbin_grid,int kbins,double dk,double abyai)
	{

		FILE *fp;
		

		double f[ntN][2],ft[ntN][2],pwr_spec[ntN],fld[ntN],pk_val,ksqr,kv,delta_pw;
		double kbincnt[kbins+1],famp2;
		double pwspctrm[kbins+1];
		int i,j,k,maxi=0;

	
		char fp_num[10];
		

		for(i=0;i<100;++i)
		{
			
			char fp_name[20]("stat_");
		
			sprintf(fp_num,"%d",i);
			strcat(fp_name,fp_num);
			fp = fopen(fp_name,"w");
			

			gen_check(k_grid,f,ft,p_k);
			printf("check i %d  %s\n",i,fp_name);
			for(j=0;j<ntN;++j)
			{
				ksqr = (k_grid[j][0]*k_grid[j][0]+k_grid[j][1]*k_grid[j][1]+k_grid[j][2]*k_grid[j][2]);
				pk_val = p_k.get_ini_spectrum(sqrt(ksqr));
				famp2 = ft[j][0]*ft[j][0]+ft[j][1]*ft[j][1];
				fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",ksqr,ft[j][0],ft[j][1],famp2,pk_val);

			}

			fclose(fp);

		}
		
		

	}




};
