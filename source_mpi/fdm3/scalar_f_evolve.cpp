
int evolve_kdk_openmp(int *n_glbl,int *n,fdm_poisson_mpi &psi,metric_potential_poisson_mpi &phi,double k_grid[][3],int kbin_grid[],
					double a_final,double a_ini,double a0,double omega_dm_0,double Xb_0,double *dx,double dk,int kbins,double da,
							int cum_lin_id,bool use_hdf_format)
{	printf("OMP Yo\n");
	double a,a_t,a_tt,t,ak,a3a03omega;
	
	double a_print;
	double potn,potn_k,potn_a,potn_a_part,potn_rhs;

	double *pwr_spec = new double[n[0]*n[1]*n[2]];
	double *dc = new double[n[0]*n[1]*n[2]];
	double *potn_val = new double[n[0]*n[1]*n[2]];

	double psi_b[2],psi_amp2;
	
	double potn_retrive;
	int i,j,k,ci,ind[3];
	int c1,c2,fail=0;
	int step_cnt;
	int num_prints = 14;
	double z_print_list[num_prints]{99.0,96.0,40.0,10.0,5.0,4.0,5.0,3.0,2.0,1.0,0.5,0.25,0.1,0.000001,0.0};
	double z_cur;
	
	int prn=0;	

	
	
	int mpi_check;
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id;

	
	
	FILE *fp_phi;

	FILE *fpwr_spec ;
	hid_t filename;
	herr_t status;
	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	
	
	a_t = a_ini*H0*sqrt(omega_dm_0*pow(a0/a_ini,3.0*(1.0+w))+ (1.0-omega_dm_0));
	
	

	double k2lin[3];
	int midk,mink,maxk;

  /*    if(cum_lin_id==0)
      {
	maxk = (int)(((double)n_glbl[0])/2.0)-1;
	mink = 1;
	//maxk = n_glbl[0]-1;
	k2lin[0] = 	sqrt(k_grid[mink][0]*k_grid[mink][0]+k_grid[mink][1]*k_grid[mink][1]+k_grid[mink][2]*k_grid[mink][2]);
	k2lin[2] = 	sqrt(k_grid[maxk][0]*k_grid[maxk][0]+k_grid[maxk][1]*k_grid[maxk][1]+k_grid[maxk][2]*k_grid[maxk][2]);
	k2lin[1] = 	0.5*k2lin[2];
	
	linear_poisson_field_mpi lin_calc(k2lin,0.001,1.0);
	printf("Running lin calc....%d\n",maxk);
	lin_calc.evolve(0.001,da,omega_dm_0,H0);
	printf("Lin_calc...Done...\n");

      }*/
	//MPI_Barrier(cart_comm);
	

	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<=a0)&&(!fail);++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);


	  a_t = a*H0*sqrt(omega_dm_0*pow(a0/a,3.0*(1.0+w))+ (1.0-omega_dm_0));
          a_tt = a*H0*H0*(-0.5*omega_dm_0*pow(a0/a,3.0*(1.0+w))*(1.0+3.0*w) + (1.0-omega_dm_0) );
	  psi_b[0] = sqrt(3.0*H0*H0*omega_dm_0*a0*a0*a0);
	  psi_b[1] = 0.0;
	 
	  
	  //printf("fb_a %lf\n",fb_a);
	
	  ak = a+da;
	  


	 if(a>=a_print)
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;
		z_cur = ((a0/a) -1.0);

		
		char fp_phi_name[20]("phi_z_");
		char fp_pwr_spec_name[20]("pwr_z_");
		char fp_hdf5_name[20]("data_z_");
		char fp_z_num[10];

		
		
		//printf("a %lf %lf  %lf\n",a,0.5*fb_a*fb_a*a_t*a_t,Xb_0*pow(a0/a,6.0/(2.0*alpha-1.0)));
		
	     

		if((z_cur<=z_print_list[prn])||(a==a_ini))
		{ 
			if(use_hdf_format)
			{
				sprintf(fp_z_num,"%.2lf",z_cur);
				strcat(fp_hdf5_name,fp_z_num); 
				//strcat(fp_pwr_spec_name,fp_z_num); 
			  	strcat(fp_hdf5_name,".hdf5"); 
			  	//strcat(fp_pwr_spec_name,".txt"); 
				

			  
			  	//fpwr_spec = fopen(fp_pwr_spec_name,"w");

		 		plist_id = H5Pcreate(H5P_FILE_ACCESS);
        			H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 		filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
				H5Pclose(plist_id);

				printf("hdf5 %s\n",fp_hdf5_name);
						
				
				evolve_hdf5_write(n_glbl,psi, phi,filename,dc,a0,a,a_t,omega_dm_0,cum_lin_id,true);
				
		
				status=H5Fclose(filename);
				//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
				//fclose(fpwr_spec);	
				//printf("pwr spec name %s\n",fp_pwr_spec_name);
				++prn;

			}

		
		}

		
	
		
		

	  }
	
	
	
	
	//printf("fb_t  %lf   fb_t_th  %lf\n",fb_t,fb_t_0*pow(a0/ak,3.0/(2.0*alpha-1.0)));
 #pragma omp parallel for private(j,k,ci,ind,c1,fa_k,fa_vel,potn,potn_k,potn_a,poisson_rhs,acc_fa,potn_der)

	for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {


			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;

			psi.update_amp2_value(ind);



		    }
		 }


	}
			


	 for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;


			
			if(step_cnt==0)
			potn_k = phi.get_value(ind);
			else
			potn_k = phi.get_potential(ci);

			if(step_cnt==0)
			potn_a = 0.0;
			else
			potn_a = (potn_k-potn_val[ci])/da;

			potn_val[ci] = potn_k;

	
			
			psi_amp2 = psi.get_value(ind);				


			phi.update_value(ind, potn_k);
		

			
			psi.update_A( ci, potn_k, a, da);

			c1 = phi.calc_vel(ind,potn_a_part,psi_amp2,potn_k,a,da,a_t,a_tt,dx,omega_dm_0);


				


			potn_rhs = potn_val[ci]+da*potn_a_part;
			


			if(isnan(potn_rhs))
			{
				//printf("Yfailed at step %d %lf %lf %lf\n",step_cnt,potn_k,potn_t,fa_k[1]);
				fail=1;
				break;

			}

			else
			{
				phi.update_4pieGpsi(ci,potn_rhs);
				

			}
			

		    }

		   }
	}

	
	

	
	
	

	//printf("Checking  %lf\n",(f_a_avg-fb_a)/fb_a);

	 a_t = ak*H0*sqrt(omega_dm_0*pow(a0/ak,3.0*(1.0+w))+ (1.0-omega_dm_0));
         a_tt = ak*H0*H0*(-0.5*omega_dm_0*pow(a0/ak,3.0*(1.0+w))*(1.0+3.0*w) + (1.0-omega_dm_0) );


	phi.solve_poisson(k_grid, a, a_t,da);
	psi.solve_poisson(k_grid,psi_b, a, a_t,da);



	//printf("Check f_a %.10lf %.10lf\n",fb_a,fcheck);
	

	a = a+da;


	
	 if(isnan(a))
	  {printf("FAILED  \n"); fail = 1;
	  }

	}

	a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;	
	z_cur = ((a0/a) -1.0);

	
	//char fp_phi_name[20]("phi_z_");
	//char fp_pwr_spec_name[20]("pwr_z_");
	char fp_hdf5_name[20]("data_z_");
	char fp_z_num[10];

		if(use_hdf_format)
		{
			sprintf(fp_z_num,"%.2lf",z_cur);
			strcat(fp_hdf5_name,fp_z_num); 
			//strcat(fp_pwr_spec_name,fp_z_num); 
		  	strcat(fp_hdf5_name,".hdf5"); 
		  	//strcat(fp_pwr_spec_name,".txt"); 
			

			  
		  	//fpwr_spec = fopen(fp_pwr_spec_name,"w");
		 	plist_id = H5Pcreate(H5P_FILE_ACCESS);
        		H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 	filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
			H5Pclose(plist_id);

			printf("hdf5 %s\n",fp_hdf5_name);
			//printf("pwr spec name %s\n",fp_pwr_spec_name);		
			
			evolve_hdf5_write(n_glbl,psi, phi,filename,dc,a0,a,a_t,omega_dm_0,cum_lin_id,true);

	
			status=H5Fclose(filename);
			//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			//fclose(fpwr_spec);

		}


	delete[] pwr_spec;
	delete[] dc;
	delete[] potn_val;

		

	return(fail);

}

	
void evolve_hdf5_write(int *ind,fdm_poisson_mpi psi,
				metric_potential_poisson_mpi phi,hid_t filename,double *dc,double a0,double a,double a_t,double omega_dm_0,
														int cum_lin_id,bool get_dc=false)
{	
	herr_t status_psi,status_phi,status;	
	hid_t file,dtype,dspace_psi,dspace_potn,dspace_dc;
	hsize_t dim[3],pdim[1];
	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];

	
	
	int tN = ind[0]*ind[1]*ind[0];
	pdim[0]= tN; 


	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);
	
	dspace_potn = H5Screate_simple(3, dim, NULL);
	dspace_dc = H5Screate_simple(1, pdim, NULL);
	dspace_psi = H5Screate_simple(3, dim, NULL);
	

	

	
	
	status_psi = psi.write_hdf5_values_mpi(filename, dtype, dspace_psi,dspace_dc,dc,a0,a,a_t,omega_dm_0,phi,cum_lin_id,get_dc);
	H5Sclose(dspace_psi);

	
	
	status_phi = phi.write_hdf5_values_mpi(filename, dtype, dspace_potn,a0,a,a_t,cum_lin_id,get_dc);
	H5Sclose(dspace_potn);	

		



	H5Sclose(dspace_dc);
	

	
	H5Tclose(dtype);


}







