

int evolve_kdk(int *n_glbl,int *n,fdm_psi_mpi &psi,metric_potential_mpi &phi,double k_grid[][3],int kbin_grid[],
					double a_final,double a_ini,double a0,double omega_dm_ini,double *dx,double dk,int kbins,double dt,bool use_hdf_format)
{	printf("Yo\n");
	double a,a_t,t,ak,a3a03omega,dti=dt;
	double a_print;
	double psi_vel[2],psi_val[n[0]*n[1]*n[2]][2],psi_upd[2],psi_k[2],potn,poisson_rhs,psi_amp,dc[n[0]*n[1]*n[2]],pwr_spec[n[0]*n[1]*n[2]];
	double psi_retrive[2];
	int i,j,k,ci,ind[3];
	int c1,c2,fail=0;
	int step_cnt;
	int num_prints = 14;
	double z_print_list[num_prints]{99.0,40.0,10.0,5.0,4.0,5.0,3.0,2.0,1.0,0.5,0.25,0.1,0.000001,0.0};
	double z_cur;
	
	int prn=0;	
	
	int mpi_check;
	
	
	FILE *fp_psi;
	FILE *fp_phi;

	FILE *fpwr_spec ;
	hid_t filename;
	herr_t status;
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id;

	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	
	

	
	
	


	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<=a0)&&(!fail);t+=dt,++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);
	 if(a>=a_print)
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0)/omega_dm_ini;
		z_cur = ((a0/a) -1.0);

		char fp_psi_name[20]("psi_z_");
		char fp_phi_name[20]("phi_z_");
		char fp_pwr_spec_name[20]("pwr_z_");
		char fp_hdf5_name[20]("data_z_");
		char fp_z_num[10];

		
		
		printf("a %lf %lf\n",a,a3a03omega);
		
	     

		if((z_cur<=z_print_list[prn])||(a==a_ini))
		{ 
			if(use_hdf_format)
			{
				sprintf(fp_z_num,"%.2lf",z_cur);
				strcat(fp_hdf5_name,fp_z_num); 
				strcat(fp_pwr_spec_name,fp_z_num); 
			  	strcat(fp_hdf5_name,".hdf5"); 
			  	strcat(fp_pwr_spec_name,".txt"); 
				

			  
			  	fpwr_spec = fopen(fp_pwr_spec_name,"w");

		 		plist_id = H5Pcreate(H5P_FILE_ACCESS);
        			H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 		filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
				H5Pclose(plist_id);
						
				
				evolve_hdf5_write(n_glbl,psi, phi,filename,dc,a3a03omega,a,true);
				status=H5Fclose(filename);

				printf("hdf5 %s\n",fp_hdf5_name);
				printf("pwr spec name %s\n",fp_pwr_spec_name);		
				
				
				//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
				fclose(fpwr_spec);

			}

			else
			{ sprintf(fp_z_num,"%.2lf",z_cur);
			  strcat(fp_psi_name,fp_z_num);
			  strcat(fp_phi_name,fp_z_num); 
			  strcat(fp_pwr_spec_name,fp_z_num); 
			  strcat(fp_psi_name,".txt"); 
			  strcat(fp_phi_name,".txt"); 
			  strcat(fp_pwr_spec_name,".txt"); 

			  fp_psi = fopen(fp_psi_name,"w");
			  fp_phi = fopen(fp_phi_name,"w");
			  fpwr_spec = fopen(fp_pwr_spec_name,"w");

			  printf("psi name %s\n",fp_psi_name);
			  printf("phi name %s\n",fp_phi_name);
			  printf("pwr spec name %s\n",fp_pwr_spec_name);

			  psi.write_psi(fp_psi,dc,dx,a3a03omega,a,true, true);
			  phi.write_potential(fp_phi,dx,a3a03omega,a);
			  cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			  printf("\nWriting at z = %lf\n",z_cur);
			  ++prn;
			  fclose(fp_psi);
	 		  fclose(fp_phi);
			  fclose(fpwr_spec);
			 }

		}
		
		

	  }
	a_t = a*sqrt(omega_dm_ini*pow(a0/a,3.0)+ (1.0-omega_dm_ini));
	  
	ak = a+a_t*dt;

	 for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;
			potn = phi.get_potential(ci);
			c1  = psi.get_psi(ind,psi_val[ci]);///Check if this as intented
			
			c1 = psi.calc_vel(ind,psi_vel, potn, a, a_t,dx);
			//printf("i %d j %d k %d psi r %10lf i  %10lf\n",i,j,k,psi_val[ci][0],psi_val[ci][1]);
			psi_k[0] = psi_val[ci][0]+ psi_vel[0]*dt;
			psi_k[1] = psi_val[ci][1]+ psi_vel[1]*dt;
			psi_amp = sqrt(psi_k[0]*psi_k[0] + psi_k[1]*psi_k[1]);
			
			psi.update(ind,psi_k[0],psi_k[1]);

			poisson_rhs = 1.5*H0*H0*ak*ak*(psi_amp*psi_amp - omega_dm_ini*pow(a0/ak,3.0));
			phi.update_4pieGpsi(ci,poisson_rhs);
			
			//fail =1;

		    }

		   }
	}
	phi.solve_poisson(psi,k_grid);
	mpi_check=psi.mpi_send_recv();
	mpi_check = MPI_Barrier(cart_comm);
	
	a_t = ak*sqrt(omega_dm_ini*pow(a0/ak,3.0)+ (1.0-omega_dm_ini));
	a = 0.5*(ak+a+a_t*dt);
	for(i=0;(i<n[0])&&(!fail);++i)
	 {
		  for(j=0;(j<n[1])&&(!fail);++j)
		  {
		    for(k=0;(k<n[2])&&(!fail);++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;
			
			c1  = psi.get_psi(ind,psi_retrive);			
			c1 = psi.calc_vel(ind,psi_vel, potn, ak, a_t,dx);
			psi_k[0] = 0.5*(psi_retrive[0]+psi_val[ci][0]+ psi_vel[0]*dt);
			psi_k[1] = 0.5*(psi_retrive[1]+psi_val[ci][1]+ psi_vel[1]*dt);
			//printf("psi_vel %lf %lf\n",psi_vel[0],psi_vel[1]);
			psi.update(ind,psi_k[0],psi_k[1]);

			psi_amp = sqrt(psi_k[0]*psi_k[0] + psi_k[1]*psi_k[1]);
			if(isnan(psi_amp))
			{	fail=1;
				//printf("FAIL %d %d %d\n",i,j,k);break;
			}

			poisson_rhs = 1.5*H0*H0*a*a*(psi_amp*psi_amp - omega_dm_ini*pow(a0/a,3.0));
			phi.update_4pieGpsi(ci,poisson_rhs);
			


		    }

		   }
	}

	phi.solve_poisson(psi,k_grid);
	mpi_check=psi.mpi_send_recv();

	
	 if(isnan(a))
	  {printf("FAILED  \n"); fail = 1;	break;
	  }

	}

	a3a03omega = pow(a/a0,3.0)/omega_dm_ini;	
	z_cur = ((a0/a) -1.0);

	char fp_psi_name[20]("psi_z_");
	char fp_phi_name[20]("phi_z_");
	char fp_pwr_spec_name[20]("pwr_z_");
	char fp_hdf5_name[20]("data_z_");
	char fp_z_num[10];

		if(use_hdf_format)
		{
			sprintf(fp_z_num,"%.2lf",z_cur);
			strcat(fp_hdf5_name,fp_z_num); 
			strcat(fp_pwr_spec_name,fp_z_num); 
		  	strcat(fp_hdf5_name,".hdf5"); 
		  	strcat(fp_pwr_spec_name,".txt"); 
			

			  
		  	fpwr_spec = fopen(fp_pwr_spec_name,"w");

		 	plist_id = H5Pcreate(H5P_FILE_ACCESS);
        		H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 	filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
			H5Pclose(plist_id);
						
				
			evolve_hdf5_write(n_glbl,psi, phi,filename,dc,a3a03omega,a,true);
			status=H5Fclose(filename);

			printf("hdf5 %s\n",fp_hdf5_name);
			printf("pwr spec name %s\n",fp_pwr_spec_name);		
			
		
			//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			fclose(fpwr_spec);

		}

		else
		{ sprintf(fp_z_num,"%.2lf",z_cur);
		  strcat(fp_psi_name,fp_z_num);
		  strcat(fp_phi_name,fp_z_num); 
		  strcat(fp_pwr_spec_name,fp_z_num); 
		  strcat(fp_psi_name,".txt"); 
		  strcat(fp_phi_name,".txt"); 
		  strcat(fp_pwr_spec_name,".txt"); 

		  fp_psi = fopen(fp_psi_name,"w");
		  fp_phi = fopen(fp_phi_name,"w");
		  fpwr_spec = fopen(fp_pwr_spec_name,"w");

		  printf("psi name %s\n",fp_psi_name);
		  printf("phi name %s\n",fp_phi_name);
		  printf("pwr spec name %s\n",fp_pwr_spec_name);

		  psi.write_psi(fp_psi,dc,dx,a3a03omega,a,true, true);
		  phi.write_potential(fp_phi,dx,a3a03omega,a);
		  cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
		  printf("\nWriting at z = %lf\n",z_cur);
		  
		  fclose(fp_psi);
	 	  fclose(fp_phi);
		  fclose(fpwr_spec);
		 }

	return(fail);

}

int evolve_kdk_openmp(int *n_glbl,int *n,fdm_psi_mpi &psi,metric_potential_mpi &phi,double k_grid[][3],int kbin_grid[],
					double a_final,double a_ini,double a0,double omega_dm_ini,double *dx,double dk,int kbins,double dt,bool use_hdf_format)
{	printf("Yo\n");
	double a,a_t,t,ak,a3a03omega,dti=dt;
	double a_print;
	double psi_vel[2],psi_val[n[0]*n[1]*n[2]][2],psi_upd[2],psi_k[2],potn,poisson_rhs,psi_amp,dc[n[0]*n[1]*n[2]],pwr_spec[n[0]*n[1]*n[2]];
	double psi_retrive[2];
	int i,j,k,ci,ind[3];
	int c1,c2,fail=0;
	int step_cnt;
	int num_prints = 14;
	double z_print_list[num_prints]{99.0,40.0,10.0,5.0,4.0,5.0,3.0,2.0,1.0,0.5,0.25,0.1,0.000001,0.0};
	double z_cur;
	
	int prn=0;	
	
	int mpi_check;
	
	FILE *fp_psi;
	FILE *fp_phi;

	FILE *fpwr_spec ;
	hid_t filename;
	herr_t status;
	
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id;

	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	
	

	
	
	


	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<=a0)&&(!fail);t+=dt,++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);
	 if(a>=a_print)
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0)/omega_dm_ini;
		z_cur = ((a0/a) -1.0);

		char fp_psi_name[20]("psi_z_");
		char fp_phi_name[20]("phi_z_");
		char fp_pwr_spec_name[20]("pwr_z_");
		char fp_hdf5_name[20]("data_z_");
		char fp_z_num[10];

		
		
		printf("a %lf %lf\n",a,a3a03omega);
		
	     

		if((z_cur<=z_print_list[prn])||(a==a_ini))
		{ 
			if(use_hdf_format)
			{
				sprintf(fp_z_num,"%.2lf",z_cur);
				strcat(fp_hdf5_name,fp_z_num); 
				strcat(fp_pwr_spec_name,fp_z_num); 
			  	strcat(fp_hdf5_name,".hdf5"); 
			  	strcat(fp_pwr_spec_name,".txt"); 
				

			  
			  	fpwr_spec = fopen(fp_pwr_spec_name,"w");

				plist_id = H5Pcreate(H5P_FILE_ACCESS);
        			H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 		filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
				H5Pclose(plist_id);
						
				
				evolve_hdf5_write(n_glbl,psi, phi,filename,dc,a3a03omega,a,true);
				status=H5Fclose(filename);
			
				printf("hdf5 %s\n",fp_hdf5_name);
				printf("pwr spec name %s\n",fp_pwr_spec_name);
	
				//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
				fclose(fpwr_spec);	
				 
				++prn;

			}

			else
			{ sprintf(fp_z_num,"%.2lf",z_cur);
			  strcat(fp_psi_name,fp_z_num);
			  strcat(fp_phi_name,fp_z_num); 
			  strcat(fp_pwr_spec_name,fp_z_num); 
			  strcat(fp_psi_name,".txt"); 
			  strcat(fp_phi_name,".txt"); 
			  strcat(fp_pwr_spec_name,".txt"); 

			  fp_psi = fopen(fp_psi_name,"w");
			  fp_phi = fopen(fp_phi_name,"w");
			  fpwr_spec = fopen(fp_pwr_spec_name,"w");

			  printf("psi name %s\n",fp_psi_name);
			  printf("phi name %s\n",fp_phi_name);
			  printf("pwr spec name %s\n",fp_pwr_spec_name);

			  psi.write_psi(fp_psi,dc,dx,a3a03omega,a,true, true);
			  phi.write_potential(fp_phi,dx,a3a03omega,a);
			  cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			  printf("\nWriting at z = %lf\n",z_cur);
			  ++prn;
			  fclose(fp_psi);
	 		  fclose(fp_phi);
			  fclose(fpwr_spec);
			 }

		}

		

		
		

	  }
	a_t = a*sqrt(omega_dm_ini*pow(a0/a,3.0)+ (1.0-omega_dm_ini));
	 
	ak = a+a_t*dt;
 #pragma omp parallel for private(j,k,ci,ind,c1,potn,psi_vel,psi_k,psi_amp,poisson_rhs)
	 for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;

			potn = phi.get_potential(ci);
			c1  = psi.get_psi(ind,psi_val[ci]);///Check if this as intented
			
			c1 = psi.calc_vel(ind,psi_vel, potn, a, a_t,dx);
			//printf("i %d j %d k %d psi r %10lf i  %10lf\n",i,j,k,psi_val[ci][0],psi_val[ci][1]);
			psi_k[0] = psi_val[ci][0]+ psi_vel[0]*dt;
			psi_k[1] = psi_val[ci][1]+ psi_vel[1]*dt;
			psi_amp = sqrt(psi_k[0]*psi_k[0] + psi_k[1]*psi_k[1]);
			
			psi.update(ind,psi_k[0],psi_k[1]);
	
			poisson_rhs = 1.5*H0*H0*ak*ak*(psi_amp*psi_amp - omega_dm_ini*pow(a0/ak,3.0));
			phi.update_4pieGpsi(ci,poisson_rhs);
			
			//fail =1;

		    }

		   }
	}
	phi.solve_poisson(psi,k_grid);
	mpi_check=psi.mpi_send_recv();
	MPI_Barrier(cart_comm);

	
	a_t = ak*sqrt(omega_dm_ini*pow(a0/ak,3.0)+ (1.0-omega_dm_ini));
	a = 0.5*(ak+a+a_t*dt);

 #pragma omp parallel for private(j,k,ci,ind,c1,potn,psi_vel,psi_k,psi_amp,poisson_rhs)
	for(i=0;(i<n[0]);++i)
	 {
		  for(j=0;(j<n[1]);++j)
		  {
		    for(k=0;(k<n[2]);++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;
			
			potn = phi.get_potential(ci);
			c1  = psi.get_psi(ind,psi_retrive);	
		
			c1 = psi.calc_vel(ind,psi_vel, potn, ak, a_t,dx);
			psi_k[0] = 0.5*(psi_retrive[0]+psi_val[ci][0]+ psi_vel[0]*dt);
			psi_k[1] = 0.5*(psi_retrive[1]+psi_val[ci][1]+ psi_vel[1]*dt);
			//printf("psi_vel %lf %lf\n",psi_vel[0],psi_vel[1]);
			psi.update(ind,psi_k[0],psi_k[1]);

			psi_amp = sqrt(psi_k[0]*psi_k[0] + psi_k[1]*psi_k[1]);
			if(isnan(psi_amp))
			{	fail=1;
				//printf("FAIL %d %d %d\n",i,j,k);break;
			}
	
			poisson_rhs = 1.5*H0*H0*a*a*(psi_amp*psi_amp - omega_dm_ini*pow(a0/a,3.0));
			phi.update_4pieGpsi(ci,poisson_rhs);
			


		    }

		}
	}

	phi.solve_poisson(psi,k_grid);
	mpi_check=psi.mpi_send_recv();
	MPI_Barrier(cart_comm);
	
	 if(isnan(a))
	  {printf("FAILED  \n"); fail = 1;	break;
	  }

	}

	a3a03omega = pow(a/a0,3.0)/omega_dm_ini;	
	z_cur = ((a0/a) -1.0);

	char fp_psi_name[20]("psi_z_");
	char fp_phi_name[20]("phi_z_");
	char fp_pwr_spec_name[20]("pwr_z_");
	char fp_hdf5_name[20]("data_z_");
	char fp_z_num[10];

		if(use_hdf_format)
		{
			sprintf(fp_z_num,"%.2lf",z_cur);
			strcat(fp_hdf5_name,fp_z_num); 
			strcat(fp_pwr_spec_name,fp_z_num); 
		  	strcat(fp_hdf5_name,".hdf5"); 
		  	strcat(fp_pwr_spec_name,".txt"); 
			

			  
		  	fpwr_spec = fopen(fp_pwr_spec_name,"w");

		 	plist_id = H5Pcreate(H5P_FILE_ACCESS);
        		H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 	filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
			H5Pclose(plist_id);
						
				
			evolve_hdf5_write(n_glbl,psi, phi,filename,dc,a3a03omega,a,true);
			status=H5Fclose(filename);
			printf("hdf5 %s\n",fp_hdf5_name);
			printf("pwr spec name %s\n",fp_pwr_spec_name);		
			
			

	
			
			//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			fclose(fpwr_spec);

		}

		else
		{ sprintf(fp_z_num,"%.2lf",z_cur);
		  strcat(fp_psi_name,fp_z_num);
		  strcat(fp_phi_name,fp_z_num); 
		  strcat(fp_pwr_spec_name,fp_z_num); 
		  strcat(fp_psi_name,".txt"); 
		  strcat(fp_phi_name,".txt"); 
		  strcat(fp_pwr_spec_name,".txt"); 

		  fp_psi = fopen(fp_psi_name,"w");
		  fp_phi = fopen(fp_phi_name,"w");
		  fpwr_spec = fopen(fp_pwr_spec_name,"w");

		  printf("psi name %s\n",fp_psi_name);
		  printf("phi name %s\n",fp_phi_name);
		  printf("pwr spec name %s\n",fp_pwr_spec_name);

		  psi.write_psi(fp_psi,dc,dx,a3a03omega,a,true, true);
		  phi.write_potential(fp_phi,dx,a3a03omega,a);
		  cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
		  printf("\nWriting at z = %lf\n",z_cur);
		  
		  fclose(fp_psi);
	 	  fclose(fp_phi);
		  fclose(fpwr_spec);
		 }

	return(fail);

}


void evolve_hdf5_write(int *ind,fdm_psi_mpi psi,metric_potential_mpi phi,hid_t filename,double *dc,double a3a03omega,double a,bool get_dc=false)
{	
	herr_t status_psi,status_phi,status;	
	hid_t file,dtype,dspace_psi,dspace_potn;
	hsize_t dim[3],pdim[2];
	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];

	

	int tN = ind[0]*ind[1]*ind[0];
	pdim[0]= tN; pdim[1] = 2;


	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    	status = H5Tset_order(dtype, H5T_ORDER_LE);
	dspace_psi = H5Screate_simple(3, dim, NULL);
	dspace_potn = H5Screate_simple(2, pdim, NULL);

	status_psi=psi.write_hdf5_psi_mpi(filename, dtype, dspace_psi,dc,a3a03omega,a,get_dc);
	
	
	status_phi=phi.write_hdf5_potn_mpi(filename, dtype,dspace_potn);

	H5Sclose(dspace_psi);
	H5Sclose(dspace_potn);
	H5Tclose(dtype);



}
