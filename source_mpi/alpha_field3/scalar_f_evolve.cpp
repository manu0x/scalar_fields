
int evolve_kdk_openmp(int *n_glbl,int *n,metric_potential_poisson_mpi &f_alpha,metric_potential_poisson_mpi &phi,double k_grid[][3],int kbin_grid[],
					double a_final,double a_ini,double a0,double omega_dm_0,double Xb_0,double *dx,double dk,int kbins,double da,
							int cum_lin_id,bool use_hdf_format)
{	printf("OMP Yo\n");
	double a,a_t,a_tt,t,ak,a3a03omega;
	
	double a_print;
	double fa_vel[2],f_a_val,potn_val[n[0]*n[1]*n[2]],fa_k[2],potn,potn_k,
					potn_a,potn_a_part,f_a_a_part,dc[n[0]*n[1]*n[2]],pwr_spec[n[0]*n[1]*n[2]],acc_fa,potn_rhs,f_a_rhs;

	double fb_a,fb_a_k,fb_acc,fb_a_0,Xb,f_a_avg;
	double fa_retrive[2],potn_der[3],potn_retrive;
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

	
	FILE *fp_falpha;
	FILE *fp_phi;

	FILE *fpwr_spec ;
	hid_t filename;
	herr_t status;
	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	
	
	a_t = a_ini*H0*sqrt(omega_dm_0*pow(a0/a_ini,3.0*(1.0+w))+ (1.0-omega_dm_0));
	fb_a = sqrt(2.0*(Xb_0*pow(a0/a_ini,6.0/(2.0*alpha-1.0))))/a_t;
	fb_a_0 = sqrt(2.0*Xb_0)/a_t;
	

	double k2lin[3];
	int midk,mink,maxk;

      if(cum_lin_id==0)
      {
	midk = (int)(((double)n_glbl[0])/2.0);
	mink = 1;
	maxk = n_glbl[0]-1;
	k2lin[0] = 	k_grid[mink][0]*k_grid[mink][0]+k_grid[mink][1]*k_grid[mink][1]+k_grid[mink][2]*k_grid[mink][2];
	k2lin[2] = 	k_grid[maxk][0]*k_grid[maxk][0]+k_grid[maxk][1]*k_grid[maxk][1]+k_grid[maxk][2]*k_grid[maxk][2];
	k2lin[1] = 	k_grid[midk][0]*k_grid[midk][0]+k_grid[midk][1]*k_grid[midk][1]+k_grid[midk][2]*k_grid[midk][2];
	
	linear_poisson_field_mpi lin_calc(k2lin,0.001,1.0);
	lin_calc.evolve(0.001,da,omega_dm_0,H0);

      }
	MPI_Barrier(cart_comm);
	

	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<=a0)&&(!fail);++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);


	  a_t = a*H0*sqrt(omega_dm_0*pow(a0/a,3.0*(1.0+w))+ (1.0-omega_dm_0));
          a_tt = a*H0*H0*(-0.5*omega_dm_0*pow(a0/a,3.0*(1.0+w))*(1.0+3.0*w) + (1.0-omega_dm_0) );
	  fb_acc = -3.0*fb_a*a_t/(a*(2.0*alpha-1.0));
	  fb_a_k = fb_a + fb_acc*da;

	  Xb = 0.5*fb_a*fb_a*a_t*a_t;
	  
	  ak = a+da;


	 if(a>=a_print)
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;
		z_cur = ((a0/a) -1.0);

		char fp_falpha_name[20]("f_alpha_z_");
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
				strcat(fp_pwr_spec_name,fp_z_num); 
			  	strcat(fp_hdf5_name,".hdf5"); 
			  	strcat(fp_pwr_spec_name,".txt"); 
				

			  
			  	fpwr_spec = fopen(fp_pwr_spec_name,"w");

		 		plist_id = H5Pcreate(H5P_FILE_ACCESS);
        			H5Pset_fapl_mpio(plist_id, cart_comm, info);

		
		 		filename = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
				H5Pclose(plist_id);

				printf("hdf5 %s\n",fp_hdf5_name);
						
				
				evolve_hdf5_write(n_glbl,f_alpha, phi,filename,dc,a0,a,a_t,0.5*f_a_avg*f_a_avg*a_t*a_t,cum_lin_id,true);
				status=H5Fclose(filename);
				//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
				fclose(fpwr_spec);	
				printf("pwr spec name %s\n",fp_pwr_spec_name);
				++prn;

			}

		
		}

		
	
		
		

	  }
	
	
	
	f_a_avg = 0.0;

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


			
			if(step_cnt==0)
			potn_k = phi.get_value(ind);
			else
			potn_k = phi.get_potential(ci);

			if(step_cnt==0)
			potn_a = 0.0;
			else
			potn_a = (potn_k-potn_val[ci])/da;

			potn_val[ci] = potn_k;
			
				

			//c1 = phi.get_potn_spt_der(ind,potn_der);
			if(step_cnt==0)
			f_a_val  = f_alpha.get_value(ind);
			else
			f_a_val = f_alpha.get_potential(ci);

			f_a_avg+= f_a_val;

			phi.update_value(ind, potn_k);
			f_alpha.update_value(ind, f_a_val);
			
			c1 = phi.calc_vel(ind,potn_a_part,f_a_val,potn_k,potn_a,a,a_t,a_tt,dx,omega_dm_0,Xb);


			c1 = f_alpha.calc_vel(ind,f_a_a_part,f_a_val,potn_k,potn_a,a,a_t,a_tt,dx,omega_dm_0,Xb);
				
			
			
			potn_rhs = potn_val[ci]+da*potn_a_part;
			f_a_rhs = f_a_val + da*f_a_a_part;
			
			
			//if(ci==10)
			//printf("fb_a  %.10lf   fb_a10  %.15lf \n",(fb_acc-acc_fa)/fb_acc,(fa_val[ci][1]/fb_a)-1.0);
			
			//printf(" %d %lf %lf %lf\n",step_cnt,potn_t,acc_fa,fa_k[1]);
			//if(ci==100)
			//printf("%lf %lf\n",fb_a,f_a_val);

			if(isnan(potn_rhs+f_a_rhs))
			{
				//printf("Yfailed at step %d %lf %lf %lf\n",step_cnt,potn_k,potn_t,fa_k[1]);
				fail=1;
				break;

			}

			else
			{
				phi.update_4pieGpsi(ci,potn_rhs);
				f_alpha.update_4pieGpsi(ci,f_a_rhs);

			}
			

		    }

		   }
	}

	
	
	
	
	
	f_a_avg=f_a_avg/((double)(n[0]*n[1]*n[2]));

	//printf("Checking  %lf\n",(f_a_avg-fb_a)/fb_a);

	 a_t = ak*H0*sqrt(omega_dm_0*pow(a0/ak,3.0*(1.0+w))+ (1.0-omega_dm_0));
         a_tt = ak*H0*H0*(-0.5*omega_dm_0*pow(a0/ak,3.0*(1.0+w))*(1.0+3.0*w) + (1.0-omega_dm_0) );
	 fb_acc = -3.0*fb_a_k*a_t/(ak*(2.0*alpha-1.0)) ;
	 fb_a = fb_a + fb_acc*da;
	 


	phi.solve_poisson(k_grid,fb_a, a, a_t,da);
	f_alpha.solve_poisson(k_grid,fb_a, a, a_t,da);
	
	double fcheck = f_alpha.get_potential(0,1);

	//printf("Check f_a %.10lf %.10lf\n",fb_a,fcheck);
	

	a = a+da;


	
	 if(isnan(a))
	  {printf("FAILED  \n"); fail = 1;	break;
	  }

	}

	a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;	
	z_cur = ((a0/a) -1.0);

	char fp_falpha_name[20]("f_alpha_z_");
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

			printf("hdf5 %s\n",fp_hdf5_name);
			printf("pwr spec name %s\n",fp_pwr_spec_name);		
			
			evolve_hdf5_write(n_glbl,f_alpha, phi,filename,dc,a0,a,a_t,0.5*f_a_avg*f_a_avg*a_t*a_t,cum_lin_id,true);

	
			status=H5Fclose(filename);
			//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			fclose(fpwr_spec);

		}

		

	return(fail);

}

	
void evolve_hdf5_write(int *ind,metric_potential_poisson_mpi f_alpha,metric_potential_poisson_mpi phi,hid_t filename,double *dc,double a0,double a,double a_t,double Xb,
														int cum_lin_id,bool get_dc=false)
{	
	herr_t status_falpha,status_phi,status;	
	hid_t file,dtype,dspace_falpha,dspace_potn,dspace_dc,dspace_dc2;
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
	dspace_falpha = H5Screate_simple(3, dim, NULL);

	

	
	
	
	status_falpha = f_alpha.write_hdf5_values_mpi(filename, dtype, dspace_falpha,dspace_dc,dc,a0,a,a_t,Xb,phi,cum_lin_id,get_dc);
	H5Sclose(dspace_falpha);

	

	status_phi = phi.write_hdf5_values_mpi(filename, dtype, dspace_potn,dspace_dc,dc,a0,a,a_t,Xb,phi,cum_lin_id,get_dc);
	H5Sclose(dspace_potn);	



	H5Sclose(dspace_dc);
	
	

	
	
	
	H5Tclose(dtype);




}








/*
int evolve_kdk(int *n_glbl,int *n,field_alpha_mpi &f_alpha,metric_potential_approx_1_t_mpi &phi,double k_grid[][3],int kbin_grid[],
					double a_final,double a_ini,double a0,double omega_dm_0,double Xb_0,double *dx,double dk,int kbins,double dt,
														int cum_lin_id,bool use_hdf_format)
{

	return(0);
}


							
{	printf("Yo\n");
	double a,a_t,t,ak,a3a03omega,dti=dt;
	double a_print;
	double fa_vel[2],fa_val[n[0]*n[1]*n[2]][2],potn_val[n[0]*n[1]*n[2]],fa_k[2],potn,potn_k,
					potn_t,dc[n[0]*n[1]*n[2]],pwr_spec[n[0]*n[1]*n[2]],acc_fa;
	double fa_retrive[2],potn_der[3],potn_retrive;
	int i,j,k,ci,ind[3];
	int c1,c2,fail=0;
	int step_cnt;
	int num_prints = 14;
	double z_print_list[num_prints]{99.0,40.0,10.0,5.0,4.0,5.0,3.0,2.0,1.0,0.5,0.25,0.1,0.000001,0.0};
	double z_cur;
	
	int prn=0;	
	
	int mpi_check;
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id;

	
	
	FILE *fp_falpha;
	FILE *fp_phi;

	FILE *fpwr_spec ;
	hid_t filename;
	herr_t status;
	//cal_spectrum(double *f,int *kbingrid,int kbins,int *s,double *pwspctrm,double dk,double abyai, FILE *fspwrite)
	
	

	
	
	


	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<=a0)&&(!fail);t+=dt,++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);
	 if(a>=a_print)
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;
		z_cur = ((a0/a) -1.0);

		char fp_falpha_name[20]("f_alpha_z_");
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

				printf("hdf5 %s\n",fp_hdf5_name);
				printf("pwr spec name %s\n",fp_pwr_spec_name);		
				
				evolve_hdf5_write(n_glbl,f_alpha, phi,filename,dc,a0,a,Xb_0,cum_lin_id,true);

				
				status=H5Fclose (filename);
				cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
				fclose(fpwr_spec);

			}

			else
			{ sprintf(fp_z_num,"%.2lf",z_cur);
			  strcat(fp_falpha_name,fp_z_num);
			  strcat(fp_phi_name,fp_z_num); 
			  strcat(fp_pwr_spec_name,fp_z_num); 
			  strcat(fp_falpha_name,".txt"); 
			  strcat(fp_phi_name,".txt"); 
			  strcat(fp_pwr_spec_name,".txt"); 

			  fp_falpha = fopen(fp_falpha_name,"w");
			  fp_phi = fopen(fp_phi_name,"w");
			  fpwr_spec = fopen(fp_pwr_spec_name,"w");

			  printf("f_alpha name %s\n",fp_falpha_name);
			  printf("phi name %s\n",fp_phi_name);
			  printf("pwr spec name %s\n",fp_pwr_spec_name);

			  f_alpha.write_f_alpha(fp_falpha,dc,dx,a3a03omega,a,false, true);
				
			  phi.write_potential(fp_phi,dx,a3a03omega,a);
			
			  cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			  printf("\nWriting at z = %lf\n",z_cur);
			  ++prn;
			  fclose(fp_falpha);
	 		  fclose(fp_phi);
			  fclose(fpwr_spec);
			 }

		}
		
		

	  }
	a_t = a*sqrt(omega_dm_0*pow(a0/a,3.0*(1.0+w))+ (1.0-omega_dm_0));
	  
	ak = a+a_t*dt;

	 for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;

			potn_val[ci] = phi.get_potential(ind);

			c1 = phi.get_potn_spt_der(ind,potn_der);
			c1  = f_alpha.get_field_alpha(ind,fa_val[ci]);///Check if this as intented

			
			//int calc_vel(int * ind,double &potn_vel,double f_t,double a,double a_t,double *dx,double omega_dm_0,double Xb_0)
			
			c1 = phi.calc_vel(ind,potn_t,fa_val[ci][1],a,a_t,dx,omega_dm_0,Xb_0);

			c1 = f_alpha.calc_acc(ind,acc_fa, potn_val[ci], potn_t,potn_der, a, a_t,dx);
				
			//printf("i %d j %d k %d psi r %10lf i  %10lf\n",i,j,k,potn_val[ci],potn_der[0]+potn_der[1]+potn_der[2]);
			fa_k[0] = fa_val[ci][0]+ fa_val[ci][1]*dt;
			fa_k[1] = fa_val[ci][1]+ acc_fa*dt;
			
			potn_k = potn_val[ci]+ potn_t*dt;
			
			
			//printf(" %d %lf %lf %lf\n",step_cnt,potn_t,acc_fa,fa_k[1]);

			if(isnan(potn_k+fa_k[0]+fa_k[1]))
			{
				printf("Yfailed at step %d %lf %lf %lf\n",step_cnt,potn_k,potn_t,fa_k[1]);
				fail=1;
				break;

			}

			else
			{

				phi.update(ind,potn_k);
				f_alpha.update(ind,fa_k[0],fa_k[1]);

			}
			

			
			
			//fail =1;

		    }

		   }
	}
	
	mpi_check = f_alpha.mpi_send_recv();
	mpi_check = phi.mpi_send_recv();
	mpi_check = MPI_Barrier(cart_comm);
	
	a_t = ak*sqrt(omega_dm_0*pow(a0/ak,3.0*(1.0+w))+ (1.0-omega_dm_0));
	a = 0.5*(ak+a+a_t*dt);
	for(i=0;(i<n[0])&&(!fail);++i)
	 {
		  for(j=0;(j<n[1])&&(!fail);++j)
		  {
		    for(k=0;(k<n[2])&&(!fail);++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;



			
			c1  = f_alpha.get_field_alpha(ind,fa_retrive);
			potn_k = phi.get_potential(ind);

				
			c1 = phi.calc_vel(ind,potn_t,fa_retrive[1],a,a_t,dx,omega_dm_0,Xb_0);	
			c1 = phi.get_potn_spt_der(ind,potn_der);	
			c1 = f_alpha.calc_acc(ind,acc_fa, potn_k, potn_t,potn_der, a, a_t,dx);
			//printf("failed at step %d %.15lf %lf %lf\n",step_cnt,potn_k,fa_retrive[1],acc_fa);

			fa_k[0] = 0.5*(fa_retrive[0]+fa_val[ci][0]+ fa_retrive[1]*dt);
			fa_k[1] = 0.5*(fa_retrive[1]+fa_val[ci][1]+ acc_fa*dt);

			potn_k = 0.5*(potn_k+potn_val[ci]+potn_t*dt);
			
			phi.update(ind,potn_k);
			f_alpha.update(ind,fa_k[0],fa_k[1]);

			if(isnan(potn_k+fa_k[0]+fa_k[1]))
			{
				//printf("failed at step %d\n",step_cnt);
				fail=1;
				break;

			}
			
			


		    }

		   }
	}

	
	mpi_check=phi.mpi_send_recv();
	mpi_check=f_alpha.mpi_send_recv();

	
	 if(isnan(a))
	  {printf("FAILED  \n"); fail = 1;	break;
	  }

	}

	a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;	
	z_cur = ((a0/a) -1.0);

	char fp_falpha_name[20]("f_alpha_z_");
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
			printf("hdf5 %s\n",fp_hdf5_name);
			printf("pwr spec name %s\n",fp_pwr_spec_name);		
			
			evolve_hdf5_write(n_glbl,f_alpha, phi,filename,dc,a0,a,Xb_0,cum_lin_id,true);
			status=H5Fclose (filename);
			cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			fclose(fpwr_spec);

		}

		else
		{ 	  sprintf(fp_z_num,"%.2lf",z_cur);
			  strcat(fp_falpha_name,fp_z_num);
			  strcat(fp_phi_name,fp_z_num); 
			  strcat(fp_pwr_spec_name,fp_z_num); 
			  strcat(fp_falpha_name,".txt"); 
			  strcat(fp_phi_name,".txt"); 
			  strcat(fp_pwr_spec_name,".txt"); 

			  fp_falpha = fopen(fp_falpha_name,"w");
			  fp_phi = fopen(fp_phi_name,"w");
			  fpwr_spec = fopen(fp_pwr_spec_name,"w");

			  printf("f_alpha name %s\n",fp_falpha_name);
			  printf("phi name %s\n",fp_phi_name);
			  printf("pwr spec name %s\n",fp_pwr_spec_name);

			  f_alpha.write_f_alpha(fp_falpha,dc,dx,a3a03omega,a,false, true);
			  phi.write_potential(fp_phi,dx,a3a03omega,a);
			  cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			  printf("\nWriting at z = %lf\n",z_cur);
			  ++prn;
			  fclose(fp_falpha);
	 		  fclose(fp_phi);
			  fclose(fpwr_spec);
		 }

	return(fail);

}



*/




