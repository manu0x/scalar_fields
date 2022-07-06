
int evolve_kdk_openmp(int *n_glbl,int *n,fdm_poisson_mpi &psi,metric_potential_poisson_mpi &phi,double **k_grid,int* kbin_grid,
					double a_final,double a_ini,double a0,double omega_dm_0,double Xb_0,double *dx,double dk,int kbins,double da,
							int cum_lin_id,bool use_hdf_format)
{	printf("OMP Yo %lf\n",3.0*H0*H0*omega_dm_0*a0*a0*a0);
	double a,a_t,a_tt,t,ak,a3a03omega;

	int krk;
	
	double a_print;
	double potn,potn_k,potn_a,potn_a_part,potn_rhs;

	//double *pwr_spec = new double[n[0]*n[1]*n[2]];
	double *dc = new double[n[0]*n[1]*n[2]];
	double *fdc = new double[n[0]*n[1]*n[2]];
	
	
	double *psi_hold=new double[n[0]*n[1]*n[2]*2];
	double *psi_sum=new double[n[0]*n[1]*n[2]*2];
	
	double *phi_hold=new double[n[0]*n[1]*n[2]];
	double *phi_sum=new double[n[0]*n[1]*n[2]];

	double psi_obt[2],phi_obt;

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

      if(cum_lin_id==0)
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

      }
	MPI_Barrier(cart_comm);
	

	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<=a0)&&(!fail);++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);


	  a_t = a*H0*sqrt(omega_dm_0*pow(a0/a,3.0*(1.0+w))+ (1.0-omega_dm_0));
          a_tt = a*H0*H0*(-0.5*omega_dm_0*pow(a0/a,3.0*(1.0+w))*(1.0+3.0*w) + (1.0-omega_dm_0) );
	  psi_b[0] = sqrt(3.0*H0*H0*omega_dm_0*a0*a0*a0);
	  psi_b[1] = 0.0;
	 
	  
	  //printf("fb_a %lf\n",fb_a);
	



	   for(i=0;i<n[0];++i)
	   {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {


			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;

			psi.update_amp2_value(ind);
			psi.get_fdm(ci,psi_obt);

			if(a==a_ini)
			 {
				phi_sum[ci] = phi.get_value(ind);	
				psi_sum[2*ci] = psi_obt[0];
				psi_sum[2*ci+1] = psi_obt[1];

			 }
			else
			 {
					

				psi_sum[2*ci] =  psi_sum[2*ci]+fsrk[3]*(psi_obt[0]-psi_hold[2*ci]);
				psi_sum[2*ci+1] =  psi_sum[2*ci+1]+fsrk[3]*(psi_obt[1]-psi_hold[2*ci+1]);
				
				phi_obt = phi.get_potential(ci);							
				phi_sum[ci] =  phi_sum[ci]+fsrk[3]*(phi_obt-phi_hold[ci]);	


				psi_b[0] = psi_sum[2*ci];
				psi_b[1] = psi_sum[2*ci+1];  

				psi.update_fdm(ci,psi_b);

				
				phi.update_4pieGpsi(ci,phi_sum[ci]);


			 }


			psi_hold[2*ci] = psi_sum[2*ci];
			psi_hold[2*ci+1] = psi_sum[2*ci+1];

			phi_hold[ci] = phi_sum[ci];



		    }
		 }


	    }
	  

	z_cur = ((a0/a) -1.0);
	 if(z_cur<=z_print_list[prn])
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0*(1.0+w))/omega_dm_0;
		

		
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
						
				
				evolve_hdf5_write(n_glbl,psi, phi,filename,dc,fdc,a0,a,a_t,omega_dm_0,cum_lin_id,true);
				
		
				status=H5Fclose(filename);
				//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
				//fclose(fpwr_spec);	
				//printf("pwr spec name %s\n",fp_pwr_spec_name);
				++prn;

			}

			psi.cal_mass(a);

		
		}

		
	
		
		

	  }



			
     ak = a;	
     for(krk = 0;krk<4;++krk)
     { 

	

	 for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;


			
			if(step_cnt==0)
			phi_obt = phi.get_value(ind);
			else
			phi_obt = phi.get_potential(ci);

			psi.get_fdm(ci,psi_obt);
			
			psi_amp2 = psi.get_amp2_value(ci);				

			//psiamp2_val[ci] = psi_amp2;
			phi.update_value(ind, phi_obt);


			if(krk>0)
			 { psi_sum[2*ci] =  psi_sum[2*ci]+fsrk[krk-1]*(psi_obt[0]-psi_hold[2*ci]);
			   psi_sum[2*ci+1] =  psi_sum[2*ci+1]+fsrk[krk-1]*(psi_obt[1]-psi_hold[2*ci+1]);
				
								
			   phi_sum[ci] =  phi_sum[ci]+fsrk[krk-1]*(phi_obt-phi_hold[ci]);
			 }
				

			
			psi.update_A( ci, phi_obt, psi_hold[2*ci], psi_hold[2*ci+1],krk, a,a_t, da);

			if(method==0)
			c1 = phi.calc_vel(ind,potn_a_part,psi_amp2,phi_obt,ak,da,a_t,a_tt,dx,omega_dm_0);

				

			  if(method==0)
			  potn_rhs = phi_hold[ci]+ rk[krk]*da*potn_a_part;
			//  if(method==1)
			//  potn_rhs = 0.5*psi_amp2-1.5*H0*H0*omega_dm_0*a0*a0*a0;


			if(isnan(potn_rhs)||isnan(psi_amp2))
			{
				//printf("Yfailed at step %d %lf %lf %lf\n",step_cnt,potn_k,potn_a,a_t);
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


	//phi.solve_poisson(k_grid, ak, a_t,da);

	psi.solve_poisson(k_grid,psi_b,krk, ak, a_t,da);
	phi.solve_poisson(k_grid,krk, ak, a_t,da);

	ak = a + rk[krk]*da;
	a_t = ak*H0*sqrt(omega_dm_0*pow(a0/ak,3.0*(1.0+w))+ (1.0-omega_dm_0));
        a_tt = ak*H0*H0*(-0.5*omega_dm_0*pow(a0/ak,3.0*(1.0+w))*(1.0+3.0*w) + (1.0-omega_dm_0) );

      }

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

			printf("final_hdf5 %s\n",fp_hdf5_name);
			//printf("pwr spec name %s\n",fp_pwr_spec_name);		
			
			evolve_hdf5_write(n_glbl,psi, phi,filename,dc,fdc,a0,a,a_t,omega_dm_0,cum_lin_id,true);

	
			status=H5Fclose(filename);
			//cal_spectrum(dc,kbin_grid, kbins,n,pwr_spec, dk,a/a_ini,fpwr_spec);
			//fclose(fpwr_spec);

		}


	//delete[] pwr_spec;
	delete[] dc;
	delete[] fdc;

	delete[] psi_hold;
	delete[] psi_sum;

	delete[] phi_hold;
	delete[] phi_sum;


		

	return(fail);

}

	
void evolve_hdf5_write(int *ind,fdm_poisson_mpi psi,
				metric_potential_poisson_mpi phi,hid_t filename,double *dc,double *fdc,double a0,double a,double a_t,double omega_dm_0,
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
	

	

	
	
	status_psi = psi.write_hdf5_values_mpi(filename, dtype, dspace_psi,dspace_dc,dc,fdc,a0,a,a_t,omega_dm_0,phi,cum_lin_id,get_dc);
	H5Sclose(dspace_psi);

	
	
	status_phi = phi.write_hdf5_values_mpi(filename, dtype, dspace_potn,a0,a,a_t,cum_lin_id,get_dc);
	H5Sclose(dspace_potn);	

		



	H5Sclose(dspace_dc);
	

	
	H5Tclose(dtype);


}







