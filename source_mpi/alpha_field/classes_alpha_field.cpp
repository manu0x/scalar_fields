


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
	   f_t = new double** [nx] ;
	   if(need_lap)
	   f_lap = new double** [n[0]] ;
	   if(need_space_grads)
	   { f_x = new double *** [3];
	     f_x[0] = new double** [n[0]];
	     f_x[1] = new double** [n[0]];
	     f_x[2] = new double** [n[0]];  
	   }

	   for(i=0;i<(nx);++i)
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
	if((need_lap)||(need_space_grads))
	  {	
		if((need_lap)&&(need_space_grads))
 		 printf("Field allocated with arrays for space der and laplacian %d %d %d\n",n[0],n[1],n[2]);
		else
		  if(need_lap)
			printf("Field allocated with array for laplacian %d %d %d %d\n",n[0],n[1],n[2],nx);	
		  else
			cout<<"Field allocated with arrays for space der\n";

	 }
	else
	   cout<<"Field allocated withOUT arrays for space der and laplacian\n";
		

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
	     
	     if(i>0)
	     {	ind_l1[i] = (n[i]+ind_lw[i]-1)%n[i];
	     	ind_l2[i] = (n[i]+ind_lw[i]-2)%n[i];

	     	ind_r1[i] = (ind_lw[i]+1)%n[i];
	     	ind_r2[i] = (ind_lw[i]+2)%n[i];
	      }
	    else

	      {	ind_l1[i] = (n[i]+ind_lw[i]-1);
	     	ind_l2[i] = (n[i]+ind_lw[i]-2);

	     	ind_r1[i] = (ind_lw[i]+1);
	     	ind_r2[i] = (ind_lw[i]+2);
	      }
		
	    
	     if(spt_grad==true)
	     {
	       m[0] = (-8.0*f[ind_l1[0]][ind_l1[1]][ind_l1[2]]+8.0*f[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	       m[1] = (f[ind_l2[0]][ind_l2[1]][ind_l2[2]]-f[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	       m[2] = m[0] + m[1] ;
		
		
	       f_x[i][ind[0]][ind[1]][ind[2]] = m[2]/(dx[i]);
	      }

	     if(laplacian==true)
	      {	m[0] = (16.0*f[ind_l1[0]][ind_l1[1]][ind_l1[2]]+16.0*f[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	      	m[1] = (-f[ind_l1[0]][ind_l1[1]][ind_l1[2]]-f[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	      	m[2] = m[0] + m[1] -30.0*f[ind_lw[0]][ind_lw[1]][ind_lw[2]];
	     	lapsum+= (m[2]/(dx[i]));

	     	
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


		mpi_check = MPI_Isend(&f[0][0][0],2, c_x_plain, left_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&f[n[0]+2][0][0], 2, c_x_plain,right_rank , recvtag,cart_comm, &recv_req);





	/*	mpi_check = MPI_Sendrecv_replace(&f[0][0][0], 2, c_x_plain,
                         left_rank, sendtag, left_rank, recvtag,
                         cart_comm, &status);

		

		mpi_check = MPI_Sendrecv_replace(&f[n[0]+2][0][0], 2, c_x_plain,
                         right_rank, sendtag, right_rank, recvtag,
                        cart_comm, &status);
*/
		printf("sending done for cart rank %d %d\n",my_cart_rank,n[0]);
		
		return(mpi_check);
	}


};



class field_alpha_mpi
{
	private:
	scalar_field_3d_mpi f_alpha;
	scalar_field_3d_mpi f_t_alpha;
	int n[3];
	
	public:  
	//int n[3];
	field_alpha_mpi(int *ind,int cum_lin_ind,bool lb=false,bool sgb=false):f_alpha(ind,cum_lin_ind,lb,sgb),f_t_alpha(ind,cum_lin_ind,lb,sgb)
	{
		n[0] = ind[0];  n[1] = ind[1];  n[2] = ind[2];
	}
	
	void test_ind()
	{
		f_alpha.test_ind();

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
						psi_r_val=f_alpha.get_field(locind,give_f);

				}
			}
		}

		printf("\nTTTTTTTTTESSST 222 Done\n");
	

	}

	int calc_acc(int * ind,double &acc,double potn,double potn_t,double potn_x[3],double a,double a_t,double *dx)
	{
		int c1;		
		double fa,fa_t,f_val,fa_t_val,fa_lap,fa_der[3],fa_t_der[3];

		//fa_val = f_alpha.get_field(ind,give_f);
		fa_t_val = f_alpha_t.get_field(ind,give_f);

		c1 = f_alpha.cal_spt_grads(ind,dx,true);
		c1 = f_alpha_t.cal_spt_grads(ind,dx,true);

		fa_lap = f_alpha.get_field(ind,give_f_lap);
		c1 = f_alpha.get_field_spt_der(ind,fa_der);
		c1 = f_alpha_t.get_field_spt_der(ind,fa_t_der);  

		
		acc = field_acc_approx(fa_t,fa_der, fa_t_der,potn,potn_t, potn_x);


		if(isnan(acc))
			return (-1);
		else
			return(1);
	}

	

	int update(int * ind,double fa,double fa_t)
	{
		int c1,c2;
		c1 = f_alpha.update_field(ind,fa);
		c2 = f_alpha_t.update_field(ind,fa_t);
		
		return (c1*c2);
	}

	int get_field_alpha(int *ind,double *psi_ret,code1 c = give_f)
	{
		int c1=1;		
		
		psi_ret[0]= f_alpha.get_field(ind,c);
		psi_ret[1] = f_alpha_t.get_field(ind,c);

		return c1;




	}


	double cal_X_4vel(int * ind,double a,double phi)
	{
		double fa_val,fa_t_val,s_der[3],X;	
		int c1;	
		
		
		fa_t_val = f_alpha_t.get_field(ind,give_f);

		c1 = get_field_spt_der(ind,s_der);

		
		X = fa_t_val*fa_t_val/(1.0+phi)  - (s_der[0]*s_der[0]+s_der[1]*s_der[1]+s_der[2]*s_der[2])/(a*a*(1.0+phi));
		X = 0.5*X;

		return(X);

	}

	int mpi_send_recv()
	{
		int mpi_check;
		mpi_check=f_alpha.mpi_send_recv();
		mpi_check=f_alpha_t.mpi_send_recv();
		return(mpi_check);


	}



	herr_t write_hdf5_f_alpha_mpi(hid_t filename,hid_t dtype,hid_t dspace_glbl,double *dc,double a3a03omega,double a,double phi,int cum_lin_ind,bool get_dc=false)
	{

		hid_t dataset,dset_glbl_r,dset_glbl_i,dspace,plist_id;
		herr_t status ;
		int tN  = n[0]*n[1]*n[3];
		int i,j,k,locind[3],ci;
		double fa_val,fa_t_val,x4val,rho_fa;

		

		

    /*
     * Create the dataset with default properties and close filespace.
     */
   		dset_glbl_r = H5Dcreate(filename, "field_alpha", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dset_glbl_i = H5Dcreate(filename, "field_alpha_t", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    		
		

		status = f_alpha.write_hdf5_mpi( filename,dtype,dset_glbl_r);
		
		
		status = f_alpha_t.write_hdf5_mpi( filename,dtype,dset_glbl_i);

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
					
						
					

					x4val = cal_X_4vel(locind,a,phi);	
							

					rho_fa = (2.0*alpha-1.0)*x4val*pow(x4val/(m*m*m*m),alpha-1.0);
					 

				}
	
			}

		  }
		


		dataset = H5Dcreate(filename, "dc_alpha_field", dtype, dspace_glbl,
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



	

	



};


class metric_potential_approx_1_t
{
	private:
	scalar_field_3d_mpi potn;
	
	int n[3];
	
	public:  
	//int n[3];
	metric_potential_approx_1_t(int *ind,int cum_lin_ind,bool lb=false,bool sgb=false):potn(ind,cum_lin_ind,lb,sgb)
	{
		n[0] = ind[0];  n[1] = ind[1];  n[2] = ind[2];
	}
	


	int calc_vel(int * ind,double &potn_vel,double f_t,double a,double a_t,double *dx)
	{
		int c1;		
		double potn_val,lap_potn;
		
		
		c1 = potn.cal_spt_grads(ind,dx,true);
		lap_potn = potn.get_field(ind,give_f_lap);
		potn_val = potn.get_field(ind,give_f);
		
		
		potn_vel = potn_approx_vel(a,a_t,potn_val,lap_potn,f_t);
		
		
		if(isnan(potn_vel))
			return (-1);
		else
			return(1);
	}


	int update(int * ind,double potn_val)
	{
		int c1;
		c1 = potn.update_field(ind,potn_val);
		
		
		return (c1);
	}

	double get_potn(int *ind,code1 c = give_f)
	{
		double potn_val;		
		
		potn_val= potn.get_field(ind,c);
		

		return potn_val;




	}


	int get_potn_spt_der(int *ind,double der[3])
	{
		int c1;
		
	
		
		c1 = potn.get_field_spt_der(ind,der); 
		
		retur c1;
	



	}

	int mpi_send_recv()
	{
		int mpi_check;

		mpi_check=potn.mpi_send_recv();
	
		return(mpi_check);


	}


	

	herr_t write_hdf5_potn_mpi(hid_t filename,hid_t dtype,hid_t dspace_glbl)
	{

		hid_t dataset,dset_glbl_potn,plist_id;
		herr_t status ;
		
		
		

		

    /*
     * Create the dataset with default properties and close filespace.
     */
   		dset_glbl_potn = H5Dcreate(filename, "potn_alpha", dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		

		status = potn.write_hdf5_mpi( filename,dtype,dset_glbl_potn);
		
		
		

		H5Dclose(dset_glbl_potn);
		
		return(status);

	}


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
		printf("	doing spline..\n");
		
		checkspl = spline(k,p,point_cnt,cload); 
		if(checkspl)
		 printf("\nALERT !!!! Spline fit failed....checkspl %d  %.10lf\n",checkspl,max_dx);
		else
		 printf("	spline done successfully..\n");
	
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
	printf("\nALERT index overrun for spline...for x %lf\n",x);

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







};
