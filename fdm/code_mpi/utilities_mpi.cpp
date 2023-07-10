

void calc_theta_zeldo(int * ind,int *ind_loc, double *k_grid,double  *dc,double *theta,double hbar_by_m,double ai,double a_t,double f_ini=1.0,double a0=1.0)
{
	double H_ini = a_t/ai;
	int n[3],n_loc[3],i,j,k,ci;
	double k2fac;

	fftw_complex *fpth;
	fftw_complex *fpth_ft;

    

	
	
	n[0]=ind[0];n[1]=ind[1];n[2]=ind[2];
	n_loc[0]=ind_loc[0];n_loc[1]=ind_loc[1];n_loc[2]=ind_loc[2];


	const ptrdiff_t n0 = n[0];
	const ptrdiff_t n1 = n[1];
	const ptrdiff_t n2 = n[2];

	ptrdiff_t alloc_local, local_n0, local_0_start;


	double dtN = (double)(n[0]*n[1]*n[2]);
	

	fftw_plan plan_f_th;
	fftw_plan plan_b_th;

	alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
                                 MPI_COMM_WORLD ,
                                 &local_n0, &local_0_start);
		
	fpth = fftw_alloc_complex(alloc_local);
	fpth_ft = fftw_alloc_complex(alloc_local);


	plan_f_th = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               fpth, fpth_ft,
                              MPI_COMM_WORLD , FFTW_FORWARD, FFTW_ESTIMATE);
	plan_b_th = fftw_mpi_plan_dft_3d(n0, n1,  n2,
                               fpth_ft, fpth,
                             MPI_COMM_WORLD , FFTW_BACKWARD, FFTW_ESTIMATE);

	for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;

			

			fpth[ci][0] = -H_ini*f_ini*dc[ci]*(ai/a0)*(ai/a0)/hbar_by_m;
			fpth[ci][1] = 0.0;	

		    }
		  }

		}

	fftw_execute(plan_f_th);

	for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;
			k2fac = (k_grid[i]*k_grid[i]+k_grid[i]*k_grid[i]+k_grid[i]*k_grid[i]);

			 if(k2fac>0.0)
			 {fpth_ft[ci][0] = -fpth_ft[ci][0]/k2fac;
			  fpth_ft[ci][1] = -fpth_ft[ci][1]/k2fac;	

			 }

			 else			
			 {fpth_ft[ci][0] = 0.0;
			  fpth_ft[ci][1] = 0.0;	

			 }


			fpth_ft[ci][0] = fpth_ft[ci][0]/dtN;
			fpth_ft[ci][1] = fpth_ft[ci][1]/dtN;

		    }
		  }

		}

	fftw_execute(plan_b_th);


	for(i=0;i<n_loc[0];++i)
		{
		  for(j=0;j<n_loc[1];++j)
		  {
		    for(k=0;k<n_loc[2];++k)
		    {
			ci = (n_loc[2]*n_loc[1])*i + n_loc[2]*j + k;
			

			theta[ci] = fpth[ci][0];

		    }
		  }

		}

	

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

	
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



void read_psi_from_hdf5(string fname,fftw_complex *psi,int *ind,int use_zeldo=1)
{

	herr_t status;	
	hid_t file,dtype,dspace_glbl,dspace,dset_glbl,plist_id,g_id;
	MPI_Info info  = MPI_INFO_NULL;

	hsize_t dim[3],odim[2];

	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];
	odim[0] = ind[0]*ind[1]*ind[2];
	odim[1] = 2;

	
	
	char fchar_name[fname.length()];
	fname.copy(fchar_name,fname.length()-1,0);
	printf("FINI_dc %s %d\n",fchar_name,fname.length());
	
	file = H5Fopen(fchar_name, H5F_ACC_RDONLY, H5P_DEFAULT);
	
	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(dtype, H5T_ORDER_LE);


	//dspace_glbl = H5Screate_simple(2, odim, NULL);

	//g_id = H5G.open(file, "/");
	if(use_zeldo)
		dset_glbl = H5Dopen(file, "/psi_zeldo",H5P_DEFAULT);
	else
		dset_glbl = H5Dopen(file, "/psi",H5P_DEFAULT);

	dspace_glbl = H5Dget_space(dset_glbl);

	hsize_t count[2],offset[2];
	count[0] = ind[0]*ind[1]*ind[2]; count[1] = odim[1]; 




	

	status = H5Dread(dset_glbl, dtype,H5S_ALL, H5S_ALL, H5P_DEFAULT, psi);

	H5Dclose(dset_glbl);
	H5Fclose(file);
	
	







}



void read_psi_from_hdf5_mpi(string fname,fftw_complex *psi,int *ind_loc, int cum_lin_ind,int use_zeldo=1)
{

	herr_t status;	
	hid_t file,dtype,dspace_glbl,dspace,dset_glbl,plist_id,g_id;
	MPI_Info info  = MPI_INFO_NULL;

	hsize_t dim[3],odim[2];

	dim[0] = ind_loc[0];
	dim[1] = ind_loc[1];
	dim[2] = ind_loc[2];
	odim[0] = ind_loc[0]*ind_loc[1]*ind_loc[2];
	odim[1] = 2;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

	
	char fchar_name[fname.length()];
	fname.copy(fchar_name,fname.length(),0);
	printf("FINI_dc %s %d\n",fchar_name,fname.length());
	
	file = H5Fopen(fchar_name, H5F_ACC_RDONLY, plist_id);
	H5Pclose(plist_id);

	
	dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(dtype, H5T_ORDER_LE);


	//dspace_glbl = H5Screate_simple(2, odim, NULL);

	//g_id = H5G.open(file, "/");
	if(use_zeldo)
		dset_glbl = H5Dopen(file, "/psi_zeldo",H5P_DEFAULT);
	else
		dset_glbl = H5Dopen(file, "/psi",H5P_DEFAULT);

	dspace_glbl = H5Dget_space(dset_glbl);

	hsize_t count[2],offset[2];
	count[0] = ind_loc[0]*ind_loc[1]*ind_loc[2]; count[1] = odim[1]; offset[0] = cum_lin_ind*ind_loc[1]*ind_loc[2]; offset[1] = 0;

	

	dspace = H5Screate_simple(2, count, NULL);

	

	H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET, offset, NULL, count, NULL);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	

	status = H5Dread(dset_glbl, dtype, dspace, dspace_glbl, plist_id, psi);

	H5Dclose(dset_glbl);
	H5Pclose(plist_id);
	H5Sclose(dspace);







}







////////////////////////////////////////////////////////////////	WRITING MODULES		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	herr_t write_hdf5_mpi(char *file,int dim,hsize_t *offset, hsize_t *hdf_dims,hsize_t *hdf_dims_loc,char *my_name,double *data)
	{
        hid_t hdf5_file;  hid_t dtype; hid_t dset_glbl; hid_t dspace_glbl,memspace_loc;
        herr_t status;
        char dset_name[20];

		MPI_Info info  = MPI_INFO_NULL;
        strcat(dset_name,my_name);


        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

		
		hdf5_file = H5Fcreate (file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

        dspace_glbl = H5Screate_simple(dim+1, hdf_dims, NULL);
        memspace_loc = H5Screate_simple(dim+1, hdf_dims_loc, NULL);

        dtype = H5Tcopy(H5T_NATIVE_DOUBLE);


        dset_glbl = H5Dcreate(hdf5_file, dset_name, dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET,offset,NULL,hdf_dims_loc,NULL);


        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_glbl, dtype, memspace_loc, dspace_glbl,
		      				plist_id, data);


        H5Dclose(dset_glbl);
        H5Pclose(plist_id);
        H5Fclose(hdf5_file);

        return status;

        

       

	}