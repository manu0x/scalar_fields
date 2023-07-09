


//////////////////////////////////////////////////////////////



void read_psi_from_hdf5(string fname,fftw_complex *psi,int *ind,int use_zeldo=1)
{

	herr_t status;	
	hid_t file,dtype,dspace_glbl,dspace,dset_glbl,plist_id,g_id;
	

	hsize_t dim[3],odim[2];

	dim[0] = ind[0];
	dim[1] = ind[1];
	dim[2] = ind[2];
	odim[0] = ind[0]*ind[1]*ind[2];
	odim[1] = 2;

	
	
	char fchar_name[32];
	fname.copy(fchar_name,fname.length(),0);
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









