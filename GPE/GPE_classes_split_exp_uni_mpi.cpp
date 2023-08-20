

///#include <omp.h>
#include <fftw3.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>
#include <hdf5.h>



class GPE_field_mpi
{
    private:
    int dim;
    double hbar_unit,c_unit,h,pc_unit;
    double vfac;

   
    
    double m_alpha;

    ptrdiff_t alloc_local, local_n0, local_0_start;

    public:
    double kppa,omega_m0,ai,ti,t0;
    double psi2_avg,fpsi2_avg;

     double nm_correction = 1.0;
     double fnm_correction = 1.0;

    fftw_complex *psi;

    fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

    fftw_complex *V_phi;
	fftw_complex *V_phi_ft;

    fftw_complex *theta;
	fftw_complex *theta_ft;



   

    //double ***im_K_psi = new double**[2];
	//double ***ex_K_psi = new double**[2];

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;

    fftw_plan plan_V_f;
    fftw_plan plan_V_b;


    fftw_plan plan_theta_f;
    fftw_plan plan_theta_b;


    int comp_j,N_tot,N,s;
    int unp_axis_N_tot;
    int myN_tot,myNx;
    int i,k;
    int cum_lin_ind,my_rank;

    double dj,dN_tot,mydN_tot;
    double energy,ini_energy,max_eng_err;
    double mass,ini_mass,max_mass_err;
    double mass_mpi_glbl;

    double vmax_mpi_glbl;

    int is_nm_correction;

    hid_t hdf5_file;  hid_t dtype; hid_t dset_glbl;hid_t dspace_glbl,memspace_loc;
    hsize_t *offset; hsize_t *hdf_dims,*hdf_dims_loc;
    
    MPI_Info info  = MPI_INFO_NULL;
    char my_name[20];
    char fp_hdf5_name[30]=("data_");



    GPE_field_mpi(int dim_p,int NN,int jj,int my_rank_p,char *field_name,int cor_fac=1,int nthreads=4)
    {
        comp_j= jj;
        N= NN;
        N_tot = N;
        my_rank=my_rank_p;

        is_nm_correction = cor_fac;

        
       dim = dim_p;

        ptrdiff_t *ptr_n = new  ptrdiff_t[dim];
      
        ptr_n[0] = N;
      
        unp_axis_N_tot = 1;
        for(i=0;i<(dim-1);++i)
        {
            myN_tot*=N;
            N_tot*=N;
            ptr_n[i+1] = N;
            unp_axis_N_tot*=N;
        

        }

        const ptrdiff_t *ptr_nc =  &ptr_n[0];
        


        alloc_local =  fftw_mpi_local_size(dim, ptr_nc,
                                     MPI_COMM_WORLD,&local_n0,&local_0_start);

        myNx = local_n0;
        myN_tot = myNx;
        cum_lin_ind = local_0_start;


        //hsize_t hdf_dims[ dim_p],hdf_dims_loc[ dim_p];
        offset = new hsize_t[dim+1];
        hdf_dims = new hsize_t[dim+1];
        hdf_dims_loc = new hsize_t[dim+1];
        offset[0] = cum_lin_ind;
        hdf_dims[0] = N;
        hdf_dims_loc[0] = myNx;


        for(i=0;i<(dim-1);++i)
        {
            myN_tot*=N;
            hdf_dims[i+1]=N;
            hdf_dims_loc[i+1]=N;
            offset[i+1] = 0;


        }

         offset[dim] = 0;
        hdf_dims[dim] = 2;
        hdf_dims_loc[dim] = 2;



       

        dj = (double)comp_j;
        dN_tot = (double)N_tot;
        mydN_tot= (double)myN_tot;

/////////////////////////   HDF5 things     ///////////////////////////////
       strcat(fp_hdf5_name,field_name);
        strcat(fp_hdf5_name,".hdf5"); 
        strcat(my_name,field_name);

       

        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

		
		hdf5_file = H5Fcreate (fp_hdf5_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

        dspace_glbl = H5Screate_simple(dim+1, hdf_dims, NULL);
        memspace_loc = H5Screate_simple(dim+1, hdf_dims_loc, NULL);

        dtype = H5Tcopy(H5T_NATIVE_DOUBLE);

  ////////////////////  FFT fields /////////////////////////////////////////////////////      

        fpGpsi = fftw_alloc_complex(alloc_local);// (fftw_complex*);fftw_malloc(sizeof(fftw_complex) * (N3));
	    fpGpsi_ft =fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        psi = fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    

	   // K = fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	   // K_ft =fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        V_phi = fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    V_phi_ft = fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        theta = fftw_alloc_complex(alloc_local);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    theta_ft = fftw_alloc_complex(alloc_local);

        


	    plan_pois_f = fftw_mpi_plan_dft(dim, ptr_nc,
                        fpGpsi, fpGpsi_ft ,MPI_COMM_WORLD,
                        FFTW_FORWARD, FFTW_ESTIMATE);//fftw_plan_dft_3d(N,N,N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	    plan_pois_b = fftw_mpi_plan_dft(dim, ptr_nc,
                        fpGpsi_ft, fpGpsi ,MPI_COMM_WORLD,
                        FFTW_BACKWARD, FFTW_ESTIMATE);// fftw_plan_dft_3d(N,N,N, fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

        plan_V_f = fftw_mpi_plan_dft(dim, ptr_nc,
                        V_phi, V_phi_ft ,MPI_COMM_WORLD,
                        FFTW_FORWARD, FFTW_ESTIMATE);//fftw_plan_dft_3d(N,N,N, V_phi, V_phi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	    plan_V_b = fftw_mpi_plan_dft(dim, ptr_nc,
                        V_phi_ft, V_phi ,MPI_COMM_WORLD,
                        FFTW_BACKWARD, FFTW_ESTIMATE);;//fftw_plan_dft_3d(N,N,N, V_phi_ft, V_phi, FFTW_BACKWARD, FFTW_ESTIMATE);

        plan_theta_b = fftw_mpi_plan_dft(dim, ptr_nc,
                        theta_ft, theta ,MPI_COMM_WORLD,
                        FFTW_BACKWARD, FFTW_ESTIMATE);;//fftw_plan_dft_3d(N,N,N, V_phi_ft, V_phi, FFTW_BACKWARD, FFTW_ESTIMATE);


        plan_theta_f = fftw_mpi_plan_dft(dim, ptr_nc,
                        theta, theta_ft ,MPI_COMM_WORLD,
                        FFTW_FORWARD, FFTW_ESTIMATE);//fftw_plan_dft_3d(N,N,N, V_phi, V_phi_ft, FFTW_FORWARD, FFTW_ESTIMATE);


      


    }


    double aoft(double t)
    {   double sinhval2;

        sinhval2 = sinh(1.5*sqrt(1.0-omega_m0)*t)*sinh(1.5*sqrt(1.0-omega_m0)*t);
        return(pow(omega_m0*sinhval2/(1.0-omega_m0),1.0/3.0));
        
    }

    double tofa(double a)
    {
        double sinarg;
        sinarg =  sqrt(a*a*a*(1.0-omega_m0)/omega_m0);

        return(2.0*asinh(sinarg/(3.0*sqrt(1.0-omega_m0))));


    }
    double HbyH0(double a,double a0=1.0)
    {

        return( sqrt(omega_m0*(a0/a)*(a0/a)*(a0/a) + (1.0 - omega_m0))    );


    }

    void calc_psi2_avg()
    {

       double psisqr,sum_ps2;  int loc_i;

        sum_ps2 = 0.0;

       
        for(loc_i=0;loc_i<myN_tot;++loc_i)
        {
            psisqr = psi[loc_i][0]*psi[loc_i][0] + psi[loc_i][1]*psi[loc_i][1];
            sum_ps2+= psisqr;


        }

       // MPI_Barrier(MPI_COMM_WORLD);

        MPI_Allreduce(&sum_ps2, &psi2_avg, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

        psi2_avg = psi2_avg/dN_tot;

        if(is_nm_correction)
        
            nm_correction = (3.0*omega_m0)/psi2_avg;

       


    }

    void calc_fpsi2_avg()
    {

       double psisqr,sum_ps2;  int loc_i;

        sum_ps2 = 0.0;

       
        for(loc_i=0;loc_i<myN_tot;++loc_i)
        {
            psisqr = fpGpsi[loc_i][0]*fpGpsi[loc_i][0] + fpGpsi[loc_i][1]*fpGpsi[loc_i][1];
            sum_ps2+= psisqr;


        }

       // MPI_Barrier(MPI_COMM_WORLD);

        MPI_Allreduce(&sum_ps2, &fpsi2_avg, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

        fpsi2_avg = fpsi2_avg/dN_tot;

        if(is_nm_correction)
        
            fnm_correction = (3.0*omega_m0)/fpsi2_avg;

       


    }



    void solve_V(double *kgrid)
    {
        double psisqr,ksqr;  int cur_i,loc_i,wii,ii,i2; 
        double dN = (double) N;
        double Npwr;

        calc_fpsi2_avg();

     //   #pragma omp parallel for private(psisqr)
        for(loc_i=0;loc_i<myN_tot;++loc_i)
        {
            psisqr = fpGpsi[loc_i][0]*fpGpsi[loc_i][0] + fpGpsi[loc_i][1]*fpGpsi[loc_i][1];
            //V_phi[loc_i][0] = (0.5/kppa)*(psisqr-3.0*omega_m0);
            V_phi[loc_i][0] = (1.5/kppa)*omega_m0*(psisqr/fpsi2_avg-1.0);
            V_phi[loc_i][1] = 0.0;


        }
        fftw_execute(plan_V_f);

        for(ii=0;ii<myN_tot;++ii)
        {
            ksqr= 0.0;
            wii = ii;
            for(i2=0;i2<dim;++i2)
            {
                Npwr = pow(dN,(double)(dim-i2-1));
                
                cur_i = (int) (((double)wii)/Npwr);
                
                wii = wii%((int)Npwr);
                if(i2==0)
                    ksqr+= (kgrid[cur_i+cum_lin_ind]*kgrid[cur_i+cum_lin_ind]);
                else
                    ksqr+= (kgrid[cur_i]*kgrid[cur_i]);

            }

            //ksqr = kgrid[loc_i]*kgrid[loc_i] + kgrid[loc_j]*kgrid[loc_j] + kgrid[loc_k]*kgrid[loc_k];

          if(ksqr>0.0)
           { V_phi_ft[ii][0] = -V_phi_ft[ii][0]/(ksqr*dN_tot);
                
             V_phi_ft[ii][1] = -V_phi_ft[ii][1]/(ksqr*dN_tot);
           }
            else
          {
            V_phi_ft[ii][0] = 0.0;
            V_phi_ft[ii][1] = 0.0;

          }


        }

        fftw_execute(plan_V_b);


        


    }



    void reset_consv_quant(int is_ini=0)
    {
        energy = 0.0;
        mass = 0.0;
        if(is_ini)
        {
            ini_energy=0.0;
            ini_mass = 0,0;

            max_eng_err=-1.0;
            max_mass_err=-1.0;

        }

    }


    void conserve_err(int is_ini=0)
    {
        double eng_err,mass_err;

        MPI_Allreduce(&mass, &mass_mpi_glbl, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

        if(is_ini)
        ini_mass = mass_mpi_glbl;
        
        eng_err = fabs(energy-ini_energy)/fabs(ini_energy);
        mass_err = fabs(mass_mpi_glbl-ini_mass)/ini_mass;

        if(eng_err>max_eng_err)
            max_eng_err=eng_err;
        
        if(mass_err>max_mass_err)
            max_mass_err=mass_err;


    }



    void split_step(double a,double da,double *kgrid,int acntr,bool cal_mass=true)
    {
        double styH,fac,V;
        double temp_psi[2];
       

        double ksqr;  int cur_i,loc_i,wii,ii,i2; 
        double dN = (double) N;
        double Npwr;

        int i;

        if(acntr==0)
            reset_consv_quant(1);

        if(cal_mass)
          mass = 0.0;

        styH = HbyH0(a);

        for(i=0;i<myN_tot;++i)
        {
            V = V_phi[i][0];
            fac = 0.5*da*V/(a*a*styH);

            if(cal_mass)
            {
                mass+= (fpGpsi[i][0]*fpGpsi[i][0]+fpGpsi[i][1]*fpGpsi[i][1]);


            }

            temp_psi[0] = fpGpsi[i][0];
            temp_psi[1] = fpGpsi[i][1];
            fpGpsi[i][0] = cos(fac)*fpGpsi[i][0] + sin(fac)*fpGpsi[i][1];
            fpGpsi[i][1] = -sin(fac)*temp_psi[0] + cos(fac)*fpGpsi[i][1];

            
        
            
        
           
            
        }

        fftw_execute(plan_pois_f);
        conserve_err(!acntr);


        for(i=0;i<myN_tot;++i)
        {
            ksqr = 0.0;
            wii = i;
            for(i2=0;i2<dim;++i2)
            {
                Npwr = pow(dN,(double)(dim-i2-1));
                
                cur_i = (int) (((double)wii)/Npwr);
                
                wii = wii%((int)Npwr);
                if(i2==0)
                    ksqr+= (kgrid[cur_i+cum_lin_ind]*kgrid[cur_i+cum_lin_ind]);
                else
                    ksqr+= (kgrid[cur_i]*kgrid[cur_i]);

            }
 
            fac = 0.5*kppa*da*ksqr/(a*a*a*styH);

            temp_psi[0] = fpGpsi_ft[i][0];
            temp_psi[1] = fpGpsi_ft[i][1];
            fpGpsi_ft[i][0] = (cos(fac)*fpGpsi_ft[i][0] + sin(fac)*fpGpsi_ft[i][1])/dN_tot;
            fpGpsi_ft[i][1] = (-sin(fac)*temp_psi[0] + cos(fac)*fpGpsi_ft[i][1])/dN_tot;

            
        }

        fftw_execute(plan_pois_b);

        solve_V(kgrid);

        

        for(i=0;i<myN_tot;++i)
        {
            V = V_phi[i][0];
            fac = 0.5*da*V/(a*a*styH);
            temp_psi[0] = fpGpsi[i][0];
            temp_psi[1] = fpGpsi[i][1];
            fpGpsi[i][0] = cos(fac)*fpGpsi[i][0] + sin(fac)*fpGpsi[i][1];
            fpGpsi[i][1] = -sin(fac)*temp_psi[0] + cos(fac)*fpGpsi[i][1];

            psi[i][0] = fpGpsi[i][0];
            psi[i][1] = fpGpsi[i][1];

            

            
        }
        
        



    }


   

    void read_param_from_file(char  *fname)
    {
        FILE *fp_param = fopen(fname,"r");

        int j,i=0;
        size_t fnd_div;
        double num_val;
        int max_len = 80;
	    char read_c_line[max_len];
	    string read_str,cur_str, read_key, read_val,cur_num,cur_den;

        while(fgets(read_c_line,max_len,fp_param)!=NULL)
        {  read_str = read_c_line; //printf("re str  %s  %d\n",read_str.c_str(),read_str.length());
            read_key = read_str.substr(0,read_str.find("="));
		  remove_if(read_key.begin(),read_key.end(),::isspace);
		  read_key.erase(remove_if(read_key.begin(),read_key.end(),::isspace),read_key.end());
		 

		 if(read_key.length()>0)
		 {  read_val = read_str.substr(read_str.find("=")+1);
		   

		   
		
		   if(read_key=="c_unit")
			{
                num_val = stod(read_val);
                c_unit = num_val;

            }
		   else if (read_key=="hbar_unit")
           {
               num_val = stod(read_val);
               hbar_unit = num_val;
           }
           else if (read_key=="pc_unit")
           {
               num_val = stod(read_val);
               pc_unit = num_val;
           }
           else if (read_key=="h")
           {
               num_val = stod(read_val);
               h = num_val;
           }
           else if (read_key=="m_alpha")
           {
               num_val = stod(read_val);
               m_alpha = num_val;
           }
           
           else if (read_key=="omega_m0")
           {
               num_val = stod(read_val);
               omega_m0= num_val;
           }
		   else if (read_key=="ai")
           {
               num_val = stod(read_val);
               ai= num_val;
           }


         }

    



        }
    }

    void print_params_set_kappa()
    {

        printf("\nPrinting params for field %d\n",comp_j+1);
        printf("[c_unit,hbar_unit,pc_unit] = [%lf,%lf,%lf]\n",c_unit,hbar_unit,pc_unit);
        printf("[h,omega_m0,m_alpha] = [%lf,%lf,%lf] and ai = %lf \n",h,omega_m0,m_alpha,ai);
        kppa = c_unit*c_unit*hbar_unit*h*(1e-5)/(m_alpha*pc_unit);

        ti = tofa(ai);
        t0 = tofa(1.0);

        printf("kappa = %lf\n",kppa);
        printf("t_ini = %lf  t_end = %lf\n\n",ti,t0);
    }

    
    

    void set_field()
    {int ii;
        for(ii=0;ii<myN_tot;++ii)
        {

            fpGpsi[ii][0] = psi[ii][0];
            fpGpsi[ii][1] = psi[ii][1];

        }

        fftw_execute(plan_pois_f);

        for(ii=0;ii<myN_tot;++ii)
        {

            fpGpsi_ft[ii][0] = fpGpsi_ft[ii][0]/(dN_tot);
            fpGpsi_ft[ii][1] = fpGpsi_ft[ii][1]/(dN_tot);

        }

        fftw_execute(plan_pois_b);

        

    }


double calc_v(double *kgrid)
{
    int dir;
    double dir_max,vmax_loc=-1.0;

    for(dir=0;dir<dim;++dir)
    {  
        dir_max = calc_v_dir(dir,kgrid);
        
        if(fabs(dir_max)>vmax_loc)
            vmax_loc = fabs(dir_max);

        
    }
    
     MPI_Allreduce(&vmax_loc, &vmax_mpi_glbl, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
   return(vmax_mpi_glbl);

}


double calc_v_dir(int dir,double *kgrid)
{
    int ii,Npwr,wii,cur_i,kcomp;
    double num,den;
    double dN = (double) N;
     
    for(ii=0;ii<myN_tot;++ii)
    {
        den = sqrt(fpGpsi[ii][0]*fpGpsi[ii][0]+fpGpsi[ii][1]*fpGpsi[ii][1]);
        if(den>0.0)
        {

            num = fpGpsi[ii][0];
            theta[ii][0] = acos(num/den);

        }
        else
        theta[ii][0]=0.0;

        theta[ii][1]=0.0;

    }

    fftw_execute(plan_theta_f);

    
       
            
    Npwr = pow(dN,(double)(dim-dir-1));
    

    printf("DGGDHD %d %d\n",dir,Npwr);

    
    for(ii=0;ii<myN_tot;++ii)
    {
                
            cur_i = (int) (((double)ii)/Npwr);
                
           
                if(dir==0)
                    kcomp= (kgrid[cur_i+cum_lin_ind]);
                else
                    kcomp= (kgrid[cur_i]);

            theta_ft[ii][0] = -kcomp*theta_ft[ii][1]/(dN_tot);
            theta_ft[ii][1] =  kcomp*theta_ft[ii][0]/(dN_tot);
    }

    fftw_execute(plan_theta_b);
    
    
    double loc_max = 0.0;
    for(ii=0;ii<myN_tot;++ii)
    {
        
        if(fabs(theta[ii][0])>loc_max)
            loc_max = fabs(theta[ii][0]);

    }


    
    return loc_max;

   



}


void read_from_initial()
{
    double a[10],fr,fi;

    FILE *fp = fopen("initial1.txt","r");
    int prev = cum_lin_ind*unp_axis_N_tot;
    int j;
    printf("prev is %d %d %d %d\n",prev,unp_axis_N_tot*my_rank*32,cum_lin_ind,my_rank);
    for(i=0;i<N_tot;++i)
    {
        j = i-prev;
        fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&a[0],&a[1],&a[2],&a[3],&fr,&fi,&a[4],&a[5],&a[6],&a[7],&a[8]);
        if(i>=prev)
        {psi[j][0] = fr/100.0;//H0=100 in given units and it is a unit conversion
         psi[j][1] = fi/100.0;
        }

    }



}


	herr_t write_hdf5_mpi(double z)
	{
        hid_t plist_id;
        herr_t status;
        char dset_name[20];
        snprintf(dset_name,20,"%.2lf",z);
        strcat(dset_name,"_");
        strcat(dset_name,my_name);


        dset_glbl = H5Dcreate(hdf5_file, dset_name, dtype, dspace_glbl,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Sselect_hyperslab(dspace_glbl, H5S_SELECT_SET,offset,NULL,hdf_dims_loc,NULL);


        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_glbl, dtype, memspace_loc, dspace_glbl,
		      				plist_id, &psi[0][0]);


        H5Dclose(dset_glbl);
        H5Pclose(plist_id);

        return status;

        

       

	}





   

};







    

/*
int main(int argc,char **argv)
{


	char fp_name[30]("imex_ft_");
	
    int stages;
	char *imex_file; 
    char *f1file;

    stages = atoi(argv[1]);
	imex_file = argv[2];
    f1file = argv[3];
	
	strcat(fp_name,imex_file);
	printf("ImEx table filename is %s and stages given by u is %d and out file %s \n",imex_file,stages,fp_name);

	FILE *fp = fopen(fp_name,"w");
	

	imex_table imx(stages);
    
 
  	imx.read_from_file(imex_file);
  	imx.print_table();


    GPE_field_2d  psi_1(0,128,3,8);
    GPE_field_2d  psi_2(1,128,3,8);

    psi_1.read_from_file(f1file);

    psi_1.print_params();

}*/
