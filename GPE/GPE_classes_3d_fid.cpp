

///#include <omp.h>
#include <fftw3.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>
#include <hdf5.h>



class GPE_field_3d
{
    private:
    
    double hbar_unit,c_unit,h,pc_unit;
    double vfac;
    
    double m_alpha;

    public:
    double kppa,omega_m0,ai,ti,t0,omg,T,S,kfsqr,kf[3];
    double err_n2,sol_n2,sol_err,vmax;

    fftw_complex *psi;

    fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

    fftw_complex *V_phi;
	fftw_complex *V_phi_ft;

	fftw_complex *K;
	fftw_complex *K_ft;

   

    double ***im_K_psi = new double**[2];
	double ***ex_K_psi = new double**[2];

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;

    fftw_plan plan_V_f;
    fftw_plan plan_V_b;



    int j,N,N2,N3,s;
    int i,k;

    double dj,dN3;
    double energy,ini_energy,max_eng_err;
    double mass,ini_mass,max_mass_err;
    double max_sol_err;



    GPE_field_3d(int jj,int NN,int imex_s,int nthreads=4)
    {
        j= jj;
        N= NN;
        N2 = N*N;
        N3 = N*N*N;
        s = imex_s;
        dj = (double)j;
        dN3 = (double)N3;

        vmax = -0.00001;

        fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    

	    K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        V_phi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    

        fftw_plan_with_nthreads(nthreads);


	    plan_pois_f = fftw_plan_dft_3d(N,N,N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	    plan_pois_b = fftw_plan_dft_3d(N,N,N, fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

        plan_V_f = fftw_plan_dft_3d(N,N,N, V_phi, V_phi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	  


	    plan_imp_b = fftw_plan_dft_3d(N,N,N, K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);
        

        for(i=0;i<2;++i)
	    {
		    im_K_psi[i] = new double*[imex_s];
		    ex_K_psi[i] = new double*[imex_s];

		    for(k=0;k<imex_s;++k)
		    {im_K_psi[i][k] =  new double[N3];
		     ex_K_psi[i][k] =  new double[N3];
		    }





	    }


    }

    double Vfunc(double x[3])
    {
       return(0.5*omg*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]/3.0));
       //return(1.0);


    }

    void sol(double t,double x[3],double solval[2])
    {

        double theta,gmma;
        theta = 2.0*M_PI*t/T;
        gmma = 2.0*M_PI*(kf[0]*x[0]+kf[1]*x[1]+kf[2]*x[2])/S;
        solval[0] = cos(theta)*cos(gmma) - sin(theta)*sin(gmma);
        solval[1] = sin(gmma)*cos(theta) + cos(gmma)*sin(theta);



    }

    double gfunc(double t,double x[3],double *gval)
    {

        double bta,solval[2];
        
        bta = (2.0*M_PI/T)   + 0.5*kppa*(2.0*M_PI/S)*(2.0*M_PI/S)*(kfsqr) + Vfunc(x);
        
        sol(t,x,solval);
        *(gval) = -bta*solval[1];
        *(gval+1) = bta*solval[0];

        return(bta);

    }

    


   
    void update_fft_fields(int ind,double ksqr,double lamda)
    {
        double temp_psi[2];
        temp_psi[0] = fpGpsi_ft[ind][0];
        temp_psi[1] = fpGpsi_ft[ind][1];
       
       	fpGpsi_ft[ind][0] = ((fpGpsi_ft[ind][0]) + lamda*fpGpsi_ft[ind][1])/(1.0+lamda*lamda);
		fpGpsi_ft[ind][1] = ((fpGpsi_ft[ind][1]) - lamda*temp_psi[0])/(1.0+lamda*lamda);
        


		fpGpsi_ft[ind][0] = fpGpsi_ft[ind][0]/((double)N3);
		fpGpsi_ft[ind][1] = fpGpsi_ft[ind][1]/((double)N3);

       
		K_ft[ind][0] = -fpGpsi_ft[ind][0]*ksqr;
		K_ft[ind][1] = -fpGpsi_ft[ind][1]*ksqr;
       
       




    }

    double ex_rhs(int ind,int stg_i,double t,double *x,double expnfac=1.0)
    {
            double gval[2]; double xv[3]; double l1,l2,bta;
            xv[0] = *(x); xv[1] = *(x+1);xv[2] = *(x+2);
            bta = gfunc(t,xv,gval);
            V_phi[ind][0] = Vfunc(xv)*(fpGpsi[ind][0]*fpGpsi[ind][0]+fpGpsi[ind][1]*fpGpsi[ind][1]);
            

            ex_K_psi[0][stg_i][ind] =  (V_phi[ind][0]-bta)*fpGpsi[ind][1];///expnfac  ;
            ex_K_psi[1][stg_i][ind] = -(V_phi[ind][0]-bta)*fpGpsi[ind][0];///expnfac ;

            if(fabs(V_phi[ind][0])>vmax)
                vmax = fabs(V_phi[ind][0]);

           // printf("%d %lf %lf %lf\n",ind,*(x),*(x+1),*(x+2));

           return(fabs(bta-V_phi[ind][0]));


            

            



    }

    void im_rhs(int ind,int stg_i,double expnfac=1.0)
    {
        im_K_psi[0][stg_i][ind]  = -K[ind][1]*0.5*kppa/expnfac; 
        im_K_psi[1][stg_i][ind]=  K[ind][0]*0.5*kppa/expnfac;


    }

    void do_forward_fft()
    {
        fftw_execute(plan_pois_f);
       // fftw_execute(plan_V_f);


    }

    void do_back_fft()
    {
        fftw_execute(plan_pois_b);
		fftw_execute(plan_imp_b);
       // fftw_execute(plan_V_b);

    }

    void read_from_file(char  *fname)
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

        printf("\nPrinting params for field %d\n",j+1);
        printf("[c_unit,hbar_unit,pc_unit] = [%lf,%lf,%lf]\n",c_unit,hbar_unit,pc_unit);
        printf("[h,omega_m0,m_alpha] = [%lf,%lf,%lf] and ai = %lf \n",h,omega_m0,m_alpha,ai);
        kppa = c_unit*c_unit*hbar_unit*h*(1e-5)/(m_alpha*pc_unit);
        omg =  -1.0/kppa;
        S = 0.5;
        T = 0.5;

        kf[0] = 4.0/(2.0*S); kf[1] = 4.0/(2.0*S); kf[2] = 4.0/(2.0*S);

        kfsqr = kf[0]*kf[0] + kf[1]*kf[1] + kf[2]*kf[2];

        ti = 0.0;//tofa(ai);
        t0 = 1.0;//tofa(1.0);

        printf("kappa = %lf\n",kppa);
        printf("t_ini = %lf  t_end = %lf\n\n",ti,t0);
    }

    

    void set_field()
    {int ii;
        for(ii=0;ii<N3;++ii)
        {

            fpGpsi[ii][0] = psi[ii][0];
            fpGpsi[ii][1] = psi[ii][1];

        }

        fftw_execute(plan_pois_f);

        for(ii=0;ii<N3;++ii)
        {

            fpGpsi_ft[ii][0] = fpGpsi_ft[ii][0]/(dN3);
            fpGpsi_ft[ii][1] = fpGpsi_ft[ii][1]/(dN3);

        }

        fftw_execute(plan_pois_b);

        

    }

    void reset_consv_quant(int is_ini=0)
    {
        energy = 0.0;
        mass = 0.0;
        err_n2 = 0.0;
        sol_n2 = 0.0;
        if(is_ini)
        {
            ini_energy=0.0;
            ini_mass = 0,0;

            max_eng_err=-1.0;
            max_mass_err=-1.0;

            max_sol_err = -1.0;

        }

    }
    void cal_conserve_at_point(int ci,double t,double xv[3],double dx,int is_ini=0)
    {
           
        
            double loc_energy,loc_mass,psi2,solv[2];
            ///////////////////////////////////////////////////
            psi2 = psi[ci][0]*psi[ci][0] + psi[ci][1]*psi[ci][1];
            
            sol(t,xv,solv);
         
      

           loc_energy = 0.0;// 0.5*der_amp2 +  V(x)*psi2 + fcntr + Rcntr;
            loc_mass = psi2;

            energy+=loc_energy;
            mass+=loc_mass;
            (err_n2)+=((solv[0]-psi[ci][0])*(solv[0]-psi[ci][0])+(solv[1]-psi[ci][1])*(solv[1]-psi[ci][1]));
			sol_n2+= (solv[0]*solv[0] + solv[1]*solv[1]);

            if(is_ini)
            {
                ini_energy+=loc_energy;
                ini_mass+=loc_mass;

            }

            
    }

    void conserve_err()
    {
        double eng_err,mass_err;
        eng_err = fabs(energy-ini_energy)/fabs(ini_energy);
        mass_err = fabs(mass-ini_mass)/ini_mass;

        sol_err = err_n2/sol_n2;

        if(sol_err>max_sol_err)
            max_sol_err = sol_err;

        if(eng_err>max_eng_err)
            max_eng_err=eng_err;
        
        if(mass_err>max_mass_err)
            max_mass_err=mass_err;


    }



void initialise_point(int ind,double x[3],double t=0.0)
{
    double solv[2];

    sol(t,x,solv);
    psi[ind][0] = solv[0];
    psi[ind][1] = solv[1];


}

void read_dc_from_hdf5(string fname,double *dc,double *theta)
{

	herr_t status;	
	hid_t file,dtype,dspace_glbl,dspace,dset_glbl,plist_id,g_id;
	

	hsize_t dim[3],odim[2];

	dim[0] = N;
	dim[1] = N;
	dim[2] = N;
	odim[0] = N3;
	odim[1] = 1;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
       

	
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
	count[0] = N3; count[1] = odim[1]; offset[0] = 0; offset[1] = 0;

	

	dspace = H5Screate_simple(2, count, NULL);

	

	

	status = H5Dread(dset_glbl, dtype, dspace, dspace_glbl, plist_id, dc);

	H5Dclose(dset_glbl);
	H5Pclose(plist_id);
	H5Sclose(dspace);


	dset_glbl = H5Dopen(file, "/theta",H5P_DEFAULT);

	dspace_glbl = H5Dget_space(dset_glbl);


	count[0] = N3; count[1] = 1; 

	

	dspace = H5Screate_simple(2, count, NULL);

	

	

	

	status = H5Dread(dset_glbl, dtype, dspace, dspace_glbl, plist_id, theta);

	H5Dclose(dset_glbl);
	H5Pclose(plist_id);
	H5Sclose(dspace);







}


   

};







    
