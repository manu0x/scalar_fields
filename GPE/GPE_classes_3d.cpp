

#include <omp.h>
#include <fftw3.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>



class GPE_field_3d
{
    private:

    double hbar_unit,c_unit,h,pc_unit,omega_m0;
    double vfac;
    
    double m_alpha;

    public:
    double kppa;

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



    GPE_field_3d(int jj,int NN,int imex_s,int nthreads=4)
    {
        j= jj;
        N= NN;
        N2 = N*N;
        N3 = N*N*N;
        s = imex_s;
        dj = (double)j;
        dN3 = (double)N3;

        fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    

	    K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        V_phi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	    V_phi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

        fftw_plan_with_nthreads(nthreads);


	    plan_pois_f = fftw_plan_dft_3d(N,N,N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	    plan_pois_b = fftw_plan_dft_3d(N,N,N, fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

        plan_V_f = fftw_plan_dft_3d(N,N,N, V_phi, V_phi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	    plan_V_b = fftw_plan_dft_3d(N,N,N, V_phi_ft, V_phi, FFTW_BACKWARD, FFTW_ESTIMATE);


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




    void solve_V(double *kgrid)
    {
        double psisqr,ksqr;  int loc_i,loc_j,loc_k,ii;

        #pragma omp parallel for private(psisqr)
        for(loc_i=0;loc_i<N3;++loc_i)
        {
            psisqr = fpGpsi[loc_i][0]*fpGpsi[loc_i][0] + fpGpsi[loc_i][1]*fpGpsi[loc_i][1];
            V_phi[loc_i][0] = (0.5/kppa)*(psisqr-3.0*omega_m0);
            V_phi[loc_i][1] = 0.0;


        }
        fftw_execute(plan_V_f);

        for(ii=0,loc_i=-1,loc_j=-1,loc_k=0;ii<N3;++ii,++loc_k)
        {
            if(ii%N ==0)
            {
                loc_k=0;
                ++loc_j;
                if(ii%N2==0)
                {
                    loc_j=0;
                    ++loc_i;


                }


            }
            ksqr = kgrid[loc_i]*kgrid[loc_i] + kgrid[loc_j]*kgrid[loc_j] + kgrid[loc_k]*kgrid[loc_k];

          if(ksqr>0.0)
           { V_phi_ft[ii][0] = -V_phi_ft[ii][0]/(ksqr*dN3);
                
             V_phi_ft[ii][1] = -V_phi_ft[ii][1]/(ksqr*dN3);
           }
            else
          {
            V_phi_ft[ii][0] = 0.0;
            V_phi_ft[ii][1] = 0.0;

          }


        }

        fftw_execute(plan_V_b);


        


    }
   
    void update_fft_fields(int ind,double ksqr,double lamda)
    {
       
       
       	fpGpsi_ft[ind][0] = ((fpGpsi_ft[ind][0]) + lamda*fpGpsi_ft[ind][1])/(1.0+lamda*lamda);
		fpGpsi_ft[ind][1] = ((fpGpsi_ft[ind][1]) - lamda*fpGpsi_ft[ind][0])/(1.0+lamda*lamda);
        


		fpGpsi_ft[ind][0] = fpGpsi_ft[ind][0]/((double)N3);
		fpGpsi_ft[ind][1] = fpGpsi_ft[ind][1]/((double)N3);

       
		K_ft[ind][0] = -fpGpsi_ft[ind][0]*ksqr;
		K_ft[ind][1] = -fpGpsi_ft[ind][1]*ksqr;
       
       




    }

    void ex_rhs(int ind,int stg_i)
    {
            

            ex_K_psi[0][stg_i][ind] =  V_phi[ind][0]*fpGpsi[ind][1];
            ex_K_psi[1][stg_i][ind] =  -V_phi[ind][0]*fpGpsi[ind][0];
            

            



    }

    void im_rhs(int ind,int stg_i)
    {
        im_K_psi[0][stg_i][ind]  = -K[ind][1]*0.5*kppa;  im_K_psi[1][stg_i][ind]=  K[ind][0]*0.5*kppa;


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
		   


         }

    



        }
    }

    void print_params_set_kappa()
    {

        printf("\nPrinting params for field %d\n",j+1);
        printf("[c_unit,hbar_unit,pc_unit] = [%lf,%lf,%lf]\n",c_unit,hbar_unit,pc_unit);
        printf("[h,omega_m0,m_alpha] = [%lf,%lf,%lf]\n",h,omega_m0,m_alpha);
        kppa = c_unit*c_unit*hbar_unit*h*(1e-5)/(m_alpha*pc_unit);

        printf("kappa = %lf\n",kppa);
    }

    void initialise_random(double *kgrid,double amp=1e-6)
    {
        int ii,loc_i,loc_j,loc_k;
        double uni_rand,theta_rnd;
        double psiamp,ksqr,loc_delta;
        double vmax;

        double rmx = (double)RAND_MAX;
        for(ii=0;ii<N3;++ii)
        {

            uni_rand = amp*(((double)rand())/rmx);

            

            fpGpsi[ii][0] = uni_rand;
            fpGpsi[ii][1] = 0.0;


        }

        fftw_execute(plan_pois_f);
        fpGpsi_ft[0][0] = 0.0;   fpGpsi_ft[0][1] = 0.0;
        fftw_execute(plan_pois_b);
        double dn3 = (double)N3;
        
        for ( ii = 0,loc_i=-1,loc_j=-1; ii < N3; ++ii,++loc_k)
        {   
            if(ii%N ==0)
            {
                loc_k=0;
                ++loc_j;
                if(ii%N2==0)
                {
                    loc_j=0;
                    ++loc_i;


                }


            }
            ksqr = kgrid[loc_i]*kgrid[loc_i] + kgrid[loc_j]*kgrid[loc_j] + kgrid[loc_k]*kgrid[loc_k];

           loc_delta = fpGpsi[ii][0]/dn3;

           psiamp = sqrt(3.0*omega_m0*(1.0+loc_delta));
           theta_rnd = (((double)rand())/rmx)*2.0*M_PI;
           psi[ii][0] = psiamp*cos(theta_rnd);
           psi[ii][1] = psiamp*sin(theta_rnd);

           V_phi_ft[ii][0] = 1.5*omega_m0*fpGpsi_ft[ii][0]/(dn3*kppa*ksqr);
           V_phi_ft[ii][1] = 1.5*omega_m0*fpGpsi_ft[ii][1]/(dn3*kppa*ksqr);

           if(isnan(psi[ii][0]+psi[ii][1]))
           {
            printf("ini break %d %lf %lf\n",ii,psi[ii][0],psi[ii][1]);
            break;

           }
           
        }
        V_phi_ft[0][0] = 0.0; V_phi_ft[0][1] = 0.0;

        fftw_execute(plan_V_b);

        vmax = -10000000000.0;
        for ( ii = 0; ii < N3; ++ii,++loc_k)
        { 

            if(fabs(V_phi[ii][0])>vmax)
                vmax = V_phi[ii][0];

        }


        printf("Initial random done with vmax %.10lf\n",vmax);



    }
    void initialise_from_file(char *f_real,char *f_img)
    {
        FILE *fp_real  = fopen(f_real,"r");
        FILE *fp_img  = fopen(f_img,"r");

       // FILE *fpcheck = fopen("ini_check.txt","w");

        int ii,jj,ci;
       
        //printf("Working...\n");
       
        for(ii=0;ii<N+1;++ii)
        {   
            for(jj=0;jj<N+1;++jj)
            {
                
                ci = ii*N+jj;
               


                fscanf(fp_real,"%lf\t",&psi[ci][0]);
                fscanf(fp_img,"%lf\t",&psi[ci][1]);
               
              //  printf("%d %d \n",ii,jj);
             //   fprintf(fpcheck,"%lf ",psi[ci][0]);

            }
          //  fprintf(fpcheck,"\n");
         


        }

      

        fclose(fp_real);
        fclose(fp_img);
        //fclose(fpcheck);

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
        if(is_ini)
        {
            ini_energy=0.0;
            ini_mass = 0,0;

            max_eng_err=-1.0;
            max_mass_err=-1.0;

        }

    }
    void cal_conserve_at_point(int ci,double dx,int is_ini=0)
    {
            
        
            double loc_energy,loc_mass,psi2;//der_x[2],der_y[2],der_amp2,fcntr,Rcntr;
         /*     int left[2],right[2],left_y,right_y,ci,cr,cl;   
           
            /////////////////////  x-derivative //////////////////////
          left[0] = ind[0]-1;
            left[1] = ind[1];
            if(ind[0]==0)
            left[0] = N-2;

            right[0] = ind[0]+1;
            right[1] = ind[1];
            if(ind[0]==(N-1))
            right[0] = 1;

            cr = right[0]*N+right[1];
            cl = left[0]*N+left[1];

            der_x[0] = (psi[cr][0]+psi[cl][0]-2.0*psi[ci][0])/dx;
            der_x[1] = (psi[cr][1]+psi[cl][1]-2.0*psi[ci][1])/dx;

            ///////////////////////////////////

            /////////////////////  y-derivative //////////////////////
            left[0] = ind[0];
            left[1] = ind[1]-1;
            if(ind[1]==0)
            left[1] = N-2;

            right[0] = ind[0];
            right[1] = ind[1]+1;
            if(ind[1]==(N-1))
            right[1] = 1;

            cr = right[0]*N+right[1];
            cl = left[0]*N+left[1];

            der_y[0] = (psi[cr][0]+psi[cl][0]-2.0*psi[ci][0])/dx;
            der_y[1] = (psi[cr][1]+psi[cl][1]-2.0*psi[ci][1])/dx;

            ///////////////////////////////////  (der_x)^2+(der_y)^2    /////////////////
            der_amp2 =    der_x[0]*der_x[0] + der_x[1]*der_x[1] + der_y[0]*der_y[0] + der_y[1]*der_y[1];
         */   ///////////////////////////////////////////////////
            psi2 = psi[ci][0]*psi[ci][0] + psi[ci][1]*psi[ci][1];

         //   fcntr = 0.5*bta[0]*psi2*psi2 + 0.5*bta[1]*psi2*psi2_other;
           // Rcntr = -(psi[ci][0]*     )
      

           loc_energy = 0.0;// 0.5*der_amp2 +  V(x)*psi2 + fcntr + Rcntr;
            loc_mass = psi2;

            energy+=loc_energy;
            mass+=loc_mass;

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

        if(eng_err>max_eng_err)
            max_eng_err=eng_err;
        
        if(mass_err>max_mass_err)
            max_mass_err=mass_err;


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
