


#include <fftw3.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>



class GPE_field_2d
{
    private:

    double gmma[2];
    double bta[2];
    double kppa;

    public:
    fftw_complex *psi;

    fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_complex *K;
	fftw_complex *K_ft;

    fftw_complex *R;
	fftw_complex *R_ft;

    double ***im_K_psi = new double**[2];
	double ***ex_K_psi = new double**[2];

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;
    fftw_plan plan_imp_R;



    int j,N,N2,s;
    int i,k;

    double dj;
    double energy,ini_energy,max_eng_err;
    double mass,ini_mass,max_mass_err;

    FILE *fp;



    GPE_field_2d(int jj,int NN,int imex_s,int nthreads=4)
    {
        j= jj;
        N= NN;
        N2 = N*N;
        s = imex_s;
        dj = (double)j;

        fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N2));
	    fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N2));

        psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N2));
	    

	    K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N2));
	    K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N2));

        R = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N2));
	    R_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N2));

        fftw_plan_with_nthreads(nthreads);


	    plan_pois_f = fftw_plan_dft_2d(N,N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	    plan_pois_b = fftw_plan_dft_2d(N,N,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);


	    plan_imp_b = fftw_plan_dft_2d(N,N,K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);
        plan_imp_R = fftw_plan_dft_2d(N,N,R_ft, R, FFTW_BACKWARD, FFTW_ESTIMATE);

        for(i=0;i<2;++i)
	    {
		    im_K_psi[i] = new double*[imex_s];
		    ex_K_psi[i] = new double*[imex_s];

		    for(k=0;k<imex_s;++k)
		    {im_K_psi[i][k] =  new double[N2];
		     ex_K_psi[i][k] =  new double[N2];
		    }





	    }

        fp = fopen("my_err.txt","w");


    }



    double V(double x[2])
    {

        return(0.5*(gmma[0]*x[0]*gmma[0]*x[0]+gmma[1]*x[1]*gmma[1]*x[1]));



    }


    double f(double psi2_me,double psi2_other)
    {
        return(bta[0]*psi2_me+bta[1]*psi2_other);

    }

    void update_fft_fields(int ind,double k[2],double lamda)
    {
        double temp_psi[2];
        temp_psi[0] = fpGpsi_ft[ind][0];
        temp_psi[1] = fpGpsi_ft[ind][1];
       
       	fpGpsi_ft[ind][0] = ((fpGpsi_ft[ind][0]) + lamda*fpGpsi_ft[ind][1])/(1.0+lamda*lamda);
		fpGpsi_ft[ind][1] = ((fpGpsi_ft[ind][1]) - lamda*temp_psi[0])/(1.0+lamda*lamda);
        


		fpGpsi_ft[ind][0] = fpGpsi_ft[ind][0]/((double)N2);
		fpGpsi_ft[ind][1] = fpGpsi_ft[ind][1]/((double)N2);

       
		K_ft[ind][0] = -fpGpsi_ft[ind][0]*(k[0]*k[0]+k[1]*k[1]);
		K_ft[ind][1] = -fpGpsi_ft[ind][1]*(k[0]*k[0]+k[1]*k[1]);
       
        R_ft[ind][0] = kppa*(k[0]*fpGpsi_ft[ind][1] - pow(-1.0,dj+1.0)*k[0]*fpGpsi_ft[ind][0]);
        R_ft[ind][1] = kppa*(-k[0]*fpGpsi_ft[ind][0]- pow(-1.0,dj+1.0)*k[1]*fpGpsi_ft[ind][1]);


    }

    void ex_rhs(int ind,int stg_i,double psi2_me,double psi2_other,double R[2],double x[2])
    {
            

            ex_K_psi[0][stg_i][ind] = R[0] + (f(psi2_me,psi2_other) + V(x))*fpGpsi[ind][1];
            ex_K_psi[1][stg_i][ind] = R[1]   -(f(psi2_me,psi2_other) + V(x))*fpGpsi[ind][0];
            //if((stg_i==0)&&ind<100)
           // printf("fcheck %d  %d %lf %lf \n",j,stg_i,fpGpsi[ind][1],psi[ind][1]);

            



    }

    void im_rhs(int ind,int stg_i)
    {
        im_K_psi[0][stg_i][ind]  = -K[ind][1]*0.5;  im_K_psi[1][stg_i][ind]=  K[ind][0]*0.5;


    }

    void do_forward_fft()
    {
        fftw_execute(plan_pois_f);


    }

    void do_back_fft()
    {
        fftw_execute(plan_pois_b);
		fftw_execute(plan_imp_b);
        fftw_execute(plan_imp_R);

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
		   

		   
		
		   if(read_key=="gamma_x")
			{
                num_val = stod(read_val);
                gmma[0] = num_val;

            }
		   else if (read_key=="gamma_y")
           {
               num_val = stod(read_val);
               gmma[1] = num_val;
           }
           else if (read_key=="beta_my")
           {
               num_val = stod(read_val);
               bta[0] = num_val;
           }
           else if (read_key=="beta_other")
           {
               num_val = stod(read_val);
               bta[1] = num_val;
           }
           else if (read_key=="kappa")
           {
               num_val = stod(read_val);
               kppa = num_val;
           }
           
		   


         }

    



        }
    }

    void print_params()
    {

        printf("\nPrinting params for field %d\n",j+1);
        printf("gamma = [%lf,%lf]\n",gmma[0],gmma[1]);
        printf("beta = [%lf,%lf]\n",bta[0],bta[1]);
        printf("kappa = %lf\n",kppa);
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
        for(ii=0;ii<N2;++ii)
        {

            fpGpsi[ii][0] = psi[ii][0];
            fpGpsi[ii][1] = psi[ii][1];

        }

        fftw_execute(plan_pois_f);

        for(ii=0;ii<N2;++ii)
        {

            fpGpsi_ft[ii][0] = fpGpsi_ft[ii][0]/((double)N2);
            fpGpsi_ft[ii][1] = fpGpsi_ft[ii][1]/((double)N2);

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
    void cal_conserve_at_point(int ind[2],double psi2_other,double R[2],double x[2],double dx,int is_ini=0)
    {
            
            int left[2],right[2],left_y,right_y,ci,cr,cl;
            double der_x[2],der_y[2],der_amp2,psi2,fcntr,Rcntr,loc_energy,loc_mass;
            
            ci = ind[0]*N+ind[1];
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
            ///////////////////////////////////////////////////
            psi2 = psi[ci][0]*psi[ci][0] + psi[ci][1]*psi[ci][1];

            fcntr = 0.5*bta[0]*psi2*psi2 + 0.5*bta[1]*psi2*psi2_other;
           // Rcntr = -(psi[ci][0]*     )
      

            loc_energy = 0.5*der_amp2 +  V(x)*psi2 + fcntr + Rcntr;
            loc_mass = psi2;

            energy+=loc_energy;
            mass+=loc_mass;

            if(is_ini)
            {
                ini_energy+=loc_energy;
                ini_mass+=loc_mass;

            }

            
    }

    void conserve_err(double t)
    {
        double eng_err,mass_err;
        eng_err = fabs(energy-ini_energy)/fabs(ini_energy);
        mass_err = fabs(mass-ini_mass)/ini_mass;

        if(eng_err>max_eng_err)
            max_eng_err=eng_err;
        
        if(mass_err>max_mass_err)
            max_mass_err=mass_err;


        fprintf(fp,"%lf\t%lf\t%lf\n",t,max_mass_err,max_eng_err);


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
