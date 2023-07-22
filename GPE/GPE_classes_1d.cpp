


#include <fftw3.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>

#define pie M_PI




class GPE_field_1d
{
    private:

   
    
    

    public:
    double bta,a,c,x0,T;
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
    



    int j,N,s;
    int i,k;

    double dj;
    double energy,ini_energy,max_eng_err;
    double mass,ini_mass,max_mass_err;
    double sol_n2,err_n2,max_sol_err;



    GPE_field_1d(int jj,int NN,int imex_s,int nthreads=4)
    {
        j= jj;
        N= NN;
        
        s = imex_s;
        dj = (double)j;

        fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	    fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N));

        psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	    

	    K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	    K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N));

      //  R = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	  //  R_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N));

      //  fftw_plan_with_nthreads(nthreads);


	    plan_pois_f = fftw_plan_dft_1d(N,fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
        plan_pois_b = fftw_plan_dft_1d(N,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);


	    plan_imp_b = fftw_plan_dft_1d(N,K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);
       //plan_imp_R = fftw_plan_dft_2d(N,N,R_ft, R, FFTW_BACKWARD, FFTW_ESTIMATE);

        for(i=0;i<2;++i)
	    {
		    im_K_psi[i] = new double*[imex_s];
		    ex_K_psi[i] = new double*[imex_s];

		    for(k=0;k<imex_s;++k)
		    {im_K_psi[i][k] =  new double[N];
		     ex_K_psi[i][k] =  new double[N];
		    }





	    }


    }



    double V(double x[2])
    {

        return(0.0);



    }


    double f(double psi2_me,double psi2_other)
    {
        return(0.0);

    }

    void update_fft_fields(int ind,double k,double lamda)
    {
       
     //  if(ind==130)
		//	printf("Lambda %.10lf  %.10lf  %.10lf  %.10lf  %.10lf\n",lamda,fpGpsi_ft[130][0],fpGpsi[130][0],fpGpsi_ft[130][1],fpGpsi[130][1]);
       	fpGpsi_ft[ind][0] = ((fpGpsi_ft[ind][0]) + lamda*fpGpsi_ft[ind][1])/(1.0+lamda*lamda); 
     //   if(ind==130)
		//	printf("Lambda %.10lf  %.10lf  %.10lf  %.10lf  %.10lf\n",lamda,fpGpsi_ft[130][0],fpGpsi[130][0],fpGpsi_ft[130][1],fpGpsi[130][1]);
		fpGpsi_ft[ind][1] = ((fpGpsi_ft[ind][1]) - lamda*fpGpsi_ft[ind][0])/(1.0+lamda*lamda);
        
        


		fpGpsi_ft[ind][0] = fpGpsi_ft[ind][0]/((double)N);
		fpGpsi_ft[ind][1] = fpGpsi_ft[ind][1]/((double)N);

       
		K_ft[ind][0] = -fpGpsi_ft[ind][0]*(k*k);
		K_ft[ind][1] = -fpGpsi_ft[ind][1]*(k*k);
       
      


    }


    void initialise(double dx,double xini,double ti=0.0)
    {


	    int i,j; double di,dk;
	    double theta,x,fs;
	    dk = 2.0*pie/(((double)(N))*dx);


	    x = xini;
	    for(i=0;i<N;++i,x+=dx)
	    {
	      di = (double)i;
	   
	     theta = 0.5*c*(x-x0)-(0.25*c*c -a)*ti;
	     fs = (2.0*a/bta)/cosh(sqrt(a)*(x-x0-c*ti));
	     
	      psi[i][0] = cos(theta)*fs;
	       psi[i][1] = sin(theta)*fs;
	     
        }


    }
    void ex_rhs(int ind,int stg_i)
    {
           

            ex_K_psi[0][stg_i][ind] =  bta*(fpGpsi[ind][0]*fpGpsi[ind][0] + fpGpsi[ind][1]*fpGpsi[ind][1])*fpGpsi[ind][1];
            ex_K_psi[1][stg_i][ind] =   -bta*(fpGpsi[ind][0]*fpGpsi[ind][0] + fpGpsi[ind][1]*fpGpsi[ind][1])*fpGpsi[ind][0];
            //if((stg_i==0)&&ind<100)
           // printf("fcheck %d  %d %lf %lf \n",j,stg_i,fpGpsi[ind][1],psi[ind][1]);

            



    }

    void im_rhs(int ind,int stg_i)
    {
        im_K_psi[0][stg_i][ind]  = -K[ind][1];  im_K_psi[1][stg_i][ind]=  K[ind][0];


    }

    void do_forward_fft()
    {
        fftw_execute(plan_pois_f);


    }

    void do_back_fft()
    {
      
       
        fftw_execute(plan_pois_b);
        
        fftw_execute(plan_imp_b);
        //fftw_execute(plan_imp_R);

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
		   

		   
            num_val = stod(read_val);
		  
		    if (read_key=="beta")
           {
               
               bta = num_val;
               a = bta*bta/16.0;
           }
           else if (read_key=="c")
           {
                c = num_val;
           }
           else if (read_key=="x0")
           {
                x0 = num_val;
           }
           
           
         
           
		   


         }

    



        }
    }

    void print_params()
    {

        printf("\nPrinting params for field %d\n",j+1);
       
        printf("beta = %lf\n",bta);
        printf("a = %lf\n",a);
        printf("c = %lf\n",c);
        printf("x0 = %lf\n",x0);
        
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
        for(ii=0;ii<N;++ii)
        {

            fpGpsi[ii][0] = psi[ii][0];
            fpGpsi[ii][1] = psi[ii][1];

        }

    }

    void reset_consv_quant(int is_ini=0)
    {
        energy = 0.0;
        sol_n2=0.0;
        err_n2=0.0;
        mass = 0.0;
        if(is_ini)
        {
            
            ini_mass = 0,0;

            max_sol_err=-1.0;
            max_mass_err=-1.0;

        }

    }
    void cal_err_at_point(int ind,double xini,double dx,double t,int is_ini=0)
    {
            
            int left,right,left_y,right_y,ci;
            double x,sol[2],theta,dbi,psi2,loc_err,loc_mass,fs;
            
            ci = ind;
       
            psi2 = psi[ci][0]*psi[ci][0] + psi[ci][1]*psi[ci][1];

            dbi = (double)(ind);
			x = xini+dbi*dx;
			theta = 0.5*c*(x-x0)-(0.25*c*c -a)*(t);
	   		fs = (2.0*a/bta)/cosh(sqrt(a)*(x-x0-c*(t)));
			sol[0] = cos(theta)*fs;
			sol[1] = sin(theta)*fs;

			(err_n2)+=((sol[0]-psi[ci][0])*(sol[0]-psi[ci][0])+(sol[1]-psi[ci][1])*(sol[1]-psi[ci][1]));
			sol_n2+= (sol[0]*sol[0] + sol[1]*sol[1]);


          
            

           
            loc_mass = psi2;

            
            mass+=loc_mass;

            if(is_ini)
            {
                
                ini_mass+=loc_mass;

            }

            
    }

    void conserve_err()
    {
        double err,mass_err;
       // eng_err = fabs(energy-ini_energy)/fabs(ini_energy);
        mass_err = fabs(mass-ini_mass)/mass;

        err_n2 = sqrt(err_n2);
		sol_n2 = sqrt(sol_n2);

		err = err_n2/sol_n2;
		if(err>(max_sol_err))
			max_sol_err = err;
        
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
