


#include <fftw3.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>
#include "../imex/imex_classes.cpp"


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



    GPE_field_2d(int jj,int NN,int imex_s,int nthreads=8)
    {
        j= jj;
        N= NN;
        N2 = N*N;
        s = imex_s;
        dj = (double)j;

        fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N2));
	    fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N2));

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


    }



    double V(double x[2])
    {

        return(gmma[0]*x[0]*gmma[0]*x[0]+gmma[1]*x[1]*gmma[1]*x[1]);



    }


    double f(double psi2_me,double psi2_other)
    {
        return(bta[0]*psi2_me+bta[1]*psi2_other);

    }

    void update_fft_fields(int ind,double k[2],double lamda)
    {
       	fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lamda*fpGpsi_ft[i][1])/(1.0+lamda*lamda);
		fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lamda*fpGpsi_ft[i][0])/(1.0+lamda*lamda);



		fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N);
		fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N);

		K_ft[i][0] = -fpGpsi_ft[i][0]*(k[0]*k[0]+k[1]*k[1]);
		K_ft[i][1] = -fpGpsi_ft[i][1]*(k[0]*k[0]+k[1]*k[1]);
       
        R_ft[ind][0] = kppa*(k[0]*fpGpsi_ft[ind][1] + pow(-1.0,dj)*k[0]*fpGpsi_ft[ind][0]);
        R_ft[ind][1] = kppa*(-k[0]*fpGpsi_ft[ind][0]+ pow(-1.0,dj)*k[1]*fpGpsi_ft[ind][1]);


    }

    void ex_rhs(int ind,int stg_i,double psi2_me,double psi2_other,double R[2],double x[2])
    {
            

            ex_K_psi[0][stg_i][ind] = R[0]  + (f(psi2_me,psi2_other) + V(x))*fpGpsi[ind][1];
            ex_K_psi[1][stg_i][ind] = R[1]  + (-f(psi2_me,psi2_other) - V(x))*fpGpsi[ind][0];

            



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

    void initialise()
    {}

    void set_field()
    {}

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
