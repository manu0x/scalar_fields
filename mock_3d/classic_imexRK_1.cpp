

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <mpi.h>
#include <limits>

using namespace std;

#define pie M_PI



//////////GLobal constants/////////////

double m,n,T;

///////////////////////////////////////


//////////// ImEx RK Butcher Tableau /////
const int imex_s = 3;
const int im_s = imex_s;
const int ex_s = imex_s;

//	/	/	/	/ SCHEME 17	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  205.0/462.0,2.0/11.0,0.0,  2033.0/4620.0,21.0/110.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,289.0/462.0,751.0/924.0};
double im_b[im_s] = {24.0/55.0,1.0/5.0,4.0/11.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  5.0/6.0,0.0,0.0,  11.0/24.0,11.0/24.0,0.0};
double ex_c[ex_s] = {0.0,5.0/6.0,11.0/12.0};
double ex_b[ex_s] = {24.0/55.0,1.0/5.0,4.0/11.0};

/*
//	/	/	/	/ SCHEME 20	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  41.0/154.0,2.0/11.0,0.0,  289.0/847.0,42.0/121.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,69.0/154.0,67.0/77.0};
double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};
double ex_c[ex_s] = {0.0,0.5,1.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};


////////////////////////////////////////////


//	/	/	/	/ SCHEME 22	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  2829.0/9317.0,2.0/11.0,0.0,  148529.0/428582.0,7.0/23.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,4523.0/9317.0,15517.0/18634.0};
double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};
double ex_c[ex_s] = {0.0,0.5,1.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};


////////////////////////////////////////////

//	/	/	/	/ SCHEME 23	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  2583.0/13310.0,2.0/11.0,0.0,  39731.0/139755.0,10.0/21.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,5003.0/13310.0,6271.0/6655.0};
double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};
double ex_c[ex_s] = {0.0,0.5,1.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

*/
////////////////////////////////////////////

double V(double psi_amp)
{

	return( 1.0 - 2.0*pie*pie*n*n/(m*m));
	


}



void vel_all(double v[2],double psi[2],double Vval,double lap_psi[2])
{

	v[0] = 2.0*pie*(Vval*psi[1] - lap_psi[1]/(2.0*m*m))/T;
	v[1] = 2.0*pie*(-Vval*psi[0] + lap_psi[0]/(2.0*m*m))/T;



}

void vel(double v[2],double psi[2],double Vval,double lap_psi[2])
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



}

void ex_vel(double v[2],double psi[2],double Vval)
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



}


void lap(double *lpsi,fftw_complex *psi, double dx,int N)
{
	int l1,l2,r1,r2,i;

	double vl1,vl2,vr1,vr2,vc;
   for(i=0;i<(N-1);++i)
   {	l1 = ((N-1)+ (i-1))%(N-1);
	l2 = ((N-1)+ (i-2))%(N-1);

	r1 = (i+1)%(N-1);
	r2 = (i+2)%(N-1);


	vl1 = psi[l1][0];
	vl2 = psi[l2][0];

	vr1 = psi[r1][0];
	vr2 = psi[r2][0];

	vc = psi[i][0];

	*(lpsi+2*i) =  (-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2)/(12.0*dx*dx);

	//if(i==0||i==1||i==2||i==(N-1))
	//{printf("\n%d\t%d\t%d\t%d\t%d\n",l2,l1,i,r1,r2);
	// printf("%lf\t%lf\t%lf\t%lf\t%lf\n",vl2,vl1,vc,vr1,vr2);
	// printf("%lf\t%lf\n\n",(-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2),*(lpsi+2*i));
	//}


	vl1 = psi[l1][1];
	vl2 = psi[l2][1];

	vr1 = psi[r1][1];
	vr2 = psi[r2][1];

	vc = psi[i][1];
	*(lpsi+2*i+1) = (-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2)/(12.0*dx*dx);

    }

	*(lpsi+2*i) = *(lpsi);
	*(lpsi+2*i+1) = *(lpsi+1);

}


void initialise(fftw_complex *psi,double *k,double dx,int N)
{

	int i,j,l,ci; double di,dj,dl,dci,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	
	for(i=0;i<N;++i)
	{
	   di = (double)i;

	   if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk;


	}

	for(i=0;i<N;++i)
	{	di = (double)i;
		for(j=0;j<N;++j)
		{	dj = (double)j;
			for(l=0;l<N;++l)
			{	dl = (double)l;
				ci = i*N*N+j*N+l;
				dci = (double)ci;
				psi[ci][0] = sin(2.0*pie*n*di*dx+2.0*pie*n*dj*dx+2.0*pie*n*dl*dx);
	    		psi[ci][1] = 0.0;
	 

			}


		}

	}




}




double run(double dt,double dx,double *abs_err,double *stb_avg,int stb_any,int printfp,int prt)
{

	int N,t_steps;
	double box_len,t_end,t_start,xval,dbi,dbj,dbk,sol[2],stb_ini;

	
///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data_ft.txt","w");
	FILE *fp2 = fopen("data2_ft.txt","w");
	FILE *fplap = fopen("lap_ft.txt","w");
	FILE *fpmass = fopen("mass_ft.txt","w");
	FILE *fptime = fopen("tm_ft.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 0.1;

	T = 2.0*pie/m;


   




/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 2.0;
	n  = 2.0/box_len;
	N = 128;
	dx = box_len/(double(N));
	//N = ((int)(box_len/dx));
	int N3 = N*N*N;
	int N2 = N*N;
	printf("N2 %d N3 %d\n",N2,N3);
////////////////////////////// Time & dt settings ////////////////////////////////////////////////

	t_start = 0.0;
	t_end = 2.0;
	if(stb_any)
	t_end = 0.01*t_end;
	t_steps = (int)((t_end-t_start)/dt);
	//dt  = (t_end-t_start)/((double)t_steps);
	if(prt)
	printf("dt %lf N %d\n",dt,N);


	printf("\nHHHH\n");
/////////////////////////////////////////RK things/////////////////////////////////////////
	int i,j;

	//double im_K_psi[2][ex_s][N*N*N],ex_K_psi[2][im_s][N*N*N];

	//double *exKswap = new double [2*ex_s*N3];
    //double *imKswap = new double [2*im_s*N3];

	double ***im_K_psi = new double**[2];
	double ***ex_K_psi = new double**[2];
	for(i=0;i<2;++i)
	{
		im_K_psi[i] = new double*[ex_s];
		ex_K_psi[i] = new double*[ex_s];

		for(j=0;j<ex_s;++j)
		{im_K_psi[i][j] =  new double[N3];
		 ex_K_psi[i][j] =  new double[N3];
		}


	}

////////////////////////////// Psi variables  /////////////////////////////////////////////////////

	printf("Ini start\n");


	double lambda;

	fftw_complex *psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	


	double k_grid[N];

	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_complex *K;
	fftw_complex *K_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;

	fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N3));

	K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));
	K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N3));


	plan_pois_f = fftw_plan_dft_3d(N,N,N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_pois_b = fftw_plan_dft_3d(N,N,N,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);


	plan_imp_b = fftw_plan_dft_3d(N,N,N,K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);

	initialise(psi,k_grid,dx,N);

	

	 	for(i=0;i<N3;++i)
		{

			fpGpsi[i][0] = psi[i][0] ;
			fpGpsi[i][1] = psi[i][1] ;

			





		}


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
    int ii,jj,kk;
	int s_cntr,tcntr,printcntr,fail=0;
	printcntr = (int)(((double) t_steps)/1000.0);
	if(t_steps<=100)
		printcntr = 1.0;	//printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp,Vval,amp,avg_amp;
	double drc = 2.0*pie*n*2.0*pie*n;
	double fdt,amp_ini;
	*stb_avg=0.0;


	ii=-1;
	jj=-1;
	kk=0;


	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&( 1);t+=dt,++tcntr)
	{

		avg_amp = 0.0;

		  fftw_execute(plan_pois_f);

	 	for(i=0,ii=-1,jj=-1,kk=0;i<(N3);++i,++kk)
		{	

			if((i%N)==0)
			{   kk=0;
				++jj;
				if(i%(N2)==0)
				{
					jj = 0;
					++ii;
				}
			}
		//	printf("i,j,k  %d %d %d\n",ii,jj,kk);
			if((tcntr%printcntr)==0||(stb_any))
			{
			  amp  = sqrt(psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);
			   avg_amp+=amp;
			   dbi = (double)(ii); dbj = (double)(jj); dbk = (double)(kk);
			   sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi+2.0*pie*n*dx*dbj+2.0*pie*n*dx*dbk);
			   sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi+2.0*pie*n*dx*dbj+2.0*pie*n*dx*dbk);
			 // if(printfp)
			 // fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*((double)i),fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			  

			  if((i==(10))&&(printfp))
			   fprintf(fptime,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}

			lambda = (k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj]+k_grid[kk]*k_grid[kk])*im_a[0][0]*dt/(2.0*m*3.0);


			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

			fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N3);
			fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N3);

			K_ft[i][0] = -fpGpsi_ft[i][0]*((k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj]+k_grid[kk]*k_grid[kk]));
			K_ft[i][1] = -fpGpsi_ft[i][1]*((k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj]+k_grid[kk]*k_grid[kk]));

			if(i==0)
			{
			 K_ft[i][0] = 0.0;
			 K_ft[i][1] = 0.0;
			}

		}
		stb_ini = avg_amp;
		if(t==t_start)
		amp_ini = avg_amp;

		fftw_execute(plan_pois_b);
		fftw_execute(plan_imp_b);

		/*
		for(i=0;i<(N);++i)
		{	dbi = (double)i;
			sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);


			fdt = (fpGpsi[i-1][0]+fpGpsi[i+1][0]-2.0*fpGpsi[i][0])/(dx*dx);
			if(i==0)
			fdt = (fpGpsi[N-1][0]+fpGpsi[i+1][0]-2.0*fpGpsi[i][0])/(dx*dx);
			if(i==(N-1))
			fdt = (fpGpsi[i-1][0]+fpGpsi[0][0]-2.0*fpGpsi[i][0])/(dx*dx);
			fprintf(fplap,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,K[i][0],K[i][1],-drc*fpGpsi[i][0],-drc*fpGpsi[i][1],fdt);

		}

		//printf("l do\n");
		*/

		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N3;++i)
			   {


				if(j==0)
				{
				/*	if(s_cntr==1)
					{
					  c_psi[0] = psi[i][0];
		    			  c_psi[1] =  psi[i][1];


					}
					else
					*/

					c_psi[0] = fpGpsi[i][0];
		    		c_psi[1] =  fpGpsi[i][1];




		    			c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    			Vval = V(c_psi_amp);
		    			ex_vel(vel_val,c_psi,Vval);

		    			ex_K_psi[0][s_cntr-1][i] = vel_val[0];  ex_K_psi[1][s_cntr-1][i] = vel_val[1];
					im_K_psi[0][s_cntr-1][i] = -K[i][1]/(2.0*m*3.0);  im_K_psi[1][s_cntr-1][i] =  K[i][0]/(2.0*m*3.0);

					fpGpsi[i][0] = psi[i][0] + dt*ex_a[s_cntr][j]*ex_K_psi[0][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[0][j][i];
					fpGpsi[i][1] = psi[i][1] + dt*ex_a[s_cntr][j]*ex_K_psi[1][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[1][j][i];

				}


				else
				{fpGpsi[i][0]+=  dt*ex_a[s_cntr][j]*ex_K_psi[0][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[0][j][i];
				 fpGpsi[i][1]+=  dt*ex_a[s_cntr][j]*ex_K_psi[1][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[1][j][i];
				}



			   }


			}

		

			fftw_execute(plan_pois_f);

			for(i=0,ii=-1,jj=-1,kk=0;i<N3;++i,++kk)
			{
			
				if((i%N)==0)
			    {   kk=0;
					++jj;
					if(i%(N2)==0)
					{
						jj = 0;
						++ii;
					}
				}
				lambda = (k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj]+k_grid[kk]*k_grid[kk])*im_a[s_cntr][s_cntr]*dt/(2.0*m*3.0);
				//if(lambda>1.0)
				//printf("ii %d jj %d kk %d lambda %lf\n",ii,jj,kk,lambda);

				fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
				fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);



				fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N3);
				fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N3);

				K_ft[i][0] = -fpGpsi_ft[i][0]*((k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj]+k_grid[kk]*k_grid[kk]));
				K_ft[i][1] = -fpGpsi_ft[i][1]*((k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj]+k_grid[kk]*k_grid[kk]));

				if(i==0)
				{
				 K_ft[i][0] = 0.0;
				 K_ft[i][1] = 0.0;
				}



			}

			fftw_execute(plan_pois_b);
			fftw_execute(plan_imp_b);







		}//RK stages concluded

		if((tcntr%printcntr)==0)
		{
			  if(prt)
			   printf("%lf\t%lf\n",t/t_end,avg_amp/((double)(N3)));
		}

		avg_amp=0.0;
		*abs_err = 0.0;


		for(i=0,ii=-1,jj=-1,kk=0;i<N3;++i,++kk)
		{

			if((i%N)==0)
			{   kk=0;
				++jj;
				if(i%(N2)==0)
				{
					jj = 0;
					++ii;
				}
			}

			c_psi[0] = fpGpsi[i][0];
		    c_psi[1] = fpGpsi[i][1];

			c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    Vval = V(c_psi_amp);
		    ex_vel(vel_val,c_psi,Vval);

		    ex_K_psi[0][imex_s-1][i] = vel_val[0];  ex_K_psi[1][imex_s-1][i] = vel_val[1];
			im_K_psi[0][imex_s-1][i] = -K[i][1]/(2.0*m*3.0);  im_K_psi[1][imex_s-1][i] =  K[i][0]/(2.0*m*3.0);


			for(j=0;j<imex_s;++j)
			{	psi[i][0]+=  dt*ex_b[j]*ex_K_psi[0][j][i]+ dt*im_b[j]*im_K_psi[0][j][i];
				psi[i][1]+=  dt*ex_b[j]*ex_K_psi[1][j][i]+ dt*im_b[j]*im_K_psi[1][j][i];




			}

			fpGpsi[i][0] = psi[i][0];
			fpGpsi[i][1] = psi[i][1];


			if(isnan((psi[i][0])+psi[i][1]))
			{
				fail =1;
				if(prt)
				printf("FAILED %d tcntr %d %lf %lf\n",i,tcntr,psi[i][0],psi[i][0]);

				break;

			}


			dbi = (double)(ii); dbj = (double)(jj); dbk = (double)(kk);
			
			sol[0] = cos(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi+2.0*pie*n*dx*dbj+2.0*pie*n*dx*dbk);
			sol[1] = -sin(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi+2.0*pie*n*dx*dbj+2.0*pie*n*dx*dbk);

			(*abs_err)+=(fabs(sol[0]-psi[i][0])+fabs(sol[1]-psi[i][1]));






			amp  = sqrt(psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);
			avg_amp+=amp;

			if((tcntr%printcntr)==0)
			{
			 // if(printfp)
			 // fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}


		   }

			*abs_err = (*abs_err)/((double)N3);
			(*stb_avg)+=(avg_amp/stb_ini);






		if((stb_any)&&(fail))
		{

			*abs_err = 1e5;
			*stb_avg = std::numeric_limits<double>::infinity();

			return(1e5);

		}
		else
		if(((100.0*fabs(avg_amp-amp_ini)/amp_ini)>=1e3)||(fail)||((*abs_err)>=1e3))
		{

			if(*abs_err>=1e3)
			printf("abs err blow %lf %d %d %d\n",*abs_err,ii,jj,kk);

			if((100.0*fabs(avg_amp-amp_ini)/amp_ini)>=1e3)
			printf("amp blow  %lf  %lf  %d %d %d\n",avg_amp,amp_ini,ii,jj,kk);

			*abs_err = -1e3;

			return(-1e3);

		}



		if((tcntr%printcntr)==0)
		{
		  if(printfp)
		  {fprintf(fp,"\n\n\n");
		   fprintf(fp2,"\n\n\n");
		   fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)(N3)));
		  }
		  if(prt)
		  printf("%lf\t%.10lf\n\n",t/t_end,avg_amp/((double)(N3)));
		}





	}///// ENd f Time Evolution /////////////////////////
	*stb_avg = (*stb_avg)/((double) tcntr);
	avg_amp = 0.0;
	for(i=0;i<N3;++i)
	{	amp  = sqrt(fpGpsi[i][0]*fpGpsi[i][0] + fpGpsi[i][1]*fpGpsi[i][1]);
		avg_amp+=amp;

	}
	if(printfp)
	fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)N3));


	if(prt)
	printf("N %d\n Run en los %lf abs err %lf\n",N,100.0*fabs(avg_amp-amp_ini)/amp_ini,*abs_err);



	return(100.0*fabs(avg_amp-amp_ini)/amp_ini);









}






int main()
{

	double dt = 3e-4;
	double dx = 4e-3;
	double abs_err,en_loss,stb_avg;

	double dx_l=2e-3,dx_u = 4e-2;
	double dt_l= 1e-5,dt_u = 1e-2;

	double ddx = (dx_u-dx_l)/(20.0);
	double ddt = (dt_u-dt_l)/(20.0);

	FILE *fp = fopen("imex_ft.txt","w");

	//for(dt=dt_l;dt<=dt_u;dt+=ddt)
	{


		//for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			dx = 2e-2;
			dt = 1e-3;
			en_loss = run(dt,dx,&abs_err,&stb_avg,0,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
		}

	}





}
