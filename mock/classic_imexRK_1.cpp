using namespace std;

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <limits>

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


void initialise(fftw_complex *psi,double *k,double dx,int N,int spd)
{

	int i,j; double di,dk;
	dk = 2.0*pie/(((double)(N+2*spd))*dx);
	int nN = N+2*spd;
	for(i=0;i<nN;++i)
	{
	   di = (double)i;
	   if(i<N)	
	  { psi[i][0] = sin(2.0*pie*n*di*dx);
	    psi[i][1] = 0.0;
	 }
	   if(i<=(nN/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) nN))*dk; 

	/*    if(i<=(N/2))
		*(k+i+spd) = di*dk;
	    else
		*(k+i+spd) = (di-((double) N))*dk; 

	//	printf("k %lf  %lf  dx %lf\n",dk,*(k+i),dx );

		if((i<spd+1)&&(i>0))
			{
			  *(k+N+spd-1+i) = *(k+i+spd) ;
			 
				

			}
		if(i>(N-spd-1))
			{
			  *(k+i-N+spd) = *(k+i+spd-1) ;
			  
				

			}
	*/

	}

	


}




double run(double dt,double dx,double *abs_err,double *stb_avg,int stb_any,int printfp,int prt)
{

	int N,t_steps;
	double box_len,t_end,t_start,xval,dbi,sol[2],stb_ini;
	

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
	//dx = box_len/(double(N-1));
	N = ((int)(box_len/dx));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 2.0;
	if(stb_any)
	t_end = 0.01*t_end;
	t_steps = (int)((t_end-t_start)/dt);
	//dt  = (t_end-t_start)/((double)t_steps);
	if(prt)
	printf("dt %lf N %d\n",dt,N);
	


/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_psi[2][ex_s][N],ex_K_psi[2][im_s][N];


////////////////////////////// Psi variables  /////////////////////////////////////////////////////




	double psi[N][2],lap_val[2*N],lambda;

	int spd =0;//N/4;
	int nN = N+ 2*spd;

	
	double k_grid[N+2*spd];
	
	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_complex *K;
	fftw_complex *K_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;

	fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+2*spd));
	fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N+2*spd));

	K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+2*spd));
	K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N+2*spd));


	plan_pois_f = fftw_plan_dft_1d(nN, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_pois_b = fftw_plan_dft_1d(nN,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

	
	plan_imp_b = fftw_plan_dft_1d(nN,K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	initialise(psi,k_grid,dx,N,spd);
	
	int i,j;
	
	 	for(i=0;i<N;++i)
		{	
			
			fpGpsi[i+spd][0] = psi[i][0] ;
			fpGpsi[i+spd][1] = psi[i][1] ;

			if((i<spd+1)&&(i>0))
			{
			  fpGpsi[N+spd-1+i][0] = psi[i][0] ;
			  fpGpsi[N+spd-1+i][1] = psi[i][1] ;	
				

			}
			if(i>(N-spd-1))
			{
			  fpGpsi[i-N+spd][0] = psi[i-1][0] ;
			  fpGpsi[i-N+spd][1] = psi[i-1][1] ;	
				

			}

					
			
		

		}


///////////////////////  Evolution ///////////////////////////////////////////////////////////////

	int s_cntr,tcntr,printcntr,fail=0;
	printcntr = (int)(((double) t_steps)/100.0);
	if(t_steps<=100)
		printcntr = 1.0;	//printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp,Vval,amp,avg_amp;
	double drc = 2.0*pie*n*2.0*pie*n;
	double fdt,amp_ini;
	*stb_avg=0.0;


			
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(tcntr<1);t+=dt,++tcntr)
	{	
		
		avg_amp = 0.0;

		  fftw_execute(plan_pois_f);

	 	for(i=0;i<(N+2*spd);++i)
		{	//dbi = (double)(i-spd);

			
			
			if((tcntr%printcntr)==0||(stb_any))  
			{
			  if(i<N)
			  {amp  = sqrt(psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);				
			   avg_amp+=amp;
			   dbi = (double)(i);
			   sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			   sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  if(printfp)
			  fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*((double)i),fpGpsi[i+spd][0],fpGpsi[i+spd][1],sol[0],sol[1],amp);
			  }
			  
			  if((i==(spd+10))&&(printfp))
			   fprintf(fptime,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}

			lambda = k_grid[i]*k_grid[i]*im_a[0][0]*dt/(2.0*m);


			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

			fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)nN);
			fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)nN);

			K_ft[i][0] = -fpGpsi_ft[i][0]*(k_grid[i]*k_grid[i]);
			K_ft[i][1] = -fpGpsi_ft[i][1]*(k_grid[i]*k_grid[i]);

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

		for(i=spd;i<(N+spd);++i)
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
	
		printf("l do\n");

		
		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N;++i)
			   {	
				
				
				if(j==0)
				{	
					if(s_cntr==1)
					{
					  c_psi[0] = psi[i][0];
		    			  c_psi[1] =  psi[i][1];
		

					}
					else
					{
					  c_psi[0] = fpGpsi[i+spd][0];
		    			  c_psi[1] =  fpGpsi[i+spd][1];
		

					}
					
		    			c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    			Vval = V(c_psi_amp);
		    			ex_vel(vel_val,c_psi,Vval);
					
		    			ex_K_psi[0][s_cntr-1][i] = vel_val[0];  ex_K_psi[1][s_cntr-1][i] = vel_val[1];
					im_K_psi[0][s_cntr-1][i] = -K[i+spd][1]/(2.0*m);  im_K_psi[1][s_cntr-1][i] =  K[i+spd][0]/(2.0*m);
					
					fpGpsi[i+spd][0] = psi[i][0] + dt*ex_a[s_cntr][j]*ex_K_psi[0][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[0][j][i];
					fpGpsi[i+spd][1] = psi[i][1] + dt*ex_a[s_cntr][j]*ex_K_psi[1][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[1][j][i];
					
				}
			
				
				else
				{fpGpsi[i+spd][0]+=  dt*ex_a[s_cntr][j]*ex_K_psi[0][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[0][j][i];
				 fpGpsi[i+spd][1]+=  dt*ex_a[s_cntr][j]*ex_K_psi[1][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[1][j][i];
				}
				
	
	
			   }

	
			}
			
			for(i=0;i<spd;++i)
			{
				fpGpsi[i][0] = fpGpsi[N-1+i][0];
				fpGpsi[N+spd+i][0] = fpGpsi[i+1+spd][0];
		
				fpGpsi[i][1] = fpGpsi[N-1+i][1];
				fpGpsi[N+spd+i][1] = fpGpsi[i+1+spd][1];
			}
	
			   fftw_execute(plan_pois_f);

			  for(i=0;i<N+2*spd;++i)
			  {	
				lambda = k_grid[i]*k_grid[i]*im_a[s_cntr][s_cntr]*dt/(2.0*m);
				
				fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
				fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

				

				fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)nN);
				fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)nN);

				K_ft[i][0] = -fpGpsi_ft[i][0]*(k_grid[i]*k_grid[i]);
				K_ft[i][1] = -fpGpsi_ft[i][1]*(k_grid[i]*k_grid[i]);

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
			   printf("%lf\t%lf\n",t/t_end,avg_amp/((double)(N)));
		}
			
		avg_amp=0.0;
		*abs_err = 0.0;
		

		for(i=0;i<N;++i)
		{	
	
			c_psi[0] = fpGpsi[i+spd][0];
		    	c_psi[1] =  fpGpsi[i+spd][1]; 
					
			c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    	Vval = V(c_psi_amp);
		    	ex_vel(vel_val,c_psi,Vval);

		    	ex_K_psi[0][imex_s-1][i] = vel_val[0];  ex_K_psi[1][imex_s-1][i] = vel_val[1];
			im_K_psi[0][imex_s-1][i] = -K[i+spd][1]/(2.0*m);  im_K_psi[1][imex_s-1][i] =  K[i+spd][0]/(2.0*m);
	
							
			   for(j=0;j<imex_s;++j)
			   {	psi[i][0]+=  dt*ex_b[j]*ex_K_psi[0][j][i]+ dt*im_b[j]*im_K_psi[0][j][i];
				psi[i][1]+=  dt*ex_b[j]*ex_K_psi[1][j][i]+ dt*im_b[j]*im_K_psi[1][j][i];
				
				

				
			    }


			fpGpsi[i+spd][0] = psi[i][0];
			fpGpsi[i+spd][1] = psi[i][1];


			if(isnan((psi[i][0])+psi[i][1]))
			{
				fail =1;
				if(prt)
				printf("FAILED %d tcntr %d %lf %lf\n",i,tcntr,psi[i][0],psi[i][0]);
			
				break;

			}

			
			dbi = (double)(i);
			sol[0] = cos(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);
			sol[1] = -sin(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);

			(*abs_err)+=(fabs(sol[0]-psi[i][0])+fabs(sol[1]-psi[i][1]));
			


			


			amp  = sqrt(psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);
			avg_amp+=amp;

			if((tcntr%printcntr)==0)  
			{ 
			  if(printfp)
			  fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}
			

		   }

			*abs_err = (*abs_err)/((double)N);
			(*stb_avg)+=(avg_amp/stb_ini);


		 for(i=0;i<spd;++i)
		 {
				fpGpsi[i][0] = fpGpsi[N-1+i][0];
				fpGpsi[N+spd+i][0] = fpGpsi[i+1+spd][0];
		
				fpGpsi[i][1] = fpGpsi[N-1+i][1];
				fpGpsi[N+spd+i][1] = fpGpsi[i+1+spd][1];
		 }

		

		if((stb_any)&&(fail))
		{

			*abs_err = 1e5;
			*stb_avg = std::numeric_limits<double>::infinity();		
			return(1e5);
			
		}
		else 
		if(((100.0*fabs(avg_amp-amp_ini)/amp_ini)>=1e3)||(fail)||((*abs_err)>=1e3))
		{

			*abs_err = -1e3;
					
			return(-1e3);
			
		}
	


		if((tcntr%printcntr)==0)  
		{ 
		  if(printfp)
		  {fprintf(fp,"\n\n\n");
		   fprintf(fp2,"\n\n\n");
		   fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)(N)));
		  }
		  if(prt)
		  printf("%lf\t%lf\n\n",t/t_end,avg_amp/((double)(N)));
		}
		


		  

	}///// ENd f Time Evolution /////////////////////////
	*stb_avg = (*stb_avg)/((double) tcntr);
	avg_amp = 0.0;
	for(i=0;i<N;++i)
	{	amp  = sqrt(fpGpsi[i+spd][0]*fpGpsi[i+spd][0] + fpGpsi[i+spd][1]*fpGpsi[i+spd][1]);
		avg_amp+=amp;

	}
	if(printfp)
	fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)N));

	
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
			dx = 1e-3;
			dt = 1e-6;
			en_loss = run(dt,dx,&abs_err,&stb_avg,0,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
		}

	}

	



}
