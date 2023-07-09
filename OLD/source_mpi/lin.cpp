using namespace std;

#define twopie 6.28
#define G 1.0

#include <math.h>
#include <stdio.h>






int main()
{
	double psi_k[3], psi_a_k[3], V_k[3],V_a_k[3];
	double kpsi_k[3], kpsi_a_k[3], kV_k[3];
	double psi_k_ini, psi_a_k_ini, V_k_ini;
	double k[3],ki[3]; 
	double delta[3],delta_ini;

	


	double da,a0,ai,H0,omega_dm_0;
	double fb_a,fb_a0,pow_arg;

	double h=0.7;
	double M4 = 1.0;
	
	double alpha = 1e5, w;


	int i;
		double a,z,a_t,ak,a_tt, acc1[3], acc2[3],kappa_lin,beta_lin,gamma_lin,omega_lin,omega,A;
		double ini_psi_k[3],ini_psi_a_k[3],ini_V_k[3];
		double ini_delta[3],U;

	da = 1e-4;

	double kiii = 0.01;


	ki[0] = kiii*h;  ki[1] = kiii*h;  ki[2] = kiii*h;

	ai = 0.001;
	a0 = 1.0;
	H0 = (h/2.99)*0.001;
	omega_dm_0 = 0.3;

	

	w = 1.0/(2.0*alpha-1.0);

	pow_arg = 3.0*(H0*H0)*(omega_dm_0)/(4.0*twopie*G*(2.0*alpha-1.0));

	a_t = ai*H0*sqrt(omega_dm_0*pow(a0/ai,3.0*(1.0+w))+ (1.0-omega_dm_0));
	
	fb_a0 = sqrt(2.0*pow_arg)/a_t ;
	fb_a = fb_a0*pow(a0/ai,3.0/(2.0*alpha-1.0));

	U = -3.0*(a_t/ai)/(2.0*alpha-1.0);

	delta_ini = ai;

	psi_k_ini = 0.0001;
	psi_a_k_ini = (-3.0*(a_t/ai)*alpha*fb_a*delta_ini/(2.0*alpha-1.0))/(U + 3.0*a_t/(ai*(2.0*alpha-1.0)) - 6.0*a_t*alpha/(ai*(2.0*alpha-1.0)));
	
	

	V_k_ini = -0.5*delta_ini + psi_a_k_ini/fb_a;
	
	printf("Initial:\ndelta_ini %lf\npsi_a_ini %lf\nV_ini %lf\nfb_a %.15lf\n",delta_ini,psi_a_k_ini,V_k_ini,fb_a);


		k[0] = ki[0]; k[1] = ki[1]; k[2] = ki[2];
		//k[0] = 0.5*h;
		psi_k[0] = psi_k_ini; psi_k[1] = psi_k_ini; psi_k[2] = psi_k_ini;
		psi_a_k[0] = psi_a_k_ini; psi_a_k[1] = psi_a_k_ini; psi_a_k[2] = psi_a_k_ini;

		V_k[0] = V_k_ini; V_k[1] = V_k_ini; V_k[2] = V_k_ini;
		
		delta[0] = 2.0*(-V_k[0] + psi_a_k[0]/fb_a);
		delta[1] = 2.0*(-V_k[1] + psi_a_k[1]/fb_a);
		delta[2] = 2.0*(-V_k[2] + psi_a_k[2]/fb_a);
	

		

		FILE *fp_lin[3];

		//da=0.01*da;
		for(i=0;i<3;++i)
		{
			ini_psi_k[i] = psi_k[i];
			ini_psi_a_k[i] = psi_a_k[i];
			ini_delta[i] = delta[i]; 
		}

		fp_lin[0] = fopen("linear_min.txt","w");
		fp_lin[1] = fopen("linear_mid.txt","w");
		fp_lin[2] = fopen("linear_max.txt","w");

		

		for(a=ai;a<=a0;a+=da)
		{
			a_t = a*H0*sqrt(omega_dm_0*pow(a0/a,3.0*(1.0+w))+ (1.0-omega_dm_0));
			omega = omega_dm_0*pow(a/a0,-3.0*(1.0+w))*H0*H0*a*a/(a_t*a_t);
			A = -1.5*(a/a0)*(1.0+w)*omega_dm_0*pow(a/a0,-3.0*(1.0+w)-1.0)/(omega_dm_0*pow(a/a0,-3.0*(1.0+w))+(1.0-omega_dm_0));
			fb_a = fb_a0*pow(a0/a,3.0/(2.0*alpha-1.0));

			a_tt = (1.0+A)*a_t*a_t/a;
			z = a0/a -1.0;
			
			kappa_lin = -3.0/(a*(2.0*alpha-1.0));
			beta_lin = 6.0*alpha*fb_a/(a*(2.0*alpha-1.0));
			gamma_lin = 2.0*fb_a*(1.0+alpha)/(2.0*alpha-1.0);
;
			
			ak = a + da;
	
			for(i=0;i<3;++i)
			{
			  
			  delta[i] = 2.0*(-V_k[i] + psi_a_k[i]/fb_a);

			  fprintf(fp_lin[i],"%lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n",a,a/ai,psi_a_k[i],delta[i],delta[i]/ini_delta[i]);

			  V_a_k[i] = -V_k[i]/(a) - k[i]*k[i]*V_k[i]/(3.0*a*a_t*a_t) - H0*H0*omega_dm_0*pow(a0/a,3.0*(1.0+w))*(-V_k[i] + psi_a_k[i]/fb_a)/a_t;
				
			  

			  acc1[i] = kappa_lin*psi_a_k[i] + beta_lin*V_k[i] + gamma_lin*V_a_k[i] - k[i]*k[i]*psi_k[i]/(a*a*a_t*a_t*(2.0*alpha-1.0))-a_tt*psi_a_k[i]/(a_t*a_t);

					  
			
			// if(a==ai)
			 //  printf("%lf for i %d k is %lf  (h/Mpc)\n",acc1[i],i,k[i]/h);
			// printf("i %d acc1 %.10lf  %.10lf  %.10lf\n",i,k[i],cs2,beta_lin);

			   V_k[i] = (V_k[i] - da*V_k[i]/a - da*H0*H0*omega_dm_0*pow(a0/a,3.0*(1.0+w))*(-V_k[i] + psi_a_k[i]/fb_a)/a_t)
														/(1.0 + k[i]*k[i]*da/(3.0*a*a_t*a_t));

			

			   psi_k[i] = ( psi_k[i] + da*(1.0+0.5*da*kappa_lin)*psi_a_k[i] + 0.5*beta_lin*da*da*V_k[i]
					+ 0.5*gamma_lin*da*da*V_a_k[i] - 0.5*a_tt*psi_a_k[i]*da*da/(a_t*a_t) )			
					/( 1.0 + 0.5*k[i]*k[i]*da*da/(a*a*a_t*a_t*(2.0*alpha-1.0))  );


		  	   if(a<(ai+5.0*da))
				printf("V%d %lf %.10lf\n",i,V_k[i],psi_a_k[i]);

		

			}




			a_t = ak*H0*sqrt(omega_dm_0*pow(a0/ak,3.0*(1.0+w))+ (1.0-omega_dm_0));
			omega = omega_dm_0*pow(a/a0,-3.0*(1.0+w))*H0*H0*ak*ak/(a_t*a_t);
			A = -1.5*(a/a0)*(1.0+w)*omega_dm_0*pow(ak/a0,-3.0*(1.0+w)-1.0)/(omega_dm_0*pow(ak/a0,-3.0*(1.0+w))+(1.0-omega_dm_0));
			a_tt = (1.0+A)*a_t*a_t/a;

			fb_a = fb_a0*pow(a0/ak,3.0/(2.0*alpha-1.0));
			
			kappa_lin = -3.0/(ak*(2.0*alpha-1.0));
			beta_lin = 6.0*alpha*fb_a/(ak*(2.0*alpha-1.0));
			gamma_lin = 2.0*fb_a*(1.0+alpha)/(2.0*alpha-1.0);
	
			for(i=0;i<3;++i)
			{

			   V_a_k[i] = -V_k[i]/(ak) - k[i]*k[i]*V_k[i]/(3.0*ak*a_t*a_t) - H0*H0*omega_dm_0*pow(a0/ak,3.0*(1.0+w))*(-V_k[i] + psi_a_k[i]/fb_a)/a_t;
				
			  
			  acc2[i] = kappa_lin*psi_a_k[i] + beta_lin*V_k[i] + gamma_lin*V_a_k[i] - k[i]*k[i]*psi_k[i]/(ak*ak*a_t*a_t*(2.0*alpha-1.0))
							-a_tt*psi_a_k[i]/(a_t*a_t);
				
			 

			 psi_a_k[i] = psi_a_k[i] + 0.5*da*(acc1[i]+acc1[i]);

				

			  

			}
			

		}


}
