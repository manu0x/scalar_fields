using namespace std;
#include "../include/include_custom.h"



//////////////////////////////Constants/////////////////////////

///////////////////////////////////////////Todo List///////////////////////////////////////////
//////////////   1) Add error check for successful memory allocation in class scalar_field_3d    /////////////////////////////////////////
//////////////   2) Check cosmo ini conditions as per fdm					 /////////////////////////////////////////
//////////////   3) Check k grid related things							 /////////////////////////////////////////




int main()
{
	int c1=1,c2=1,fail=0;
	int ind[3]{128,128,128};	
	int tN = ind[0]*ind[1]*ind[2];
	fdm_psi psi(ind,true);
	int use_omp{1};
	metric_potential phi(ind,true);

	const char name[] = "lcdm_00_pk.dat";
	printf("\n Building pk...\n");
	ini_power_generator pk(name);
	pk.check(0.0);
	gauss_rand_field_gen grf(ind);
	
	//gen.stats_check(3.1);

	

	double a0,ai,omega_dm_ini;
	double k_grid[tN][3],dx[3],dk;
	int kbins,kbin_grid[tN];
	
	set_back_cosmo(a0,ai,Hi,omega_dm_ini);
	printf("Hi %lf\nOmega_dm_ini %lf\nai %lf\n",Hi,omega_dm_ini,ai);

	initialise(ind,psi,phi,k_grid,kbin_grid,a0,ai,Hi,omega_dm_ini,dx,dk,kbins,pk,grf);
	//printf("\ndk is %lf\n",dk);
	
	//if(use_omp)
	// fail = evolve_kdk_openmp(ind,psi,phi,k_grid,kbin_grid,a0,ai,a0,omega_dm_ini,dx,dk,kbins,0.4e-4);
	//else
	// fail = evolve_kdk(ind,psi,phi,k_grid,kbin_grid,a0,ai,a0,omega_dm_ini,dx,dk,kbins,0.4e-4);
	//*/
	printf("fail is %d\n",fail);
	
	
}






