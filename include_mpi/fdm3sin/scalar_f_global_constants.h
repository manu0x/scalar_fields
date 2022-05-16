double hbar_by_m,h,H0,Hi,lenfac;
double c_box,pc_box,hbar_box;
double w,cs2;
double space_mpc_to_dimless;
double alpha;
int method;

double nwave[3], nwave_amp;
double T;
double tot_steps;

FILE *fp_sim_info;
////////////////////// MPI Globals....///////////////////
MPI_Datatype c_x_plain,c_y_plain,c_z_plain;
MPI_Comm cart_comm;
/////////////////////////////////////////////////////////

#define twopie  2.0*M_PI
#define mass 1.0
#define G 1.0

