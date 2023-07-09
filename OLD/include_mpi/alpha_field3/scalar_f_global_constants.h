double h,H0,Hi,space_mpc_to_dimless,lenfac;
double c_box,pc_box,hbar_box;
double alpha,w,cs2;

int X_POWER, binomial_n;

FILE *fp_sim_info;

////////////////////// MPI Globals....///////////////////
MPI_Datatype c_x_plain,c_y_plain,c_z_plain;
MPI_Comm cart_comm;
/////////////////////////////////////////////////////////

#define twopie  2.0*M_PI
#define Mfield 1.0
#define G 1.0

