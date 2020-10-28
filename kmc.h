#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>


int write_mdp(int mdtime,float Temp,int imp);

int random_number(float *x);

int kmc(int argc,char *argv[]);

int mdprod(int argc,char *argv[],int mdtime,float Temp,int imp,char *grofile);

int noe_sep(int argc,char *argv[],int range_checker_A,int range_checker_B,
int eventmax);

//int read_xvg_files(char *xvgfile,int *sizeof_x,float x[100000],float y[100000]);

int dih_trans(int argc,char *argv[],float stddev_dih[10000],int *sizeof_stddev);

int dihedralbreak(int argc,char *argv[],int scansteps,int aminoresnum,float rot_change,float energy_scan2[1000][7][50],int num_list1,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event);

int energy_shift(int argc,char *argv[],float *energy_val,float sim_temp,int int_ex);

int replex(int argc,char *argv[],float T1,float T2,float T3,float T4,float T5,float *sim_temp,int mdtime,float Temp_orig,int imp);

int break_trans(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float sim_temp,float expconst,float OMEGA_B[7],int event_or_scan,
		float eigenval2_x,float eigenval2_y,float eigenval2_z, float *eig_x, float *eig_y, float *eig_z,int num_event);

int form_trans(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float sim_temp,float expconst,float OMEGA_B[7],int event_or_scan,
                float eigenval2_x,float eigenval2_y,float eigenval2_z,float *eig_x2,float *eig_y2,float *eig_z2,int num_event);

int phi_trans(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float delta_phi,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event);

int psi_trans(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float delta_phi,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event);

int phi_trans_neg(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float delta_phi,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event);

int psi_trans_neg(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float delta_phi,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event);

int delta_g_search(int icounter,int scansteps,float energy_scan2[1000][7][50],int *execsteps,int *selectamino,
		   int *selectevent,float sim_temp,float OMEGA_B[7],float OMEGA_ZERO[7],float *delta_time,int number_of_event);

int dih_processor(int k,int aminoresnum,float rot_change,int scansteps,int num_event);

int kramer_search(int argc,char *argv[],float *OMEGA_B_VAL,float sim_temp,int aminoresnum);

int kramer_zero(int argc,char *argv[],float *OMEGA_B_VAL,float sim_temp,int aminoresnum);

int solvent_config_separator(int *num_na,int *num_cl,int imp);

int solvent(int argc,char *argv[],int imp);
 
int solvent2(int argc,char *argv[],int imp);
 
int minimize_struct(int argc,char *argv[],int int_ex);

int read_xvg_files_int(char *xvgfile,int *sizeof_x,int x[100000],float y[100000]);

int barint(float *free_energy);

int repeated_stretch(int argc,char *argv[],int scansteps,int range_checker_A,int range_checker_B,float energy_scan2[1000][7][50],int num_list1,float sim_temp,float expconst,float OMEGA_B[7],int event_or_scan,
		float eigenval2_x,float eigenval2_y,float eigenval2_z, float *eig_x, float *eig_y, float *eig_z,int num_event); 
