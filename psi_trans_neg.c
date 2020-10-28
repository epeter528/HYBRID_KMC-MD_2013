#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kmc.h"


int psi_trans_neg(int argc,char *argv[],int scansteps,int aminoresnum,float energy_scan2[1000][7][50],int num_list1,float delta_phi,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event)
{
  
  float ntrans;
  
  int a,i,k,l,p,m;

  char char1[50000][100];
  
  char charA[50000][30];  
  
  char res_type[50000][5],atom_type[50000][5];
  
  int atom_num[50000];
  
  float atom_x[50000],atom_y[50000],atom_z[50000];
  
  float atom_x2[50000],atom_y2[50000],atom_z2[50000];  
  
  ntrans = (float) scansteps;

  int res_number[50000];

  FILE *fp;
  
  float PI = 3.141592654;
  
  float energy_val;  

  char chari;
  
  delta_phi = delta_phi*(-1);
  
  float  atom_1_x,atom_1_y,atom_1_z;
  float  atom_2_x,atom_2_y,atom_2_z;  
  float  atom_3_x,atom_3_y,atom_3_z;
  float  atom_4_x,atom_4_y,atom_4_z;
  float  diff_x1,diff_y1,diff_z1;
  float  diff_tot;
  float  n1,n2,n3;
  
  float OMEGA_B_VAL, OMEGA_ZERO_VAL;  
  float box_x,box_y,box_z;

  int int_ex;
  
for(k=1;k<=scansteps;k++){ 

  if(k==1){
  
  fp = fopen("minimized3.gro","r");
      
      a = 0;
      
      while(fgets(char1[a],sizeof(char1[a]),fp) != NULL){
	
	a ++;
	
      };
      
 //     atom_x = malloc(a*sizeof(float));
 //     atom_y = malloc(a*sizeof(float));
 //     atom_z = malloc(a*sizeof(float));
 //   res_number = malloc(a*sizeof(int));
 //     atom_num   = malloc(a*sizeof(int));
  //    res_type = malloc(1000);
  //    atom_type = malloc(1000);
            
      for(i=2;i<=a;i++){
    
	
        if(i<=a-2)  { 

            sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x[i],&atom_y[i],&atom_z[i]);	  
	  
//            printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",res_number[i],res_type[i],atom_type[i],atom_num[i],atom_x[i],atom_y[i],atom_z[i]);
	};
	
	   
        if(i < 10000) {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                      };	
	    
	if(i==a-1)  {

	            sscanf(char1[i],"%f%f%f%s[^\n]",&box_x,&box_y,&box_z,&chari);
	//            printf("%f%f%f\n",box_x,box_y,box_z);
      };
 	
      };   

  fclose(fp);

  };
  
  if(k>=2){

  
  fp = fopen("minimized.gro","r");
      
      a = 0;
      
      while(fgets(char1[a],sizeof(char1[a]),fp) != NULL){
	
	a ++;
	
      };
      
 //     atom_x = malloc(a*sizeof(float));
 //     atom_y = malloc(a*sizeof(float));
 //     atom_z = malloc(a*sizeof(float));
 //   res_number = malloc(a*sizeof(int));
 //     atom_num   = malloc(a*sizeof(int));
  //    res_type = malloc(1000);
  //    atom_type = malloc(1000);
            
      for(i=2;i<=a;i++){
    
	
        if(i<=a-2)  { 

            sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x[i],&atom_y[i],&atom_z[i]);	  
	  
//            printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",res_number[i],res_type[i],atom_type[i],atom_num[i],atom_x[i],atom_y[i],atom_z[i]);
	};
	
	   
        if(i < 10000) {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                      };	
	    
	if(i==a-1)  {

	            sscanf(char1[i],"%f%f%f%s[^\n]",&box_x,&box_y,&box_z,&chari);
	//            printf("%f%f%f\n",box_x,box_y,box_z);
      };
 	
      };   

  fclose(fp);
    
    
  };

        for(m=2;m<=a-2;m++){
	    
	   
	//    printf("%s\n",&atom_type[m]);
	  if(aminoresnum+1 == res_number[m]) {
	    
	        if(strcmp(atom_type[m],"C")==0) {
	    
	           atom_3_x = atom_x[m];
		   atom_3_y = atom_y[m];
		   atom_3_z = atom_z[m];
		
	 	};
		
	  };
		
	if(aminoresnum+1 == res_number[m]){
	  
	        if(strcmp(atom_type[m],"C")==0){
		  
		   atom_1_x = atom_x[m];
		   atom_1_y = atom_y[m];
		   atom_1_z = atom_z[m];		  
		}; 
		
		if(strcmp(atom_type[m],"CA")==0){

		   atom_2_x = atom_x[m];
		   atom_2_y = atom_y[m];
		   atom_2_z = atom_z[m];		  
		  
		}; 
		  
	      }; 
	      
	 };
	 
	 
diff_x1   = atom_1_x - atom_2_x;
diff_y1   = atom_1_y - atom_2_y;
diff_z1   = atom_1_z - atom_2_z;

diff_tot    = sqrt(pow(diff_x1,2)+pow(diff_y1,2)+pow(diff_z1,2));

n1      = diff_x1 / diff_tot;
n2      = diff_y1 / diff_tot;
n3      = diff_z1 / diff_tot;	    
	   
    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum) {
	
          atom_x[i] = atom_x[i] - atom_3_x;
          atom_y[i] = atom_y[i] - atom_3_y;
          atom_z[i] = atom_z[i] - atom_3_z;	  
	
      };
      
    };
    
    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum) {      
      
	 atom_x2[i] = atom_x[i]*(cos(delta_phi*PI/180) + pow(n1,2)*(1 - cos(delta_phi*PI/180))) + atom_y[i]*((n2*n1)*(1 - cos(delta_phi*PI/180))
	   	- n3 * sin(delta_phi*PI/180)) + atom_z[i]*((n3*n1)*(1 - cos(delta_phi*PI/180))+ n2*sin(delta_phi*PI/180));

	 atom_y2[i] = atom_x[i]*(n1*n2*(1 - cos(delta_phi*PI/180))+n3*(sin(delta_phi*PI/180))) + atom_y[i]*( cos(delta_phi*PI/180) + pow(n2,2)*	
	        (1 - cos(delta_phi*PI/180))) + atom_z[i]*(n3*n1*(1 - cos(delta_phi*PI/180))-n1*sin(delta_phi*PI/180));
		
	 atom_z2[i] = atom_x[i]*(n1*n3*(1 - cos(delta_phi*PI/180))-n2*sin(delta_phi*PI/180)) + atom_y[i]*( n2*n3*	(1 - cos(delta_phi*PI/180)) +
	        n1*sin(delta_phi*PI/180)) + atom_z[i]*( cos(delta_phi*PI/180) + pow(n3,2) *(1 - cos(delta_phi*PI/180)));
      };
      
    };  

    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum) {
	
          atom_x2[i] = atom_x2[i] + atom_3_x;
          atom_y2[i] = atom_y2[i] + atom_3_y;
          atom_z2[i] = atom_z2[i] + atom_3_z;	  
	
      };
      
    };    
    
     fp = fopen("minimized.gro","w");
    
    fprintf(fp,"%s\n","intermediate-struct");
    fprintf(fp,"%d\n", a-3);
    
    for(i=2;i<=a-2;i++){
      
      if(res_number[i] == aminoresnum) {
	
	 atom_x[i] = atom_x2[i];
	 atom_y[i] = atom_y2[i];
	 atom_z[i] = atom_z2[i];
	 
      };

       fprintf(fp,"%s%8.3f%8.3f%8.3f\n",charA[i],atom_x[i],atom_y[i],atom_z[i]);      
      
    };
    
    fprintf(fp,"%s",char1[a-1]);

    fclose(fp);
    
/*    fp = fopen("t.ndx","w");
	
	fprintf(fp,"%s\n","[ group ]");
	
	for(m=2;m<=a-2;m++){
	  
	  if(res_number[m] == aminoresnum+1) {

	   if(strcmp(atom_type[m],"N")==0){	    
	    
	    if(atom_num[m] > 0) fprintf(fp,"%d\n",atom_num[m]);
	    
	   };
	   
	   if(strcmp(atom_type[m],"CA")==0){	    
	    
	    if(atom_num[m] > 0) fprintf(fp,"%d\n",atom_num[m]);
	    
	   };		    
	    
	   if(strcmp(atom_type[m],"C")==0){	    
	    
	    if(atom_num[m] > 0) fprintf(fp,"%d\n",atom_num[m]);
	    
	   };
	    
	  };
	  
	  if(res_number[m] == aminoresnum+2) {
	    
	   if(strcmp(atom_type[m],"N")==0){	    
	    
	    if(atom_num[m] > 0) fprintf(fp,"%d\n",atom_num[m]);
	    
	   };	   
	    
	  };	  

	};
	
	fclose(fp);
	
	system("cat t.ndx index.ndx > 2.ndx"); */
	
    
    i = energy_shift(argc,argv,&energy_val,sim_temp,int_ex);


   if(int_ex == 1) return 0;
    
    energy_scan2[num_list1][num_event][k] = energy_val;
    
    printf("%f%s\n",energy_scan2[num_list1][num_event][k],"energy");    
    
    
};


  	fp = fopen("1.ndx","w");
	
	fprintf(fp,"%s\n","[ group ]");
	
	for(m=2;m<=a-2;m++){
	  
//	  if(res_number[m] == aminoresnum) {
	    
	    fprintf(fp,"%d\n",atom_num[m]);
	    
//	  };
	  
	};
	
	fclose(fp);
	  
    if(event_or_scan == 0){
	
    i = kramer_search(argc,argv,&OMEGA_B_VAL,sim_temp,aminoresnum);

    };
    
    OMEGA_B[num_event] = OMEGA_B_VAL;


    
}
