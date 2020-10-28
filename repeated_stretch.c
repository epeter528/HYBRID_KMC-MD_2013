#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kmc.h"

int repeated_stretch(int argc,char *argv[],int scansteps,int range_checker_A,int range_checker_B,float energy_scan2[1000][7][50],int num_list1,float sim_temp,float expconst,float OMEGA_B[7],int event_or_scan,
  float eigenval2_x,float eigenval2_y,float eigenval2_z, float *eig_x, float *eig_y, float *eig_z,int num_event){
  
  int a,i,k,l,p,m;

  char char1[50000][100]; 

  char char5[200];
  
  int o;
  
  char a1[3];
  
  float eig1,eig2,eig3;

  char charA[50000][30];
   
  char res_type[20000][5],atom_type[20000][5];
  
  int atom_num[20000];
  
  float atom_x[50000],atom_y[50000],atom_z[50000];   

  float box_x,box_y,box_z;
  
  char *trajfile, *ndxfile;
  
  int range_1or2;

  int res_number[20000];
  
  char chari;  
  
  int don,acc,hyd;
  
  int num_hbond2;
   
  FILE *fp;

  int res_atom_number;
  
  int res_id2;
  
  float don_x,don_y,don_z,acc_x,acc_y,acc_z;
  
  float D_x2,D_y2,D_z2,D_A2;
  
  float atom_x2[50000],atom_y2[50000],atom_z2[50000];
  
  float scale,scalekey;
  
  float energy_val,scansteps2;

  scansteps2 = (float) scansteps;
  
  float OMEGA_B_VAL, OMEGA_ZERO_VAL;

  float x;

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
    
	
        if(i<=a-2)     { 

            sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x[i],&atom_y[i],&atom_z[i]);	    
	    
                	};
	   
        if(i < 10000)  {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                        };	    
	if(i==a-1)      {

 	              sscanf(char1[i],"%f%f%f%s[^\n]",&box_x,&box_y,&box_z,&chari);

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
            
          for(i=2;i<=a;i++){
    
	
           if(i<=a-2) sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x[i],&atom_y[i],&atom_z[i]);
	     // sscanf(char1[i],"%5d%-5s%5s%5d%8.3f%8.3f%8.3f",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i],&atom_x[i],&atom_y[i],&atom_z[i]);
	   
           if(i < 10000) {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                         };	   
	   
           if(i==a-1)    {

	            sscanf(char1[i],"%f%f%f%s[^\n]",&box_x,&box_y,&box_z,&chari);

                         };
 	
         };   

         fclose(fp);
    
    
  };
  
  fp = fopen("t.ndx","w");
  
  
      fprintf(fp,"%s\n","[ 1 ]");  
      fprintf(fp,"%s\n","[ group ]");
  
      for(i=2;i<=a-2;i++){

	if(res_number[i] == range_checker_A) {
	  
	   if(atom_num[i] > 0) fprintf(fp,"%d\n",atom_num[i]);
	  
	};
	
      };

  fclose(fp);     

  system("cat t.ndx index.ndx > 2.ndx");
  
     if(k == 1) { 
      
 //    gmx_rmsf3(argc,argv,range_1or2,trajfile,ndxfile,&eigenval2_x,&eigenval2_y,&eigenval2_z);
 
   fp = fopen("rmsf.sh","w");
   
   fprintf(fp,"%s\n","#!/bin/bash");
   fprintf(fp,"%s\n","~/GRO453/bin/g_rmsf -f md_traj.trr -n 2.ndx -s run.tpr -dir eigrmsf.log << EOF");
   fprintf(fp,"%s\n","1");
   fprintf(fp,"%s\n","1");   
   fprintf(fp,"%s\n","<< EOF");
   
   fclose(fp);
   
   system("chmod 744 rmsf.sh");

   system("./rmsf.sh"); 

 //    fp = fopen("eigrmsf.log","r");
 
     if(event_or_scan == 0){


    o = 0; 
     
    fp = fopen ("eigrmsf.log","r");
  
  while(fgets(char5,sizeof(char5),fp) != NULL){
    
        o++;
    
	if(o >= 10){
	  
	  sscanf(char5,"%s%f%f%f",&a1,&eig1,&eig2,&eig3);
	  
//	  printf("%s%f%f%f\n",a,eig1,eig2,eig3);
	  
	  if(o == 10) eigenval2_x = eig1;
	  if(o == 11) eigenval2_y = eig1;
	  if(o == 12) eigenval2_z = eig1;
	  
	};
    
  };
  
  fclose(fp); 
     
     *eig_z = eigenval2_z;
     *eig_y = eigenval2_y;
     *eig_x = eigenval2_x;     
     
     };
     
    
     printf("%f\n",x);     
     printf("%f%f%f\n",eigenval2_x,eigenval2_y,eigenval2_z);
     
     range_1or2 = 2;
     
//     gmx_hbond2(argc,argv,range_1or2,trajfile,ndxfile,&don,&hyd,&acc,&num_hbond2);
  
     don = 1;
     acc = 1;
     
     };
     
     m = 1;
     
     if(don == 1 && acc == 1){
       
      for(i=2;i<=a-2;i++){       

	if(res_number[i] == range_checker_A && m ==1 ){
	  
	      don = atom_num[i];
	      hyd = atom_num[i+1];
	  if(don != 1)    acc = 1;
	  if(don == 1)    acc = atom_num[a-2];
	      
	      m++;
       
       
	};
	
      };
//       for(i=1;i<=scansteps;i++) {
	
//        energy_scan2[num_list1][num_event][i] = 10000;	 

//	printf("%s\n","h11");
	
//      };    
       
//      return 0;
      
    };

    
      D_x2 = eigenval2_x;
      D_y2 = eigenval2_y;
      D_z2 = eigenval2_z;
      
      D_A2 = sqrt(pow(D_x2,2)+pow(D_y2,2)+pow(D_z2,2));
      
      scalekey = expconst/(scansteps*D_A2);
//      scalekey = 0.0011;
  
      printf("%f\n",scalekey);
             
      fp = fopen("2.ndx","w");
      
      fprintf(fp,"%s\n","[ group ]");
      
      fprintf(fp,"%d\n%d\n%d\n",don,hyd,acc);
  
//      for(i=2;i<=a-2;i++) {
	
//	  if(res_number[i] == aminoresnum) fprintf(fp,"%d\n",atom_num[i]);
	
//      };      
      
      fprintf(fp,"%s\n","[ protein ]");
      
      for(i=2;i<=a-2;i++) {
	
	  fprintf(fp,"%d\n",atom_num[i]);
	
      };
      
       fprintf(fp,"%s\n","[ system ]");
      
      for(i=2;i<=a-2;i++) {
	
	  fprintf(fp,"%d\n",atom_num[i]);
	
      };     
    
      fclose(fp);
      res_atom_number = p-1;
      
      scale = 0;      

           for(i=2;i<=a-2;i++) {
	
        	if(atom_num[i] == don) res_id2 = res_number[i];
	
            };
      
            p = 1;
      
            for(i=2;i<=a-2;i++) {
	 
           	if(res_number[i] == res_id2) {
	  
	  
	           atom_x2[p] = scale*scalekey*(eigenval2_x)+atom_x[i];
	           atom_y2[p] = scale*scalekey*(eigenval2_y)+atom_y[i];
	           atom_z2[p] = scale*scalekey*(eigenval2_z)+atom_z[i];	  
	  
		   p++;
		   
		};   
		   
	  };
	  
	  
          res_atom_number = p-1;
	  
	
	scale = scale + 1/scansteps2;
	
        m = 1;           
           
          for(i=1;i<=a-2;i++){
	    
	    if(res_number[i] == res_id2 && m==1) {
	      
	       for(p=1;p<=res_atom_number;p++){
		
		 atom_x[i+p-1] = atom_x2[p];
		 atom_y[i+p-1] = atom_y2[p];
		 atom_z[i+p-1] = atom_z2[p];
		 
		 m++;
		 
	      };
	      
	    };
	    
	  }; 

	  
	  fp = fopen("minimized.gro","w");
	  
	  fprintf(fp,"%s\n","");
	  fprintf(fp,"%d\n", a - 3);
  
	  for(i=2;i<=a-2;i++){
	    
	    fprintf(fp,"%s%8.3f%8.3f%8.3f\n",charA[i],atom_x[i],atom_y[i],atom_z[i]);

	  };
	  
	  fprintf(fp,"%s",char1[a-1]);
	  
	  fclose(fp);
	  
	  	  
	  i = energy_shift(argc,argv,&energy_val,sim_temp,int_ex);

    
          if(int_ex == 1) return 0;

          energy_scan2[num_list1][num_event][k] = energy_val;          
	  
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
	
    i = kramer_search(argc,argv,&OMEGA_B_VAL,sim_temp,range_checker_A);

    };
    
    OMEGA_B[num_event] = OMEGA_B_VAL;
    
    //   free(atom_x);
   //   free(atom_y);
   //   free(atom_z);
   //   free(res_number);
   //   free(atom_num);  

     
}


