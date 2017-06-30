#include <stdio.h>
#include <math.h>
#include "stdlib.h"
#include "time.h"
#include "ejercicios.h"
void rk4(void deri(int , double [], double [], double ), double h[], int n, double t, double dt);
void ecuaciones1a(int n, double v[], double dv[], double t);
int ej_1a(double v_x,double v_y,double v_z,double t,double  paso,double t_max);


int main(int argc, char **argv) {
  int Programa;
  if(argc>1){
    sscanf(argv[1], "%x", &Programa);
	printf("%d \n",Programa);
  }else{
    Programa =0;
  }
	if(Programa == 26){
    	 printf("Ejecutando simulacion ejercicio 1.a)\n");
	    double v_x,v_y,v_z, t,t_max,paso;
	    sscanf(argv[2], "%lg", &v_x);
	    sscanf(argv[3], "%lg", &v_y);
	    sscanf(argv[4], "%lg", &v_z);
	    sscanf(argv[5], "%lg", &t);
	    sscanf(argv[6], "%lg", &paso);
	    sscanf(argv[7], "%lg", &t_max);
   	 ej_1a(v_x, v_y,v_z,t, paso, t_max);
  	}

	if(Programa == 42){
    	 printf("Ejecutando simulacion ejercicio 2.a)\n");
	    double v_x,v_y,v_z, t,t_max,paso,c	;
	    sscanf(argv[2], "%lg", &v_x); //6
	    sscanf(argv[3], "%lg", &v_y);//-4
	    sscanf(argv[4], "%lg", &v_z);//2
	    sscanf(argv[5], "%lg", &t);//0
	    sscanf(argv[6], "%lg", &paso);//0.001
	    sscanf(argv[7], "%lg", &t_max);//300
	    sscanf(argv[8], "%lg", &c);//5
   	 ej_2a(v_x, v_y,v_z,t, paso, t_max,c);
  	}

  if(Programa == 43){
    	 printf("Ejecutando simulacion ejercicio 2.b)\n");
	    double v_x,v_y,v_z, t,t_max,paso,c	;
	    sscanf(argv[2], "%lg", &v_x); //-1 //6
	    sscanf(argv[3], "%lg", &v_y); //-5 //-4 
	    sscanf(argv[4], "%lg", &v_z); //2 //2
	    sscanf(argv[5], "%lg", &t); //0
	    sscanf(argv[6], "%lg", &paso); //0.001
	    sscanf(argv[7], "%lg", &t_max); //300
	    sscanf(argv[8], "%lg", &c); //5
   	 ej2_b(v_x, v_y,v_z,t, paso, t_max,c);
  	}

  if(Programa == 44){
    	 printf("Ejecutando simulacion ejercicio 2.c)\n");
	    double v_x,v_y,v_z, t,t_max,paso,c	;
	    sscanf(argv[2], "%lg", &v_x); //-1 
	    sscanf(argv[3], "%lg", &v_y); //-5
	    sscanf(argv[4], "%lg", &v_z); //2
	    sscanf(argv[5], "%lg", &t); //0
	    sscanf(argv[6], "%lg", &paso); //0.001
	    sscanf(argv[7], "%lg", &t_max); //300
	    sscanf(argv[8], "%lg", &c); //5
   	 ej2_c_2(v_x, v_y,v_z,t, paso, t_max,c);
  	}
  return 0;


}
