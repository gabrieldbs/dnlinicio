#include <stdio.h> 
#include <math.h>
#include "stdlib.h"
#include "ejerciciosmapas.h"

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
	    double v_x,v_y,v_z, t,t_max;
		 float a_max,paso;
		 int iter,iteracion_deseada;
	    sscanf(argv[2], "%d", &iter);
	    sscanf(argv[3], "%lg", &v_x);
	    sscanf(argv[4], "%f", &a_max);
	    sscanf(argv[5], "%f", &paso);
	    sscanf(argv[6], "%d", &iteracion_deseada);
		 mapalog(iter,v_x,a_max,paso,iteracion_deseada);
  	}
	
  return 0;

}


