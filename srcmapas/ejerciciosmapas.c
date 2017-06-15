// ejerciciomapas
#include <stdio.h> 
#include <math.h>
#include "stdlib.h"
#include "ejerciciosmapas.h"

//	mapa logistico

int mapalog(int iter, double v_x,float a_min,float a_max,float paso,int iteracion_deseada){
	double *x;
	FILE *pt;	
	pt=fopen("mapaslog.dat","w");
	FILE *ptt;
	ptt=fopen("mapaslogej.dat","w");
	fprintf(pt,"Iteraciones de X \n");
	fprintf(ptt,"Corrida en a \n");
	x =(double *) malloc((iter)*sizeof(double));
	fprintf(ptt," Valores de la iteracion %d  de x \n", iteracion_deseada);
	for(float l=a_min;l<a_max;l=l+paso){	
		for(int i=0;i<iter;i++){
			x[i]=0;}
		x[0]=v_x;
		for(int i=1;i<iter;i++){
			x[i]=l*x[i-1]*(1-x[i-1]);}
		for(int j=0;j<iter;j++){
			fprintf(pt,"%lg  ,",x[j]);}
		fprintf(pt," \n");		
	printf("Terminado %f \n",l);
	fprintf(ptt,"%lg  ,", x[iteracion_deseada]);
	free(x);
	}
	fprintf(ptt,"\n");
	fprintf(ptt,"Valores de a \n");
	for(float l=a_min;l<a_max;l=l+paso){
		fprintf(ptt," %f  ,",l);}
	
	fprintf(pt,"\n");
	fprintf(ptt,"\n");
	free(x);
	fclose(pt);
	fclose(ptt);
	return(0);
}
