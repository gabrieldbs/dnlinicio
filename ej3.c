#include <stdio.h> 
#include <math.h>
#include "rk4.c"
#define a -1
// Ecuaciones del sistema
void ecuaciones(int n, double v[], double dv[],double t);

//Programa principal
int main(){
  int j;
  FILE *ptr;
  double i, v[2],t,dt,t_pre,t_max;
  //Archivo de salida
  ptr=fopen("ej3.dat","w");
  dt=0.01;
  t_max=5;
  //Condiciones Iniciales 
  for(i=-10;i<10;i+=0.1){
    v[0]=i;
    v[1]=0;			
    t=0.;					
    while(t<t_max){
      //Integra las ecuaciones utilizando el metodo de Runge Kutta
      rk4(ecuaciones,v,2,t,dt);
      //Imprime la integracion
      fprintf(ptr,"%lg\t%lg\t%lg\n",t,v[0],v[1]);
      t+=dt;
      }
    fprintf(ptr,"\n");
   }
  // printf("%g \n",mu);
  fclose(ptr);
  return(0);
}

void ecuaciones(int n, double v[], double dv[],double t){
  double x,y;
  double mu=4;
  x=v[0];
  y=v[1];
  dv[0]= y-x*x*x+x;
  dv[1]= -0.1*x;
  return;
}

