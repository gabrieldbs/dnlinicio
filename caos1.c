#include <stdio.h> 
#include <math.h>
#include "rk4.c"
void ecuaciones(int n, double v[], double dv[], double t)
{
  double x,y,z;
  x=v[0];
  y=v[1];
  z=v[2];
  dv[0]= 10*(y-x);
  dv[1]= 28*x-y-x*z;
  dv[2]=x*y-8*z/3.; 
  return;
}


int main(){
	
	FILE *ptr;
	double v[3],t,dt,t_pre,t_max;
	ptr=fopen("caos1.dat","w");

	v[0]=9.9;
	v[1]=10.;
	v[2]=40.;
	
	dt=0.001;
		
	t=0.;
	t_pre=0;
	t_max=100.;
			
	while(t<t_max){
		rk4(ecuaciones,v,3,t,dt);
		if(t>t_pre)
		  fprintf(ptr,"%lg\t%lg\t%lg\t%lg\n",t,v[0],v[1],v[2]);
		  t+=dt;
	         }	

	fclose(ptr);
	return(0);
}	
