#include <stdio.h> 
#include <math.h>

#define c 5

void ecuaciones(int n, double v[], double dv[], double t,double A, double w )
{
  double x,y,z;
  x=v[0];
  y=v[1];
  z=v[2];
  dv[0]= -y-z;
  dv[1]= x+0.2*y;
  dv[2]= 0.2+z*(x-c); 

  return;
}


int main(){
	FILE *ptr;
	double v[3],t,dt,t_pre,t_max;
	ptr=fopen("caos2a.dat","w");
	v[0]=-1.;
	v[1]=-5.;
	v[2]=2.;
        dt=0.001;
	t=0.;
	t_pre=0;
	t_max=300;	
	while(t<t_max){
		rk4(ecuaciones,v,3,t,dt);
		if(t>t_pre) fprintf(ptr,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
		t+=dt;
	}	

	fclose(ptr);
	return(0);

}


