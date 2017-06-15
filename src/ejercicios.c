#include <stdio.h>
#include <math.h>
#include "stdlib.h"
#include "time.h"
#include "ejercicios.h"
void ecuaciones1a(int n, double v[], double dv[], double t);
void ecuaciones2a(int n, double v[], double dv[], double t);

void ecuaciones1a(int n, double v[], double dv[], double t){
  double x,y,z;
  x=v[0];
  y=v[1];
  z=v[2];
  dv[0]= 10*(y-x);
  dv[1]= 28*x-y-x*z;
  dv[2]=x*y-8*z/3.;
  return;
}


int ej_1a(double v_x,double v_y,double v_z,double t,double paso,double  t_max){
	double v[3],dt,t_pre,T_max ,T;
	FILE *ptr;
	ptr=fopen("caos1.dat","w");
	FILE *pt;
	pt=fopen("caos1t.dat","a");
	FILE *ptt;
	ptt=fopen("caos1x.dat","a");
	fprintf(pt,"Barrido temporal /n");
	fprintf(ptt,"Barrido x");
	v[0]=v_x;	//10
	v[1]=v_y;	//10
	v[2]=v_z; 	//40
	dt=paso;		//0.001
	T = t;		//0.
	t_pre=0;
	T_max=100.;
	while(t<t_max){
		rk4(ecuaciones1a,v,3,t,dt);
		if(t>t_pre){
			fprintf(ptr,"%lg\t%lg\t%lg\t%lg\n",t,v[0],v[1],v[2]);
			fprintf(pt,"%lg ,",t);
			fprintf(ptt,"%lg ,",v[0]);
		}
		  t+=dt;
	}
	fprintf(pt,"\n");
	fprintf(ptt,"\n");

	fclose(ptr);
	fclose(pt);
	fclose(ptt);
	return(0);
}


void ecuaciones2a(int n, double v[], double dv[], double t)
{
  double x,y,z;
  x=v[0];
  y=v[1];
  z=v[2];
  dv[0]= -y-z;
  dv[1]= x+0.2*y;
  dv[2]= 0.2+z*(x-5);  // c va  haber que cambiarlo a mano

  return;
}
int ej_2a(double v_x,double v_y,double v_z,double T,double paso,double  t_max, double c){
	FILE *ptr;
	double v[3],t,dt,t_pre;
	ptr=fopen("caos2a.dat","w");
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//2
   dt=paso;	//0.001
	t= T;//0
	t_pre=0;
	// t_ max = 300
	while(t<t_max){
		rk4(ecuaciones2a,v,3,t,dt);
		if(t>t_pre) {
			fprintf(ptr,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
		t+=dt;
	}

	fclose(ptr);
	return(0);
}

int ej2_b(double v_x,double v_y,double v_z,double T,double paso,double  t_max, double c){
	FILE *ptr;
	double v[3],t,dt,t_pre,*x, y_pre;
  double a;
  ptr=fopen("caos2b.dat","a");
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//3
   dt=paso;	//0.001
	t= T;//0
	t_pre=0;
	// t_ max = 300
  while(t<t_max){
		rk4(ecuaciones2a,v,3,t,dt);
    if(t>t_pre & v[0]>0  & v[1]*y_pre <0) {
			fprintf(ptr,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
    y_pre=v[1];
		t+=dt;
	}
  fprintf(ptr,"\n");
	fclose(ptr);
	return(0);
}
/*
planos
a*x+b*y+c*z =d;
c=0
en el caso puesto es a =0
hacer para variios
*/


/*
ejercicios c
es ver un punto luego si el siguietne es parecido si eso pasa es de perioso 1
ir al graf grande y listo

 [2,2,5]
[2,2,5.1] periodo 1
 si pàsa en 2 es de periodo dos
 [2,2,5]
[2.3.3]
 [2,2,5.1] periodo 1


*/
