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
	ptr=fopen("caos2a_2.5.dat","w");
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
	double v[3],t,dt,t_pre,*x, y_pre,x_pre;
  double a;
	FILE *pt;
  pt=fopen("caos2b_2.dat","w");
	FILE *ptt;	
  ptt=fopen("caos2b_3.dat","w");
  ptr=fopen("caos2b_1.dat","w");
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//2
   dt=paso;	//0.001
	t= T;//0
	t_pre=0;
	// t_ max = 300
  while(t<t_max){
		rk4(ecuaciones2a,v,3,t,dt);
    if(t>t_pre & v[0]>0  & v[1]*y_pre <0) {
			fprintf(ptr,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
	if(t>t_pre  & v[1] >0& v[0]*x_pre <0) {
			fprintf(pt,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
	if(t>t_pre  & v[0] >0& sqrt((v[0]-v[1])*(v[0]-v[1])) <0.01) {
			fprintf(ptt,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
    y_pre=v[1];
	 x_pre=v[0];
		t+=dt;
	}
   fprintf(ptr,"\n");
	fclose(ptr);   
	fprintf(ptt,"\n");
	fclose(ptt);
   fprintf(pt,"\n");
	fclose(pt);
	return(0);
}
/*
planos
a*x+b*y =d;
en el caso puesto es a =0
hacer para variios
*/



int ej2_c(double v_x,double v_y,double v_z,double T,double paso,double  t_max, double c){
	double v[3],t,dt,t_pre,*x, y_pre,x_pre, z_pre;
  double a;
	FILE *pt;
  pt=fopen("caos2c_1.dat","w");
	
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//2
   dt=paso;	//0.001
	t= T;//0
	t_pre=0;
	// t_ max = 300
  while(t<t_max){
		rk4(ecuaciones2a,v,3,t,dt);
    if(t>t_pre & sqrt((v[0]-x_pre)*(v[0]-x_pre)+(v[1]-y_pre)*(v[1]-y_pre)+(v[2]-z_pre)*(v[2]-z_pre))<0.1) { //0.00295
			fprintf(pt,"%lg\t%lg\t%lg\n%lg\t%lg\t%lg\n",v[0],v[1],v[2],x_pre,y_pre,z_pre);
			}
	 z_pre=v[2];
    y_pre=v[1];
	 x_pre=v[0];
		t=t+dt;
		t_pre=0.01;	
	}
   
   fprintf(pt,"\n");
	fclose(pt);


	return(0);
}


int ej2_c_2(double v_x,double v_y,double v_z,double T,double paso,double  t_max, double c){
	double v[3],t,dt,t_pre,*x, y_pre,x_pre, z_pre;
  double a;
	int i=0;
	FILE *pt;
  pt=fopen("caos2c_1.dat","w");	
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//2
   dt=paso;	//0.001
	t= T;//0
	double *vecx=malloc((t_max/dt)*sizeof(double));
	double *vecy=malloc((t_max/dt)*sizeof(double));
	double *vecz=malloc((t_max/dt)*sizeof(double));
	double *vect=malloc((t_max/dt)*sizeof(double));
	for(int i=0;i<t_max/dt;i++){
		t=i*dt;
		rk4(ecuaciones2a,v,3,t,dt);
		vecx[i]=v[0];	
		vecy[i]=v[1];
		vecz[i]=v[2];
		vect[i]=t;}
		i=0;


	while(i<t_max/dt){
		i=i+1;
	 z_pre=vecx[i];
    y_pre=vecy[i];
	 x_pre=vecz[i];
	 t_pre=vect[i];
		for(int j =i;j<t_max/dt;j++){
			if( vect[j]-t_pre>0.1 & sqrt((vecx[j]-x_pre)*(vecx[j]-x_pre)+(vecy[j]-y_pre)*(vecy[j]-y_pre)+(vecz[j]-z_pre)*(vecz[j]-z_pre))<1){
				fprintf(pt,"%lg\t%lg\t%lg\n%lg\t%lg\t%lg\n",vecx[j],vecy[j],vecz[j],x_pre,y_pre,z_pre);
				printf("\n t=%lg tpre=%lg",vect[j],t_pre);}			
		}
	}	
   fprintf(pt,"\n");
	fclose(pt);
	free(vecx);
	free(vecy);
	free(vecz);	
	free(vect);
	return(0);
}

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
