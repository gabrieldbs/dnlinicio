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
//	FILE *pt;
//	pt=fopen("caos1t.dat","a");
	//FILE *ptt;
	//ptt=fopen("caos1x.dat","a");
	//fprintf(pt,"Barrido temporal /n");
//	fprintf(ptt,"Barrido x");
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
			//fprintf(pt,"%lg ,",t);
			//fprintf(ptt,"%lg ,",v[0]);
		}
		  t+=dt;
	}
	//fprintf(pt,"\n");
//	fprintf(ptt,"\n");

	fclose(ptr);
	//fclose(pt);
	//fclose(ptt);
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
	ptr=fopen("caos2d_2.dat","w");
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
   ptr=fopen("caos2b_0.dat","w");
	double v[3],t,dt,t_pre,*x, y_pre,x_pre;
   double a,theta;
	FILE *pt;
   pt=fopen("caos2b_1_6.dat","w");
	//FILE *ptt;	
  //  ptt=fopen("caos2b_3.dat","w");
	FILE *ptt;
	ptt=fopen("caos2b_2_6.dat","w");
	FILE *pttr;
	pttr=fopen("caos2b_3_6.dat","w");
	FILE *pttrr;
	pttrr=fopen("caos2b_4_6.dat","w");
	FILE *pttt;
	pttt=fopen("caos2b_5_6.dat","w");
	FILE *ptm;
	ptm=fopen("caos2b_1_2.dat","w");

	FILE *ptmm;
	ptmm=fopen("caos2b_1.dat","w");
	
  for (int j=0;j<6;j++){
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//2
   dt=paso;	//0.001
	t= T;//0
	t_pre=0;
	theta=j*3.1415/3;
	double l=tan(theta);
	printf("theta=%lg,tang(theta)=%lg \n",theta,l);
	// t_ max = 300
  while(t<t_max){
		y_pre=v[1];
		x_pre=v[0];
		rk4(ecuaciones2a,v,3,t,dt);
    if(t>t_pre &j==0& y_pre<tan(theta)*x_pre && v[1]>tan(theta)*v[0]) {
			fprintf(ptr,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);
			}
    if(t>t_pre &j==1& y_pre<tan(-theta)*x_pre && v[1]>tan(-theta)*v[0]) {
			fprintf(pt,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);			
			}
    if(t>t_pre &j==2& y_pre<(tan(-theta))*x_pre && v[1]>(tan(-theta))*v[0]) {
			fprintf(ptt,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);			
			}
    if(t>t_pre &j==3& y_pre<(tan(theta))*x_pre && v[1]>(tan(theta))*v[0]) {
			fprintf(pttr,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);			
			}
    if(t>t_pre &j==4& y_pre<tan(theta)*x_pre && v[1]>tan(theta)*v[0]) {
			fprintf(pttrr,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);			
			}
    if(t>t_pre &j==5& y_pre<tan(theta)*x_pre && v[1]>tan(theta)*v[0]) {
			fprintf(pttt,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);	
			}
		if(t>t_pre  & v[1] >0& v[0]*x_pre <0) {
			fprintf(ptm,"%lg\t%lg\n",sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);
			}
		if(t>t_pre  & v[0] >0& v[1]*y_pre <0) {
			fprintf(ptmm,"%lg\t%lg\n", sqrt(v[0]*v[0]+v[1]*v[1]),v[2]);
			}

/*
	if(t>t_pre  & v[1] >0& v[0]*x_pre <0) {
			fprintf(pt,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
*/

	/*if(t>t_pre  & v[0] >0& sqrt((v[0]-v[1])*(v[0]-v[1])) <0.01) {
			fprintf(ptt,"%lg\t%lg\t%lg\n",v[0],v[1],v[2]);
			}
    */

		t+=dt;
	}
  }
   fprintf(ptm,"\n");
	fclose(ptm);   

   fprintf(ptmm,"\n");
	fclose(ptmm);   
   fprintf(ptr,"\n");
	fclose(ptr);   
	fprintf(ptt,"\n");
	fclose(ptt);
   fprintf(pt,"\n");
	fclose(pt);
   fprintf(pttt,"\n");
	fclose(pttt);
   fprintf(pttr,"\n");
	fclose(pttr);
   fprintf(pttrr,"\n");
	fclose(pttrr);
	return(0);
}



int ej2_c(double v_x,double v_y,double v_z,double T,double paso,double  t_max, double c){
	double v[3],t,dt,t_pre,*x, y_pre,x_pre, z_pre;
  double a,DT;
	FILE *pt;
  pt=fopen("caos2e_1.dat","w");
	
	v[0]=v_x; 	//-1
	v[1]=v_y;	// -5
	v[2]=v_z;		//2
   dt=paso;	//0.001
	t= T;//0
	t_pre=10;
	z_pre=v_z;
 	y_pre=v_y;
	x_pre=v_x;
	
	// t_ max = 300
  while(t<t_max){
		rk4(ecuaciones1a,v,3,t,dt);
		
		if(t<2 & 3<t){
			z_pre=v[2];
   	 	y_pre=v[1];
	 		x_pre=v[0];
			t_pre = t;

		}
    if(t>t_pre+1 & sqrt((v[0]-x_pre)*(v[0]-x_pre)+(v[1]-y_pre)*(v[1]-y_pre)+(v[2]-z_pre)*(v[2]-z_pre))<0.4) { //0.00295
			DT=t-t_pre;
			fprintf(pt,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",x_pre,y_pre,z_pre,v[0],v[1],v[2],t_pre,t,DT);
			t_pre = t;
			z_pre=v[2];
   	 	y_pre=v[1];
	 		x_pre=v[0];
			//printf("l"); 
			}
	 	t=t+dt;
		//t_pre= t + 0.01;	
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
  pt=fopen("caos2c_1_2.dat","w");	
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
		vect[i]=t;
  //printf("%lg \n",t);
		}
		i=0;


	while(i<t_max/dt){
		i=i+1;
	 z_pre=vecx[i];
    y_pre=vecy[i];
	 x_pre=vecz[i];
	 t_pre=vect[i];
		for(int j =i;j<t_max/dt;j++){
			if( vect[j]-t_pre>0.03 && sqrt((vecx[j]-x_pre)*(vecx[j]-x_pre)+(vecy[j]-y_pre)*(vecy[j]-y_pre)+(vecz[j]-z_pre)*(vecz[j]-z_pre))<0.01){
				for (int k=i;k<j;k++){
					fprintf(pt,"%lg\t%lg\t%lg\n",vecx[k],vecy[k],vecz[k]);
				printf("\n t_pre=%lg ",t_pre);}
				break;
			}			
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
