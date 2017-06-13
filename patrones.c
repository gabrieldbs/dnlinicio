#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
//------------------------------------------------------------------------
//Defino constantes

#define N 200
#define Nh N
#define Na N
#define Ntot Na + Nh
double gap=1;
//-----------------------------------------------------------------------------
//defino la funcion para llamar al gnuplot
void gnuplot(const char *gnucommand)
{
  char syscommand[4000];
  sprintf(syscommand, "echo \"%s\" | gnuplot -persist", gnucommand);
  system(syscommand);
}
//-----------------------------------------------------------------------------
//Nombrador de files
void nombrar_par1_par2(char *name_space,char* nom,int par1, int par2){
  sprintf(name_space,"temp/salida-%s.i%d.j%d.dat",nom,par1,par2);
return;
}
//------------------------------------------------------------------------
//Parámetros
struct Par {
	double tita;
	double sigma;
	double mu; 
	double s;
	double ba;
	double rb;
	double ra;
	double bb;
	double alfa;
	double sigmaa;
	double sigmab;
	} par;
//------------------------------------------------------------------------
//Runge Kutta 4 pasos
void rk4(void deri(int , double [], double [], double, double [], double [] ), \
double h[], int n, double t, double dt, double par1[], double par2[]){
#define naux 2000 

int i;
double k1[naux],k2[naux],k3[naux],k4[naux],h0[naux];
double dt2, dt6;

dt2=dt/2.;
dt6=dt/6.;

for (i = 0 ; i<n; i++)
	h0[i] = h[i];

deri(n,h0,k1,t, par1, par2);
for (i =0 ; i<n ; i++)
	h0[i]=h[i]+dt2*k1[i];

deri(n,h0,k2,t+dt2, par1, par2);
for (i =0 ; i<n ; i++)
	h0[i]=h[i]+dt2*k2[i];

deri(n,h0,k3,t+dt2, par1, par2);
for (i =0 ; i<n ; i++)
	h0[i]=h[i]+dt*k3[i];

deri(n,h0,k4,t+dt, par1, par2);
for (i = 0 ; i<n ; i++)
	{h0[i]=h[i]+dt*k4[i];};

for (i =0; i<n ; i++)
	h[i]=h[i]+dt6*(2.*(k2[i]+k3[i])+k1[i]+k4[i]);

return;
}
//------------------------------------------------------------------------
//Ecuaciones a integrar (N*2), v[0 a N-1] población a, v[N a 2N-1] población b
void takens(int n, double v[], double dv[], double t, double par1[],double par2[]){
  int ii,jj;
  double laplace1,laplace2, laplace0,laplaceN,laplaceN1,laplace2N; //laplacianos discretos
    
    //---------------------------------------------------------------------    
    //Condiciones de contorno, derivada nula
           
    laplace0 = 2*(v[1]-v[0])/(gap*gap);
    laplaceN1 = 2*(v[N-2]-v[N-1])/(gap*gap);
    laplaceN =	2*(v[N+1]-v[N])/(gap*gap);
    laplace2N =	2*(v[N-2+N]-v[N+N-1])/(gap*gap);												
                                                               
    dv[0] = par.s*(v[0]*v[0]/v[0+N]+par.ba)-par.ra*v[0]+par.sigmaa*laplace0;
    dv[N-1] = par.s*(v[N-1]*v[N-1]/v[N-1+N]+par.ba)-par.ra*v[N-1]+par.sigmaa*laplaceN1;
    
    
    dv[N] =par.s*v[0]*v[0]-par.rb*v[N]+par.sigmab*laplaceN+par.bb;
    dv[N+N-1] =par.s*v[N-1]*v[N-1]-par.rb*v[N+N-1]+par.sigmab*laplace2N+par.bb;
    //---------------------------------------------------------------------
                                                                                     
for(ii=1;ii<N-1;ii++){
	laplace1 = (v[ii-1]+v[ii+1]-2*v[ii])/(gap*gap);
	laplace2 = (v[ii-1+N]+v[ii+1+N]-2*v[ii+N])/(gap*gap);
  	
  	dv[ii]   =  par.s*(v[ii]*v[ii]/v[ii+N]+par.ba)-par.ra*v[ii]+par.sigmaa*laplace1;
  	 	
  	dv[ii+N] =  par.s*v[ii]*v[ii]-par.rb*v[ii+N]+par.sigmab*laplace2+par.bb;
  	
  	
  	//-v[ii]*v[ii]/(1+v[ii+N])-v[ii]+par.sigma*laplace1;
  	//par.mu*(v[ii]*v[ii]-v[ii+N])+laplace2; 
 	}                                                                                          
return;
}
//------------------------------------------------------------------------
//Programa principal
int main(){

  FILE *ptr1,*ptr2,*ptr3,*ptr4;
  
  int i,j,ii,jj,ni,nj,count, gnuplotcount,k;
  
  double v[Ntot],t,dt,par1[N],par2[N];
  
  double vari, varj, parj0,pari0,ventanash, ventanasv,concent0, t_trans;
  
  char nombre_archivo1[30], nombre_archivo2[30], nombre_archivo3[30], nombre_archivo4[30];
 
//------------------------------------------------------------------------ 
//carpeta_archivos temporarios
mode_t process_mask = umask(0);
mkdir("./temp",S_IRWXU | S_IRWXG | S_IRWXO);
umask(process_mask); 
//------------------------------------------------------------------------
//Parámetros
par.sigma = 10;
par. mu = 1;

par.s = 1;
par.ba = 3;
par.bb = 0.1;
par.ra = 1;
par.rb = 8;
par.bb= 1;
par.sigmaa = 0.0001;
par.sigmab=50;

dt = 0.01;  
concent0=1; 			//inicio de las concentraciones


char *nombreaux="";
char nombre[8000];

//---------------------------------------------------
		//------------------------------------------------------------------------
		//------------------------------------------------------------------------
		//Condiciones iniciales de todas 
		for(i=0;i<N;i++){
			v[i]=concent0;
			v[i+N]=concent0*0.1;
			}
		//------------------------------------------------------------------------
		// La perturbación, o función inicial del activador	
		v[88] = 1;
		v[89] = 1;
		v[90] = 1.001;
		v[91] = 1;
		v[92] =1;
		//------------------------------------------------------------------------
		//------------------------------------------------------------------------	
	//Abro files 

ii=0;
jj=0;	
	
	nombreaux="concA";	
	nombrar_par1_par2(nombre_archivo1,nombreaux,ii,jj);
	ptr1=fopen(nombre_archivo1,"w");
	
	nombreaux="concB";	
	nombrar_par1_par2(nombre_archivo1,nombreaux,ii,jj);
	ptr2=fopen(nombre_archivo1,"w");
	
		 
    	//----------------------------------------------------------------------------
    	//Comienza integración temporal
 		count=0;									//reinicio tiempos
 		while(count<25000){
 			t=count*dt;    //calculo el tiempo
      		
      		if(count % 500 ==0){//Guardo data cada 500 pasos temporales
		    			
			for(k=0;k<N;k++){
				fprintf(ptr1,"%d \t %lg \t %lg \n",k,t,v[k]);
				fprintf(ptr2,"%d \t %lg \t %lg \n",k,t,v[k+N]);
			}
			fprintf(ptr1,"\n");
			fprintf(ptr2,"\n");
			
			}
      		
      		rk4(takens, v, Ntot, t, dt, par1, par2);
						
	  		//------------------------------------------------------------------------
		    					
			count++;							
		}
			
			//Termino integración
			//------------------------------------------------------------------------
			fclose(ptr1);
			fclose(ptr2);
	
	
//Para el gnuplot ---------------------
					 	
                    
                     sprintf(nombre,"load './contorno.txt'\n");
					 //---------------------------------------

					//Imprimo en un multiplot
					gnuplot(nombre);
					printf("%s",nombre);   //imprimo en pantalla el comando para gnuplot

return(0);
}
