#ifndef EJERCICIOS_H
#define EJERCICIOS_H
void rk4(void deri(int , double [], double [], double ), double h[], int n, double t, double dt);
void ecuaciones1a(int n, double v[], double dv[], double t);
void ecuaciones2a(int n, double v[], double dv[], double t);
int ej_1a(double v_x,double v_y,double v_z,double t,double  paso,double t_max);
int ej_2a(double v_x,double v_y,double v_z,double t,double  paso,double t_max, double c);
#endif
