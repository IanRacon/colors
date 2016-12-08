#pragma once
void licz_krok_fct(double dt,double dx, int n, int nx,int ny, double *vxn,  double *vyn, double *un, double *un1);

double rplus(int l, int ny,double dx, double *ax, double *ay,double *qmax, double*qtd);
double rminus(int l, int ny,double dx, double *ax, double *ay,double *qmin, double*qtd);

void licz_predkosci(int nx, int ny, double dx, double x0, double y0,double rmax,\
			  double *zeta, double *psi, double *cos_psi, double *vxn1,double *vyn1,double *dt_cfl);

void bicgstab(int n, double *b, int *jkb, int *iwb, double *xn, double *vecb);

void macierz_A_B(int n, int nx, int ny, double dt, double dx,\
	           double *vxn, double *vyn, double *vxn1, double *vyn1,\
		     double *a, int *jka, int *iwa, double *b, int *jkb, int *iwb);

void licz_axb(int n,double *b, int *jkb,int *iwb, double *un, double *vecb);
double skalar(int n,double *a,double *b);
