#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include "funkcje.h"

#include "imnmath.hpp"
			
using namespace std;
int main() {
    
    double pi = 4.*atan(1.0);    
    int i,j,k,l;	
    int it,itmax;
    int nx,ny,n,ntot,norma;	
    double dt_cfl;
    double dt,dx,xp,yp,x0,y0,rmax;		
    double x,y,r,r2;
    double c;
    double *zeta,*psi,*cos_psi;
    double *urn,*ugn,*ubn; //red green blue
    double *urn1,*ugn1,*ubn1; //red green blue
    double vmax,vv;
    double *a,*b,*vecb;
    double *vxn,*vyn,*vxn1,*vyn1;
    int *jka,*jkb,*iwa,*iwb;
    
    double *ff,*sk1,*sk2,*sk3,*sk4;
    double *utemp1,*utemp2;
    
    double dkx,dky,skx,sky,sk;
    int ir0,jr0,ig0,jg0,ib0,jb0;
    double sigma,div;
    FILE *fp;

//==== parametry symulacji============================    
    dt=0.01;
    dx=0.1;
    nx=150;
    ny=150;
    n=nx*ny;
    itmax=10000;
    
    xp=dx*nx/2.;
    yp=dx*ny/2.;
    rmax=xp/10.;
    
    ntot=5*n-2*nx-2*ny;
//===============================
    zeta=new double[n];
    psi=new double[n];
    cos_psi=new double[n];
    
    urn=new double[n];
    ugn=new double[n];
    //ubn=new double[n];
    
    urn1=new double[n];
    ugn1=new double[n];
    //ubn1=new double[n];

    vecb=new double[n];
    
    vxn=new double[n];
    vyn=new double[n];
    vxn1=new double[n];
    vyn1=new double[n];
    
    
    iwa=new int[n+1];
	jka=new int[ntot];
	a=new double[ntot];
	
	iwb=new int[n+1];
	jkb=new int[ntot];
	b=new double[ntot];
	vecb=new double[ntot];
	
	
	ff=new double[ntot];
	sk1=new double[ntot];
	sk2=new double[ntot];
	sk3=new double[ntot];
	sk4=new double[ntot];
	utemp1=new double[ntot];
	utemp2=new double[ntot];
    
    printf("ntot= %d \n",ntot);
    
	//===============================    
	//pakiet startowy     
    fp=fopen("u0.dat","w");
    ir0=100;
    jr0=100;
    sigma=dx*5;
    c=0.;
    for(i=0;i<nx;i++)
	{
      	for(j=0;j<ny;j++)
		{
	  		l=j+i*ny;
	  		x=dx*(i-ir0);
	  		y=dx*(j-jr0);
	  		r2=x*x+y*y;
			urn[l] = exp(-r2/2./pow(sigma,2))/2./pi/sigma/sigma;
			//urn[l] = exp(-r2/2./pow(sigma,2))/((2./pi)*sigma/sigma);
			c=c+urn[l]*dx*dx;
			fprintf(fp,"%d  %d  %f  \n",i,j,urn[l]);
		}
		fprintf(fp,"\n");
	}
    fclose(fp); 
//################################################################################################ 
    for(it=1;it<=itmax;it++)
	{
		//schemat calkujacy: RK4		    
		//procedura do wyznaczania predkosci: zalezy od aktualnego polozenia srodka wiru (x0,y0)	    
		x0=xp;
		y0=xp;
		licz_predkosci(nx, ny, dx, x0, y0, rmax, zeta, psi, cos_psi, vxn1, vyn1, &dt_cfl);
		dt=dt_cfl/4;
		
		//k1
		for(l=0;l<n;l++)
			utemp1[l]=urn[l];

		licz_krok_fct(dt,dx,n,nx,ny,vxn1,vyn1,utemp1,utemp2);
		  for(l=0;l<n;l++) 
			  sk1[l]=(utemp2[l]-utemp1[l])/dt;
		//k2
			for(l=0;l<n;l++)
				utemp1[l]=urn[l]+dt/2*sk1[l];
		  //UWAGA:trzeba przeliczyc predkosci dla tn+dt/2 i wstawic!!!!!!!!!!
		licz_krok_fct(dt,dx,n,nx,ny,vxn1,vyn1,utemp1,utemp2);
			for(l=0;l<n;l++) 
				sk2[l]=(utemp2[l]-utemp1[l])/dt;
		//k3
			for(l=0;l<n;l++)
				utemp1[l]=urn[l]+dt/2*sk2[l];
		  //uwaga: uzywamy prekosci dla tn+dt/2 i wstawic!!!!!!!!!!
		licz_krok_fct(dt,dx,n,nx,ny,vxn1,vyn1,utemp1,utemp2);
			for(l=0;l<n;l++) 
				sk3[l]=(utemp2[l]-utemp1[l])/dt;
		
		//k4
			for(l=0;l<n;l++)
				utemp1[l]=urn[l]+dt*sk3[l];
		//UWAGA:trzeba przeliczyc predkosci dla tn+dt i wstawic!!!!!!!!!!
		licz_krok_fct(dt,dx,n,nx,ny,vxn1,vyn1,utemp1,utemp2);
		for(l=0;l<n;l++) 
			sk4[l]=(utemp2[l]-utemp1[l])/dt;
		//skladamy rozwiazania dla RK4
		for(l=0;l<n;l++)
			urn1[l]=urn[l]+dt/6*(sk1[l]+2*sk2[l]+2*sk3[l]+sk4[l]); 	    
		c=0.;
		for(l=0;l<n;l++)
		{
			urn[l]=urn1[l]; 
			vxn[l]=vxn1[l]; 
			vyn[l]=vyn1[l]; 
			c=c+urn[l]*dx*dx;
		}
		
		if(it%20 == 0)
			printf("it= %d  dt=  %g  c= %g\n",it,dt_cfl,c);
		
		if(it%100 == 0)
		{
			fp=fopen("psi.dat","w");
			for(i=0;i<nx;i++)
			{
				for(j=0;j<ny;j++)
				{
					l=j+i*ny;
					div=0.;
					if( (i>0) && (i< (nx-1)))
						div=div+(vxn1[l+ny]-vxn1[l-ny])/2/dx;
					if( (j>0) && (j< (ny-1)) )
						div=div+(vyn1[l+1]-vyn1[l-1])/2/dx;
					fprintf(fp,"%d  %d  %f  %f   %f   %f   %f    %f\n",i,j,zeta[l],psi[l],vxn1[l],vyn1[l],div,urn[l]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp); 
			imnd::push_data2D(urn, nx, ny);
		}
		//=====================================================================================================	        
		/*  
		 //nie dziala   
		//......liczymy Bi-CG-stab...........................
			macierz_A_B(n,nx,ny,dt,dx,vxn,vyn,vxn1,vyn1,a,jka,iwa,b,jkb,iwb);
		//......B_n*u_n......................................	
			licz_axb(n,b,jkb,iwb,urn,vecb);
			for(l=0;l<n;l++)urn1[i]=0;//urn[i];//przyblizenie startowe
			bicgstab(n,a,jka,iwa,urn1,vecb);
			c=0.;
			for(l=0;l<n;l++){
				urn[l]=urn1[l];
				c=c+urn[l]*dx*dx;
			}		
			printf("it= %d   cn= %f  \n",it,c);
		*/
	}
	imnd::write_data2D("density.dat", nx, ny, dx, dx);
	imnd::free_data2D(nx, ny);
	//it
     //trzeba wykasowac pozostale tablice!!!!!!!!!!!!!
    delete [] psi;
    delete [] zeta;
    delete [] vxn;
    delete [] vyn;
    delete [] vxn1;
    delete [] vyn1;
    delete [] vecb;

    return 0;
}
