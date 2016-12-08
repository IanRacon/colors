#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "funkcje.h"


//======== Liczymy jeden krok dla FCT =====================================================================	   
// uwaga: wartosci un1 sa wyznaczane do odleglosci 2dx od brzegu - na brzegu i w jego poblizu un=un1=0

void licz_krok_fct(double dt,double dx, int n, int nx,int ny, double *vxn,  double *vyn, double *un, double *un1){
	
	double q1,q2,q3,q4,q5,c;
	double vxp,vyp;
	double *flx, *gly, *fhx, *ghy;
	double *qlx,*qly;
	double qhx,qhy;
	
	double *ax,*ay;
	double *cx,*cy;
	double *qmin,*qmax;
	double *qtd;
	
	int i,j,l;
	
	flx=new double [n];
	gly=new double [n];
	fhx=new double [n];
	ghy=new double [n];
	
	qlx=new double [n]; //q^L_(i+1/2,j)
	qly=new double [n]; //q^L_(i,j+1/2)
		
	ax=new double [n]; //a_(i+1/2,j)
	ay=new double [n]; //a_(i,j+1/2)
	
	cx=new double [n]; //c_(i+1/2,j)
	cy=new double [n]; //c_(i,j+1/2)
	
	
	qtd=new double [n]; //q_td
	qmin=new double [n]; //q_min
	qmax=new double [n]; //q_max
	
	for(l=0;l<n;l++){
		flx[l]=0.;
		gly[l]=0.;
		fhx[l]=0.;
		ghy[l]=0.;
		qlx[l]=0.;
		qly[l]=0.;
		ax[l]=0.;
		ay[l]=0.;
		qtd[l]=0.;
		cx[l]=0.;
		cy[l]=0.;
		qmin[l]=0.;
		qmax[l]=0.;
	}
//lower order fluxes

	for(i=0;i<(nx-1);i++){
		for(j=0;j<(ny-1);j++){
		l=j+i*ny;
		
		vxp=(vxn[l]+vxn[l+ny])/2.;
		if(vxp>0){
			qlx[l]=un[l];	
		}else{
			qlx[l]=un[l+ny];	
		}
		
		vyp=(vyn[l]+vyn[l+1])/2.;
		if(vyp>0){
			qly[l]=un[l];	
		}else{
			qly[l]=un[l+1];	
		}
		
		flx[l]=qlx[l]*vxp*dx*dt;
		gly[l]=qly[l]*vyp*dx*dt;
		
		}	
	}
	
//higher order fluxes

	for(i=1;i<(nx-2);i++){
		for(j=1;j<(ny-2);j++){
		l=j+i*ny;
		
		qhx=7./12.*(un[l+ny]+un[l])-1./12.*(un[l+2*ny]+un[l-ny]);
		qhy=7./12.*(un[l+1]+un[l])-1./12.*(un[l+2]+un[l-1]);
		
		vxp=(vxn[l]+vxn[l+ny])/2.;
		vyp=(vyn[l]+vyn[l+1])/2.;
		
		fhx[l]=qhx*vxp*dx*dt;
		ghy[l]=qhy*vyp*dx*dt;
		
		}	
	}	
	
//antidiffusive fluxes

	for(i=1;i<(nx-2);i++){
		for(j=1;j<(ny-2);j++){
		l=j+i*ny;
		ax[l]=fhx[l]-flx[l];
		ay[l]=ghy[l]-gly[l];
		}	
	}	
	
//low-order transported-diffusedsolution

	for(i=1;i<(nx-1);i++){
		for(j=1;j<(ny-1);j++){
		l=j+i*ny;	
		qtd[l]=un[l]-(flx[l]-flx[l-ny]+gly[l]-gly[l-1])/dx/dx;
		}	
	}	
	
//qmin,qmax

	for(i=1;i<(nx-1);i++){
		for(j=1;j<(ny-1);j++){
		l=j+i*ny;
		
		q1=fmax(un[l],qtd[l]);
		q2=fmax(un[l-ny],qtd[l-ny]);
		q3=fmax(un[l+ny],qtd[l+ny]);
		q4=fmax(un[l-1],qtd[l-1]);
		q5=fmax(un[l+1],qtd[l+1]);
		qmax[l]=fmax(q1,fmax(q2,fmax(q3,fmax(q4,q5))));
	
		q1=fmin(un[l],qtd[l]);
		q2=fmin(un[l-ny],qtd[l-ny]);
		q3=fmin(un[l+ny],qtd[l+ny]);
		q4=fmin(un[l-1],qtd[l-1]);
		q5=fmin(un[l+1],qtd[l+1]);
		qmin[l]=fmin(q1,fmin(q2,fmin(q3,fmin(q4,q5))));
	
		}	
	}		
	
	
//Cx,Cy

	for(i=1;i<(nx-2);i++){
		for(j=1;j<(ny-2);j++){
		l=j+i*ny;
		
		if(ax[l]>0){		
			cx[l]=fmin(rplus(l+ny,ny,dx,ax,ay,qmax,qtd),rminus(l,ny,dx,ax,ay,qmin,qtd));	
		}else{
			cx[l]=fmin(rplus(l,ny,dx,ax,ay,qmax,qtd),rminus(l+ny,ny,dx,ax,ay,qmin,qtd));
		}
		
		if(ay[l]>0){
			cy[l]=fmin(rplus(l+1,ny,dx,ax,ay,qmax,qtd),rminus(l,ny,dx,ax,ay,qmin,qtd));	
		}else{
			cy[l]=fmin(rplus(l,ny,dx,ax,ay,qmax,qtd),rminus(l+1,ny,dx,ax,ay,qmin,qtd));
		}
	   }	
	}	
	
//modyfikacja: Ax,Ay
	for(i=1;i<(nx-2);i++){
		for(j=1;j<(ny-2);j++){
		l=j+i*ny;
	
		c=ax[l]*(qtd[l+ny]-qtd[l]);
		if(c<0)ax[l]=0.;
		
		c=ay[l]*(qtd[l+1]-qtd[l]);
		if(c<0)ay[l]=0.;
		
		
	   }	
	}	
	
	
	
	
	
//Acx,Acy
	for(i=1;i<(nx-2);i++){
		for(j=1;j<(ny-2);j++){
		l=j+i*ny;
	
		if(cx[l]>0 && cx[l]<=1.)ax[l]=ax[l]*cx[l];
		if(cy[l]>0 && cy[l]<=1.)ay[l]=ay[l]*cy[l];
	   }	
	}	
	
	
	
//un1, prawa strona rowania rozniczkowego (ff)
	for(i=1;i<(nx-1);i++){
		for(j=1;j<(ny-1);j++){
		l=j+i*ny;
		un1[l]=qtd[l]-(ax[l]-ax[l-ny]+ay[l]-ay[l-1])/dx/dx;
		}	
	}	
	
	
	
	delete [] flx;
	delete [] gly;
	delete [] fhx;
	delete [] ghy;
	delete [] qlx;
	delete [] qly;
	
	delete [] ax;
	delete [] ay;
	delete [] cx;
	delete [] cy;
	delete [] qmin;
	delete [] qmax;
	delete [] qtd;
	
	return;
}

double rplus(int l, int ny,double dx, double *ax, double *ay,double *qmax, double*qtd){
			double pplus,qplus,rp;
			pplus=fmax(ax[l-ny],0)-fmin(ax[l],0)+fmax(ay[l-1],0)-fmin(ay[l],0);
			qplus=(qmax[l]-qtd[l])*dx*dx;
			if(pplus>0){
				rp=fmin(1.,qplus/pplus);
			}else{
				rp=0.;
			}
			return rp;
}
			
double rminus(int l, int ny,double dx, double *ax, double *ay,double *qmin, double*qtd){
	double pminus,qminus,rm;
	pminus=fmax(ax[l],0)-fmin(ax[l-ny],0)+fmax(ay[l],0)-fmin(ay[l-1],0);
	qminus=(qtd[l]-qmin[l])*dx*dx;
	if(pminus>0){
		rm=fmin(1.,qminus/pminus);
	}else{
		rm=0.;
	}
	return rm;
}

//======== Liczymy: zeta -> psi -> predkosci ================================================	    
void licz_predkosci(int nx, int ny, double dx, double x0, double y0,double rmax,\
			  double *zeta, double *psi, double *cos_psi, double *vxn1,double *vyn1, double *dt_cfl){    
	double y,x,dkx,dky,pi,skx,sky,sk,norma,r,vmax,vv;
	double gamma;
	fftw_plan plan_1,plan_2;	    
	pi=4.*atan(1.0);

//========  zeta_n ================================================
	for(int i=0;i<nx;i++)
	{
      	for(int j=0;j<ny;j++)
		{
	  		int l=j+i*ny;
	  		x=dx*i-x0;
	  		y=dx*j-y0;
	  		r=sqrt(x*x+y*y);
			//zeta[l]=exp(-pow(r/rmax,2));
			
			gamma=1.1; //zmiana gammy pozwala wyzerowac predkosci w poblizu scianek	
			zeta[l]= -(gamma*pow(r/rmax,2)-1)*exp(-pow(r/rmax,2));		
		}
	}	    
//================= liczymy poissona: nabla^2*psi=zeta ================================================ 
//================= rozwiazujemy rownanie uzywajac DCT ===============================================
	
	plan_1 = fftw_plan_r2r_2d(nx, ny, zeta, cos_psi, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(plan_1);
	
	dkx=2.*pi/nx;
	dky=2.*pi/ny;
	norma=nx*ny;
	for(int i=0;i<nx;i++)
	{
      	for(int j=0;j<ny;j++)
		{
	  		int l=j+i*ny;
		 	skx=dkx*i;
			sky=dky*j;
			sk=skx*skx+sky*sky;
			if(i==0 && j==0){
				cos_psi[l]=0.;
			}else{
				cos_psi[l]=cos_psi[l]/sk/norma;
			}	
		}
	}
	plan_2 = fftw_plan_r2r_2d(nx,ny,cos_psi,psi, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
      fftw_execute(plan_2);
   //nowe predkosci: vxn1, vyn1
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
	  		int l=j+i*ny;
			vxn1[l]=0.;
			vyn1[l]=0.;
	  		if( (j!=0) && (j!=(ny-1)) )vxn1[l]=(psi[l+1]-psi[l-1])/2./dx;
			if((i!=0) && (i!=(nx-1)) )vyn1[l]=-(psi[l+ny]-psi[l-ny])/2./dx;
		}
	}
    
    /*
    	i=0; j=0; l=j+i*ny;
    	vxn1[l]=(vxn1[l+1]+vxn1[l+ny])/2;
	vyn1[l]=(vyn1[l+1]+vyn1[l+ny])/2;

	i=0; j=ny-1; l=j+i*ny;
	vxn1[l]=(vxn1[l-1]+vxn1[l+ny])/2;
	vyn1[l]=(vyn1[l-1]+vyn1[l+ny])/2;

	i=nx-1; j=ny-1; l=j+i*ny;
	vxn1[l]=(vxn1[l-1]+vxn1[l-ny])/2;
	vyn1[l]=(vyn1[l-1]+vyn1[l-ny])/2;
    
    	i=nx-1; j=0; l=j+i*ny;
	vxn1[l]=(vxn1[l+1]+vxn1[l-ny])/2;
	vyn1[l]=(vyn1[l+1]+vyn1[l-ny])/2;
    */
 
	fftw_destroy_plan(plan_1);       
	fftw_destroy_plan(plan_2);       
	fftw_cleanup();
    
    // maksymalny krok czasowy -> CFL
	vmax=0.;
	for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){
	  		int l=j+i*ny;		 	
			vv=sqrt(pow(vxn1[l],2)+pow(vyn1[l],2));
			if(vv>vmax)vmax=vv;
		}
	}
    
	*dt_cfl=dx/vmax;   
}

//=====================================================================================
//.....wypelniamy macierze A,B dla schematu CN.........................................
//=====================================================================================     
/*
void macierz_A_B(int n, int nx, int ny, double dt, double dx,
	           double *vxn, double *vyn, double *vxn1, double *vyn1,
			   */
//======================================================================================

void bicgstab(int n, double *b, int *jkb, int *iwb, double *xn, double *vecb){	
	int i,k;
	double rr,alfa,beta,rj0,omega;
	double *pn,*rn, *temp, *r0, *sn, *tems;
	pn=new double[n];
	rn=new double[n];
	temp=new double[n];
	r0=new double[n];
	sn=new double[n];
	tems=new double[n];
 	//inicjalizacja	
	licz_axb(n,b,jkb,iwb,xn,temp);
	for(i=0;i<n;i++){
		rn[i]=vecb[i]-temp[i];
		r0[i]=(double)i/n;
		pn[i]=r0[i];
	}
	//iteracja
	k=0;
	rr=skalar(n,rn,rn);
	
	while(rr > (1.E-16) && k<=2000){
		licz_axb(n,b,jkb,iwb,pn,temp);
		rj0=skalar(n,rn,r0);
		alfa=rj0/skalar(n,temp,r0);
		for(i=0;i<n;i++){
			sn[i]=rn[i]-alfa*temp[i];
		}
		licz_axb(n,b,jkb,iwb,sn,tems);
		omega=skalar(n,tems,sn)/skalar(n,tems,tems);
		for(i=0;i<n;i++){
	 		xn[i]=xn[i]+alfa*pn[i]+omega*sn[i];
 	 		rn[i]=sn[i]-omega*tems[i];
		}
		beta=skalar(n,rn,r0)/rj0*alfa/omega;
		for(i=0;i<n;i++){
	 		pn[i]=rn[i]+beta*(pn[i]-omega*temp[i]);
		}
		rr=skalar(n,rn,rn);	
		k=k+1;

//		if((k % 100) == 0)printf("%d  %g\n",k,rr);
	}
	if(rr > 1.E-6)printf("UWAGA:      rr=  %g \n",rr);	
	return;	
}

//=====================================================
	void licz_axb(int n,double *b, int *jkb,int *iwb, double *un, double *vecb){
		int l,l1,l2;
		int i,j,kolumna;
	
		for(l=0;l<n;l++){
			l1=iwb[l];
			l2=iwb[l+1]-1;
			vecb[l]=0.;
			for(i=l1;i<=l2;i++){
				kolumna=jkb[i];
				vecb[l]=vecb[l]+b[i]*un[kolumna];

			}
		}
		return;
	}
//=====================================================
	double skalar(int n,double *a,double *b){
		double c=0.;
		for(int l=0;l<n;l++)
			c=c+a[l]*b[l];
		return c;
	}
