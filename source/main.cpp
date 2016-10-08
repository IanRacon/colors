#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "imnmath.hpp"
#include "easylogging++.h"
#include "INIReader.h"
#include "Configuration.h"

INITIALIZE_EASYLOGGINGPP

const int xsize = 201;
const int ysize = 61;
const double dx = 0.02, dy = 0.02;
const int OBSTACLE=1;
const int NONE=0;
const int SEQUENCES=100;

void fill(int x1, int y1, int x2, int y2, double **array, double val){
		for(int i=x1;i<x2;++i)
			for(int j=y1;j<y2;++j)
				array[i][j] = val;	
}
void printArray(double **array){
	for(int i=0;i<xsize;++i){
		for(int j=0;j<ysize;++j)
			cout <<	array[i][j]; 
		cout << endl;
	}
}
void fillValue(double **array, double val){
	for(int i=0;i<xsize;++i)
		for(int j=0;j<ysize;++j)
			array[i][j] = val; 
}
double **initArray(){
	double **array = new double*[xsize];
	for(int i=0;i<xsize;++i)
		array[i] = new double[ysize];
	fillValue(array, 0);
	return array;
}
void freeArray(double** array){
	for(int i=0;i<xsize;++i)	
		delete[] array[i];
	delete[] array;
}
//tutaj zadajemy gestosc plamy w danym punkcie, czyli "objetosc" w danym punkcie
double p0(double x, double y){
	double part1;
	double part2;
	//part1 = pow(x-0.2,2);
	part1 = pow(x-0.4,2);
	part2 = pow(y-0.6,2);
	return exp(-30*(part1+part2));
}
//do P_old wpisujemy poczatkowa gestosc plamy w danym punkcie
void initP_old(double **array){
	for(int i=1;i<xsize-1;++i)	
		for(int j=1;j<ysize-1;++j){
			array[i][j] = p0(i*dx, j*dy);
		}	
}
//inicjalizujemy plame w srednim punkcie czasowym juz na podstawie predkosci
void initP(double **p_old, double **p, double **u, double **v, double dt){
	for(int i=0;i<xsize;++i)	
		for(int j=0;j<ysize;++j){
			p[i][j] = p0(i*dx-u[i][j]*dt, j*dy-v[i][j]*dt);
		}	
}
void leapfrog(double **p_old, double **p, double **p_new, double **v, double **u, double dt){
	for(int i=0;i<xsize;++i){
		if(i==0){
			for(int j=1;j<ysize-1;++j){
				p_new[i][j] = p_old[i][j] - dt*(u[i][j]*(p[i+1][j]-p[xsize-1][j])/dx + v[i][j]*(p[i][j+1]-p[i][j-1])/dy);
			}	
		}else if(i==xsize-1){
			for(int j=1;j<ysize-1;++j){
				p_new[i][j] = p_old[i][j] - dt*(u[i][j]*(p[0][j]-p[i-1][j])/dx + v[i][j]*(p[i][j+1]-p[i][j-1])/dy);
			}	
		}else{
			for(int j=1;j<ysize-1;++j){
				p_new[i][j] = p_old[i][j] - dt*(u[i][j]*(p[i+1][j]-p[i-1][j])/dx + v[i][j]*(p[i][j+1]-p[i][j-1])/dy);
			}	
		}
	}	
}
void laxfriedrichs(double **flag, double **p, double **p_new, double **u, double **v, double dt){
	double part1;
	double part2;
	for(int i=1;i<xsize-1;++i)	
		for(int j=1;j<ysize-1;++j){
			if(i==0){
				if(flag[i][j]!=OBSTACLE){
					part1=p[i+1][j]+p[xsize-1][j]+p[i][j-1]+p[i][j+1];
					part2=dt*(u[i][j]*(p[i+1][j]-p[xsize-1][j])/(2*dx) + v[i][j]*(p[i][j+1]-p[i][j-1])/(2*dy));
					p_new[i][j] = part1/4.0 - part2; 
				}
			}else if(i==xsize-1){
				if(flag[i][j]!=OBSTACLE){
					part1=p[0][j]+p[i-1][j]+p[i][j-1]+p[i][j+1];
					part2=dt*(u[i][j]*(p[0][j]-p[i-1][j])/(2*dx) + v[i][j]*(p[i][j+1]-p[i][j-1])/(2*dy));
					p_new[i][j] = part1/4.0 - part2; 
				}
			}else{
				if(flag[i][j]!=OBSTACLE){
					part1=p[i+1][j]+p[i-1][j]+p[i][j-1]+p[i][j+1];
					part2=dt*(u[i][j]*(p[i+1][j]-p[i-1][j])/(2*dx) + v[i][j]*(p[i][j+1]-p[i][j-1])/(2*dy));
					p_new[i][j] = part1/4.0 - part2; 
				}
			}
		}
}
double u(double x, double y){
	double mi=2, ymin=0, ymax=1.2, Q=-10;
	return Q/(2*mi)*(y-ymin)*(y-ymax);
}
void initU(double **u_arr){
	for(int i=0;i<xsize;++i)	
		for(int j=0;j<ysize;++j)
			u_arr[i][j] = u(i*dx, j*dy);	
}
double max(double **u, double **v){
	double umax=u[0][0];
	double vmax=v[0][0];
	double temp;
	double max=u[0][0];
	for(int i=0;i<xsize;++i)	
		for(int j=0;j<ysize;++j){
			//umax=u[i][j]>umax?u[i][j]:umax;
			//vmax=v[i][j]>vmax?v[i][j]:vmax;
			temp=sqrt(pow(u[i][j],2) + pow(v[i][j],2));
			max=temp>max?temp:max;
		}
	return max;
}
double densityI(double **p){	
	double sum=0;
	double temp=0;
	for(int i=0;i<xsize;++i){
		for(int j=0;j<ysize;++j){
			temp+=p[i][j]*dy;	
		}
		sum+=temp*dx;
		temp=0;
	}
	return sum;
}
double centerX(double **p, double I){
	double sum=0;
	double temp=0;
	for(int i=0;i<xsize;++i)	
		for(int j=0;j<ysize;++j)
			sum+=p[i][j]*dx*dy*i*dx;
	return sum/I;
}
void zad1(){
	double **p_new, **p, **p_old;
	//predkosc w kierunku (poziomym)
	double **u;
	//predkosc w kierunku (pionowym)
	double **v;
	//tymczasowa tablica 2D
	double **temp;
	// krok czasowy delta time
	double dt;
	// czas
	double t;
	//czas maksymalny
	double tmax=15; 
	//wartosc calki
	double I; 
	//centrum 
	double center_x;
	//inicjalizujemy strumien do zapisu
	ofstream file;
	// licznik
	int counter=0;
	// inicjalizacja tablic p_new, na poczatku wypelniona zerami
	p_new = initArray();
	// inicjalizacja tablicy p
	p = initArray();
	// inicjalizacja tablicy p_old
	p_old = initArray();
	// inicjalizacja tablicy predkosci poziomej
	u = initArray();
	// inicjalizacja tablicy predkosci pionowej
	v = initArray();	
	
	// inicjalizacja predkosci poziomych w ksztalt poziomej paraboli
	initU(u);
	// wypelnienie predkosci pionowych zerami
	fillValue(v, 0);

	// dt rowna sie dx podzielone przez 4*maksymalna predkosc zlozonych predkosci poziomej i pionowej
	//mamy predkosc czastek i krok odleglosciowy - mozemy wyliczyc krok czasowy zeby wszystko sie zgadzalo
	dt = dx/(4*max(u, v));	
	
	// inicjalizcja tablicy zawierajaca plame, to tutaj ustalamy gdzie ta plama jest!
	initP_old(p_old);

	//inicjalizacja tablicy plamy w pierwszym kroku czasowym na podstawie predkosci
	initP(p_old, p, u, v, dt);

	//otwieramy plik do zapisu
	file.open("zad1.dat");

	//dopoki nie przekroczymy calego czasu
	while(t<tmax){
		//schemat leapfrog
		leapfrog(p_old, p, p_new, v, u, dt);
		//kopiujemy wskaznik na tablice p_old do temp
		temp = p_old;
		//p_old wskazuje teraz na p
		p_old = p;
		//p wskazuje teraz na p_new
		p = p_new;
		//p_new zawiera teraz p_old
		p_new = temp;
		//przechodzimy do kolejnego kroku czasowego
		t+=dt;
		//obliczamy calke
		I = densityI(p);
		//jakis punkt centralny
		center_x = centerX(p, I);
		//zapisujemy do pliku calke i punkt centralny
		file << t << " " << I << " " << center_x << endl;
		//zapisujemy do pliku wyglad plamy oleju
		if(counter++%int(tmax/dt/SEQUENCES)==0)
			imnd::push_data2D(p,xsize,ysize);
	}
	file.close();
	imnd::write_data2D("zad1_1.dat", xsize, ysize, dx, dy);
	imnd::free_data2D(xsize, ysize);

	freeArray(p_new);
	freeArray(p);
	freeArray(p_old);
	freeArray(u);
	freeArray(v);
}
void zad2(){
	double **p_new, **p, **p_old, **u, **v, **temp;
	double dt=0, t=0, tmax=15, I=0, center_x=0;
	ofstream file;
	p_new = initArray();
	p = initArray();
	p_old = initArray();
	u = initArray();
	v = initArray();	

	//wczytanie predkosci z pliku	
	//----------------------------------------
	ifstream fuv("predkosc.txt");
	string line;
	double x, y;
	int counter=0;
	if(fuv.is_open()){
		for(int i=0;i<xsize;++i){
			for(int j=0;j<ysize;++j){
				fuv >> x >> y >> u[i][j] >> v[i][j];
		}
		getline(fuv, line);
	}
		fuv.close();
	}else{
		cout << "error" << endl;
	}
	//----------------------------------------
	initP_old(p_old);
	initP(p_old, p, u, v, dt);
	
	//printArray(u);
	//printArray(v);
	dt = dx/(4*max(u, v));	
	file.open("zad2.dat");
	while(t<tmax){
		leapfrog(p_old, p, p_new, v, u, dt);
		temp = p_old;
		p_old = p; 	
		p = p_new;
		p_new = temp;
		//printArray(p);
		t+=dt;
		I = densityI(p);
		center_x = centerX(p, I);
		file << t << " " << I << " " << center_x << endl;
		if(counter++%int(tmax/dt/SEQUENCES)==0)
			imnd::push_data2D(p,xsize,ysize);
	}
	file.close();
	imnd::write_data2D("zad2_1.dat", xsize, ysize, dx, dy);
	imnd::free_data2D(xsize, ysize);

	freeArray(p_new);
	freeArray(p);
	freeArray(p_old);
	freeArray(u);
	freeArray(v);
}
void zad3(){
	double **p_new, **flag, **p, **u, **v, **temp;
	double dt=0, t=0, tmax=15, I, center_x;
	ofstream file;
	p_new = initArray();
	p = initArray();
	u = initArray();
	v = initArray();	
	flag = initArray();

	//wczytanie predkosci z pliku	
	//----------------------------------------
	ifstream fuv("predkosc.txt");
	string line;
	double x, y;
	if(fuv.is_open()){
		for(int i=0;i<xsize;++i){
			for(int j=0;j<ysize;++j){
				fuv >> x >> y >> u[i][j] >> v[i][j];
		}
		getline(fuv, line);
	}
		fuv.close();
	}else{
		cout << "error" << endl;
	}
	//----------------------------------------

	dt = dx/(4*max(u, v));	

	initP_old(p);	
	dt=dt/2.0;
	fill(50, 40, 56, 61, flag, OBSTACLE);
	
	int counter = 0;
	file.open("zad3.dat");
	while(t<tmax){
		laxfriedrichs(flag, p, p_new, u, v, dt);
		temp = p;
		p = p_new; 	
		p_new = temp;
		t+=dt;
		I = densityI(p);
		center_x = centerX(p, I);
		file << t << " " << I << " " << center_x << endl;
		if(counter++%int(tmax/dt/SEQUENCES)==0)
			imnd::push_data2D(p,xsize,ysize);
	}
	file.close();
	
	imnd::write_data2D("zad3_1.dat", xsize, ysize, dx, dy);
	imnd::free_data2D(xsize, ysize);
	
	freeArray(flag);
	freeArray(p_new);
	freeArray(p);
	freeArray(u);
	freeArray(v);
}

int main(void){

	// LOG(INFO) << reader.GetInteger("test", "test", 0);
	zad1();
	zad2();
	zad3();
	return 0;
}
