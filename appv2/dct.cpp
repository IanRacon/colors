#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include "funkcje.h"
#include "utils.h"
#include "imnmath.hpp"

void rk4(double* densityMap, double* newDensityMap, double* vMapX, double* vMapY, double dt, double dx, int rows, int cols)
{
    int size = rows * cols;
    double* temp = new double[size];
    double* temp2 = new double[size];
    double* k1 = new double[size];
    double* k2 = new double[size];
    double* k3 = new double[size];
    double* k4 = new double[size];

    for (int l = 0; l < size; l++)
        temp[l] = densityMap[l];
    licz_krok_fct(dt, dx, size, cols, rows, vMapX, vMapY, temp, temp2);
    for (int l = 0; l < size; l++)
        k1[l] = (temp2[l] - temp[l]) / dt;

    for (int l = 0; l < size; l++)
        temp[l] = densityMap[l] + dt / 2 * k1[l];
    // UWAGA:trzeba przeliczyc predkosci dla tn+dt/2 i wstawic!!!!!!!!!!
    licz_krok_fct(dt, dx, size, cols, rows, vMapX, vMapY, temp, temp2);
    for (int l = 0; l < size; l++)
        k2[l] = (temp2[l] - temp[l]) / dt;

    for (int l = 0; l < size; l++)
        temp[l] = densityMap[l] + dt / 2 * k2[l];
    // uwaga: uzywamy prekosci dla tn+dt/2 i wstawic!!!!!!!!!!
    licz_krok_fct(dt, dx, size, cols, rows, vMapX, vMapY, temp, temp2);
    for (int l = 0; l < size; l++)
        k3[l] = (temp2[l] - temp[l]) / dt;

    for (int l = 0; l < size; l++)
        temp[l] = densityMap[l] + dt * k3[l];
    // UWAGA:trzeba przeliczyc predkosci dla tn+dt i wstawic!!!!!!!!!!
    licz_krok_fct(dt, dx, size, cols, rows, vMapX, vMapY, temp, temp2);
    for (int l = 0; l < size; l++)
        k4[l] = (temp2[l] - temp[l]) / dt;

    for (int l = 0; l < size; l++)
        newDensityMap[l] = densityMap[l] + dt / 6 * (k1[l] + 2 * k2[l] + 2 * k3[l] + k4[l]);

    delete[] temp;
    delete[] temp2;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
}
using namespace std;
int main()
{
    double pi = 4. * atan(1.0);
    int i, j, k, l;
    int it, itmax;
    int nx, ny, n, ntot, norma;
    double dt_cfl;
    double dt, dx, xp, yp, x0, y0, rmax;
    double x, y, r, r2;
    double c, cr, cg, cb;
    double *zeta, *psi, *cos_psi;
    double *urn, *ugn, *ubn; // red green blue
    double *urn1, *ugn1, *ubn1; // red green blue
    double vmax, vv;
    double *a, *b, *vecb;
    double *vxn, *vyn, *vxn1, *vyn1;
    int *jka, *jkb, *iwa, *iwb;
    double *ff, *sk1, *sk2, *sk3, *sk4;
    double *utemp1, *utemp2;
    double dkx, dky, skx, sky, sk;
    int ir0, jr0, ig0, jg0, ib0, jb0;
    double sigma, div;
    FILE* fp;

    //==== parametry symulacji============================
    dt = 0.02;
    dx = 0.2;
    nx = 250;
    ny = 250;
    n = nx * ny;
    itmax = 60000;

    xp = dx * nx / 2.;
    yp = dx * ny / 2.;
    rmax = xp / 10.;

    ntot = 5 * n - 2 * nx - 2 * ny;
    //===============================
    zeta = new double[n];
    psi = new double[n];
    cos_psi = new double[n];
    urn = new double[n];
    ugn = new double[n];
    ubn = new double[n];
    urn1 = new double[n];
    ugn1 = new double[n];
    ubn1 = new double[n];
    vecb = new double[n];
    vxn = new double[n];
    vyn = new double[n];
    vxn1 = new double[n];
    vyn1 = new double[n];
    iwa = new int[n + 1];
    jka = new int[ntot];
    a = new double[ntot];
    iwb = new int[n + 1];
    jkb = new int[ntot];
    b = new double[ntot];
    vecb = new double[ntot];
    ff = new double[ntot];
    sk1 = new double[ntot];
    sk2 = new double[ntot];
    sk3 = new double[ntot];
    sk4 = new double[ntot];
    utemp1 = new double[ntot];
    utemp2 = new double[ntot];

    printf("ntot= %d \n", ntot);

    //===============================
    // pakiet startowy
    fp = fopen("u0.dat", "w");
    sigma = dx * 12;
    c = 0.;
    cr = 0.;
    cg = 0.;
    cb = 0.;
    ir0 = 125;
    jr0 = 178;
    // r
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            l = j + i * ny;
            x = dx * (i - ir0);
            y = dx * (j - jr0);
            r2 = x * x + y * y;
            urn[l] = exp(-r2 / 2. / pow(sigma, 2)) / 2. / pi / sigma / sigma;
            // urn[l] = exp(-r2/2./pow(sigma,2))/((2./pi)*sigma/sigma);
            c = c + urn[l] * dx * dx;
            fprintf(fp, "%d  %d  %f  \n", i, j, urn[l]);
        }
        fprintf(fp, "\n");
    }

    sigma = dx * 7;
    ir0 = 151;
    jr0 = 109;
    // g
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            l = j + i * ny;
            x = dx * (i - ir0);
            y = dx * (j - jr0);
            r2 = x * x + y * y;
            ugn[l] = exp(-r2 / 2. / pow(sigma, 2)) / 2. / pi / sigma / sigma;
            // urn[l] = exp(-r2/2./pow(sigma,2))/((2./pi)*sigma/sigma);
            c = c + ugn[l] * dx * dx;
            fprintf(fp, "%d  %d  %f  \n", i, j, ugn[l]);
        }
        fprintf(fp, "\n");
    }
    sigma = dx * 12;
    ir0 = 98;
    jr0 = 109;
    // b
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            l = j + i * ny;
            x = dx * (i - ir0);
            y = dx * (j - jr0);
            r2 = x * x + y * y;
            ubn[l] = exp(-r2 / 2. / pow(sigma, 2)) / 2. / pi / sigma / sigma;
            // urn[l] = exp(-r2/2./pow(sigma,2))/((2./pi)*sigma/sigma);
            c = c + ubn[l] * dx * dx;
            fprintf(fp, "%d  %d  %f  \n", i, j, ubn[l]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    x0 = xp;
    y0 = xp;
    licz_predkosci(nx, ny, dx, x0, y0, rmax, zeta, psi, cos_psi, vxn1, vyn1,
        &dt_cfl);
    //################################################################################################
    int snapshotCounter = 0;
    for (it = 0; it <= itmax; it++) {
        // schemat calkujacy: RK4
        // procedura do wyznaczania predkosci: zalezy od aktualnego polozenia srodka
        // wiru (x0,y0)
        // x0 = xp;
        // y0 = xp;
        // licz_predkosci(nx, ny, dx, x0, y0, rmax, zeta, psi, cos_psi, vxn1, vyn1,
        // &dt_cfl);
        dt = dt_cfl / 4;
        // r
        ////////////////////////////////////////////////////////////////////////////
        // k1
        for (l = 0; l < n; l++)
            utemp1[l] = urn[l];

        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk1[l] = (utemp2[l] - utemp1[l]) / dt;
        // k2
        for (l = 0; l < n; l++)
            utemp1[l] = urn[l] + dt / 2 * sk1[l];
        // UWAGA:trzeba przeliczyc predkosci dla tn+dt/2 i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk2[l] = (utemp2[l] - utemp1[l]) / dt;
        // k3
        for (l = 0; l < n; l++)
            utemp1[l] = urn[l] + dt / 2 * sk2[l];
        // uwaga: uzywamy prekosci dla tn+dt/2 i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk3[l] = (utemp2[l] - utemp1[l]) / dt;

        // k4
        for (l = 0; l < n; l++)
            utemp1[l] = urn[l] + dt * sk3[l];
        // UWAGA:trzeba przeliczyc predkosci dla tn+dt i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk4[l] = (utemp2[l] - utemp1[l]) / dt;
        // skladamy rozwiazania dla RK4
        for (l = 0; l < n; l++)
            urn1[l] = urn[l] + dt / 6 * (sk1[l] + 2 * sk2[l] + 2 * sk3[l] + sk4[l]);
        ////////////////////////////////////////////////////////////////////////////

        // g
        ////////////////////////////////////////////////////////////////////////////
        // k1
        for (l = 0; l < n; l++)
            utemp1[l] = ugn[l];

        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk1[l] = (utemp2[l] - utemp1[l]) / dt;
        // k2
        for (l = 0; l < n; l++)
            utemp1[l] = ugn[l] + dt / 2 * sk1[l];
        // UWAGA:trzeba przeliczyc predkosci dla tn+dt/2 i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk2[l] = (utemp2[l] - utemp1[l]) / dt;
        // k3
        for (l = 0; l < n; l++)
            utemp1[l] = ugn[l] + dt / 2 * sk2[l];
        // uwaga: uzywamy prekosci dla tn+dt/2 i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk3[l] = (utemp2[l] - utemp1[l]) / dt;

        // k4
        for (l = 0; l < n; l++)
            utemp1[l] = ugn[l] + dt * sk3[l];
        // UWAGA:trzeba przeliczyc predkosci dla tn+dt i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk4[l] = (utemp2[l] - utemp1[l]) / dt;
        // skladamy rozwiazania dla RK4
        for (l = 0; l < n; l++)
            ugn1[l] = ugn[l] + dt / 6 * (sk1[l] + 2 * sk2[l] + 2 * sk3[l] + sk4[l]);
        ////////////////////////////////////////////////////////////////////////////

        // b
        ////////////////////////////////////////////////////////////////////////////
        // k1
        for (l = 0; l < n; l++)
            utemp1[l] = ubn[l];

        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk1[l] = (utemp2[l] - utemp1[l]) / dt;
        // k2
        for (l = 0; l < n; l++)
            utemp1[l] = ubn[l] + dt / 2 * sk1[l];
        // UWAGA:trzeba przeliczyc predkosci dla tn+dt/2 i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk2[l] = (utemp2[l] - utemp1[l]) / dt;
        // k3
        for (l = 0; l < n; l++)
            utemp1[l] = ubn[l] + dt / 2 * sk2[l];
        // uwaga: uzywamy prekosci dla tn+dt/2 i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk3[l] = (utemp2[l] - utemp1[l]) / dt;

        // k4
        for (l = 0; l < n; l++)
            utemp1[l] = ubn[l] + dt * sk3[l];
        // UWAGA:trzeba przeliczyc predkosci dla tn+dt i wstawic!!!!!!!!!!
        licz_krok_fct(dt, dx, n, nx, ny, vxn1, vyn1, utemp1, utemp2);
        for (l = 0; l < n; l++)
            sk4[l] = (utemp2[l] - utemp1[l]) / dt;
        // skladamy rozwiazania dla RK4
        for (l = 0; l < n; l++)
            ubn1[l] = ubn[l] + dt / 6 * (sk1[l] + 2 * sk2[l] + 2 * sk3[l] + sk4[l]);
        ////////////////////////////////////////////////////////////////////////////
        cr = 0.;
        cg = 0.;
        cb = 0.;
        // przepisywanie wartości
        for (l = 0; l < n; l++) {
            urn[l] = urn1[l];
            ugn[l] = ugn1[l];
            ubn[l] = ubn1[l];
            vxn[l] = vxn1[l];
            vyn[l] = vyn1[l];
            cr = cr + urn[l] * dx * dx;
            cg = cg + ugn[l] * dx * dx;
            cb = cb + ubn[l] * dx * dx;
        }
        if (it % int(itmax / 10) == 0) {
            ++snapshotCounter;
            std::string filename = std::string("FCT") + "_" + std::to_string(snapshotCounter) + ".dat";
            std::cout << "Zapisywanie danych do pliku: " << filename << std::endl;
            utils::saveData2DRGB(filename, urn, ugn, ubn, ny, nx);
        }
        if (it % 20 == 0)
            printf("it= %d  dt=  %g  cr= %g cg= %g cb= %g\n", it, dt_cfl, cr, cg, cb);
        // zapis do pliku i obliczanie dywergencji
        if (it % 200 == 0) {
            fp = fopen("psi.dat", "w");
            for (i = 0; i < nx; i++) {
                for (j = 0; j < ny; j++) {
                    l = j + i * ny;
                    div = 0.;
                    if ((i > 0) && (i < (nx - 1)))
                        div = div + (vxn1[l + ny] - vxn1[l - ny]) / 2 / dx;
                    if ((j > 0) && (j < (ny - 1)))
                        div = div + (vyn1[l + 1] - vyn1[l - 1]) / 2 / dx;
                    fprintf(fp, "%d  %d  %f  %f   %f   %f   %f    %f\n", i, j, zeta[l],
                        psi[l], vxn1[l], vyn1[l], div, urn[l]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            imnd::push_data2D(urn, nx, ny);
        }
        //=====================================================================================================

        // //nie dziala
        // //......liczymy Bi-CG-stab...........................
        // macierz_A_B(n, nx, ny, dt, dx, vxn, vyn, vxn1, vyn1, a, jka, iwa, b, jkb,
        // iwb);
        // //......B_n*u_n......................................
        // licz_axb(n, b, jkb, iwb, urn, vecb);
        // for (l = 0; l < n; l++)
        //     urn1[i] = 0; //urn[i];//przyblizenie startowe
        // bicgstab(n, a, jka, iwa, urn1, vecb);
        // c = 0.;
        // for (l = 0; l < n; l++)
        // {
        //     urn[l] = urn1[l];
        //     c = c + urn[l] * dx * dx;
        // }

        // if (it % 40 == 0)
        // {
        //     printf("it= %d   cn= %f  \n", it, c);
        //     imnd::push_data2D(urn, nx, ny);
        // }
    }
    imnd::plot_2d_system("velocityX.png", vxn, nx, ny, dx, dx);
    imnd::write_data2D("density.dat", nx, ny, dx, dx);
    imnd::free_data2D(nx, ny);
    // it
    // trzeba wykasowac pozostale tablice!!!!!!!!!!!!!
    delete[] urn;
    delete[] ugn;
    delete[] ubn;
    delete[] psi;
    delete[] zeta;
    delete[] vxn;
    delete[] vyn;
    delete[] vxn1;
    delete[] vyn1;
    delete[] vecb;

    return 0;
}
